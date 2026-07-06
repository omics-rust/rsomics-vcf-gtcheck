use std::io::{BufRead, Write};
use std::path::Path;

use rsomics_common::{Result, RsomicsError};

use crate::{open_vcf_reader, write_all};

/// Format a float in C `printf("%e")` style: 6 decimal places, sign on exponent, minimum 2-digit exponent.
/// E.g. 4.605170e+01, 0.000000e+00
fn fmt_e(v: f64) -> String {
    if v == 0.0 || !v.is_finite() {
        return format!("{v:.6e}")
            .replace("e", "e+0")
            .replace("e+0+", "e+")
            .replace("e+0-", "e-");
    }
    let exp = v.abs().log10().floor() as i32;
    let mantissa = v / 10f64.powi(exp);
    let exp_sign = if exp >= 0 { '+' } else { '-' };
    format!("{mantissa:.6}e{exp_sign}{:02}", exp.unsigned_abs())
}

/// Dosage bitmask matching bcftools vcfgtcheck.c gt_to_dsg()/pl_to_dsg():
/// 0 = missing, 1<<0 = hom-ref, 1<<1 = het, 1<<2 = hom-alt. Two dosages are concordant
/// when their bitmasks share a bit; an uncertain PL (e.g. 0,0,10) sets several bits.
type Dosage = u8;

const DSG_MISSING: Dosage = 0;
const DSG_HOM_REF: Dosage = 1;
const DSG_HET: Dosage = 2;
const DSG_HOM_ALT: Dosage = 4;

pub struct GtcheckArgs {
    pub error_prob: u32,
    pub no_hwe_prob: bool,
    pub homs_only: bool,
    pub use_gt: bool,
}

pub enum GtcheckMode<'a> {
    CrossCheck(&'a Path),
    QueryVsGt { query: &'a Path, gt: &'a Path },
}

/// Which tag a record was read with, or why it was unusable. Mirrors the per-record decision in
/// bcftools set_data(): a site whose chosen tag is not diploid is skipped whole, counted separately
/// for GT vs PL.
#[derive(Clone, Copy, PartialEq)]
enum RowStatus {
    Ok,
    NotDiploid,
    NoData,
}

#[derive(Clone)]
struct Row {
    chrom: String,
    pos: u64,
    is_variant: bool,
    used_gt: bool,
    status: RowStatus,
    /// Alternate allele frequency from the record's GT (bcf_calc_ac), independent of the
    /// comparison tag; falls back to the comparison dosage when no GT is present.
    af: f64,
}

/// Flat-column dosage table. `dsg` and `prob` are row-major: record i, sample j at i*ncols + j.
/// `prob` is populated only under the error model (-E > 0).
#[derive(Clone)]
struct VcfTable {
    samples: Vec<String>,
    dsg: Vec<Dosage>,
    prob: Vec<[f64; 3]>,
    rows: Vec<Row>,
}

impl VcfTable {
    fn nrows(&self) -> usize {
        self.rows.len()
    }
    fn ncols(&self) -> usize {
        self.samples.len()
    }
    fn dsg_row(&self, i: usize) -> &[Dosage] {
        let nc = self.ncols();
        &self.dsg[i * nc..(i + 1) * nc]
    }
    fn prob_row(&self, i: usize) -> &[[f64; 3]] {
        let nc = self.ncols();
        &self.prob[i * nc..(i + 1) * nc]
    }
}

/// GT subfield → (dosage, ploidy). Ploidy is the allele count (separators + 1); a sample is only
/// scored when it is diploid with both alleles present.
#[inline]
fn gt_to_dsg_ploidy(field: &str) -> (Dosage, usize) {
    let ploidy = 1 + field.bytes().filter(|&b| b == b'/' || b == b'|').count();
    if ploidy != 2 {
        return (DSG_MISSING, ploidy);
    }
    let mut it = field.split(['/', '|']);
    let a = it.next().unwrap_or(".");
    let b = it.next().unwrap_or(".");
    if a == "." || b == "." {
        return (DSG_MISSING, ploidy);
    }
    let alt = (a != "0") as usize + (b != "0") as usize;
    (1 << alt, ploidy)
}

/// PL subfield → (dosage, value-count, raw values). Values are needed for the soft likelihood score;
/// a biallelic diploid PL carries exactly three.
#[inline]
fn pl_parse(field: &str) -> (Dosage, usize, Option<[i32; 3]>) {
    let mut vals = [0i32; 3];
    let mut count = 0usize;
    let mut ok = true;
    for (i, s) in field.split(',').enumerate() {
        count += 1;
        if i < 3 {
            match s.parse::<i32>() {
                Ok(v) if v >= 0 => vals[i] = v,
                _ => ok = false,
            }
        }
    }
    if count >= 3 && ok {
        let min = vals[0].min(vals[1]).min(vals[2]);
        let mut d = 0;
        if vals[0] == min {
            d |= DSG_HOM_REF;
        }
        if vals[1] == min {
            d |= DSG_HET;
        }
        if vals[2] == min {
            d |= DSG_HOM_ALT;
        }
        (d, count, Some(vals))
    } else {
        (DSG_MISSING, count, None)
    }
}

/// PL likelihoods → -log of the normalized genotype probabilities, matching bcftools pl_to_prob().
/// pl2prob[i] = 10^(-0.1*i) with PL clamped at 255.
#[inline]
fn pl_to_prob(vals: [i32; 3]) -> [f64; 3] {
    let l = |v: i32| 10f64.powf(-0.1 * (if v >= 255 { 255 } else { v }) as f64);
    let p = [l(vals[0]), l(vals[1]), l(vals[2])];
    let sum = p[0] + p[1] + p[2];
    [-(p[0] / sum).ln(), -(p[1] / sum).ln(), -(p[2] / sum).ln()]
}

/// GT dosage → -log P(geno|dosage) under the phred error model, matching bcftools dsg2prob.
/// `neg_log_e` = -ln(eprob), eprob = 10^(-0.1*E).
#[inline]
fn gt_to_prob(dsg: Dosage, neg_log_e: f64) -> [f64; 3] {
    match dsg {
        DSG_HOM_REF => [0.0, neg_log_e, 2.0 * neg_log_e],
        DSG_HET => [neg_log_e, 0.0, neg_log_e],
        DSG_HOM_ALT => [2.0 * neg_log_e, neg_log_e, 0.0],
        _ => [f64::INFINITY; 3],
    }
}

/// Alternate-allele frequency from a record's GT subfields (bcf_calc_ac over FORMAT/GT), counting
/// every present allele regardless of ploidy.
fn af_from_gt(sample_data: &str, gt_idx: usize) -> f64 {
    let mut alt = 0usize;
    let mut tot = 0usize;
    for sample_str in sample_data.split('\t') {
        let field = sample_str.split(':').nth(gt_idx).unwrap_or(".");
        for allele in field.split(['/', '|']) {
            if allele == "." {
                continue;
            }
            tot += 1;
            if allele != "0" {
                alt += 1;
            }
        }
    }
    if tot == 0 {
        1e-6
    } else {
        alt as f64 / tot as f64
    }
}

fn af_from_dsg(dsg_row: &[Dosage]) -> f64 {
    let mut alt = 0usize;
    let mut tot = 0usize;
    for &d in dsg_row {
        let (a, t) = match d {
            DSG_HOM_REF => (0, 2),
            DSG_HET => (1, 2),
            DSG_HOM_ALT => (2, 2),
            _ => continue,
        };
        alt += a;
        tot += t;
    }
    if tot == 0 {
        1e-6
    } else {
        alt as f64 / tot as f64
    }
}

/// Tag preference: the comparison tag is the first present in a record's FORMAT.
struct ReadOpts {
    /// Preference order; the first tag found in FORMAT is used (["GT","PL"] or ["PL","GT"]).
    order: [&'static str; 2],
    /// -ln(eprob) when the error model is active; None disables soft-probability storage.
    neg_log_e: Option<f64>,
}

/// Reads a VCF into a dosage table. Returns the table and the count of multiallelic records skipped.
fn read_vcf(reader: &mut dyn BufRead, opts: &ReadOpts) -> Result<(VcfTable, usize)> {
    let mut samples: Vec<String> = Vec::new();
    let mut dsg: Vec<Dosage> = Vec::new();
    let mut prob: Vec<[f64; 3]> = Vec::new();
    let mut rows: Vec<Row> = Vec::new();
    let mut n_skipped_multiallelic = 0usize;

    let mut line = String::new();
    loop {
        line.clear();
        let n = reader.read_line(&mut line).map_err(RsomicsError::Io)?;
        if n == 0 {
            break;
        }
        let line = line.trim_end_matches('\n').trim_end_matches('\r');

        if line.starts_with("##") {
            continue;
        }
        if line.starts_with('#') {
            let cols: Vec<&str> = line.split('\t').collect();
            if cols.len() > 9 {
                samples = cols[9..].iter().map(|s| s.to_string()).collect();
            }
            continue;
        }

        let mut col_iter = line.splitn(10, '\t');
        let chrom = col_iter.next().unwrap_or("");
        let pos_str = col_iter.next().unwrap_or("");
        let _ = col_iter.next(); // ID
        let _ = col_iter.next(); // REF
        let alt_allele = col_iter.next().unwrap_or(".");
        let _ = col_iter.next(); // QUAL
        let _ = col_iter.next(); // FILTER
        let _ = col_iter.next(); // INFO
        let fmt_str = col_iter.next().unwrap_or("");
        let sample_data = col_iter.next().unwrap_or("");

        let pos: u64 = pos_str
            .parse()
            .map_err(|_| RsomicsError::InvalidInput(format!("bad POS: {pos_str}")))?;

        let n_alleles = if alt_allele == "." {
            1usize
        } else {
            alt_allele.split(',').count() + 1
        };
        if n_alleles > 2 {
            n_skipped_multiallelic += 1;
            continue;
        }
        let is_variant = alt_allele != "." && !alt_allele.is_empty();

        let ncols = samples.len();
        let fmt_cols: Vec<&str> = fmt_str.split(':').collect();
        let gt_idx = fmt_cols.iter().position(|&f| f == "GT");

        let chosen = opts
            .order
            .iter()
            .find_map(|&t| fmt_cols.iter().position(|&f| f == t).map(|i| (t, i)));

        let row_start = dsg.len();
        dsg.resize(row_start + ncols, DSG_MISSING);
        if opts.neg_log_e.is_some() {
            prob.resize(row_start + ncols, [0.0; 3]);
        }

        let (used_gt, status) = match chosen {
            None => (false, RowStatus::NoData),
            Some((tag, idx)) => {
                let used_gt = tag == "GT";
                let expected = if used_gt { 2 } else { 3 };
                let mut max_width = 0usize;
                for (si, sample_str) in sample_data.split('\t').enumerate() {
                    if si >= ncols {
                        break;
                    }
                    let field = sample_str.split(':').nth(idx).unwrap_or(".");
                    if used_gt {
                        let (d, ploidy) = gt_to_dsg_ploidy(field);
                        max_width = max_width.max(ploidy);
                        dsg[row_start + si] = d;
                        if let Some(nle) = opts.neg_log_e {
                            prob[row_start + si] = gt_to_prob(d, nle);
                        }
                    } else {
                        let (d, count, vals) = pl_parse(field);
                        max_width = max_width.max(count);
                        dsg[row_start + si] = d;
                        if opts.neg_log_e.is_some() {
                            prob[row_start + si] = match vals {
                                Some(v) => pl_to_prob(v),
                                None => [f64::INFINITY; 3],
                            };
                        }
                    }
                }
                let status = if max_width == expected {
                    RowStatus::Ok
                } else {
                    RowStatus::NotDiploid
                };
                (used_gt, status)
            }
        };

        let af = match gt_idx {
            Some(gi) => af_from_gt(sample_data, gi),
            None => af_from_dsg(&dsg[row_start..row_start + ncols]),
        };

        rows.push(Row {
            chrom: chrom.to_string(),
            pos,
            is_variant,
            used_gt,
            status,
            af,
        });
    }

    Ok((
        VcfTable {
            samples,
            dsg,
            prob,
            rows,
        },
        n_skipped_multiallelic,
    ))
}

/// -log P(HWE) for a concordant dosage combination; takes the min over the set bits of the shared
/// mask, matching bcftools hwe_dsg[].
fn hwe_neg_log(af: f64, match_mask: Dosage) -> f64 {
    let hwe = [
        -((1.0 - af) * (1.0 - af)).ln(),
        -(2.0 * af * (1.0 - af)).ln(),
        -(af * af).ln(),
    ];
    let mut best = f64::INFINITY;
    for (k, &h) in hwe.iter().enumerate() {
        if match_mask & (1 << k) != 0 && h < best {
            best = h;
        }
    }
    best
}

struct PairStats {
    ndiff: u32,
    pdiff: f64,
    ncnt: u32,
    nmatch: u32,
    hwe_prob: f64,
}

impl PairStats {
    fn new() -> Self {
        Self {
            ndiff: 0,
            pdiff: 0.0,
            ncnt: 0,
            nmatch: 0,
            hwe_prob: 0.0,
        }
    }
}

pub fn run_gtcheck(mode: &GtcheckMode<'_>, args: &GtcheckArgs, out: &mut impl Write) -> Result<()> {
    let use_error_model = args.error_prob > 0;
    let neg_log_e = if use_error_model {
        // -ln(eprob), eprob = 10^(-0.1*E)
        Some(0.1 * args.error_prob as f64 * std::f64::consts::LN_10)
    } else {
        None
    };

    // Query defaults to PL then GT; the -g panel defaults to GT then PL. --use-gt forces GT first.
    let qry_order: [&str; 2] = if args.use_gt {
        ["GT", "PL"]
    } else {
        ["PL", "GT"]
    };
    let gt_order: [&str; 2] = ["GT", "PL"];

    let (qry_tbl, skip_multi_qry, gt_tbl, skip_multi_gt, cross_check) = match mode {
        GtcheckMode::CrossCheck(path) => {
            let mut rdr = open_vcf_reader(path)?;
            let (tbl, skip) = read_vcf(
                &mut *rdr,
                &ReadOpts {
                    order: qry_order,
                    neg_log_e,
                },
            )?;
            let tbl2 = tbl.clone();
            (tbl, skip, tbl2, 0, true)
        }
        GtcheckMode::QueryVsGt { query, gt } => {
            let mut qrdr = open_vcf_reader(query)?;
            let (qt, qskip) = read_vcf(
                &mut *qrdr,
                &ReadOpts {
                    order: qry_order,
                    neg_log_e,
                },
            )?;
            let mut grdr = open_vcf_reader(gt)?;
            let (gt_t, gskip) = read_vcf(
                &mut *grdr,
                &ReadOpts {
                    order: gt_order,
                    neg_log_e,
                },
            )?;
            (qt, qskip, gt_t, gskip, false)
        }
    };

    let nqry = qry_tbl.ncols();
    let ngt = gt_tbl.ncols();

    if nqry == 0 {
        return Err(RsomicsError::InvalidInput("no samples in query VCF".into()));
    }

    let npairs = if cross_check {
        nqry * (nqry + 1) / 2
    } else {
        nqry * ngt
    };
    let mut pairs_stats: Vec<PairStats> = (0..npairs).map(|_| PairStats::new()).collect();

    let mut ncmp: u32 = 0;
    let nskip_no_match: u32 = 0;
    let nskip_not_ba: u32 = skip_multi_qry as u32 + skip_multi_gt as u32;
    let mut nskip_mono: u32 = 0;
    let mut nskip_no_data: u32 = 0;
    // Indexed by used_gt: [0] = PL-not-diploid, [1] = GT-not-diploid.
    let mut nskip_dip = [0u32; 2];
    // nused[qry_used_gt][gt_used_gt], mirroring bcftools nused[2][2].
    let mut nused = [[0u32; 2]; 2];

    // Two-file mode maps each query row to the matching -g row by (chrom, pos).
    let paired_gt_rows: Vec<Option<usize>> = if cross_check {
        (0..qry_tbl.nrows()).map(Some).collect()
    } else {
        let mut out_idx = vec![None; qry_tbl.nrows()];
        let mut gi = 0usize;
        for (qi, slot) in out_idx.iter_mut().enumerate() {
            let qchrom = &qry_tbl.rows[qi].chrom;
            let qpos = qry_tbl.rows[qi].pos;
            while gi < gt_tbl.nrows() {
                let gchrom = &gt_tbl.rows[gi].chrom;
                let gpos = gt_tbl.rows[gi].pos;
                if gchrom == qchrom && gpos == qpos {
                    *slot = Some(gi);
                    gi += 1;
                    break;
                }
                if (gchrom.as_str(), gpos) < (qchrom.as_str(), qpos) {
                    gi += 1;
                } else {
                    break;
                }
            }
        }
        out_idx
    };

    for (rec_idx, paired_gi) in paired_gt_rows.iter().enumerate() {
        let Some(gi) = paired_gi else {
            continue;
        };
        let qrow = &qry_tbl.rows[rec_idx];
        let grow = &gt_tbl.rows[*gi];

        // Skip order follows bcftools: monoallelic (both sides REF), then per-file no-data / not-diploid.
        if !qrow.is_variant && !grow.is_variant {
            nskip_mono += 1;
            continue;
        }
        match qrow.status {
            RowStatus::NoData => {
                nskip_no_data += 1;
                continue;
            }
            RowStatus::NotDiploid => {
                nskip_dip[qrow.used_gt as usize] += 1;
                continue;
            }
            RowStatus::Ok => {}
        }
        if !cross_check {
            match grow.status {
                RowStatus::NoData => {
                    nskip_no_data += 1;
                    continue;
                }
                RowStatus::NotDiploid => {
                    nskip_dip[grow.used_gt as usize] += 1;
                    continue;
                }
                RowStatus::Ok => {}
            }
        }

        ncmp += 1;
        nused[qrow.used_gt as usize][grow.used_gt as usize] += 1;

        // HWE allele frequency comes from the -g panel record when present, else the query record.
        let af = if cross_check { qrow.af } else { grow.af };

        let qdsg = qry_tbl.dsg_row(rec_idx);
        let gdsg = gt_tbl.dsg_row(*gi);
        let qprob = if use_error_model {
            qry_tbl.prob_row(rec_idx)
        } else {
            &[]
        };
        let gprob = if use_error_model {
            gt_tbl.prob_row(*gi)
        } else {
            &[]
        };

        if cross_check {
            // Lower triangle: pair (i>j) packed at idx = i*(i-1)/2 + j.
            for (i, &dsg_i) in qdsg.iter().enumerate() {
                if dsg_i == DSG_MISSING {
                    continue;
                }
                let base_idx = i * i.wrapping_sub(1) / 2;
                for (j, &dsg_j) in qdsg[..i].iter().enumerate() {
                    if dsg_j == DSG_MISSING {
                        continue;
                    }
                    let ps = &mut pairs_stats[base_idx + j];
                    let match_mask = dsg_i & dsg_j;
                    if use_error_model {
                        let pi = &qprob[i];
                        let pj = &qprob[j];
                        let min = (pi[0] + pj[0]).min(pi[1] + pj[1]).min(pi[2] + pj[2]);
                        ps.pdiff += min;
                        if !args.no_hwe_prob && match_mask != 0 {
                            ps.hwe_prob += hwe_neg_log(af, match_mask);
                            ps.nmatch += 1;
                        }
                    } else if match_mask == 0 {
                        ps.ndiff += 1;
                    } else if !args.no_hwe_prob {
                        ps.hwe_prob += hwe_neg_log(af, match_mask);
                        ps.nmatch += 1;
                    }
                    ps.ncnt += 1;
                }
            }
        } else {
            // Two-file mode: all qry × gt pairs at idx = i*ngt + j.
            for (i, &dsg_i) in qdsg.iter().enumerate() {
                if dsg_i == DSG_MISSING {
                    continue;
                }
                let base_idx = i * ngt;
                for (j, &dsg_j) in gdsg.iter().enumerate() {
                    if dsg_j == DSG_MISSING {
                        continue;
                    }
                    if args.homs_only && dsg_j & (DSG_HOM_REF | DSG_HOM_ALT) == 0 {
                        continue;
                    }
                    let ps = &mut pairs_stats[base_idx + j];
                    let match_mask = dsg_i & dsg_j;
                    if use_error_model {
                        let pi = &qprob[i];
                        let pj = &gprob[j];
                        let min = (pi[0] + pj[0]).min(pi[1] + pj[1]).min(pi[2] + pj[2]);
                        ps.pdiff += min;
                        if !args.no_hwe_prob && match_mask != 0 {
                            ps.hwe_prob += hwe_neg_log(af, match_mask);
                            ps.nmatch += 1;
                        }
                    } else if match_mask == 0 {
                        ps.ndiff += 1;
                    } else if !args.no_hwe_prob {
                        ps.hwe_prob += hwe_neg_log(af, match_mask);
                        ps.nmatch += 1;
                    }
                    ps.ncnt += 1;
                }
            }
        }
    }

    write_all(out, format!("INFO\tsites-compared\t{ncmp}\n").as_bytes())?;
    write_all(
        out,
        format!("INFO\tsites-skipped-no-match\t{nskip_no_match}\n").as_bytes(),
    )?;
    write_all(
        out,
        format!("INFO\tsites-skipped-multiallelic\t{nskip_not_ba}\n").as_bytes(),
    )?;
    write_all(
        out,
        format!("INFO\tsites-skipped-monoallelic\t{nskip_mono}\n").as_bytes(),
    )?;
    write_all(
        out,
        format!("INFO\tsites-skipped-no-data\t{nskip_no_data}\n").as_bytes(),
    )?;
    write_all(
        out,
        format!("INFO\tsites-skipped-GT-not-diploid\t{}\n", nskip_dip[1]).as_bytes(),
    )?;
    write_all(
        out,
        format!("INFO\tsites-skipped-PL-not-diploid\t{}\n", nskip_dip[0]).as_bytes(),
    )?;
    write_all(out, b"INFO\tsites-skipped-filtering-expression\t0\n")?;
    write_all(
        out,
        format!("INFO\tsites-used-PL-vs-PL\t{}\n", nused[0][0]).as_bytes(),
    )?;
    write_all(
        out,
        format!("INFO\tsites-used-PL-vs-GT\t{}\n", nused[0][1]).as_bytes(),
    )?;
    write_all(
        out,
        format!("INFO\tsites-used-GT-vs-PL\t{}\n", nused[1][0]).as_bytes(),
    )?;
    write_all(
        out,
        format!("INFO\tsites-used-GT-vs-GT\t{}\n", nused[1][1]).as_bytes(),
    )?;
    write_all(
        out,
        b"# DCv2, discordance version 2:\n\
          #     - Query sample\n\
          #     - Genotyped sample\n\
          #     - Discordance, given either as an abstract score or number of mismatches, see the options -E/-u\n\
          #       in man page for details. Note that samples with high missingness have fewer sites compared,\n\
          #       which results in lower overall discordance. Therefore it is advisable to use the average score\n\
          #       per site rather than the absolute value, i.e. divide the value by the number of sites compared\n\
          #       (smaller value = better match)\n\
          #     - Average negative log of HWE probability at matching sites, attempts to quantify the following\n\
          #       intuition: rare genotype matches are more informative than common genotype matches, hence two\n\
          #       samples with similar discordance can be further stratified by the HWE score (bigger value = better\n\
          #       match, the observed concordance was less likely to occur by chance)\n\
          #     - Number of sites compared for this pair of samples (bigger = more informative)\n\
          #     - Number of matching genotypes\n\
          #DCv2\t[2]Query Sample\t[3]Genotyped Sample\t[4]Discordance\t[5]Average -log P(HWE)\t[6]Number of sites compared\t[7]Number of matching genotypes\n",
    )?;

    let qry_samples = &qry_tbl.samples;
    let gt_samples = &gt_tbl.samples;

    let write_row = |out: &mut dyn Write, qname: &str, gname: &str, ps: &PairStats| -> Result<()> {
        let hwe_avg = if !args.no_hwe_prob && ps.nmatch > 0 {
            ps.hwe_prob / ps.nmatch as f64
        } else {
            0.0
        };
        let nmatch_col = if args.no_hwe_prob { 0 } else { ps.nmatch };
        let disc = if use_error_model {
            fmt_e(ps.pdiff)
        } else {
            ps.ndiff.to_string()
        };
        write_all(
            out,
            format!(
                "DCv2\t{qname}\t{gname}\t{disc}\t{}\t{}\t{nmatch_col}\n",
                fmt_e(hwe_avg),
                ps.ncnt,
            )
            .as_bytes(),
        )
    };

    if cross_check {
        let mut idx = 0usize;
        for (i, sname_i) in qry_samples.iter().enumerate().skip(1) {
            for sname_j in &qry_samples[..i] {
                write_row(out, sname_i, sname_j, &pairs_stats[idx])?;
                idx += 1;
            }
        }
    } else {
        for (i, qsname) in qry_samples.iter().enumerate() {
            for (j, gtsname) in gt_samples.iter().enumerate() {
                write_row(out, qsname, gtsname, &pairs_stats[i * ngt + j])?;
            }
        }
    }

    Ok(())
}
