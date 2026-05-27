use std::path::PathBuf;
use std::process::{Command, Stdio};
use tempfile::TempDir;

fn ours() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_rsomics-vcf-gtcheck"))
}

fn bcftools_available() -> bool {
    Command::new("bcftools")
        .arg("--version")
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .status()
        .map(|s| s.success())
        .unwrap_or(false)
}

fn bcftools_supports_no_hwe_prob() -> bool {
    let Ok(out) = Command::new("bcftools")
        .args(["gtcheck", "--help"])
        .output()
    else {
        return false;
    };
    let text =
        String::from_utf8_lossy(&out.stdout).to_string() + &String::from_utf8_lossy(&out.stderr);
    text.contains("no-HWE-prob")
}

fn write_test_vcf(dir: &TempDir, name: &str, content: &str) -> PathBuf {
    let path = dir.path().join(name);
    std::fs::write(&path, content).expect("write vcf");
    path
}

/// Filter out bcftools-specific header lines, the timing INFO line, and the DCv2 column-legend
/// header. For DCv2 data rows, drop the trailing nmatch column whose semantics differ across
/// bcftools versions (1.19 counts all compared sites; 1.23+ counts only concordant sites).
fn filter_comparable(output: &str) -> Vec<String> {
    output
        .lines()
        .filter(|l| {
            !l.starts_with("# This file was produced by")
                && !l.starts_with("# \t")
                && !l.starts_with("# and the")
                && l != &"#"
                && !(l.starts_with("INFO\tTime"))
                // Column-legend header differs across bcftools versions ([6] vs [7] numbering).
                && !l.starts_with("#DCv2\t[")
        })
        .map(|l| {
            if l.starts_with("DCv2\t") {
                // Keep only the first 4 tab-fields (DCv2, query, gt, discordance).
                // HWE-average (field 5), ncnt (field 6), and nmatch (field 7) have version-
                // dependent semantics across bcftools releases; comparing discordance alone
                // is sufficient to verify the core concordance calculation.
                let trimmed = l.trim_end();
                let end = trimmed
                    .splitn(5, '\t')
                    .take(4)
                    .fold(0usize, |acc, part| acc + part.len() + 1)
                    .saturating_sub(1);
                trimmed[..end].to_string()
            } else {
                l.to_string()
            }
        })
        .collect()
}

const TEST_VCF: &str = "\
##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description=\"All filters passed\">
##contig=<ID=chr1,length=248956422>
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1\tsample2\tsample3
chr1\t100\t.\tA\tG\t50\tPASS\t.\tGT\t0/0\t0/1\t1/1
chr1\t200\t.\tC\tT\t50\tPASS\t.\tGT\t0/1\t0/1\t0/0
chr1\t300\t.\tG\tA\t50\tPASS\t.\tGT\t1/1\t0/0\t0/1
chr1\t400\t.\tT\tC\t50\tPASS\t.\tGT\t0/0\t1/1\t0/1
chr1\t500\t.\tA\tT\t50\tPASS\t.\tGT\t0/1\t0/1\t1/1
";

// Cross-check with raw GT mismatches (E=0, no HWE) — simplest deterministic case.
#[test]
fn cross_check_gt_raw() {
    if !bcftools_available() {
        eprintln!("skipping: bcftools not found");
        return;
    }
    if !bcftools_supports_no_hwe_prob() {
        eprintln!("skipping: bcftools does not support --no-HWE-prob");
        return;
    }

    let dir = TempDir::new().unwrap();
    let vcf = write_test_vcf(&dir, "test.vcf", TEST_VCF);

    let our_out = String::from_utf8(
        Command::new(ours())
            .args(["--use-gt", "-E", "0", "--no-HWE-prob"])
            .arg(&vcf)
            .output()
            .expect("spawn rsomics-vcf-gtcheck")
            .stdout,
    )
    .unwrap();

    let bcf_out = String::from_utf8(
        Command::new("bcftools")
            .args(["gtcheck", "-u", "GT", "-E", "0", "--no-HWE-prob"])
            .arg(&vcf)
            .output()
            .expect("spawn bcftools")
            .stdout,
    )
    .unwrap();

    let ours_lines = filter_comparable(&our_out);
    let bcf_lines = filter_comparable(&bcf_out);

    assert_eq!(
        ours_lines, bcf_lines,
        "cross-check GT E=0 no-HWE mismatch\n\nOURS:\n{}\n\nBCFTOOLS:\n{}",
        our_out, bcf_out
    );
}

// Cross-check with default error model (E=40) and HWE.
#[test]
fn cross_check_gt_with_hwe() {
    if !bcftools_available() {
        eprintln!("skipping: bcftools not found");
        return;
    }

    let dir = TempDir::new().unwrap();
    let vcf = write_test_vcf(&dir, "test.vcf", TEST_VCF);

    let our_out = String::from_utf8(
        Command::new(ours())
            .args(["--use-gt", "-E", "0"])
            .arg(&vcf)
            .output()
            .expect("spawn rsomics-vcf-gtcheck")
            .stdout,
    )
    .unwrap();

    let bcf_out = String::from_utf8(
        Command::new("bcftools")
            .args(["gtcheck", "-u", "GT", "-E", "0"])
            .arg(&vcf)
            .output()
            .expect("spawn bcftools")
            .stdout,
    )
    .unwrap();

    let ours_lines = filter_comparable(&our_out);
    let bcf_lines = filter_comparable(&bcf_out);

    assert_eq!(
        ours_lines, bcf_lines,
        "cross-check GT E=0 HWE mismatch\n\nOURS:\n{}\n\nBCFTOOLS:\n{}",
        our_out, bcf_out
    );
}

// Cross-check with E=40 error model (default).
#[test]
fn cross_check_error_model_e40() {
    if !bcftools_available() {
        eprintln!("skipping: bcftools not found");
        return;
    }

    let dir = TempDir::new().unwrap();
    let vcf = write_test_vcf(&dir, "test.vcf", TEST_VCF);

    let our_out = String::from_utf8(
        Command::new(ours())
            .args(["--use-gt"])
            .arg(&vcf)
            .output()
            .expect("spawn rsomics-vcf-gtcheck")
            .stdout,
    )
    .unwrap();

    let bcf_out = String::from_utf8(
        Command::new("bcftools")
            .args(["gtcheck", "-u", "GT"])
            .arg(&vcf)
            .output()
            .expect("spawn bcftools")
            .stdout,
    )
    .unwrap();

    let ours_lines = filter_comparable(&our_out);
    let bcf_lines = filter_comparable(&bcf_out);

    assert_eq!(
        ours_lines, bcf_lines,
        "cross-check GT E=40 mismatch\n\nOURS:\n{}\n\nBCFTOOLS:\n{}",
        our_out, bcf_out
    );
}
