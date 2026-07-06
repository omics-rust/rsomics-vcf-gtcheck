use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

fn ours() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_rsomics-vcf-gtcheck"))
}

fn golden(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/golden")
        .join(name)
}

/// Keep only the value-bearing lines shared with bcftools: the INFO stat counters (minus the
/// nondeterministic timing line) and the DCv2 data rows. Drops the produced-by banner, the DCv2
/// column legend, and blank comment lines that differ across tools/versions.
fn comparable(output: &str) -> String {
    let mut lines: Vec<&str> = output
        .lines()
        .filter(|l| {
            (l.starts_with("INFO\t") && !l.starts_with("INFO\tTime")) || l.starts_with("DCv2\t")
        })
        .collect();
    // The golden files end without a trailing blank; normalise by joining with '\n'.
    lines.retain(|l| !l.is_empty());
    lines.join("\n")
}

fn run_ours(args: &[&str]) -> String {
    let out = Command::new(ours())
        .args(args)
        .output()
        .expect("spawn rsomics-vcf-gtcheck");
    assert!(
        out.status.success(),
        "rsomics-vcf-gtcheck exited non-zero for args {args:?}\nstderr: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    String::from_utf8(out.stdout).expect("utf8 stdout")
}

/// Run ours against a committed input and assert its comparable output equals the golden captured
/// from bcftools 1.23.1. No live oracle at test time — CI needs no bcftools.
fn check(golden_name: &str, args: &[&str]) {
    let expected = std::fs::read_to_string(golden(golden_name))
        .unwrap_or_else(|e| panic!("read golden {golden_name}: {e}"));
    let got = run_ours(args);
    assert_eq!(
        comparable(&got),
        comparable(&expected),
        "\n=== {golden_name} mismatch ===\nOURS:\n{got}\nEXPECTED (bcftools golden):\n{expected}"
    );
}

fn p(name: &str) -> String {
    golden(name).to_str().unwrap().to_string()
}

// Default query tag is PL: uncertain PLs must drive a soft likelihood-weighted discordance
// (2.809369e+00), not a hard argmax (0). Provenance must read PL-vs-PL, not GT-vs-GT.
#[test]
fn pl_soft_crosscheck_default() {
    check("pl_default.expected", &[&p("pl.vcf")]);
}

#[test]
fn pl_hard_crosscheck_e0() {
    check("pl_e0.expected", &["-E", "0", &p("pl.vcf")]);
}

// -g panel: the HWE column uses the panel allele frequency (7.253123e-01), not the query's.
#[test]
fn panel_hwe_uses_panel_af_e0() {
    check(
        "panel_e0.expected",
        &["-g", &p("panel.vcf"), "-E", "0", &p("query.vcf")],
    );
}

#[test]
fn panel_default_soft() {
    check(
        "panel_default.expected",
        &["-g", &p("panel.vcf"), &p("query.vcf")],
    );
}

// PL query vs GT panel: provenance must report PL-vs-GT and mix the soft models correctly.
#[test]
fn mixed_pl_query_gt_panel() {
    check(
        "mixed_default.expected",
        &["-g", &p("panel.vcf"), &p("qpl.vcf")],
    );
}

// Haploid GT sites are skipped whole: sites-compared excludes them and sites-skipped-GT-not-diploid
// counts them.
#[test]
fn haploid_sites_skipped_e0() {
    check("hap_e0.expected", &["--use-gt", "-E", "0", &p("hap.vcf")]);
}

#[test]
fn haploid_sites_skipped_default() {
    check("hap_default.expected", &["--use-gt", &p("hap.vcf")]);
}

// The valid GT-vs-GT cross-check stays value-exact in both the raw-count and the soft-score modes.
#[test]
fn gt_crosscheck_e0_no_hwe() {
    check(
        "gt3_e0_nohwe.expected",
        &["--use-gt", "-E", "0", "--no-HWE-prob", &p("gt3.vcf")],
    );
}

#[test]
fn gt_crosscheck_default_soft() {
    check("gt3_default.expected", &["--use-gt", &p("gt3.vcf")]);
}

fn expect_fail(args: &[&str]) {
    let out = Command::new(ours())
        .args(args)
        .stdout(Stdio::null())
        .output()
        .expect("spawn");
    assert!(
        !out.status.success(),
        "expected non-zero exit for args {args:?}"
    );
    assert!(
        !out.stderr.is_empty(),
        "expected non-empty stderr for args {args:?}"
    );
}

#[test]
fn fails_loud_on_missing_file() {
    expect_fail(&["/nonexistent/does-not-exist.vcf"]);
}

#[test]
fn fails_loud_on_no_samples() {
    let dir = tempfile::tempdir().unwrap();
    let vcf = dir.path().join("nosamp.vcf");
    std::fs::write(
        &vcf,
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
    )
    .unwrap();
    expect_fail(&[vcf.to_str().unwrap()]);
}

/// Optional live differential against bcftools when it is installed (e.g. on the oracle host);
/// re-derives every golden and asserts ours still matches the freshly produced answer.
#[test]
fn live_oracle_matches_goldens() {
    let have_bcftools = Command::new("bcftools")
        .arg("--version")
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .status()
        .map(|s| s.success())
        .unwrap_or(false);
    if !have_bcftools {
        eprintln!("skipping live oracle: bcftools not found");
        return;
    }
    // Cross-check reads plain VCF directly; the -g cases need bgzip+index for bcftools.
    let dir = tempfile::tempdir().unwrap();
    let bgzip_ok = Command::new("bgzip").arg("--version").output().is_ok();
    let mut gz = std::collections::HashMap::new();
    if bgzip_ok {
        for name in ["query", "panel", "qpl"] {
            let src = golden(&format!("{name}.vcf"));
            let dst = dir.path().join(format!("{name}.vcf"));
            std::fs::copy(&src, &dst).unwrap();
            let z = Command::new("bgzip").arg("-f").arg(&dst).status().unwrap();
            assert!(z.success());
            let gzp = dir.path().join(format!("{name}.vcf.gz"));
            let t = Command::new("tabix")
                .args(["-f", "-p", "vcf"])
                .arg(&gzp)
                .status()
                .unwrap();
            assert!(t.success());
            gz.insert(name, gzp.to_str().unwrap().to_string());
        }
    }

    let norm = |s: &str| comparable(s);
    let bcf = |args: &[&str]| -> String {
        let out = Command::new("bcftools")
            .arg("gtcheck")
            .args(args)
            .output()
            .expect("spawn bcftools");
        String::from_utf8(out.stdout).unwrap()
    };

    // Cross-check cases (plain VCF on both sides).
    let cases: &[(&[&str], &[&str])] = &[
        (&[&p("pl.vcf")], &[&p("pl.vcf")]),
        (&["-E", "0", &p("pl.vcf")], &["-E", "0", &p("pl.vcf")]),
        (
            &["--use-gt", "-E", "0", "--no-HWE-prob", &p("gt3.vcf")],
            &["-u", "GT", "-E", "0", "--no-HWE-prob", &p("gt3.vcf")],
        ),
        (&["--use-gt", &p("gt3.vcf")], &["-u", "GT", &p("gt3.vcf")]),
        (
            &["--use-gt", "-E", "0", &p("hap.vcf")],
            &["-u", "GT", "-E", "0", &p("hap.vcf")],
        ),
    ];
    for (our_args, bcf_args) in cases {
        assert_eq!(
            norm(&run_ours(our_args)),
            norm(&bcf(bcf_args)),
            "live cross-check mismatch for {our_args:?}"
        );
    }

    if bgzip_ok {
        let panel = &gz["panel"];
        let query = &gz["query"];
        let qpl = &gz["qpl"];
        let g_cases: &[(&[&str], &[&str])] = &[
            (
                &["-g", panel, "-E", "0", query],
                &["-g", panel, "-E", "0", query],
            ),
            (&["-g", panel, query], &["-g", panel, query]),
            (&["-g", panel, qpl], &["-g", panel, qpl]),
        ];
        for (our_args, bcf_args) in g_cases {
            assert_eq!(
                norm(&run_ours(our_args)),
                norm(&bcf(bcf_args)),
                "live -g mismatch for {our_args:?}"
            );
        }
    }
}
