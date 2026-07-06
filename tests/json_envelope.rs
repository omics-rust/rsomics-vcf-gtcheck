//! `--json` must emit exactly ONE enveloped document with a populated `result`. A trailing
//! `result: null` (or a second doc) makes `serde_json::from_str` fail with trailing data. No live
//! oracle: spawns the own binary and parses with serde_json.

use std::io::Write;
use std::path::PathBuf;
use std::process::Command;

fn ours() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_rsomics-vcf-gtcheck"))
}

const VCF: &str = "\
##fileformat=VCFv4.2
##contig=<ID=chr1,length=1000>
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Phred likelihoods\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tX\tY
chr1\t100\t.\tA\tG\t50\tPASS\t.\tGT:PL\t0/0:0,3,20\t0/0:0,2,15
chr1\t200\t.\tC\tT\t50\tPASS\t.\tGT:PL\t0/1:5,0,18\t0/1:6,0,20
chr1\t300\t.\tG\tA\t50\tPASS\t.\tGT:PL\t1/1:20,5,0\t0/0:0,4,22
";

fn write_vcf() -> tempfile::TempDir {
    let dir = tempfile::tempdir().expect("tempdir");
    let mut f = std::fs::File::create(dir.path().join("in.vcf")).unwrap();
    f.write_all(VCF.as_bytes()).unwrap();
    dir
}

#[test]
fn json_is_single_populated_envelope() {
    let dir = write_vcf();
    let vcf = dir.path().join("in.vcf");
    let out = Command::new(ours())
        .arg(&vcf)
        .arg("--json")
        .output()
        .expect("spawn rsomics-vcf-gtcheck");
    assert!(
        out.status.success(),
        "non-zero exit\nstderr: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    let stdout = String::from_utf8(out.stdout).expect("utf8");

    // A second trailing document (or text before the envelope) makes this fail with trailing data.
    let env: serde_json::Value =
        serde_json::from_str(&stdout).expect("stdout must be a single JSON document");
    assert_eq!(env["status"], "ok");

    let result = &env["result"];
    assert!(!result.is_null(), "result must be populated, got: {env}");
    assert_eq!(result["sites_compared"], 3);
    let pairs = result["pairs"].as_array().expect("pairs array");
    assert_eq!(pairs.len(), 1);
    assert_eq!(pairs[0]["query_sample"], "Y");
    assert_eq!(pairs[0]["genotyped_sample"], "X");
    assert!(pairs[0]["n_sites_compared"].as_u64().unwrap() >= 1);
}
