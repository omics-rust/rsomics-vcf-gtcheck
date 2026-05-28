use criterion::{Criterion, criterion_group, criterion_main};
use std::hint::black_box;
use std::path::PathBuf;
use std::process::Command;

fn bench_vcf_gtcheck(c: &mut Criterion) {
    let bin = env!("CARGO_BIN_EXE_rsomics-vcf-gtcheck");
    // Use small.vcf from vcf-filter golden (neighbour in the formats tree)
    let manifest = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let vcf = manifest
        .parent()
        .unwrap()
        .join("rsomics-vcf-filter/tests/golden/two.vcf");
    c.bench_function("rsomics-vcf-gtcheck golden", |b| {
        b.iter(|| {
            let out = Command::new(black_box(bin))
                .args(["--use-gt", "-E", "0"])
                .arg(vcf.to_str().unwrap())
                .output()
                .unwrap();
            assert!(out.status.success());
        });
    });
}

criterion_group!(benches, bench_vcf_gtcheck);
criterion_main!(benches);
