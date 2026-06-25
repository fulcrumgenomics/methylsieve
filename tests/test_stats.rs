//! `--stats` multi-row TSV: per-context tallies, conversion rates, sample
//! resolution, and empty cells where a context has no monitored sites.

mod helpers;
use helpers::*;

// C@0,2,4,6,8 → CpA, CpG, CpA, CpG, CpA = 3 CpA + 2 CpG. No CpC/CpT sites.
const REF: &str = "CACGCACGCA";

#[test]
fn per_context_counts_and_conv_rates() {
    let env = TestEnv::new();
    let reference = RefBuilder::new().contig("chr1", REF);
    let stats = env.metrics_prefix.to_str().unwrap().to_string();
    let sam = SamBuilder::new().rg("rg1", "sampleX").sq("chr1", REF.len()).record(
        "r",
        0,
        "chr1",
        1,
        "10M",
        REF,
        &q40(10),
    ); // fully unconverted
    run_ok(&sam, &reference, &env, &["--metrics-prefix", &stats]);

    let g = genome_stats(&env.stats);
    assert_eq!(g["sample"], "sampleX", "sample resolved from @RG SM:");
    assert_eq!(g["CA_total"], "3");
    assert_eq!(g["CA_unconv"], "3");
    assert_eq!(g["CG_total"], "2");
    assert_eq!(g["CG_unconv"], "2");
    // Fully unconverted → conversion rate 0 for contexts with sites.
    assert_eq!(g["conv_rate_CpA"], "0.000000");
    assert_eq!(g["conv_rate_CpG"], "0.000000");
    // No CpC or CpT sites → empty conversion-rate cells.
    assert_eq!(g["CC_total"], "0");
    assert_eq!(g["conv_rate_CpC"], "");
    assert_eq!(g["conv_rate_CpT"], "");
}

#[test]
fn sample_override_wins_over_read_group() {
    let env = TestEnv::new();
    let reference = RefBuilder::new().contig("chr1", REF);
    let stats = env.metrics_prefix.to_str().unwrap().to_string();
    let sam = SamBuilder::new().rg("rg1", "sampleX").sq("chr1", REF.len()).record(
        "r",
        0,
        "chr1",
        1,
        "10M",
        REF,
        &q40(10),
    );
    run_ok(&sam, &reference, &env, &["--sample", "forced", "--metrics-prefix", &stats]);
    assert_eq!(genome_stats(&env.stats)["sample"], "forced");
}

#[test]
fn metrics_prefix_writes_all_files() {
    // `--metrics-prefix` produces the summary plus the M-bias / matrix TSVs.
    let env = TestEnv::new();
    let reference = RefBuilder::new().contig("chr1", REF);
    let sam =
        SamBuilder::new().sq("chr1", REF.len()).record("r", 0, "chr1", 1, "10M", REF, &q40(10));
    let r =
        run_methylsieve(&sam, &reference, &env, &["--metrics-prefix", &env.metrics_prefix_arg()]);
    assert!(r.status_ok, "stderr: {}", r.stderr);
    let pfx = env.metrics_prefix.to_str().unwrap();
    for suffix in ["summary.tsv", "mbias.tsv", "mbias_bounds.tsv", "conversion_matrix.tsv"] {
        let p = std::path::PathBuf::from(format!("{pfx}.{suffix}"));
        assert!(p.exists(), "expected metrics file {}", p.display());
    }
}
