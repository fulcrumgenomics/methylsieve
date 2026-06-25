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
    assert_eq!(g["CpA_obs"], "3");
    assert_eq!(ctx_unconv(&g, "CpA"), 3);
    assert_eq!(g["CpG_obs"], "2");
    assert_eq!(ctx_unconv(&g, "CpG"), 2);
    // Fully unconverted → conversion rate 0 for contexts with sites, and CpG
    // methylation rate 1.
    assert_eq!(g["CpA_conv_rate"], "0.000000");
    assert_eq!(g["CpG_conv_rate"], "0.000000");
    assert_eq!(g["CpG_meth_rate"], "1.000000");
    // No CpC or CpT sites → empty conversion-rate cells.
    assert_eq!(g["CpC_obs"], "0");
    assert_eq!(g["CpC_conv_rate"], "");
    assert_eq!(g["CpT_conv_rate"], "");
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
    for suffix in ["summary.tsv", "mbias.tsv", "conversion-matrix.tsv"] {
        let p = std::path::PathBuf::from(format!("{pfx}.{suffix}"));
        assert!(p.exists(), "expected metrics file {}", p.display());
    }
    // The separate mbias-bounds file was folded into the summary.
    assert!(
        !std::path::PathBuf::from(format!("{pfx}.mbias_bounds.tsv")).exists(),
        "mbias_bounds.tsv should no longer be written"
    );
}
