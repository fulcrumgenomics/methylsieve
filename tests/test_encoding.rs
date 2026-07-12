//! The reference loads to the same 2-bit store whether or not a `.fai` index is
//! present: the index-driven bulk read and the sequential fallback must produce
//! identical end-to-end results. (The unit tests in `reference.rs` cover the
//! packed bytes directly; this guards the whole run.)

mod helpers;
use helpers::*;

// Mixed-context reference exercising both monitored strands.
const REF: &str = "CACGCACGCATTACGTCACA";

#[test]
fn indexed_and_unindexed_reference_agree() {
    let sam = SamBuilder::new()
        .sq("chr1", REF.len())
        .record("fwd", 0, "chr1", 1, "20M", REF, &q40(20))
        .record("rev", FLAG_REVERSE, "chr1", 1, "20M", REF, &q40(20));

    let indexed = {
        let env = TestEnv::new();
        let stats = env.metrics_prefix_arg();
        run_ok(&sam, &RefBuilder::new().contig("chr1", REF), &env, &["--metrics-prefix", &stats]);
        genome_stats(&env.stats)
    };

    let unindexed = {
        let env = TestEnv::new();
        let stats = env.metrics_prefix_arg();
        run_ok(
            &sam,
            &RefBuilder::new().contig("chr1", REF).without_fai(),
            &env,
            &["--metrics-prefix", &stats],
        );
        genome_stats(&env.stats)
    };

    assert_eq!(indexed, unindexed, "indexed and sequential loads disagree");
    // Sanity: the input actually produced monitored sites.
    assert!(ctx_obs(&indexed, "CpA") + ctx_obs(&indexed, "CpG") > 0, "expected monitored sites");
}
