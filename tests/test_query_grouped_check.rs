//! Non-query-grouped input is detected and rejected: a QNAME block containing
//! only a secondary/supplementary record (no primary) is the tell-tale sign.

mod helpers;
use helpers::*;

const REF: &str = "CACACACACA";

#[test]
fn block_without_primary_bails() {
    let env = TestEnv::new();
    let reference = RefBuilder::new().contig("chr1", REF);
    // A lone secondary alignment for QNAME "x" — its primary is supposedly
    // elsewhere in the stream (i.e. the input isn't query-grouped).
    let sam = SamBuilder::new().sq("chr1", REF.len()).record(
        "x",
        FLAG_SECONDARY,
        "chr1",
        1,
        "10M",
        "CACACACACA",
        &q40(10),
    );
    let r = run_methylsieve(&sam, &reference, &env, &[]);
    assert!(!r.status_ok, "should fail on non-query-grouped input");
    assert!(
        r.stderr.contains("query-grouped"),
        "error should mention query-grouping; got: {}",
        r.stderr
    );
}
