//! Reference / `@SQ` cross-checks: FASTA-superset is allowed; a missing contig
//! or a length mismatch is fatal.

mod helpers;
use helpers::*;

const SEQ10: &str = "CACACACACA";

#[test]
fn fasta_superset_of_bam_contigs_is_allowed() {
    // FASTA has chr1 + chr2; the BAM only references chr1. Extra FASTA contigs
    // are ignored.
    let env = TestEnv::new();
    let reference = RefBuilder::new().contig("chr1", SEQ10).contig("chr2", "GGGGGGGGGG");
    let sam =
        SamBuilder::new().sq("chr1", SEQ10.len()).record("r", 0, "chr1", 1, "10M", SEQ10, &q40(10));
    let r = run_methylsieve(&sam, &reference, &env, &[]);
    assert!(r.status_ok, "FASTA superset should be fine: {}", r.stderr);
}

#[test]
fn missing_contig_is_fatal() {
    // BAM references chr2 but the FASTA lacks it.
    let env = TestEnv::new();
    let reference = RefBuilder::new().contig("chr1", SEQ10);
    let sam = SamBuilder::new().sq("chr1", SEQ10.len()).sq("chr2", SEQ10.len()).record(
        "r",
        0,
        "chr2",
        1,
        "10M",
        SEQ10,
        &q40(10),
    );
    let r = run_methylsieve(&sam, &reference, &env, &[]);
    assert!(!r.status_ok);
    assert!(r.stderr.contains("not present"), "expected missing-contig error; got: {}", r.stderr);
}

#[test]
fn length_mismatch_is_fatal() {
    // BAM @SQ says chr1 is 20 bp; the FASTA contig is 10 bp.
    let env = TestEnv::new();
    let reference = RefBuilder::new().contig("chr1", SEQ10); // 10 bp
    let sam = SamBuilder::new().sq("chr1", 20).record("r", 0, "chr1", 1, "10M", SEQ10, &q40(10));
    let r = run_methylsieve(&sam, &reference, &env, &[]);
    assert!(!r.status_ok);
    assert!(
        r.stderr.contains("Length mismatch"),
        "expected length-mismatch error; got: {}",
        r.stderr
    );
}
