//! The output BAM ends with the BGZF EOF marker only when the run succeeds, so a
//! run that fails partway never leaves output that passes as complete.

mod helpers;
use helpers::*;

const REF: &str = "CACACACACA";

/// The 28-byte empty BGZF block that marks the end of a complete BAM.
const BGZF_EOF: [u8; 28] = [
    0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02, 0x00,
    0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
];

fn ends_with_eof_marker(path: &std::path::Path) -> bool {
    std::fs::read(path).map(|bytes| bytes.ends_with(&BGZF_EOF)).unwrap_or(false)
}

#[test]
fn successful_run_ends_output_with_eof_marker() {
    let env = TestEnv::new();
    let sam = SamBuilder::new().sq("chr1", REF.len()).record(
        "r",
        0,
        "chr1",
        1,
        "10M",
        "TATATATATA",
        &q40(10),
    );
    let r = run_methylsieve(&sam, &RefBuilder::new().contig("chr1", REF), &env, &[]);
    assert!(r.status_ok, "{}", r.stderr);
    assert!(ends_with_eof_marker(&r.output));
}

#[test]
fn failed_run_leaves_output_without_eof_marker() {
    // The second template has no primary record, so the run fails after the
    // first template has already been written.
    let env = TestEnv::new();
    let sam = SamBuilder::new()
        .sq("chr1", REF.len())
        .record("ok", 0, "chr1", 1, "10M", "TATATATATA", &q40(10))
        .record("orphan", FLAG_SUPPLEMENTARY, "chr1", 1, "10M", "TATATATATA", &q40(10));
    let r = run_methylsieve(&sam, &RefBuilder::new().contig("chr1", REF), &env, &[]);
    assert!(!r.status_ok);
    assert!(r.stderr.contains("no primary alignment"), "{}", r.stderr);
    assert!(!ends_with_eof_marker(&r.output));
}
