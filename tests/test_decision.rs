//! Count-threshold decision: with the default `--max-unconverted-count 3`,
//! templates with 0/2 unconverted CpH pass; 3/4 are tagged + QC-failed.

mod helpers;
use helpers::*;

/// Reference: alternating C/A so every monitored top-strand C is in CpA
/// context: positions 0,2,4,...,18 are C (10 sites), each followed by A.
const REF: &str = "CACACACACACACACACACA";

/// A forward read identical to the reference except the first `unconv` of the
/// ten C-positions keep their C (unconverted); the rest read T (converted).
fn read_with_unconverted(unconv: usize) -> String {
    let mut s = String::with_capacity(20);
    for i in 0..20 {
        if i % 2 == 1 {
            s.push('A');
        } else if (i / 2) < unconv {
            s.push('C');
        } else {
            s.push('T');
        }
    }
    s
}

#[test]
fn count_threshold_tags_at_three_and_four_not_zero_or_two() {
    let env = TestEnv::new();
    let reference = RefBuilder::new().contig("chr1", REF);
    let mut sam = SamBuilder::new().sq("chr1", REF.len());
    for unconv in [0usize, 2, 3, 4] {
        let seq = read_with_unconverted(unconv);
        sam = sam.record(&format!("r{unconv}"), 0, "chr1", 1, "20M", &seq, &q40(20));
    }

    let recs = run_ok(&sam, &reference, &env, &[]);
    assert_eq!(recs.len(), 4);

    for rec in &recs {
        let name = rec.name().unwrap().to_string();
        let tagged = has_tag(rec, *b"XX");
        let qc_failed = u16::from(rec.flags()) & FLAG_QC_FAIL != 0;
        match name.as_str() {
            "r0" | "r2" => {
                assert!(!tagged, "{name} should not be tagged");
                assert!(!qc_failed, "{name} should not be QC-failed");
            }
            "r3" | "r4" => {
                assert!(tagged, "{name} should be tagged unconverted");
                assert_eq!(tag_string(rec, *b"XX").as_deref(), Some("UC"));
                assert!(qc_failed, "{name} should be QC-failed");
            }
            other => panic!("unexpected record {other}"),
        }
    }
}
