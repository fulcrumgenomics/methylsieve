//! `--count-tag`: stamp every record with `<TAG>:Z:u/n`, where u is the
//! unconverted count and n the total monitored sites in the `--contexts` subset
//! (the decision's numerator/denominator) for that record's template.

mod helpers;
use helpers::*;

/// Alternating C/A: positions 0,2,..,18 are CpA top-strand C — 10 monitored
/// sites. These short reads sit below the default `--min-sites`, so the count
/// test decides; the tag's denominator is the full 10 either way.
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
fn count_tag_records_unconverted_over_total_on_every_template() {
    let env = TestEnv::new();
    let reference = RefBuilder::new().contig("chr1", REF);
    let mut sam = SamBuilder::new().sq("chr1", REF.len());
    for unconv in [0usize, 3] {
        let seq = read_with_unconverted(unconv);
        sam = sam.record(&format!("r{unconv}"), 0, "chr1", 1, "20M", &seq, &q40(20));
    }

    // On by default — no flag needed.
    let recs = run_ok(&sam, &reference, &env, &[]);
    for rec in &recs {
        let name = rec.name().unwrap().to_string();
        let expected = match name.as_str() {
            "r0" => "0/10", // present even though this template is NOT flagged
            "r3" => "3/10",
            other => panic!("unexpected record {other}"),
        };
        assert_eq!(tag_string(rec, [b'c', b'h']).as_deref(), Some(expected), "{name} ch tag");
    }
}

#[test]
fn count_tag_name_is_configurable() {
    let env = TestEnv::new();
    let reference = RefBuilder::new().contig("chr1", REF);
    let sam = SamBuilder::new().sq("chr1", REF.len()).record(
        "r",
        0,
        "chr1",
        1,
        "20M",
        &read_with_unconverted(2),
        &q40(20),
    );

    let recs = run_ok(&sam, &reference, &env, &["--count-tag", "xy"]);
    assert_eq!(tag_string(&recs[0], [b'x', b'y']).as_deref(), Some("2/10"));
    assert!(!has_tag(&recs[0], [b'c', b'h']), "default name not used when overridden");
}

#[test]
fn count_tag_disabled_with_no_count_tag() {
    let env = TestEnv::new();
    let reference = RefBuilder::new().contig("chr1", REF);
    let sam = SamBuilder::new().sq("chr1", REF.len()).record(
        "r",
        0,
        "chr1",
        1,
        "20M",
        &read_with_unconverted(3),
        &q40(20),
    );

    let recs = run_ok(&sam, &reference, &env, &["--no-count-tag"]);
    assert!(!has_tag(&recs[0], [b'c', b'h']), "--no-count-tag suppresses the count tag");
}
