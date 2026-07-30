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
        assert_eq!(tag_string(rec, *b"ch").as_deref(), Some(expected), "{name} ch tag");
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
    assert_eq!(tag_string(&recs[0], *b"xy").as_deref(), Some("2/10"));
    assert!(!has_tag(&recs[0], *b"ch"), "default name not used when overridden");
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
    assert!(!has_tag(&recs[0], *b"ch"), "--no-count-tag suppresses the count tag");
}

#[test]
fn count_tag_is_appended_alongside_pre_existing_aux() {
    // Every other test here stamps a record with no aux at all. This one already
    // carries fields, so it covers appending to a populated aux block: the new tag
    // lands and the existing fields survive intact.
    let env = TestEnv::new();
    let reference = RefBuilder::new().contig("chr1", REF);
    let sam = SamBuilder::new().sq("chr1", REF.len()).record_with_aux(
        "r",
        0,
        "chr1",
        1,
        "20M",
        &read_with_unconverted(3),
        &q40(20),
        &["XA:Z:chr1,+100,20M,0;", "MD:Z:20"],
    );

    let recs = run_ok(&sam, &reference, &env, &[]);
    assert_eq!(tag_string(&recs[0], *b"ch").as_deref(), Some("3/10"));
    assert_eq!(tag_string(&recs[0], *b"XA").as_deref(), Some("chr1,+100,20M,0;"));
    assert_eq!(tag_string(&recs[0], *b"MD").as_deref(), Some("20"));
}

#[test]
fn stale_count_tag_from_a_previous_run_is_overwritten() {
    // A tag may appear only once per record, so an inherited value has to be
    // replaced rather than appended-beside or left alone. Leaving it alone would
    // make the count disagree with the decision stamped next to it — the counts
    // legitimately change between runs, since masking drops bases below the
    // base-quality gate.
    let env = TestEnv::new();
    let reference = RefBuilder::new().contig("chr1", REF);
    let sam = SamBuilder::new().sq("chr1", REF.len()).record_with_aux(
        "r",
        0,
        "chr1",
        1,
        "20M",
        &read_with_unconverted(3),
        &q40(20),
        &["ch:Z:99/99"],
    );

    let recs = run_ok(&sam, &reference, &env, &[]);
    assert_eq!(tag_string(&recs[0], *b"ch").as_deref(), Some("3/10"), "stale value replaced");
}

#[test]
fn count_tag_name_already_used_with_another_type_is_replaced() {
    // Nothing stops another tool from defining `ch` — and picking a different
    // type. Appending ours beside it would emit two `ch` fields, which is invalid.
    let env = TestEnv::new();
    let reference = RefBuilder::new().contig("chr1", REF);
    let sam = SamBuilder::new().sq("chr1", REF.len()).record_with_aux(
        "r",
        0,
        "chr1",
        1,
        "20M",
        &read_with_unconverted(3),
        &q40(20),
        &["ch:i:5", "MD:Z:20"],
    );

    let recs = run_ok(&sam, &reference, &env, &[]);
    assert_eq!(physical_tag_counts(&env.output, *b"ch"), vec![1], "exactly one ch field");
    assert_eq!(tag_string(&recs[0], *b"ch").as_deref(), Some("3/10"));
    assert_eq!(tag_string(&recs[0], *b"MD").as_deref(), Some("20"), "later tags survive");
}

/// 10 bp of the same alternating reference: 5 monitored sites, so the count
/// renders as three characters. That width matters — see the test below.
const REF10: &str = "CACACACACA";

/// C at the first `unconv` monitored positions of [`REF10`], T at the rest.
fn read10_with_unconverted(unconv: usize) -> String {
    (0..10)
        .map(|i| {
            if i % 2 == 1 {
                'A'
            } else if i / 2 < unconv {
                'C'
            } else {
                'T'
            }
        })
        .collect()
}

#[test]
fn count_tag_replacing_a_four_byte_foreign_tag_keeps_its_string_type() {
    // The dangerous width coincidence: a four-byte `i`/`I`/`f` field occupies the
    // same bytes a three-character `Z` value would, and a replacement that only
    // compares those lengths will overwrite the payload while leaving the type byte
    // alone — destroying the value instead of replacing it. Three characters is not
    // exotic: it is what `0/0` renders, on every zero-site template.
    let env = TestEnv::new();
    let reference = RefBuilder::new().contig("chr1", REF10);
    let sam = SamBuilder::new().sq("chr1", REF10.len()).record_with_aux(
        "r",
        0,
        "chr1",
        1,
        "10M",
        &read10_with_unconverted(2),
        &q40(10),
        &["ch:i:100000", "MD:Z:10"], // 100000 needs the full four bytes
    );

    let recs = run_ok(&sam, &reference, &env, &[]);
    assert_eq!(physical_tag_counts(&env.output, *b"ch"), vec![1], "exactly one ch field");
    assert_eq!(
        tag_string(&recs[0], *b"ch").as_deref(),
        Some("2/5"),
        "count tag is a Z string, not the old field's type with new bytes in it"
    );
}
