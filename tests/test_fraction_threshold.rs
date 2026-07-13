//! Fraction threshold: gated by `--min-sites`, and OR-combined with the count
//! threshold (count can trigger even when the fraction path is gated out).

mod helpers;
use helpers::*;

const REF: &str = "CACACACACACA"; // 6 CpA C's at 0,2,4,6,8,10

#[test]
fn fraction_is_gated_below_min_sites_and_fires_above() {
    let env = TestEnv::new();
    let reference = RefBuilder::new().contig("chr1", REF);
    let sam = SamBuilder::new()
        .sq("chr1", REF.len())
        // few: 6M covering C@0,2,4 → 2 unconv of 3 sites (frac 0.67, sites < 5).
        .record("few", 0, "chr1", 1, "6M", "CACATA", &q40(6))
        // many: 12M → C@0,2,4,6 unconv, 8,10 converted → 4 of 6 sites.
        .record("many", 0, "chr1", 1, "12M", "CACACACATATA", &q40(12));

    // Count threshold set absurdly high so only the fraction path can fire.
    let recs = run_ok(
        &sam,
        &reference,
        &env,
        &[
            "--max-unconverted-fraction",
            "0.5",
            "--min-sites",
            "5",
            "--max-unconverted-count",
            "1000",
        ],
    );
    for rec in &recs {
        let name = rec.name().unwrap().to_string();
        let tagged = has_tag(rec, *b"XX");
        match name.as_str() {
            "few" => assert!(!tagged, "3 sites < min-sites 5 → fraction gated → not tagged"),
            "many" => assert!(tagged, "6 sites ≥ 5 and 0.67 > 0.5 → tagged"),
            other => panic!("unexpected {other}"),
        }
    }
}

#[test]
fn count_threshold_fires_even_when_fraction_is_gated() {
    // 6M covering C@0,2,4, all unconverted → count 3 ≥ 3 fires even though only
    // 3 sites (< min-sites 5) gate out the fraction path.
    let env = TestEnv::new();
    let reference = RefBuilder::new().contig("chr1", REF);
    let sam =
        SamBuilder::new().sq("chr1", REF.len()).record("r", 0, "chr1", 1, "6M", "CACACA", &q40(6));
    let recs =
        run_ok(&sam, &reference, &env, &["--max-unconverted-fraction", "0.5", "--min-sites", "5"]);
    assert!(has_tag(&recs[0], *b"XX"), "count path (3 ≥ 3) tags regardless of gating");
}
