//! Control contigs are excluded from the main decision, never tagged, and
//! reported in their own stats scope. A main template with a supplementary on
//! a control contig still counts toward the genome and bumps the chimeric
//! diagnostic.

mod helpers;
use helpers::*;

const SEQ10: &str = "CACACACACA"; // 5 CpA top-strand C's

#[test]
fn control_reads_are_separated_and_chimeric_is_counted() {
    let env = TestEnv::new();
    let reference = RefBuilder::new().contig("chr1", SEQ10).contig("lambda", SEQ10);
    let stats = env.stats.to_str().unwrap().to_string();
    let sam = SamBuilder::new()
        .sq("chr1", SEQ10.len())
        .sq("lambda", SEQ10.len())
        // main genome read, fully unconverted → tagged.
        .record("main", 0, "chr1", 1, "10M", "CACACACACA", &q40(10))
        // control read, fully unconverted → NOT tagged (control scope).
        .record("ctrl", 0, "lambda", 1, "10M", "CACACACACA", &q40(10))
        // chimeric: primary on chr1 (1 unconv) + supplementary on lambda (2 unconv).
        .record("chim", 0, "chr1", 1, "10M", "CATATATATA", &q40(10))
        .record("chim", FLAG_SUPPLEMENTARY, "lambda", 1, "10M", "CACATATATA", &q40(10));

    let recs = run_ok(&sam, &reference, &env, &["--control-contig", "lambda", "--stats", &stats]);

    for rec in &recs {
        let name = rec.name().unwrap().to_string();
        let tagged = has_tag(rec, [b'X', b'X']);
        match name.as_str() {
            "main" | "chim" => assert!(tagged, "{name} should be tagged"),
            "ctrl" => assert!(!tagged, "control read must never be tagged"),
            other => panic!("unexpected {other}"),
        }
    }

    let rows = read_stats_rows(&env.stats);
    assert_eq!(rows.len(), 2, "genome + one control row");
    let genome = &rows[0];
    assert_eq!(genome["scope"], "genome");
    assert_eq!(genome["chimeric_to_control_templates"], "1");
    let lambda = &rows[1];
    assert_eq!(lambda["scope"], "lambda");
    assert_eq!(lambda["n_templates"], "1");
    assert_eq!(lambda["CA_total"], "5", "only the ctrl read tallies into lambda");
    // CA_total counts every monitored CpA site (converted + unconverted). Each
    // 10 bp read covers all 5 ref-C positions, so the chimeric supplementary's
    // sites count toward the genome (not lambda): main(5) + chim primary(5) +
    // chim supp(5) = 15, with 5 + 1 + 2 = 8 of them unconverted.
    assert_eq!(genome["CA_total"], "15");
    assert_eq!(genome["CA_unconv"], "8");
}
