//! End-to-end `--mbias-mask`: the two-phase buffer → learn → drain → stream
//! pipeline, with Q2 actually written over the learned 5' window in the emitted
//! BAM. Exercises both the single-phase path (whole file buffered, drained at
//! EOF) and the two-phase path (file larger than the buffer).
//!
//! The fixture uses a CpG-dense reference and a designed M-bias ramp: a
//! forward read monitors the top-strand C of every CpG; the first few cycles
//! read `T` (unmethylated) and the rest read `C` (methylated), so the per-cycle
//! CpG methylation rises from 0 to a plateau — exactly the shape the detector
//! masks. The mask-length detector requires ≥100 observations per cycle, so the
//! fixtures generate enough single-end templates to clear that floor.

mod helpers;
use helpers::*;

// CpG at every even reference position (C followed by G); 40 bp / 20 CpG.
const REF: &str = "CGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCG";

/// Forward-read SEQ over `REF`: at the monitored CpG C's (even positions) the
/// first 8 cycles read `T` (unmethylated) and the rest read `C` (methylated);
/// odd positions (ref G) are unmonitored, filled with `A`.
fn ramp_seq() -> String {
    (0..REF.len())
        .map(|i| {
            if i % 2 == 1 {
                'A'
            } else if i < 8 {
                'T'
            } else {
                'C'
            }
        })
        .collect()
}

/// A SAM of `n` single-end forward reads over `REF`, one template each.
fn ramp_sam(n: usize) -> SamBuilder {
    let mut b = SamBuilder::new().sq("chr1", REF.len());
    let seq = ramp_seq();
    let qual = q40(REF.len());
    for i in 0..n {
        b = b.record(&format!("r{i}"), 0, "chr1", 1, "40M", &seq, &qual);
    }
    b
}

#[test]
fn masks_learned_five_prime_window_in_output() {
    // Whole file fits the (default) buffer → learn over all reads, drain at EOF.
    let env = TestEnv::new();
    let reference = RefBuilder::new().contig("chr1", REF);
    let recs = run_ok(&ramp_sam(150), &reference, &env, &["--mbias-mask"]);
    assert_eq!(recs.len(), 150);

    // The ramp clears the plateau around cycle 8, so ~8 leading cycles are Q2.
    let masked = leading_quality_run(&recs[0], 2);
    assert!((5..=15).contains(&masked), "expected ~8 leading Q2 cycles, got {masked}");
    let q = quality_scores(&recs[0]);
    assert_eq!(q[q.len() - 1], 40, "3' end (plateau) is left untouched");
}

#[test]
fn two_phase_masks_both_buffered_and_streamed_reads() {
    // Buffer smaller than the file → learn on the first 120, then stream the
    // rest applying the same frozen masks. 120 reads ≥ the per-cycle floor.
    let env = TestEnv::new();
    let reference = RefBuilder::new().contig("chr1", REF);
    let recs = run_ok(
        &ramp_sam(250),
        &reference,
        &env,
        &["--mbias-mask", "--mbias-buffer-templates", "120"],
    );
    assert_eq!(recs.len(), 250);

    let buffered = leading_quality_run(&recs[0], 2); // from the learn/drain phase
    let streamed = leading_quality_run(&recs[249], 2); // from the stream phase
    assert!((5..=15).contains(&buffered), "buffered read masked, got {buffered}");
    assert_eq!(buffered, streamed, "streamed reads use the same frozen mask length");
}

#[test]
fn masking_off_leaves_qualities_untouched() {
    // Without --mbias-mask, qualities pass through unchanged (no Q2).
    let env = TestEnv::new();
    let reference = RefBuilder::new().contig("chr1", REF);
    let recs = run_ok(&ramp_sam(150), &reference, &env, &[]);
    assert!(
        quality_scores(&recs[0]).iter().all(|&b| b == 40),
        "no masking when --mbias-mask is off"
    );
}
