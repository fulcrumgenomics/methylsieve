//! Per-read-cycle M-bias accumulation and 5'/3' mask-length detection.
//!
//! M-bias is the methylation rate as a function of sequencing cycle (distance
//! from a read's 5' start). End-repair fill-in in bisulfite / EM-seq libraries
//! skews the first cycles — especially read 2 — so we measure the per-cycle CpG
//! methylation rate, find where it first reaches the plateau, and mask the
//! biased 5' (and, for single-end reads, 3') cycles.
//!
//! Counts are bucketed by `(read role, read end, context, cycle)`. Each site is
//! classified with methylsieve's own reference-based call (the same one the
//! tally uses), so a learned mask stays consistent with what the tally excludes.
//!
//! The accumulator is fed only at *matched* monitored sites — never on the
//! per-base reference scan — and only when M-bias output or masking is enabled,
//! so the default fast path pays nothing (see the `SiteSink` plumbing in
//! `classify`/`sieve`).

use crate::reference::Context;

/// Number of `(role, end, context)` cells: 3 roles × 2 ends × 4 contexts.
const N_CELLS: usize = 3 * 2 * 4;

/// Defensive cap on the cycle index, so corrupt/extreme read lengths can't grow
/// a per-cycle vector without bound. Far above any real read length.
const MAX_CYCLE: usize = 100_000;

/// Which read of a template a record belongs to. Orphans (one mate unmapped)
/// keep their `R1`/`R2` role; only genuinely unpaired reads are [`ReadRole::Se`].
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum ReadRole {
    R1,
    R2,
    Se,
}

impl ReadRole {
    /// Dense index 0..3 for array storage.
    #[inline]
    fn index(self) -> usize {
        match self {
            ReadRole::R1 => 0,
            ReadRole::R2 => 1,
            ReadRole::Se => 2,
        }
    }

    /// All roles, in report order.
    pub(crate) const ALL: [ReadRole; 3] = [ReadRole::R1, ReadRole::R2, ReadRole::Se];

    /// The paired-mate role: R1↔R2. `Se` has no mate and maps to itself (callers
    /// only use this for paired reads).
    pub(crate) fn mate(self) -> ReadRole {
        match self {
            ReadRole::R1 => ReadRole::R2,
            ReadRole::R2 => ReadRole::R1,
            ReadRole::Se => ReadRole::Se,
        }
    }

    /// Lowercase short label used in TSV output.
    pub(crate) fn label(self) -> &'static str {
        match self {
            ReadRole::R1 => "R1",
            ReadRole::R2 => "R2",
            ReadRole::Se => "SE",
        }
    }
}

/// Which end of the read a cycle is measured from. Paired/orphan reads use only
/// [`ReadEnd::FivePrime`]; single-end reads track both ends (their far template
/// terminus is unknown, so both ends may need masking).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum ReadEnd {
    FivePrime,
    ThreePrime,
}

impl ReadEnd {
    #[inline]
    fn index(self) -> usize {
        match self {
            ReadEnd::FivePrime => 0,
            ReadEnd::ThreePrime => 1,
        }
    }

    pub(crate) const ALL: [ReadEnd; 2] = [ReadEnd::FivePrime, ReadEnd::ThreePrime];

    pub(crate) fn label(self) -> &'static str {
        match self {
            ReadEnd::FivePrime => "5p",
            ReadEnd::ThreePrime => "3p",
        }
    }
}

/// Methylated / total cytosine counts at one cycle of one (role, end, context).
#[derive(Clone, Copy, Default, Debug)]
pub(crate) struct CycleCounts {
    /// Methylated (unconverted) calls at this cycle.
    meth: u64,
    /// Total monitored & called cytosines at this cycle.
    total: u64,
}

impl CycleCounts {
    /// Methylated (unconverted) call count.
    #[must_use]
    pub(crate) fn meth(self) -> u64 {
        self.meth
    }

    /// Total monitored & called count.
    #[must_use]
    pub(crate) fn total(self) -> u64 {
        self.total
    }

    /// Methylation fraction, or `None` when no calls were seen.
    #[must_use]
    pub(crate) fn frac(self) -> Option<f64> {
        (self.total > 0).then(|| self.meth as f64 / self.total as f64)
    }
}

/// Per-cycle M-bias counts, growable along the cycle axis (so any read length is
/// captured exactly) but bounded by [`MAX_CYCLE`]. Indexed densely by
/// `(role, end, context)`.
///
/// Storage is a single contiguous buffer of `N_CELLS × stride` cells laid out
/// `slot-major` (cell `(slot, cycle)` at `slot * stride + cycle`), rather than a
/// `Vec` per slot. The hot [`record`](Self::record) path then touches one buffer
/// with a single bounds-checked index instead of chasing a per-slot `Vec`
/// pointer and bounds-checking it twice — the per-site accumulation is the
/// dominant cost of the `--metrics-prefix` M-bias walk, so removing that
/// indirection matters. `stride` starts at zero and grows (rarely — typically
/// once, on the first record, to cover the read length) to fit the largest cycle
/// seen; `lens[slot]` tracks each slot's highest recorded cycle so reports see
/// exactly the cycles that were observed.
pub(crate) struct MbiasAccumulator {
    /// Flat `N_CELLS × stride` buffer, `slot`-major.
    data: Vec<CycleCounts>,
    /// Cycles allocated per slot. All slots share one stride for O(1) indexing.
    stride: usize,
    /// Highest recorded cycle + 1, per slot (the logical length of each slice).
    lens: [usize; N_CELLS],
}

impl MbiasAccumulator {
    /// An empty accumulator.
    #[must_use]
    pub(crate) fn new() -> Self {
        Self { data: Vec::new(), stride: 0, lens: [0; N_CELLS] }
    }

    /// Dense slot for a `(role, end, context)` triple.
    #[inline]
    fn slot(role: ReadRole, end: ReadEnd, ctx: Context) -> usize {
        (role.index() * 2 + end.index()) * 4 + ctx.index()
    }

    /// Grow the shared stride to fit `min_stride` cycles, re-laying every slot's
    /// existing counts at the new (larger) stride. Geometric growth keeps this
    /// amortized O(1): for fixed-length reads it runs once.
    #[cold]
    fn grow_stride(&mut self, min_stride: usize) {
        let new_stride = min_stride.next_power_of_two().max(64);
        let mut new_data = vec![CycleCounts::default(); N_CELLS * new_stride];
        for slot in 0..N_CELLS {
            let used = self.lens[slot];
            if used > 0 {
                let src = slot * self.stride;
                let dst = slot * new_stride;
                new_data[dst..dst + used].copy_from_slice(&self.data[src..src + used]);
            }
        }
        self.data = new_data;
        self.stride = new_stride;
    }

    /// Record one classified site at `cycle` (distance from `end`).
    #[inline]
    pub(crate) fn record(
        &mut self,
        role: ReadRole,
        end: ReadEnd,
        ctx: Context,
        cycle: usize,
        unconverted: bool,
    ) {
        if cycle > MAX_CYCLE {
            return;
        }
        if cycle >= self.stride {
            self.grow_stride(cycle + 1);
        }
        let slot = Self::slot(role, end, ctx);
        let cell = &mut self.data[slot * self.stride + cycle];
        cell.total += 1;
        cell.meth += u64::from(unconverted);
        if cycle >= self.lens[slot] {
            self.lens[slot] = cycle + 1;
        }
    }

    /// Per-cycle counts for a `(role, end, context)`, oldest→newest cycle.
    #[must_use]
    pub(crate) fn cycles(&self, role: ReadRole, end: ReadEnd, ctx: Context) -> &[CycleCounts] {
        let slot = Self::slot(role, end, ctx);
        let start = slot * self.stride;
        &self.data[start..start + self.lens[slot]]
    }

    /// Whether any site was recorded for this `(role, end)` in any context — used
    /// to skip emitting empty rows / deciding masks for absent read classes.
    #[must_use]
    pub(crate) fn has_data(&self, role: ReadRole, end: ReadEnd) -> bool {
        Context::ALL.iter().any(|&c| self.lens[Self::slot(role, end, c)] > 0)
    }
}

impl Default for MbiasAccumulator {
    fn default() -> Self {
        Self::new()
    }
}

/// Tunables for [`detect_mask_length`].
#[derive(Clone, Copy, Debug)]
pub(crate) struct DetectParams {
    /// Keep from the first cycle whose smoothed rate reaches this fraction of the
    /// plateau (`--mbias-plateau-fraction`, e.g. 0.90).
    pub(crate) plateau_fraction: f64,
    /// Never mask more than this many leading cycles (`--mbias-max-mask`).
    pub(crate) max_mask: usize,
    /// Ignore cycles with fewer than this many calls when estimating the plateau
    /// and scanning, so sparse cycles don't drive the decision.
    pub(crate) min_cycle_cov: u64,
}

impl Default for DetectParams {
    fn default() -> Self {
        Self { plateau_fraction: 0.90, max_mask: 30, min_cycle_cov: 100 }
    }
}

/// Decide how many cycles to mask from `end` of `role`, from the CpG curve.
///
/// The plateau is the median of well-covered rates *beyond* `max_mask` (past the
/// ramp); the mask length is the first cycle whose lightly-smoothed rate reaches
/// `plateau_fraction × plateau`, minus one — i.e. we cut just before the read
/// becomes trustworthy. Scanning for the *first* crossing (rather than the last
/// deviation) means a noisy interior never inflates the mask. Returns 0 when
/// there is no usable signal (no data, or no estimable plateau).
#[must_use]
pub(crate) fn detect_mask_length(
    acc: &MbiasAccumulator,
    role: ReadRole,
    end: ReadEnd,
    p: DetectParams,
) -> usize {
    let cycles = acc.cycles(role, end, Context::CpG);
    if cycles.is_empty() {
        return 0;
    }

    // Per-cycle CpG methylation rate where coverage is adequate (else NaN).
    let rates: Vec<f64> =
        cycles
            .iter()
            .map(|cc| {
                if cc.total >= p.min_cycle_cov {
                    cc.meth as f64 / cc.total as f64
                } else {
                    f64::NAN
                }
            })
            .collect();

    let Some(plateau) = estimate_plateau(&rates, p.max_mask) else {
        return 0;
    };
    let threshold = plateau * p.plateau_fraction;

    let smoothed = moving_average(&rates, 3);
    let scan_to = smoothed.len().min(p.max_mask + 1);
    for (cycle, &r) in smoothed.iter().enumerate().take(scan_to) {
        if r.is_finite() && r >= threshold {
            return cycle; // first trustworthy cycle → mask the `cycle` cycles before it
        }
    }
    // No cycle reached the plateau within the cap → mask up to the cap.
    p.max_mask.min(cycles.len())
}

/// Median of finite per-cycle rates from cycle `max_mask` onward (the region
/// past the typical ramp). Falls back to all finite rates when too few cycles
/// lie beyond `max_mask` (short reads). `None` if there are no covered cycles.
fn estimate_plateau(rates: &[f64], max_mask: usize) -> Option<f64> {
    let interior: Vec<f64> =
        rates.iter().skip(max_mask).copied().filter(|r| r.is_finite()).collect();
    let pool = if interior.len() >= 3 {
        interior
    } else {
        rates.iter().copied().filter(|r| r.is_finite()).collect()
    };
    median(pool)
}

/// Median of a set of values, or `None` when empty.
fn median(mut v: Vec<f64>) -> Option<f64> {
    if v.is_empty() {
        return None;
    }
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let n = v.len();
    Some(if n % 2 == 1 { v[n / 2] } else { (v[n / 2 - 1] + v[n / 2]) / 2.0 })
}

/// Centered moving average of half-width `(w-1)/2`, skipping NaN inputs. A
/// window whose neighbors are all NaN stays NaN. Tames single-cycle spikes
/// before the crossing scan without shifting the ramp.
fn moving_average(rates: &[f64], w: usize) -> Vec<f64> {
    let half = w / 2;
    (0..rates.len())
        .map(|i| {
            let lo = i.saturating_sub(half);
            let hi = (i + half + 1).min(rates.len());
            let (sum, n) = rates[lo..hi]
                .iter()
                .filter(|r| r.is_finite())
                .fold((0.0, 0u32), |(s, n), r| (s + r, n + 1));
            if n == 0 { f64::NAN } else { sum / n as f64 }
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Fill a CpG 5' curve for `role` from an explicit per-cycle (meth, total)
    /// list, starting at cycle 0.
    fn fill_cpg_5p(acc: &mut MbiasAccumulator, role: ReadRole, per_cycle: &[(u64, u64)]) {
        for (cycle, &(meth, total)) in per_cycle.iter().enumerate() {
            for _ in 0..meth {
                acc.record(role, ReadEnd::FivePrime, Context::CpG, cycle, true);
            }
            for _ in 0..(total - meth) {
                acc.record(role, ReadEnd::FivePrime, Context::CpG, cycle, false);
            }
        }
    }

    /// A flat plateau of `n` cycles at rate `p` (cov 1000/cycle).
    fn plateau_curve(n: usize, p: f64) -> Vec<(u64, u64)> {
        (0..n).map(|_| ((p * 1000.0).round() as u64, 1000)).collect()
    }

    #[test]
    fn clean_ramp_masks_through_last_low_cycle() {
        // Cycles 0..3 are depressed (end-repair fill-in), then a flat 0.75
        // plateau. With fraction 0.90, threshold = 0.675; cycle 4 (0.75) is the
        // first to clear it → mask the 4 leading cycles.
        let mut acc = MbiasAccumulator::new();
        let mut curve = vec![(100, 1000), (300, 1000), (500, 1000), (640, 1000)];
        curve.extend(plateau_curve(40, 0.75));
        fill_cpg_5p(&mut acc, ReadRole::R2, &curve);
        let k = detect_mask_length(&acc, ReadRole::R2, ReadEnd::FivePrime, DetectParams::default());
        assert_eq!(k, 4, "first cycle ≥ 0.675 is cycle 4");
    }

    #[test]
    fn already_at_plateau_masks_nothing() {
        let mut acc = MbiasAccumulator::new();
        fill_cpg_5p(&mut acc, ReadRole::R1, &plateau_curve(40, 0.72));
        let k = detect_mask_length(&acc, ReadRole::R1, ReadEnd::FivePrime, DetectParams::default());
        assert_eq!(k, 0);
    }

    #[test]
    fn noisy_interior_does_not_inflate_mask() {
        // Ramp ends by cycle 3, but cycle 20 has a single-cycle dip well below
        // threshold. First-crossing semantics ignore it (a last-deviation rule
        // would wrongly mask 20+).
        let mut acc = MbiasAccumulator::new();
        let mut curve = vec![(200, 1000), (400, 1000), (600, 1000)];
        curve.extend(plateau_curve(40, 0.75));
        curve[20] = (300, 1000); // interior dip to 0.30
        fill_cpg_5p(&mut acc, ReadRole::R2, &curve);
        let k = detect_mask_length(&acc, ReadRole::R2, ReadEnd::FivePrime, DetectParams::default());
        assert!(k <= 3, "interior dip must not inflate the mask, got {k}");
    }

    #[test]
    fn no_data_masks_nothing() {
        let acc = MbiasAccumulator::new();
        assert_eq!(
            detect_mask_length(&acc, ReadRole::Se, ReadEnd::ThreePrime, DetectParams::default()),
            0
        );
    }

    #[test]
    fn mask_capped_at_max_mask() {
        // A curve that rises so slowly it never reaches 90% of the plateau within
        // the cap → mask exactly max_mask. Cycles 0..35 ramp 0.00→0.34 (still well
        // below the ~0.70 plateau at 0.90× = 0.63), so no cycle ≤ cap clears it.
        let mut acc = MbiasAccumulator::new();
        let curve: Vec<(u64, u64)> =
            (0..50).map(|c| if c < 35 { (10 * c as u64, 1000) } else { (700, 1000) }).collect();
        fill_cpg_5p(&mut acc, ReadRole::R2, &curve);
        let p = DetectParams { plateau_fraction: 0.90, max_mask: 10, min_cycle_cov: 100 };
        let k = detect_mask_length(&acc, ReadRole::R2, ReadEnd::FivePrime, p);
        assert_eq!(k, 10, "no recovery within cap → mask exactly max_mask");
    }

    #[test]
    fn record_and_cycles_roundtrip() {
        let mut acc = MbiasAccumulator::new();
        acc.record(ReadRole::R1, ReadEnd::FivePrime, Context::CpG, 5, true);
        acc.record(ReadRole::R1, ReadEnd::FivePrime, Context::CpG, 5, false);
        let cy = acc.cycles(ReadRole::R1, ReadEnd::FivePrime, Context::CpG);
        assert_eq!(cy.len(), 6);
        assert_eq!(cy[5].meth, 1);
        assert_eq!(cy[5].total, 2);
        assert_eq!(cy[5].frac(), Some(0.5));
        assert!(acc.has_data(ReadRole::R1, ReadEnd::FivePrime));
        assert!(!acc.has_data(ReadRole::R2, ReadEnd::FivePrime));
    }
}
