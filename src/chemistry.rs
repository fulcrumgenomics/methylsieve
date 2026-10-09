//! Library chemistry: which cytosines the conversion turns into T, and a check
//! that the data's conversion pattern agrees with the chemistry the user named.
//!
//! Bisulfite and EM-seq convert unmethylated cytosines; TAPS and Illumina 5-base
//! convert methylated ones. Everywhere else, tallies are kept in
//! chemistry-neutral retained/converted terms, and [`Readout`] is the one place
//! that turns them into methylation.

use anyhow::{Result, bail};
use clap::ValueEnum;

use crate::reference::Context;

/// TAPS-like data converts little CpH, because CpH is nearly all unmethylated.
/// Illumina 5-base runs near 1% from false conversion, so the ceiling leaves room.
const TAPS_LIKE_MAX_CPH_CONVERSION: f64 = 0.05;
/// ...while converting a visible share of CpG (the methylated ones). Kept above
/// the CpH ceiling, because a failed EM-seq library's CpG conversion never exceeds
/// its CpH conversion by more than the C/T SNP rate at CpG (about 1%).
const TAPS_LIKE_MIN_CPG_CONVERSION: f64 = 0.10;
/// EM-seq-like data converts most CpH, because most CpH is unmethylated.
const EMSEQ_LIKE_MIN_CPH_CONVERSION: f64 = 0.80;
/// Genome CpH sites needed before a mismatch stops the run.
const CHECK_MIN_CPH_SITES: u64 = 100_000;
/// Genome CpG sites needed before a mismatch stops the run.
const CHECK_MIN_CPG_SITES: u64 = 10_000;

/// The library chemistry (`--chemistry`).
#[derive(Copy, Clone, Debug, PartialEq, Eq, ValueEnum)]
pub(crate) enum Chemistry {
    /// Enzymatic methyl-seq: unmethylated cytosines read as T.
    #[value(name = "em-seq")]
    EmSeq,
    /// Bisulfite sequencing: unmethylated cytosines read as T.
    #[value(name = "bisulfite")]
    Bisulfite,
    /// TAPS and TAPS+: methylated cytosines read as T.
    #[value(name = "taps")]
    Taps,
    /// Illumina 5-base: methylated cytosines read as T.
    #[value(name = "illumina-5base")]
    Illumina5Base,
}

impl Chemistry {
    /// Which cytosines this chemistry converts.
    #[must_use]
    pub(crate) fn readout(self) -> Readout {
        match self {
            Chemistry::EmSeq | Chemistry::Bisulfite => Readout::UnmethylatedConverts,
            Chemistry::Taps | Chemistry::Illumina5Base => Readout::MethylatedConverts,
        }
    }

    /// The `--chemistry` spelling, used in metric output and messages.
    #[must_use]
    pub(crate) fn name(self) -> &'static str {
        match self {
            Chemistry::EmSeq => "em-seq",
            Chemistry::Bisulfite => "bisulfite",
            Chemistry::Taps => "taps",
            Chemistry::Illumina5Base => "illumina-5base",
        }
    }
}

/// Which cytosines a chemistry converts to T. Every [`Chemistry`] reduces to one
/// of these two.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub(crate) enum Readout {
    /// Bisulfite and EM-seq: a cytosine that still reads as C was methylated.
    UnmethylatedConverts,
    /// TAPS and Illumina 5-base: a cytosine that reads as T was methylated.
    MethylatedConverts,
}

impl Readout {
    /// Methylated calls among `total` monitored cytosines, of which `retained`
    /// still read as C.
    #[must_use]
    pub(crate) fn methylated(self, retained: u64, total: u64) -> u64 {
        match self {
            Readout::UnmethylatedConverts => retained,
            Readout::MethylatedConverts => total - retained,
        }
    }

    /// Whether a template can be called unconverted. That needs a per-site
    /// expectation: when unmethylated cytosines convert, nearly every CpH should
    /// read as T, and a failed molecule visibly keeps them as C. When methylated
    /// cytosines convert, a C that stayed C looks the same as an unmethylated one.
    #[must_use]
    pub(crate) fn calls_unconverted(self) -> bool {
        self == Readout::UnmethylatedConverts
    }
}

/// Converted and total cytosine counts at CpG and at CpH (CpA, CpC and CpT
/// pooled): the evidence the chemistry check classifies.
#[derive(Copy, Clone, Debug, Default, PartialEq, Eq)]
pub(crate) struct ConversionCounts {
    cpg_converted: u64,
    cpg_total: u64,
    cph_converted: u64,
    cph_total: u64,
}

impl ConversionCounts {
    /// Add `total` monitored cytosines in `ctx`, of which `retained` still read as C.
    pub(crate) fn add(&mut self, ctx: Context, retained: u64, total: u64) {
        let converted = total - retained;
        if ctx == Context::CpG {
            self.cpg_converted += converted;
            self.cpg_total += total;
        } else {
            self.cph_converted += converted;
            self.cph_total += total;
        }
    }

    /// The conversion pattern these counts show.
    fn classify(&self) -> DataClass {
        let cph = self.cph_rate();
        if cph.is_some_and(|r| r >= EMSEQ_LIKE_MIN_CPH_CONVERSION) {
            DataClass::EmSeqLike
        } else if cph.is_some_and(|r| r < TAPS_LIKE_MAX_CPH_CONVERSION)
            && self.cpg_rate().is_some_and(|r| r >= TAPS_LIKE_MIN_CPG_CONVERSION)
        {
            DataClass::TapsLike
        } else {
            DataClass::Ambiguous
        }
    }

    /// Whether there are enough sites for a mismatch to stop the run.
    fn meets_floors(&self) -> bool {
        self.cph_total >= CHECK_MIN_CPH_SITES && self.cpg_total >= CHECK_MIN_CPG_SITES
    }

    fn cpg_rate(&self) -> Option<f64> {
        (self.cpg_total > 0).then(|| self.cpg_converted as f64 / self.cpg_total as f64)
    }

    fn cph_rate(&self) -> Option<f64> {
        (self.cph_total > 0).then(|| self.cph_converted as f64 / self.cph_total as f64)
    }

    /// Both rates for messages, with `n/a` where no sites were seen.
    fn describe(&self) -> String {
        let fmt = |r: Option<f64>| r.map_or_else(|| "n/a".to_string(), |r| format!("{r:.4}"));
        format!(
            "genome CpG conversion {} over {} sites, CpH conversion {} over {} sites",
            fmt(self.cpg_rate()),
            self.cpg_total,
            fmt(self.cph_rate()),
            self.cph_total
        )
    }
}

/// The conversion pattern of a set of counts, independent of what was declared.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
enum DataClass {
    /// Methylated cytosines converted: CpH barely converted, CpG visibly so.
    TapsLike,
    /// Unmethylated cytosines converted: most CpH converted.
    EmSeqLike,
    /// Neither pattern. An unconverted or failed library under either chemistry
    /// lands here, as does TAPS data with very low CpG or high CpH methylation.
    Ambiguous,
}

impl DataClass {
    fn label(self) -> &'static str {
        match self {
            DataClass::TapsLike => "TAPS-like",
            DataClass::EmSeqLike => "EM-seq-like",
            DataClass::Ambiguous => "ambiguous",
        }
    }
}

/// Checks once per run that the genome's conversion pattern agrees with the
/// declared chemistry. Only a positive match to the *other* chemistry's pattern
/// counts as a mismatch, so a failed library is never mistaken for the wrong
/// chemistry: under EM-seq, CpG conversion can never exceed CpH conversion,
/// however many molecules failed.
pub(crate) struct ChemistryCheck {
    chemistry: Chemistry,
    done: bool,
}

impl ChemistryCheck {
    /// A check that has not run yet.
    #[must_use]
    pub(crate) fn new(chemistry: Chemistry) -> Self {
        Self { chemistry, done: false }
    }

    /// Whether the check has yet to run.
    #[must_use]
    pub(crate) fn is_pending(&self) -> bool {
        !self.done
    }

    /// Run the check once `counts` clears the site floors; until then, do nothing.
    ///
    /// # Errors
    /// Returns an error when the data shows the other chemistry's pattern.
    pub(crate) fn check_if_ready(&mut self, counts: &ConversionCounts) -> Result<()> {
        if self.done || !counts.meets_floors() {
            return Ok(());
        }
        self.done = true;
        let class = counts.classify();
        if self.contradicts(class) {
            bail!("{}", self.mismatch_message(counts));
        }
        log::info!(
            "Chemistry check: {}: {}, {} with --chemistry {}.",
            counts.describe(),
            class.label(),
            if class == DataClass::Ambiguous { "so not verified" } else { "consistent" },
            self.chemistry.name()
        );
        Ok(())
    }

    /// At end of run, if the floors were never reached, classify what was seen
    /// and warn on a mismatch instead of failing: too few sites to be sure.
    pub(crate) fn finish(&mut self, counts: &ConversionCounts) {
        if self.done {
            return;
        }
        self.done = true;
        if self.contradicts(counts.classify()) {
            log::warn!(
                "{} Too few sites to be certain ({CHECK_MIN_CPG_SITES} CpG and \
                 {CHECK_MIN_CPH_SITES} CpH needed), so the run was not stopped.",
                self.mismatch_message(counts)
            );
        }
    }

    /// Whether `class` is the pattern of the other readout.
    fn contradicts(&self, class: DataClass) -> bool {
        matches!(
            (self.chemistry.readout(), class),
            (Readout::UnmethylatedConverts, DataClass::TapsLike)
                | (Readout::MethylatedConverts, DataClass::EmSeqLike)
        )
    }

    fn mismatch_message(&self, counts: &ConversionCounts) -> String {
        let name = self.chemistry.name();
        let (expected, observed, suggestion) = match self.chemistry.readout() {
            Readout::UnmethylatedConverts => (
                "unmethylated",
                "TAPS or Illumina 5-base data, where methylated cytosines convert",
                "--chemistry taps or --chemistry illumina-5base",
            ),
            Readout::MethylatedConverts => (
                "methylated",
                "bisulfite or EM-seq data, where unmethylated cytosines convert",
                "--chemistry em-seq or --chemistry bisulfite",
            ),
        };
        format!(
            "--chemistry {name} expects {expected} cytosines to convert, but the reads look \
             like {observed} ({}). Re-run with {suggestion}.",
            counts.describe()
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Counts with the given converted fractions over `n` CpG and `n` CpH sites.
    fn counts(cpg_conv: f64, cph_conv: f64, n: u64) -> ConversionCounts {
        let mut c = ConversionCounts::default();
        c.add(Context::CpG, n - (cpg_conv * n as f64).round() as u64, n);
        c.add(Context::CpA, n - (cph_conv * n as f64).round() as u64, n);
        c
    }

    #[test]
    fn taps_and_illumina_5base_convert_methylated_cytosines() {
        assert_eq!(Chemistry::Taps.readout(), Readout::MethylatedConverts);
        assert_eq!(Chemistry::Illumina5Base.readout(), Readout::MethylatedConverts);
    }

    #[test]
    fn emseq_and_bisulfite_convert_unmethylated_cytosines() {
        assert_eq!(Chemistry::EmSeq.readout(), Readout::UnmethylatedConverts);
        assert_eq!(Chemistry::Bisulfite.readout(), Readout::UnmethylatedConverts);
    }

    #[test]
    fn methylated_counts_retained_cytosines_when_unmethylated_convert() {
        assert_eq!(Readout::UnmethylatedConverts.methylated(30, 100), 30);
    }

    #[test]
    fn methylated_counts_converted_cytosines_when_methylated_convert() {
        assert_eq!(Readout::MethylatedConverts.methylated(30, 100), 70);
    }

    #[test]
    fn only_unmethylated_converting_chemistries_call_unconverted_templates() {
        assert!(Readout::UnmethylatedConverts.calls_unconverted());
        assert!(!Readout::MethylatedConverts.calls_unconverted());
    }

    #[test]
    fn mammalian_emseq_pattern_is_emseq_like() {
        assert_eq!(counts(0.25, 0.99, 1000).classify(), DataClass::EmSeqLike);
    }

    #[test]
    fn mammalian_taps_pattern_is_taps_like() {
        assert_eq!(counts(0.47, 0.006, 1000).classify(), DataClass::TapsLike);
    }

    #[test]
    fn hypomethylated_taps_sample_is_still_taps_like() {
        // The ZymoUnmethylated TAPS+ library measured 0.127 / 0.006.
        assert_eq!(counts(0.127, 0.006, 1000).classify(), DataClass::TapsLike);
    }

    #[test]
    fn taps_like_requires_cpg_conversion_of_at_least_10_percent() {
        assert_eq!(counts(0.099, 0.005, 1000).classify(), DataClass::Ambiguous);
        assert_eq!(counts(0.100, 0.005, 1000).classify(), DataClass::TapsLike);
    }

    #[test]
    fn taps_like_requires_cph_conversion_below_5_percent() {
        assert_eq!(counts(0.47, 0.050, 1000).classify(), DataClass::Ambiguous);
        assert_eq!(counts(0.47, 0.049, 1000).classify(), DataClass::TapsLike);
    }

    #[test]
    fn illumina_5base_pattern_is_taps_like() {
        // The HG002 5-base demo library measured CpG 0.649, CpH 0.0106.
        assert_eq!(counts(0.649, 0.0106, 10_000).classify(), DataClass::TapsLike);
    }

    #[test]
    fn emseq_library_that_barely_converted_is_ambiguous() {
        // 5% efficiency: CpH conversion just over the TAPS ceiling, CpG below the
        // floor plus a 1% C/T SNP allowance, so never mistaken for TAPS.
        assert_eq!(
            counts(0.05 * 0.25 + 0.01, 0.05 * 0.99, 10_000).classify(),
            DataClass::Ambiguous
        );
    }

    #[test]
    fn emseq_like_starts_at_80_percent_cph_conversion() {
        assert_eq!(counts(0.25, 0.799, 1000).classify(), DataClass::Ambiguous);
        assert_eq!(counts(0.25, 0.800, 1000).classify(), DataClass::EmSeqLike);
    }

    #[test]
    fn entirely_unconverted_data_is_ambiguous() {
        assert_eq!(counts(0.0, 0.0, 1000).classify(), DataClass::Ambiguous);
    }

    #[test]
    fn partly_failed_emseq_library_is_ambiguous() {
        // 30% of molecules failed: CpG and CpH conversion both shrink by 0.7.
        assert_eq!(counts(0.7 * 0.25, 0.7 * 0.99, 1000).classify(), DataClass::Ambiguous);
    }

    #[test]
    fn no_sites_is_ambiguous() {
        assert_eq!(ConversionCounts::default().classify(), DataClass::Ambiguous);
    }

    #[test]
    fn floors_need_both_cph_and_cpg_sites() {
        let mut c = ConversionCounts::default();
        c.add(Context::CpA, 0, CHECK_MIN_CPH_SITES);
        assert!(!c.meets_floors());
        c.add(Context::CpG, 0, CHECK_MIN_CPG_SITES);
        assert!(c.meets_floors());
    }

    #[test]
    fn check_errors_on_taps_data_declared_emseq() {
        let mut check = ChemistryCheck::new(Chemistry::EmSeq);
        let err = check.check_if_ready(&counts(0.47, 0.006, 200_000)).unwrap_err();
        assert!(err.to_string().contains("--chemistry taps"), "{err}");
    }

    #[test]
    fn check_errors_on_emseq_data_declared_taps() {
        let mut check = ChemistryCheck::new(Chemistry::Taps);
        let err = check.check_if_ready(&counts(0.25, 0.99, 200_000)).unwrap_err();
        assert!(err.to_string().contains("--chemistry em-seq"), "{err}");
    }

    #[test]
    fn check_passes_matching_data() {
        let mut check = ChemistryCheck::new(Chemistry::Illumina5Base);
        assert!(check.check_if_ready(&counts(0.47, 0.006, 200_000)).is_ok());
        assert!(!check.is_pending());
    }

    #[test]
    fn check_passes_ambiguous_data_under_either_readout() {
        for chemistry in [Chemistry::EmSeq, Chemistry::Taps] {
            let mut check = ChemistryCheck::new(chemistry);
            assert!(check.check_if_ready(&counts(0.0, 0.0, 200_000)).is_ok());
        }
    }

    #[test]
    fn check_waits_below_the_site_floors() {
        let mut check = ChemistryCheck::new(Chemistry::EmSeq);
        assert!(check.check_if_ready(&counts(0.47, 0.006, 1000)).is_ok());
        assert!(check.is_pending());
    }

    #[test]
    fn finish_never_errors_below_the_site_floors() {
        let mut check = ChemistryCheck::new(Chemistry::EmSeq);
        check.finish(&counts(0.47, 0.006, 1000));
        assert!(!check.is_pending());
    }
}
