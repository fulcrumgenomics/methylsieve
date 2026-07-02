//! methylsieve — the command-line application.
//!
//! This binary is the whole tool: CLI parsing ([`Args`]), end-to-end run
//! orchestration ([`run`]), and the end-of-run resource-usage footer. The
//! per-template tagging engine (`sieve`), M-bias learning/masking (`mbias`,
//! `mask`, `buffer`), metric output (`metrics`), shared record geometry
//! (`record`), the reference, and IO live in sibling modules — this file wires
//! them together and talks to the user.
//!
//! Long-form flags follow GNU style (`--kebab-case`). Short flags mirror the
//! conventions of sibling tools: `-i` input, `-o` output, `-r` reference,
//! `-q` quiet.

mod buffer;
mod io_threading;
mod mask;
mod mbias;
mod metrics;
mod plots;
mod raw_reader;
mod raw_writer;
mod record;
mod reference;
mod sam_reader;
mod sieve;

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, Write};
use std::path::PathBuf;
use std::process::ExitCode;
use std::time::Instant;

use anyhow::{Context as _, Result, bail};
use bgzf::CompressionLevel;
use clap::{Parser, ValueEnum};
use fgumi_raw_bam::RawRecord;
use noodles_sam::Header;
use noodles_sam::header::record::value::Map;
use noodles_sam::header::record::value::map::Program;
use noodles_sam::header::record::value::map::program::tag as program_tag;

use crate::buffer::TemplateArena;
use crate::mask::{MaskPlan, MaskWindows, apply_windows, compute_mask_windows};
use crate::mbias::{DetectParams, MbiasAccumulator};
use crate::raw_reader::RawBamReader;
use crate::raw_writer::RawBamWriter;
use crate::reference::{Context, Reference};
use crate::sam_reader::SamReader;
use crate::sieve::{DecisionMode, Disposition, ProcessorOptions, RecordProcessor, Stats};

/// Crate build identifier shown in `--version` and the `@PG VN:` tag.
const METHYLSIEVE_BUILD: &str = env!("CARGO_PKG_VERSION");

/// Global allocator — mimalloc keeps the per-block temporary allocations that
/// dominate the streaming worker cheap.
#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

// ── Command-line interface ──────────────────────────────────────────────────

/// CLI mirror of [`DecisionMode`] so the kebab-case spellings and per-variant
/// help live with the CLI.
#[derive(Copy, Clone, Debug, PartialEq, Eq, ValueEnum)]
pub(crate) enum DecisionModeCli {
    /// Count only: flag when the unconverted count reaches
    /// `--max-unconverted-count`. `--min-sites` and the fraction are ignored.
    #[value(name = "count")]
    Count,
    /// Proportion only: flag when a template has at least `--min-sites` sites
    /// AND its unconverted fraction exceeds `--max-unconverted-fraction`.
    /// Templates with fewer than `--min-sites` sites are NEVER flagged — the
    /// proportion can't be estimated, so they pass through unevaluated.
    #[value(name = "proportion")]
    Proportion,
    /// Either: flag when the count OR the proportion test fires.
    #[value(name = "either")]
    Either,
    /// (Default) Adaptive: use the proportion test at/above `--min-sites` and
    /// the count test below it. Low-site templates are still evaluated (by
    /// count), while high-site templates are judged on rate — avoiding an
    /// absolute count that over-penalizes long reads / read pairs.
    #[value(name = "adaptive")]
    Adaptive,
}

impl From<DecisionModeCli> for DecisionMode {
    fn from(c: DecisionModeCli) -> Self {
        match c {
            DecisionModeCli::Count => DecisionMode::Count,
            DecisionModeCli::Proportion => DecisionMode::Proportion,
            DecisionModeCli::Either => DecisionMode::Either,
            DecisionModeCli::Adaptive => DecisionMode::Adaptive,
        }
    }
}

/// Parse a `--compression-level` argument, delegating range validation to
/// `bgzf::CompressionLevel::new` (0 = stored, 1-9 zlib, 10-12 libdeflate).
fn parse_compression_level(s: &str) -> Result<CompressionLevel, String> {
    let n: u8 = s.parse().map_err(|e| format!("not a u8: {e}"))?;
    CompressionLevel::new(n).map_err(|e| format!("{e}"))
}

/// The subset of the four contexts (CpA/CpC/CpG/CpT) that count toward the
/// unconverted-decision threshold. Stored as a presence flag per context,
/// indexed by [`Context::index`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct ContextMask([bool; 4]);

impl ContextMask {
    /// Whether `ctx` participates in the threshold.
    #[inline]
    #[must_use]
    pub(crate) fn contains(self, ctx: Context) -> bool {
        self.0[ctx.index()]
    }

    /// Short human label for the selected contexts: `"CpH"` for the CpA/CpC/CpT
    /// set (the conversion default), otherwise the selected context names joined
    /// with `+` (e.g. `"CpG"`, `"CpA+CpG"`). Used in plot titles/axes.
    #[must_use]
    pub(crate) fn label(self) -> String {
        let is_cph = self.contains(Context::CpA)
            && self.contains(Context::CpC)
            && self.contains(Context::CpT)
            && !self.contains(Context::CpG);
        if is_cph {
            return "CpH".to_string();
        }
        Context::ALL
            .iter()
            .filter(|&&c| self.contains(c))
            .map(|c| c.label())
            .collect::<Vec<_>>()
            .join("+")
    }
}

/// Parse a comma-separated context list (e.g. `CpA,CpC,CpT`). Token matching is
/// case-insensitive and accepts both `CpA` and `CA` spellings.
pub(crate) fn parse_contexts(s: &str) -> Result<ContextMask, String> {
    let mut mask = [false; 4];
    let mut any = false;
    for tok in s.split(',') {
        let t = tok.trim().to_ascii_uppercase();
        let ctx = match t.as_str() {
            "CPA" | "CA" => Context::CpA,
            "CPC" | "CC" => Context::CpC,
            "CPG" | "CG" => Context::CpG,
            "CPT" | "CT" => Context::CpT,
            "" => continue,
            other => {
                return Err(format!(
                    "unknown context '{other}'; valid contexts are CpA, CpC, CpG, CpT"
                ));
            }
        };
        mask[ctx.index()] = true;
        any = true;
    }
    if !any {
        return Err("at least one context must be specified".to_string());
    }
    Ok(ContextMask(mask))
}

/// A parsed aux-tag specification (`--tag`). Only string (`Z`) tags are
/// supported.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct TagSpec {
    /// The two-byte aux tag (e.g. `XX`).
    pub(crate) tag: [u8; 2],
    /// The string value to set (e.g. `UC`).
    pub(crate) value: Vec<u8>,
}

/// Parse a two-character aux-tag name (e.g. `ch`) for `--count-tag`.
fn parse_tag_name(s: &str) -> Result<[u8; 2], String> {
    let b = s.as_bytes();
    if b.len() != 2 || !b.iter().all(u8::is_ascii_alphanumeric) {
        return Err(format!("tag name '{s}' must be exactly two alphanumeric characters"));
    }
    Ok([b[0], b[1]])
}

/// Parse a SAM-style `TAG:Z:VALUE` spec (e.g. `XX:Z:UC`). Only the `Z`
/// (string) type is accepted.
fn parse_tag_spec(s: &str) -> Result<TagSpec, String> {
    let parts: Vec<&str> = s.splitn(3, ':').collect();
    if parts.len() != 3 {
        return Err(format!("tag spec '{s}' must be of the form TAG:Z:VALUE (e.g. XX:Z:UC)"));
    }
    let tag_bytes = parts[0].as_bytes();
    if tag_bytes.len() != 2 || !tag_bytes.iter().all(u8::is_ascii_alphanumeric) {
        return Err(format!("tag '{}' must be exactly two alphanumeric characters", parts[0]));
    }
    if !parts[1].eq_ignore_ascii_case("Z") {
        return Err(format!(
            "only string (Z) tags are supported; got type '{}' in '{s}'",
            parts[1]
        ));
    }
    if parts[2].is_empty() {
        return Err(format!("tag spec '{s}' has an empty value"));
    }
    Ok(TagSpec { tag: [tag_bytes[0], tag_bytes[1]], value: parts[2].as_bytes().to_vec() })
}

/// Parsed command-line arguments — see `--help` for per-flag descriptions.
#[derive(Parser, Debug, Clone)]
#[command(name = "methylsieve", version = METHYLSIEVE_BUILD, about = SHORT_ABOUT, long_about = LONG_ABOUT)]
pub struct Args {
    /// Input SAM/BAM file (defaults to stdin). Must be query-grouped (all records
    /// for one QNAME adjacent, typically straight from the aligner).
    #[arg(short = 'i', long = "input", help_heading = "Inputs / outputs")]
    pub(crate) input: Option<PathBuf>,

    /// Output BAM file (defaults to stdout). The path must end in `.bam` — output
    /// is always BAM, never SAM. Use `-` for stdout (no extension check).
    #[arg(short = 'o', long = "output", help_heading = "Inputs / outputs")]
    pub(crate) output: Option<PathBuf>,

    /// BGZF compression level for the output BAM (0-9 typical, up to 12 for
    /// libdeflate's strongest tier). Level 0 (the default) emits "stored" BGZF
    /// blocks — a downstream sort will recompress.
    #[arg(long = "compression-level", default_value = "0", value_parser = parse_compression_level,
          help_heading = "Inputs / outputs")]
    pub(crate) compression_level: CompressionLevel,

    /// Reference FASTA (a sibling `.fai` index is required). The FASTA may be a
    /// superset of the BAM's contigs; every BAM `@SQ` must be present with a
    /// matching length. NOTE: every referenced contig is loaded into memory
    /// (~3.1 GB for a human genome), since query-grouped input touches any
    /// contig at any time.
    #[arg(short = 'r', long = "reference", help_heading = "Inputs / outputs")]
    pub(crate) reference: PathBuf,

    /// Contexts (comma-separated CpA,CpC,CpT,CpG) counted toward the
    /// unconverted-decision threshold. The default is CpH (CpA,CpC,CpT); CpG is
    /// excluded because genuine methylation lives there. All four contexts are
    /// reported in the `--metrics-prefix` summary regardless of this subset.
    #[arg(long = "contexts", default_value = "CpA,CpC,CpT", value_parser = parse_contexts,
          help_heading = "Decision")]
    pub(crate) contexts: ContextMask,

    /// How to combine the count and proportion tests. `adaptive` (the default)
    /// uses the proportion at/above `--min-sites` and the count below it; `count`
    /// and `proportion` use only that one test; `either` flags when either fires.
    /// In `proportion` mode, templates with fewer than `--min-sites` sites are
    /// never flagged (the proportion is unestimable) — they pass through.
    #[arg(long = "mode", value_enum, default_value_t = DecisionModeCli::Adaptive,
          help_heading = "Decision")]
    pub(crate) mode: DecisionModeCli,

    /// Count test: a template is called unconverted when its count of
    /// unconverted cytosines (summed over the `--contexts` subset, across all
    /// evaluated records) is at least this value. Used by `count`, `either`, and
    /// (below `--min-sites`) `adaptive`.
    #[arg(long = "max-unconverted-count", default_value_t = 3, help_heading = "Decision")]
    pub(crate) max_unconverted_count: u32,

    /// Proportion test: a template is called unconverted when it has at least
    /// `--min-sites` sites AND its unconverted fraction exceeds this value. Used
    /// by `proportion`, `either`, and (at/above `--min-sites`) `adaptive`;
    /// ignored by `count`. Must be in (0, 1].
    #[arg(long = "max-unconverted-fraction", default_value_t = 0.05, help_heading = "Decision")]
    pub(crate) max_unconverted_fraction: f64,

    /// Floor for the proportion test (in `--contexts`-subset sites): below it the
    /// proportion is unestimable and abstains. In `adaptive` this is also the
    /// switch point — proportion at/above it, count below. The default (40) is
    /// the smallest value that makes the `adaptive` switch continuous given
    /// count=3 and fraction=0.05 (so a read doesn't flip its call as it crosses
    /// the threshold).
    #[arg(long = "min-sites", default_value_t = 40, help_heading = "Decision")]
    pub(crate) min_sites: u32,

    /// Skip read bases with base quality below this value when tallying. Higher
    /// values reduce sequencing-error-driven false unconverted calls but shrink
    /// the site count (so fewer templates clear `--min-sites`).
    #[arg(long = "min-base-quality", default_value_t = 20, help_heading = "Decision")]
    pub(crate) min_base_quality: u8,

    /// Ignore the outermost N bases at each end of the *template* (the original
    /// DNA fragment) when tallying — the bases prone to end-repair fill-in and
    /// A-tailing artifacts. For a mapped pair these are the 5' sequenced ends of
    /// R1 and R2 (the two fragment termini), whose reference positions are
    /// skipped in *both* mates (so an overlapped terminus is trimmed once, in
    /// each read that covers it). For single-end or orphan reads the far end is
    /// unknown, so both ends of the read are trimmed instead. Counted over the
    /// stored SEQ in sequencing order, so 5'-end soft-clips count toward N;
    /// hard-clipped bases are absent from SEQ and do not. Default 0 (off).
    /// Superseded by `--mbias-mask`: when masking is enabled this trim is forced
    /// to 0 (both target the same fragment-end bias) and a warning is logged if a
    /// non-zero value was set explicitly.
    #[arg(long = "ignore-template-ends", default_value_t = 0, help_heading = "Decision")]
    pub(crate) ignore_template_ends: u32,

    /// Exclude supplementary alignments from conversion tallying (they still
    /// receive the tag / qc-fail flag like every other record of an unconverted
    /// template). By default supplementaries contribute evidence.
    #[arg(long = "ignore-supplementary-evidence", help_heading = "Decision")]
    pub(crate) ignore_supplementary_evidence: bool,

    /// Aux tag to set on every record of an unconverted template, as a
    /// SAM-style `TAG:Z:VALUE` spec.
    #[arg(long = "tag", default_value = "XX:Z:UC", value_parser = parse_tag_spec,
          help_heading = "Actions on unconverted templates")]
    pub(crate) tag: TagSpec,

    /// Do not OR the QC-fail flag (0x200) into unconverted records' FLAG.
    /// (QC-fail marking is on by default.)
    #[arg(long = "no-qc-fail", help_heading = "Actions on unconverted templates")]
    pub(crate) no_qc_fail: bool,

    /// Drop every record of an unconverted template from the output entirely,
    /// instead of just tagging/flagging it.
    #[arg(long = "remove-unconverted", help_heading = "Actions on unconverted templates")]
    pub(crate) remove_unconverted: bool,

    /// Spike-in control contig (repeatable). Reads whose primary R1 maps here
    /// are excluded from the main decision, never tagged, and tallied into a
    /// separate metrics-summary row.
    #[arg(long = "control-contig", help_heading = "Spike-in controls")]
    pub(crate) control_contig: Vec<String>,

    /// Enable M-bias-aware masking: buffer the first reads to learn the per-cycle
    /// CpG methylation curves, freeze 5' (and, for single-end, 3') mask lengths,
    /// then set the biased bases' qualities to `--mbias-mask-quality` so
    /// base-quality-aware callers ignore them. The alignment is otherwise
    /// unchanged (no clip, no coordinate/CIGAR/tag/mate rewrite). Effective only
    /// for downstream tools that honor base quality.
    #[arg(long = "mbias-mask", help_heading = "M-bias masking")]
    pub(crate) mbias_mask: bool,

    /// Templates to buffer while learning M-bias before masking begins. Memory
    /// scales with this (≈0.5–1.3 KB/template). A pathological input stops
    /// buffering early and decides on what it has.
    #[arg(
        long = "mbias-buffer-templates",
        default_value_t = 500_000,
        help_heading = "M-bias masking"
    )]
    pub(crate) mbias_buffer_templates: usize,

    /// Keep from the first 5' cycle whose smoothed CpG methylation reaches this
    /// fraction of the plateau (so masking stops as soon as the read is
    /// trustworthy). Must be in (0, 1].
    #[arg(
        long = "mbias-plateau-fraction",
        default_value_t = 0.90,
        help_heading = "M-bias masking"
    )]
    pub(crate) mbias_plateau_fraction: f64,

    /// Never mask more than this many leading (or, for single-end, trailing)
    /// cycles, regardless of the curve.
    #[arg(long = "mbias-max-mask", default_value_t = 30, help_heading = "M-bias masking")]
    pub(crate) mbias_max_mask: u32,

    /// Quality value assigned to masked bases (default 2; keep below
    /// `--min-base-quality` so the masked bases also drop from this tool's tally).
    #[arg(long = "mbias-mask-quality", default_value_t = 2, help_heading = "M-bias masking")]
    pub(crate) mbias_mask_quality: u8,

    /// Write metric files under this path prefix: `PREFIX.summary.tsv` (one
    /// per-context conversion row per scope — the genome, then each control contig
    /// — with the applied 5'/3' mask lengths as `r1_mask_5p`/`r2_mask_5p`/
    /// `se_mask_5p`/`se_mask_3p` columns), `PREFIX.mbias.tsv` (per-read-cycle
    /// methylation), and `PREFIX.conversion-matrix.tsv` (the per-template decision
    /// histogram), plus PDF plots `PREFIX.mbias.pdf` (M-bias curves) and
    /// `PREFIX.conversion-matrix.pdf` (the decision hexbin). Computing these is a
    /// single streaming pass; the output BAM is unchanged.
    #[arg(long = "metrics-prefix", help_heading = "Stats & misc")]
    pub(crate) metrics_prefix: Option<PathBuf>,

    /// Sample name for the `sample` column of the metric TSVs. If omitted: the
    /// unique `@RG SM:` values from the header, else the input file stem, else
    /// `unknown`.
    #[arg(long = "sample", help_heading = "Stats & misc")]
    pub(crate) sample: Option<String>,

    /// Aux tag (2 characters) for the per-record count annotation: `<TAG>:Z:u/n`,
    /// where u is the unconverted count and n the total monitored sites in the
    /// `--contexts` subset — the exact numerator/denominator of the decision. It
    /// is a per-TEMPLATE aggregate stamped on every record of the template (not a
    /// per-read count), so a user can see why any read was (or wasn't) flagged.
    /// On by default with tag `ch`; pass a name to rename, or `--no-count-tag`
    /// to disable.
    #[arg(long = "count-tag", value_parser = parse_tag_name, default_value = "ch",
          conflicts_with = "no_count_tag", help_heading = "Stats & misc")]
    pub(crate) count_tag: [u8; 2],

    /// Disable the per-record count annotation (see `--count-tag`).
    #[arg(long = "no-count-tag", help_heading = "Stats & misc")]
    pub(crate) no_count_tag: bool,

    /// Output fewer statistics.
    #[arg(short = 'q', long = "quiet", help_heading = "Stats & misc")]
    pub(crate) quiet: bool,

    /// Verify BGZF CRC32 on input. Default: on for file input, off for stdin.
    #[arg(long = "check-crc", conflicts_with = "no_check_crc", help_heading = "Stats & misc")]
    pub(crate) check_crc: bool,

    /// Skip BGZF CRC32 verification on input regardless of source.
    #[arg(long = "no-check-crc", conflicts_with = "check_crc", help_heading = "Stats & misc")]
    pub(crate) no_check_crc: bool,

    /// Size (MB) of the ring buffer between the input IO thread and the worker.
    #[arg(long = "read-buffer-mb", default_value_t = 16, value_parser = clap::value_parser!(u32).range(1..=4096),
          help_heading = "Stats & misc")]
    pub(crate) read_buffer_mb: u32,

    /// Size (MB) of the ring buffer between the worker and the output IO thread.
    #[arg(long = "write-buffer-mb", default_value_t = 64, value_parser = clap::value_parser!(u32).range(1..=4096),
          help_heading = "Stats & misc")]
    pub(crate) write_buffer_mb: u32,
}

const SHORT_ABOUT: &str = "Tag or filter unconverted reads in bisulfite / EM-seq SAM/BAM files.";
const LONG_ABOUT: &str = "Tag or filter incompletely-converted reads in directional bisulfite or \
EM-seq data.\n\
Makes one per-template decision using all of a QNAME's primary and supplementary records \
and propagates it to every record; with --metrics-prefix it also writes per-context / \
per-spike-in conversion and M-bias metrics. \
Input must be query-grouped and should be adapter-trimmed first: untrimmed adapter \
read-through on short inserts can be force-aligned and read as spurious unconverted CpH. \
Output is always BAM.";

impl Args {
    /// Resolved CRC-verify setting: explicit flags win; otherwise on for file
    /// input and off for stdin (trusted producer, e.g. piped from bwa-meth).
    #[must_use]
    pub(crate) fn effective_check_crc(&self) -> bool {
        if self.check_crc {
            return true;
        }
        if self.no_check_crc {
            return false;
        }
        matches!(self.input.as_deref(), Some(p) if p.to_string_lossy() != "-")
    }

    /// Whether unconverted templates should be QC-fail flagged (0x200).
    #[must_use]
    pub(crate) fn qc_fail(&self) -> bool {
        !self.no_qc_fail
    }

    /// Resolved per-record count tag: `None` if disabled, else the tag name.
    #[must_use]
    pub(crate) fn effective_count_tag(&self) -> Option<[u8; 2]> {
        if self.no_count_tag { None } else { Some(self.count_tag) }
    }

    /// Validate invariants clap can't express directly.
    ///
    /// # Errors
    /// Returns an error if `-o` is not `-`/`.bam` or `--max-unconverted-fraction`
    /// is out of `(0, 1]`.
    pub(crate) fn validate(&self) -> Result<()> {
        if let Some(p) = &self.output {
            let s = p.to_string_lossy();
            if s != "-" && !s.ends_with(".bam") {
                bail!(
                    "output path {} must end in `.bam` (methylsieve only writes BAM); \
                     use `-` to send BAM to stdout",
                    p.display()
                );
            }
        }
        // Metric TSVs are always written to files under `--metrics-prefix`, so
        // only the BAM can contend for stdout — no inter-stream collision check
        // is needed.
        if !(self.max_unconverted_fraction > 0.0 && self.max_unconverted_fraction <= 1.0) {
            bail!(
                "--max-unconverted-fraction must be in (0, 1]; got {}",
                self.max_unconverted_fraction
            );
        }
        if self.max_unconverted_count == 0 {
            bail!(
                "--max-unconverted-count must be >= 1 (a threshold of 0 would tag every \
                 template). Use --max-unconverted-fraction for fraction-based filtering."
            );
        }
        if self.mbias_mask
            && !(self.mbias_plateau_fraction > 0.0 && self.mbias_plateau_fraction <= 1.0)
        {
            bail!(
                "--mbias-plateau-fraction must be in (0, 1]; got {}",
                self.mbias_plateau_fraction
            );
        }
        Ok(())
    }
}

// ── Binary entry point ──────────────────────────────────────────────────────

fn main() -> ExitCode {
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("info"))
        .format_timestamp(None)
        .init();

    let args = Args::parse();
    match run(args) {
        Ok(()) => ExitCode::SUCCESS,
        Err(err) => {
            eprintln!("methylsieve: {err:#}");
            eprintln!("methylsieve: Premature exit (return code 1).");
            ExitCode::FAILURE
        }
    }
}

// ── Run orchestration ───────────────────────────────────────────────────────

/// Input format auto-detected from the first byte.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Format {
    Sam,
    Bam,
}

/// Peek the first byte to determine SAM (`@`) vs BAM (BGZF magic `0x1f`).
fn detect_format<R: BufRead>(reader: &mut R) -> Result<Format> {
    let head = reader.fill_buf().context("reading input to detect format")?;
    if head.is_empty() {
        bail!("Empty input: no SAM/BAM header detected");
    }
    match head[0] {
        0x1f => Ok(Format::Bam),
        b'@' => Ok(Format::Sam),
        b => bail!(
            "Input doesn't look like SAM or BAM (first byte 0x{:02x}); expected '@' or 0x1f",
            b
        ),
    }
}

/// Unified reader enum so the main loop is format-agnostic.
enum Reader {
    Bam(RawBamReader<Box<dyn BufRead>>),
    Sam(SamReader<Box<dyn BufRead>>),
}

impl Reader {
    fn read_header(&mut self) -> Result<Header> {
        match self {
            Reader::Bam(r) => r.read_header().context("reading BAM header"),
            Reader::Sam(r) => r.read_header().context("reading SAM header"),
        }
    }

    fn read_record(&mut self, rec: &mut RawRecord) -> std::io::Result<bool> {
        match self {
            Reader::Bam(r) => r.read_record(rec),
            Reader::Sam(r) => r.read_record(rec),
        }
    }
}

/// Main entry point used by both the binary and the test harness.
///
/// # Errors
/// Returns an error on IO failure, malformed input, a missing/mismatched
/// reference contig, or an unknown `--control-contig`.
fn run(args: Args) -> Result<()> {
    log::info!("methylsieve by Fulcrum Genomics - https://github.com/fulcrumgenomics/methylsieve");
    args.validate()?;
    let started = StartedRun::now();
    if !args.quiet {
        eprintln!("methylsieve: Version {METHYLSIEVE_BUILD}");
    }

    // Input goes through a dedicated IO read thread + ring buffer so the worker
    // never blocks on the kernel pipe.
    let raw_source: Box<dyn std::io::Read + Send> = match args.input.as_deref() {
        Some(p) if p.to_string_lossy() != "-" => {
            let f = File::open(p).with_context(|| format!("opening {} for read", p.display()))?;
            Box::new(f)
        }
        _ => Box::new(std::io::stdin()),
    };
    let read_buf_bytes = (args.read_buffer_mb as usize).saturating_mul(1024 * 1024);
    let mut reader_box: Box<dyn BufRead> =
        Box::new(crate::io_threading::ThreadedReader::new(raw_source, read_buf_bytes));
    let input_name =
        args.input.as_deref().map(|p| p.display().to_string()).unwrap_or_else(|| "stdin".into());
    let input_format = detect_format(&mut reader_box)?;
    let check_crc = args.effective_check_crc();
    let mut reader: Reader = match input_format {
        Format::Bam => {
            if !args.quiet {
                eprintln!(
                    "methylsieve: Reading BAM from {input_name} (CRC verify: {}).",
                    if check_crc { "on" } else { "off" }
                );
            }
            Reader::Bam(RawBamReader::new(reader_box, check_crc))
        }
        Format::Sam => {
            if !args.quiet {
                eprintln!("methylsieve: Reading SAM from {input_name}.");
            }
            Reader::Sam(SamReader::new(reader_box))
        }
    };

    let mut header = reader.read_header()?;
    if header.reference_sequences().is_empty() {
        bail!("Input has no @SQ reference sequences. Exiting.");
    }
    if !args.quiet {
        eprintln!(
            "methylsieve: Loaded {} header sequence entries.",
            header.reference_sequences().len()
        );
    }

    // Load the reference (every @SQ contig, 2-bit packed) and resolve control
    // contigs to a per-tid scope map. Both cross-check against the header.
    let reference = Reference::load(&args.reference, &header)
        .with_context(|| format!("loading reference {}", args.reference.display()))?;
    if !args.quiet {
        eprintln!("methylsieve: Loaded reference {}.", args.reference.display());
    }
    let (scope_of_tid, control_names) = resolve_control_scopes(&header, &args.control_contig)?;

    append_methylsieve_pg(&mut header)?;

    let write_buf_bytes = (args.write_buffer_mb as usize).saturating_mul(1024 * 1024);
    let mut out = RawBamWriter::open(
        args.output.as_deref(),
        &header,
        write_buf_bytes,
        args.compression_level,
    )?;

    // `--mbias-mask` automates and supersedes the manual `--ignore-template-ends`
    // trim: both neutralize the same physical fragment-end bias (an FR pair's two
    // outer ends are exactly R1's and R2's 5' ends), so running them together
    // would double-trim. When masking is on, the manual trim is forced off — and
    // a warning is logged if the user set it explicitly.
    let ignore_template_ends = if args.mbias_mask {
        if args.ignore_template_ends > 0 {
            log::warn!(
                "--ignore-template-ends={} is ignored because --mbias-mask is set; masking learns \
                 and removes fragment-end bias automatically",
                args.ignore_template_ends
            );
        }
        0
    } else {
        args.ignore_template_ends
    };

    let opts = ProcessorOptions {
        contexts: args.contexts,
        mode: args.mode.into(),
        max_unconverted_count: args.max_unconverted_count,
        max_unconverted_fraction: args.max_unconverted_fraction,
        min_sites: args.min_sites,
        min_base_quality: args.min_base_quality,
        ignore_template_ends,
        ignore_supplementary_evidence: args.ignore_supplementary_evidence,
        tag: args.tag.clone(),
        count_tag: args.effective_count_tag(),
        qc_fail: args.qc_fail(),
        remove_unconverted: args.remove_unconverted,
        scope_of_tid,
        record_matrix: args.metrics_prefix.is_some(),
    };
    let processor = RecordProcessor::new(reference, opts);
    let mut stats = Stats::new(&control_names);

    // M-bias accumulation is enabled only when metric output or masking is
    // requested; otherwise it stays `None` and the decision path is untouched —
    // no per-cycle bookkeeping, no measurable overhead.
    let want_mbias = args.metrics_prefix.is_some() || args.mbias_mask;
    let mut mbias_acc = want_mbias.then(MbiasAccumulator::new);
    let detect = DetectParams {
        plateau_fraction: args.mbias_plateau_fraction,
        max_mask: args.mbias_max_mask as usize,
        ..DetectParams::default()
    };

    // methylsieve requires query-grouped input: records are read in maximal
    // runs sharing a QNAME ("blocks") and each block is processed as a unit.
    let mut pool: Vec<RawRecord> = Vec::with_capacity(8);
    // The mask plan actually applied to the data (frozen from the learn-phase
    // subset). `None` when masking is off — the summary then reports no mask
    // lengths rather than a misleading full-file recomputation.
    let applied_plan = if args.mbias_mask {
        run_masking(
            &mut reader,
            &mut pool,
            &processor,
            &mut stats,
            mbias_acc.as_mut().expect("masking enables the accumulator"),
            &mut out,
            MaskingRun {
                detect,
                buffer_templates: args.mbias_buffer_templates,
                mask_quality: args.mbias_mask_quality,
                min_base_quality: args.min_base_quality,
                want_metrics: args.metrics_prefix.is_some(),
                quiet: args.quiet,
            },
        )?
    } else {
        for_each_block(&mut reader, &mut pool, |block| {
            // No masking: no stored exclusions, so pass an empty window set.
            if processor.process_block(block, &mut stats, mbias_acc.as_mut(), &[])?
                == Disposition::Keep
            {
                write_block(block, &mut out)?;
            }
            Ok(())
        })
        .context("processing record block")?;
        None
    };

    print_run_stats(&stats, &args);
    warn_proportion_blind_spot(&stats, &args);

    if let Some(prefix) = args.metrics_prefix.as_deref() {
        let mbias = mbias_acc.as_ref().expect("mbias accumulator present when metrics requested");
        crate::metrics::write_all(
            prefix,
            &stats,
            mbias,
            applied_plan.as_ref(),
            &header,
            args.sample.as_deref(),
            args.input.as_deref(),
            &args.contexts.label(),
            |unconv, monitored| processor.classify(unconv, monitored),
        )
        .context("writing metric TSVs and plots")?;
    }

    out.finish().context("finishing main output")?;
    report(&started, stats.total_templates, args.quiet);
    Ok(())
}

/// Resolve `--control-contig` names against the header into a per-tid scope map
/// (`None` → genome, `Some(i)` → `controls[i]`) plus the ordered control names.
///
/// # Errors
/// Returns an error if a named control contig is absent from the header.
fn resolve_control_scopes(
    header: &Header,
    control_contigs: &[String],
) -> Result<(Vec<Option<usize>>, Vec<String>)> {
    let n = header.reference_sequences().len();
    let name_to_tid: HashMap<String, usize> = header
        .reference_sequences()
        .keys()
        .enumerate()
        .map(|(tid, name)| (String::from_utf8_lossy(name.as_ref()).into_owned(), tid))
        .collect();

    let mut scope_of_tid = vec![None; n];
    let mut control_names = Vec::with_capacity(control_contigs.len());
    for (ci, name) in control_contigs.iter().enumerate() {
        let tid = *name_to_tid.get(name).ok_or_else(|| {
            anyhow::anyhow!(
                "--control-contig '{name}' is not present in the input @SQ header lines."
            )
        })?;
        scope_of_tid[tid] = Some(ci);
        control_names.push(name.clone());
    }
    Ok((scope_of_tid, control_names))
}

/// Drive the QNAME-block grouping loop over `reader`, invoking `on_block` for
/// each maximal run of records sharing a QNAME. `pool` allocations are reused
/// across blocks.
fn for_each_block(
    reader: &mut Reader,
    pool: &mut Vec<RawRecord>,
    mut on_block: impl FnMut(&mut [RawRecord]) -> Result<()>,
) -> Result<()> {
    if pool.is_empty() {
        pool.push(RawRecord::new());
    }
    let mut block_len: usize = 0;
    let mut current_qname: Vec<u8> = Vec::new();
    loop {
        if pool.len() == block_len {
            pool.push(RawRecord::new());
        }
        let read_idx = block_len;
        let got = reader.read_record(&mut pool[read_idx]).context("reading record")?;
        if !got {
            break;
        }
        let new_qname_differs =
            block_len > 0 && pool[read_idx].read_name() != current_qname.as_slice();
        if new_qname_differs {
            on_block(&mut pool[..block_len])?;
            if read_idx != 0 {
                pool.swap(0, read_idx);
            }
            block_len = 1;
            current_qname.clear();
            current_qname.extend_from_slice(pool[0].read_name());
        } else {
            if block_len == 0 {
                current_qname.clear();
                current_qname.extend_from_slice(pool[0].read_name());
            }
            block_len += 1;
        }
    }
    if block_len > 0 {
        on_block(&mut pool[..block_len])?;
    }
    Ok(())
}

/// Write every record of a block to `out`.
fn write_block(block: &[RawRecord], out: &mut RawBamWriter) -> Result<()> {
    for rec in block {
        out.write_record(rec).context("writing record")?;
    }
    Ok(())
}

/// Resolved configuration for a [`run_masking`] pass.
struct MaskingRun {
    /// Mask-length detection tunables.
    detect: DetectParams,
    /// Templates to buffer in the learn phase.
    buffer_templates: usize,
    /// Quality value masked bases are set to.
    mask_quality: u8,
    /// The decision's base-quality gate. A mask window is excluded from the tally
    /// only when `mask_quality < min_base_quality` — i.e. when masking actually
    /// drops the base below the gate a downstream caller would apply.
    min_base_quality: u8,
    /// Keep accumulating M-bias after the plan freezes. Only needed to write the
    /// whole-file curves under `--metrics-prefix`; when false, the post-freeze
    /// stream skips the (now-useless) second M-bias walk entirely.
    want_metrics: bool,
    /// Suppress the run-summary line.
    quiet: bool,
}

/// Mask one template and tally/decide/stamp it on the **resulting** data, so the
/// reported metrics and the unconverted call both describe what a downstream
/// caller will see. The mask geometry (`compute_mask_windows`) drives both the
/// tally exclusion and the Q2 write, so the two never disagree; the exclusion is
/// applied only when `mask_quality < min_base_quality` (otherwise a post-mask
/// tally would still count those bases). M-bias is *not* accumulated here — the
/// curve is a pre-mask measurement taken before this point. Returns whether the
/// template should be emitted; the Q2 write is skipped for dropped templates.
fn mask_and_process(
    processor: &RecordProcessor,
    stats: &mut Stats,
    mbias: Option<&mut MbiasAccumulator>,
    plan: &MaskPlan,
    mask_quality: u8,
    min_base_quality: u8,
    block: &mut [RawRecord],
) -> Result<Disposition> {
    let windows = compute_mask_windows(plan, block);
    let tally_windows: &[MaskWindows] =
        if mask_quality < min_base_quality { &windows } else { &[] };
    let disp = processor.process_block(block, stats, mbias, tally_windows)?;
    if disp == Disposition::Keep {
        apply_windows(block, &windows, mask_quality);
    }
    Ok(disp)
}

/// Two-phase M-bias masking run. **Learn:** buffer complete templates into the
/// arena while accumulating the per-cycle M-bias curve (pre-mask), *deferring*
/// the tally/decision/stamp so they can later run on the masked data. **Drain +
/// stream:** freeze the mask lengths, then for every template — buffered first,
/// then the streamed remainder — mask it, tally/decide/stamp on the masked
/// result, and emit. If the file ends before the target, the whole file was
/// buffered — drain it.
fn run_masking(
    reader: &mut Reader,
    pool: &mut Vec<RawRecord>,
    processor: &RecordProcessor,
    stats: &mut Stats,
    mbias: &mut MbiasAccumulator,
    out: &mut RawBamWriter,
    cfg: MaskingRun,
) -> Result<Option<MaskPlan>> {
    let mut arena = TemplateArena::with_target(cfg.buffer_templates);
    let mut plan: Option<MaskPlan> = None;

    for_each_block(reader, pool, |block| {
        match &plan {
            // Stream phase: tally/decide/stamp on the masked geometry, then emit.
            // M-bias keeps accumulating (pre-mask) only when the whole-file curves
            // are needed for metrics; otherwise the frozen plan never changes.
            Some(p) => {
                let feed = cfg.want_metrics.then_some(&mut *mbias);
                let disp = mask_and_process(
                    processor,
                    stats,
                    feed,
                    p,
                    cfg.mask_quality,
                    cfg.min_base_quality,
                    block,
                )?;
                if disp == Disposition::Keep {
                    write_block(block, out)?;
                }
                Ok(())
            }
            // Learn phase: accumulate M-bias (pre-mask) only and buffer the raw
            // template; the tally/decision is deferred to the drain pass below so
            // it sees the masked data. Once the buffer fills, freeze and drain.
            None => {
                processor.accumulate_mbias(block, mbias);
                let buffered = arena.push_template(block);
                if buffered && !arena.is_full() {
                    return Ok(());
                }
                let frozen = MaskPlan::learn(mbias, cfg.detect, cfg.mask_quality);
                drain_masked(&arena, processor, stats, &frozen, &cfg, out)?;
                if !buffered {
                    // This block didn't fit the (now-full) arena → process + emit it.
                    let disp = mask_and_process(
                        processor,
                        stats,
                        None,
                        &frozen,
                        cfg.mask_quality,
                        cfg.min_base_quality,
                        block,
                    )?;
                    if disp == Disposition::Keep {
                        write_block(block, out)?;
                    }
                }
                plan = Some(frozen);
                Ok(())
            }
        }
    })
    .context("processing record block")?;

    // The file ended while still learning (fewer than the target): the entire
    // input is buffered and undrained — freeze on what we have and emit it.
    if plan.is_none() {
        let frozen = MaskPlan::learn(mbias, cfg.detect, cfg.mask_quality);
        drain_masked(&arena, processor, stats, &frozen, &cfg, out)?;
        plan = Some(frozen);
    }
    if let Some(p) = &plan
        && !cfg.quiet
    {
        eprintln!(
            "methylsieve: M-bias masking applied — {} (learned from {} buffered templates)",
            p.summary(),
            arena.template_count()
        );
    }
    Ok(plan)
}

/// Drain the buffered learn-phase templates: mask each, tally/decide/stamp it on
/// the masked result (M-bias already captured during learn, so none here), and
/// write it. Preserves arrival order.
fn drain_masked(
    arena: &TemplateArena,
    processor: &RecordProcessor,
    stats: &mut Stats,
    plan: &MaskPlan,
    cfg: &MaskingRun,
    out: &mut RawBamWriter,
) -> Result<()> {
    arena.drain(|recs| {
        let disp = mask_and_process(
            processor,
            stats,
            None,
            plan,
            cfg.mask_quality,
            cfg.min_base_quality,
            recs,
        )?;
        if disp == Disposition::Keep {
            write_block(recs, out)?;
        }
        Ok(())
    })
}

/// Append methylsieve's `@PG` line to the header. noodles auto-chains via `PP:`.
fn append_methylsieve_pg(header: &mut Header) -> Result<()> {
    // Defensive: a broken PP chain makes noodles' `programs.add` panic; convert
    // that into a clean error.
    let programs = header.programs().as_ref();
    let known_ids: std::collections::HashSet<&[u8]> =
        programs.keys().map(|k| k.as_slice()).collect();
    for (id, map) in programs.iter() {
        if let Some(pp) = map.other_fields().get(&program_tag::PREVIOUS_PROGRAM_ID) {
            let pp_bytes = pp.as_ref();
            if !known_ids.contains(pp_bytes) {
                bail!(
                    "input header @PG ID:{} has PP:{} but no @PG with that ID exists. \
                     Strip the broken PP tag (e.g. via `samtools reheader`) before re-running.",
                    String::from_utf8_lossy(id.as_slice()),
                    String::from_utf8_lossy(pp_bytes),
                );
            }
        }
    }

    let cl = command_line_for_pg();
    let mut map = Map::<Program>::default();
    map.other_fields_mut().insert(program_tag::NAME, "methylsieve".into());
    map.other_fields_mut().insert(program_tag::VERSION, METHYLSIEVE_BUILD.into());
    map.other_fields_mut().insert(program_tag::COMMAND_LINE, cl.into());
    header.programs_mut().add("methylsieve", map).context("appending @PG methylsieve record")?;
    Ok(())
}

/// Build the `@PG CL:` string from `std::env::args`, with `argv[0]` reduced to
/// its basename.
fn command_line_for_pg() -> String {
    let mut args = std::env::args();
    let prog = args
        .next()
        .map(|a| {
            std::path::Path::new(&a)
                .file_name()
                .map(|s| s.to_string_lossy().into_owned())
                .unwrap_or(a)
        })
        .unwrap_or_else(|| "methylsieve".to_string());
    let rest: Vec<String> = args.collect();
    if rest.is_empty() { prog } else { format!("{prog} {}", rest.join(" ")) }
}

fn print_run_stats(stats: &Stats, args: &Args) {
    if stats.total_templates == 0 {
        eprintln!("methylsieve: No reads processed.");
        return;
    }
    let g = &stats.genome;
    let verb = if args.remove_unconverted { "Removed" } else { "Tagged " };
    let frac = if g.n_evaluated > 0 { g.n_unconverted as f64 / g.n_evaluated as f64 } else { 0.0 };
    eprintln!(
        "methylsieve: {} {:>10} of {:>10} (frac {:.5}) evaluated genome templates as unconverted.",
        verb, g.n_unconverted, g.n_evaluated, frac
    );
    if stats.unmapped_templates > 0 || stats.zero_site_templates > 0 {
        eprintln!(
            "methylsieve: {} templates had no mapped primary; {} produced no monitored sites.",
            stats.unmapped_templates, stats.zero_site_templates
        );
    }
}

/// In `proportion` mode, warn about the blind spot: templates below `--min-sites`
/// can't be evaluated by the proportion test and pass through unflagged. Other
/// modes cover those templates with the count test, so no warning is needed.
fn warn_proportion_blind_spot(stats: &Stats, args: &Args) {
    if args.quiet || args.mode != DecisionModeCli::Proportion {
        return;
    }
    if stats.below_min_sites_templates > 0 {
        eprintln!(
            "methylsieve: WARNING proportion mode — {} templates had fewer than {} sites and \
             were NOT evaluated (passed through). Use --mode adaptive or either to catch them \
             with the count test.",
            stats.below_min_sites_templates, args.min_sites
        );
    }
}

// ── End-of-run resource-usage footer ────────────────────────────────────────
//
// Mirrors the C++ samblaster footer that prints memory + timing after the
// summary. Suppressed by `--quiet`. CPU and max-RSS reads use
// `getrusage(RUSAGE_SELF)` and are Unix-only — on other platforms only wall
// time is reported.

/// Captures a starting wall-clock timestamp; pair with [`report`] at end.
struct StartedRun {
    /// Snapshot of `Instant::now()` at process start.
    wall_start: Instant,
}

impl StartedRun {
    /// Snapshot the current wall-clock timestamp.
    fn now() -> Self {
        Self { wall_start: Instant::now() }
    }
}

/// Print a single-line resource-usage footer to stderr. No-op if `quiet`.
fn report(started: &StartedRun, n_templates: u64, quiet: bool) {
    if quiet {
        return;
    }
    let wall = started.wall_start.elapsed().as_secs_f64();
    let stderr = std::io::stderr();
    let mut stderr = stderr.lock();

    #[cfg(unix)]
    if let Some(ru) = read_rusage() {
        let user = ru.user_secs;
        let sys = ru.sys_secs;
        let rss_mb = ru.max_rss_bytes as f64 / (1024.0 * 1024.0);
        let _ = writeln!(
            stderr,
            "methylsieve: Processed {n_templates} templates in {wall:.2}s wall, \
             {user:.2}s user CPU, {sys:.2}s system CPU, max RSS {rss_mb:.1} MB.",
        );
        return;
    }

    let _ = writeln!(stderr, "methylsieve: Processed {n_templates} templates in {wall:.2}s wall.");
}

#[cfg(unix)]
struct Rusage {
    user_secs: f64,
    sys_secs: f64,
    max_rss_bytes: u64,
}

/// Snapshot the current process's user/sys CPU and max RSS via
/// `getrusage(RUSAGE_SELF)`. Returns `None` if the syscall fails.
///
/// `ru_maxrss` is in **bytes** on macOS and **kilobytes** on Linux — we
/// normalize to bytes here so the caller doesn't need to care.
#[cfg(unix)]
fn read_rusage() -> Option<Rusage> {
    // SAFETY: `getrusage` writes a `rusage` struct that we just allocated.
    // Both the syscall number and the struct layout come from `libc`.
    let mut ru: libc::rusage = unsafe { std::mem::zeroed() };
    let rc = unsafe { libc::getrusage(libc::RUSAGE_SELF, &mut ru) };
    if rc != 0 {
        return None;
    }
    let user_secs = ru.ru_utime.tv_sec as f64 + ru.ru_utime.tv_usec as f64 * 1e-6;
    let sys_secs = ru.ru_stime.tv_sec as f64 + ru.ru_stime.tv_usec as f64 * 1e-6;
    let max_rss = ru.ru_maxrss as u64;
    // macOS reports bytes; Linux + BSD report kilobytes.
    #[cfg(target_os = "macos")]
    let max_rss_bytes = max_rss;
    #[cfg(not(target_os = "macos"))]
    let max_rss_bytes = max_rss.saturating_mul(1024);
    Some(Rusage { user_secs, sys_secs, max_rss_bytes })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_contexts_default_is_cph() {
        let m = parse_contexts("CpA,CpC,CpT").unwrap();
        assert!(m.contains(Context::CpA));
        assert!(m.contains(Context::CpC));
        assert!(m.contains(Context::CpT));
        assert!(!m.contains(Context::CpG));
    }

    #[test]
    fn parse_contexts_accepts_short_spelling_and_case() {
        let m = parse_contexts("ca,CG").unwrap();
        assert!(m.contains(Context::CpA));
        assert!(m.contains(Context::CpG));
        assert!(!m.contains(Context::CpC));
    }

    #[test]
    fn parse_contexts_rejects_unknown() {
        assert!(parse_contexts("CpA,CpZ").is_err());
        assert!(parse_contexts("").is_err());
    }

    #[test]
    fn parse_tag_spec_default() {
        let t = parse_tag_spec("XX:Z:UC").unwrap();
        assert_eq!(t.tag, [b'X', b'X']);
        assert_eq!(t.value, b"UC");
    }

    #[test]
    fn parse_tag_spec_rejects_non_z_and_malformed() {
        assert!(parse_tag_spec("XX:i:3").is_err());
        assert!(parse_tag_spec("X:Z:UC").is_err());
        assert!(parse_tag_spec("XX:Z:").is_err());
        assert!(parse_tag_spec("garbage").is_err());
    }

    #[test]
    fn validate_rejects_non_bam_output() {
        let mut a = minimal_args();
        a.output = Some(PathBuf::from("out.sam"));
        assert!(a.validate().is_err());
    }

    #[test]
    fn validate_allows_bam_stdout_with_metrics_prefix() {
        // Metric TSVs are file-only (prefix), so the BAM may take stdout freely.
        let mut a = minimal_args();
        a.output = None;
        a.metrics_prefix = Some(PathBuf::from("run_metrics"));
        assert!(a.validate().is_ok());
    }

    #[test]
    fn validate_rejects_bad_fraction() {
        let mut a = minimal_args();
        a.max_unconverted_fraction = 1.5;
        assert!(a.validate().is_err());
        a.max_unconverted_fraction = 0.0;
        assert!(a.validate().is_err());
        a.max_unconverted_fraction = 0.5;
        assert!(a.validate().is_ok());
    }

    fn minimal_args() -> Args {
        Args {
            input: None,
            output: None,
            compression_level: CompressionLevel::new(0).unwrap(),
            reference: PathBuf::from("ref.fa"),
            contexts: parse_contexts("CpA,CpC,CpT").unwrap(),
            mode: DecisionModeCli::Adaptive,
            max_unconverted_count: 3,
            max_unconverted_fraction: 0.05,
            min_sites: 40,
            min_base_quality: 20,
            ignore_template_ends: 0,
            ignore_supplementary_evidence: false,
            tag: parse_tag_spec("XX:Z:UC").unwrap(),
            count_tag: [b'c', b'h'],
            no_count_tag: false,
            no_qc_fail: false,
            remove_unconverted: false,
            control_contig: vec![],
            mbias_mask: false,
            mbias_buffer_templates: 500_000,
            mbias_plateau_fraction: 0.90,
            mbias_max_mask: 30,
            mbias_mask_quality: 2,
            metrics_prefix: None,
            sample: None,
            quiet: true,
            check_crc: false,
            no_check_crc: false,
            read_buffer_mb: 16,
            write_buffer_mb: 64,
        }
    }
}
