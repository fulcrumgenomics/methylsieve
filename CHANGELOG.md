# Changelog

All notable changes to methylsieve are documented here.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- **M-bias-aware masking** (`--mbias-mask`). A two-phase mode that buffers the
  first `--mbias-buffer-templates` templates (default 500,000) to learn the
  per-cycle CpG methylation curve for read 1 / read 2 / single-end reads, freezes
  5' (and, for single-end, 3') mask lengths, then sets the biased bases'
  qualities to `--mbias-mask-quality` (default 2) on every primary mapped record.
  Mask length is the first 5' cycle whose smoothed methylation reaches
  `--mbias-plateau-fraction` of the plateau (default 0.90), minus one, capped at
  `--mbias-max-mask` (default 30). Beyond the own-5' mask: single-end reads mask
  the 3' end too; orphans (mate unmapped) mirror the mate role's length onto the
  3' end; proper pairs mask any 3' bases extending past the mate's post-mask 5'
  end. Masking only lowers base qualities — no clip, no coordinate/CIGAR/tag/mate
  rewrite — so it is effective for base-quality-aware callers and is not
  idempotent. Masking off is performance-neutral with the prior release.
- **Metric TSVs** under `--metrics-prefix PREFIX`, computed in a single streaming
  pass: `PREFIX.mbias.tsv` (per-read-cycle methylation with Agresti–Coull CIs),
  `PREFIX.mbias_bounds.tsv` (suggested mask lengths), and the conversion summary
  folded over a per-read (`R1`/`R2`/`SE`) dimension.

### Changed
- **Breaking:** `--stats` and `--conversion-matrix` are replaced by a single
  `--metrics-prefix PREFIX`, which writes `PREFIX.summary.tsv` (the former
  `--stats`, now with a `read` column) and `PREFIX.conversion_matrix.tsv`. All
  metric rates are emitted as fractions in `[0, 1]`.

## [0.1.0] - 2026-06-13

### Added
- Initial implementation: streaming, query-grouped SAM/BAM in, BAM out.
- Per-template unconverted-read decision using all primary + supplementary
  records of a QNAME, with the decision (an `XX:Z:UC` aux tag and/or the
  `0x200` QC-fail flag) propagated to every record of the template.
- Per-record strand determination (`monitor_C = (R1 or unpaired) XOR reverse`),
  matching MethylDackel's `getStrand()` and correctly handling reverse-mapped
  supplementaries.
- Reference-based context determination for CpA/CpC/CpT/CpG.
- `--mode` for combining the count and proportion tests: `count`, `proportion`,
  `either`, or `adaptive` (the default — proportion at/above `--min-sites`,
  count below it, so low-site templates are still evaluated while high-site
  templates are judged on rate rather than an absolute count that over-penalizes
  long reads). `--min-sites` is the proportion floor (below it the proportion is
  unestimable and abstains) and the count↔proportion switch point in `adaptive`.
  Default thresholds: count 3, fraction 0.05, min-sites 40 (the smallest floor
  that keeps the adaptive switch continuous at those values). In `proportion`
  mode a stderr warning reports how many templates fell below the floor and
  passed through unevaluated; the `below_min_sites_templates` diagnostic in
  `--stats` exposes that population in every mode.
- `--ignore-template-ends N`: ignore the outermost N bases at each end of the
  template (fragment) when tallying — the end-repair fill-in / A-tailing–prone
  positions. Trimmed by genomic position at the 5' sequenced ends of R1 and R2
  (the fragment termini), so an overlapped end is dropped in both mates while
  interior read ends are kept; single-end / orphan reads fall back to trimming
  both ends of the read. Replaces the per-read `--ignore-5p` / `--ignore-3p`.
- Per-record count annotation `ch:Z:x/y` (on by default): x is the unconverted
  count and y the total monitored sites in the `--contexts` subset — the exact
  numerator/denominator of the decision — as a per-template aggregate stamped on
  every record, so any read carries the evidence behind its call. Rename with
  `--count-tag <NAME>` or disable with `--no-count-tag`.
- `--min-base-quality` filtering (default 20), and
  `--ignore-supplementary-evidence`.
- Spike-in `--control-contig` scopes and a multi-row per-context conversion-rate
  `--stats` TSV.
- Verified concordant with NEB `mark-nonconverted-reads` on the shared
  (proper-pair) code path.
- `--ref-encoding {bytes,nibble,twobit}`: pack the in-memory reference to trade a
  little throughput for a lot of memory. **Default is `twobit`** (2-bit, ~¼ the
  resident genome) — chosen because in an input-rate-limited pipe its small CPU
  cost is hidden while the memory saving multiplies across parallel sample
  pipelines. `twobit` folds non-ACGT bases to A, which never changes a conversion
  call and only relabels the CpH/CpG context of a monitored C/G adjacent to a
  former-N (gap edges; below measurement noise). `nibble` (4-bit, ~½ RAM) is
  bit-identical; `bytes` (1 byte/base) is fastest for a single non-rate-limited
  stream. The tally hot path is generic over a `RefCodes` trait, monomorphized
  per encoding.
- PE-overlap deduplication: reference positions covered by both mates of an
  overlapping proper pair are counted once. The overlap is split at its midpoint
  and each mate keeps the half nearer its own 5' end (where base quality is
  highest), so neither read's calls dominate the whole overlap. Improves accuracy
  for short-insert / cfDNA libraries and avoids redundant work.

[Unreleased]: https://github.com/fulcrumgenomics/methylsieve/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/fulcrumgenomics/methylsieve/releases/tag/v0.1.0
