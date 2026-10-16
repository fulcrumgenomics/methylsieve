# Changelog

All notable changes to methylsieve are documented here.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- **TAPS and Illumina 5-base support** (`--chemistry`). `--chemistry` takes `em-seq` (the default), `bisulfite`, `taps` or `illumina-5base`. Under `taps` and `illumina-5base`, which convert methylated cytosines, M-bias masking and metrics work as for EM-seq, with every methylation figure read in the right polarity. There is no unconverted-read filtering, because a read from those chemistries carries nothing to test it against. The decision and unconverted-template options are still accepted so a pipeline can switch chemistry with one flag. Any given on the command line draws a warning that they have no effect. On COLO829 TAPS+ data, masking learns a 15-cycle read 2 mask. Read as retention, as an EM-seq run would, the same curve gives no mask at all.
- **Chemistry check.** A run stops with an error when the genome's conversion pattern plainly belongs to the other kind of chemistry. TAPS-like means CpH conversion below 5% with CpG conversion of at least 10%. EM-seq-like means CpH conversion of at least 80%. Anything else is ambiguous and never stops the run, so an unconverted or failed library is never mistaken for the wrong chemistry. The check waits for 100,000 CpH and 10,000 CpG genome sites and stops the run as soon as they are met. With `--mbias-mask` it first runs when the mask is learned, so a mismatch stops the run before any record is written whenever the buffered templates clear those floors. Smaller inputs get a warning instead.
- **Conversion matrices for spike-in controls.** With `--metrics-prefix`, every `--control-contig` now gets its own rows in `conversion-matrix.tsv` and its own page in `conversion-matrix.pdf`, which becomes a multi-page PDF. Each control counts the sites that should convert in it: the decision contexts under `em-seq` / `bisulfite`, as in unmethylated lambda, and CpG under `taps` / `illumina-5base`, as in a CpG-methylated control. Under `em-seq` / `bisulfite` the converted/unconverted boundary is drawn on every page. Under `taps` / `illumina-5base` there is no boundary, because their under-conversion is spread across most molecules rather than concentrated in failed ones; those pages draw the `x = y` line of complete conversion instead.

### Changed
- `summary.tsv` has a new `chemistry` column after `methylsieve_version`, so scripts that read its columns by position need updating.
- `conversion-matrix.tsv` has new `scope` and `context` columns after `sample`, so scripts that read its columns by position need updating. Its `decision` / `decided_by` columns are blank under `taps` / `illumina-5base`.
- Under `em-seq` / `bisulfite`, a control row's `n_unconverted` and `frac_unconverted` in `summary.tsv` now count the control templates past the boundary, instead of always 0. Control templates are still never tagged or removed.
- The M-bias plot's y-axis is labelled "Apparent Methylation" and its title names the chemistry.
- In `--help`, `--min-base-quality`, `--ignore-template-ends` and `--ignore-supplementary-evidence` move to a new "Tallying" section, and `--count-tag` / `--no-count-tag` move under "Decision".

### Fixed
- A run that fails partway no longer leaves an output BAM that passes as complete. The BGZF end-of-file marker is now written only when the run succeeds, so `samtools quickcheck` and downstream readers see a failed run's partial output as truncated.
- The conversion-matrix plot no longer leaves white holes in dense regions or draws two neighbouring cells summed into one brighter hex. The counts sit on a whole-number grid that the hexagonal grid can't line up with; each cell holding at least 4 templates is now spread evenly over its unit square before binning, so every hex collects templates in proportion to the area it covers. Sparser cells stay single points, so a lone template still draws as one hex, and cells on a zero edge are not spread below it. The axes also start just below 0, so the row of templates with no converted sites draws whole.
- The conversion-matrix boundary is traced from the decision at every site count rather than through the observed cells, so it no longer breaks or runs flat across site counts no template has, and it still shows on a page where every template is past it.

## [0.2.1] - 2026-07-30

### Changed
- **Substantially faster.** Reading uncompressed BAM and writing uncompressed BAM with `--mbias-mask --metrics-prefix`, both wall time and CPU time drop by roughly **40–45%** — a throughput increase of about **70–80%** — for byte-identical output. Three changes account for it: each record's CIGAR / packed-SEQ / QUAL layout is decoded **once per record** instead of being re-derived from the record bytes on every monitored cytosine; the reference scan tests **32 bases per 64-bit word**, which removes a poorly-predicted branch from every aligned base rather than paying one per base to reject the ~79% that are not the monitored C/G; and the `--count-tag` value is rendered into a fixed buffer instead of an allocated `String` per template. Nothing about the output, the CLI, or the metric files changes. The gain is dominated by work removed per base and per record, so the proportion should hold across architectures, but the exact figure will vary with hardware, read length, and how much of each record is auxiliary data.

### Fixed
- The conversion tag (`--tag`) and count tag (`--count-tag`) are now **updated** rather than appended-only-if-absent, fixing two defects in how methylsieve handled a record that already carried the tag name.
  - A count tag inherited from an earlier run was left in place, even though this run's counts differ — masking lowers base qualities below `--min-base-quality`, so a second pass tallies fewer sites. The tag documents the numerator and denominator of the decision, so a re-run left it describing a different decision than the flag stamped beside it.
  - When another tool had already used the tag name with a **non-string** type the check did not fire at all, because a string lookup does not match a non-`Z` value. An input carrying `ch:i:5` produced output carrying both `ch:i:5` and `ch:Z:u/n` — two fields sharing one tag, which the SAM specification does not permit.

## [0.2.0] - 2026-07-28

### Added
- **M-bias-aware masking** (`--mbias-mask`). A two-phase mode that buffers the
  first `--mbias-buffer-templates` templates (default 500,000) to learn the
  per-cycle CpG methylation curve for read 1 / read 2 / single-end reads, freezes
  5' (and, for single-end, 3') mask lengths, then sets the biased bases'
  qualities to `--mbias-mask-quality` (default 2) on every record carrying SEQ
  and qualities except secondary alignments (primary, supplementary, and unmapped
  records are all masked). Mask length is the first 5' cycle whose smoothed
  methylation reaches `--mbias-plateau-fraction` of the plateau (default 0.90),
  minus one, capped at `--mbias-max-mask` (default 30). Beyond the own-5' mask:
  single-end reads mask the 3' end too; reads whose mate is unmapped or absent
  mirror the mate role's 5' length onto their 3' end; proper pairs mask any 3'
  bases extending past the mate's post-mask 5' end. Masking only lowers base
  qualities — no clip, no coordinate/CIGAR/tag/mate rewrite — so it is effective
  for base-quality-aware callers and is not idempotent. Masking off is
  performance-neutral with the prior release.
- **Metric files** under `--metrics-prefix PREFIX`, computed in a single
  streaming pass: `PREFIX.summary.tsv` (one decision-basis per-context conversion
  row per scope — the genome, then each `--control-contig` — with the applied
  mask lengths as `r1_mask_5p` / `r2_mask_5p` / `se_mask_5p` / `se_mask_3p`
  columns), `PREFIX.conversion-matrix.tsv` (the per-template decision histogram),
  and `PREFIX.mbias.tsv` (per-read-cycle methylation with Agresti–Coull CIs),
  plus PDF plots `PREFIX.mbias.pdf` and `PREFIX.conversion-matrix.pdf`.

### Changed
- **Breaking:** `--stats` and `--conversion-matrix` are replaced by a single
  `--metrics-prefix PREFIX`, which writes `PREFIX.summary.tsv` (one per-context
  conversion row per scope, with the applied mask-length columns),
  `PREFIX.conversion-matrix.tsv`, and `PREFIX.mbias.tsv`. All metric rates are
  emitted as fractions in `[0, 1]`.
- `--mbias-mask` supersedes `--ignore-template-ends`: when masking is enabled the
  manual end-trim is forced to 0 (both target the same fragment-end bias), with a
  warning logged if a non-zero value was set explicitly.
- Reference loading is now index-driven and single-pass: with a `.fai` beside the
  FASTA, each contig is read by its byte span in one bulk read that strips
  newlines (by the index line geometry) and packs to 2-bit in a single sweep,
  building each packed byte from its four bases in a register, instead of the
  previous read → encode → pack three-pass loop with a per-base read-modify-write.
  An unindexed FASTA falls back to a sequential read (logged). On a human genome
  this cuts startup user-CPU ~4.6× (≈8.5 s → ≈1.8 s) and peak load memory ~22%
  (≈1.19 GB → ≈0.93 GB), and removes the direct `noodles-core` dependency.

### Removed
- **Breaking:** `--ref-encoding` and the `bytes` (1 byte/base) and `nibble`
  (4-bit) reference layouts. methylsieve now always uses the 2-bit packing that
  was the default, dropping the unused layout options and the code generic over
  them. The 2-bit fold of non-ACGT bases to A is unchanged (see "Reference
  memory" in the README).

### Fixed
- The threaded BAM writer no longer hangs if the downstream consumer closes the
  pipe (EPIPE) or the disk fills mid-write: the IO error is now surfaced from the
  next `write`/`finish` instead of blocking forever on a full ring.
- M-bias learning and `mbias.tsv` exclude `--control-contig` reads, so spike-in
  methylation (e.g. methylated pUC19 / unmethylated lambda) no longer skews the
  learned 5'/3' mask lengths or the reported per-cycle curve. Control reads are
  still tallied into their own summary rows.
- `--mbias-mask`: a record with more than eight combined exclusion intervals
  (many propagated mate windows on a chimeric template) no longer has the excess
  intervals silently counted in the decision tally while masked in the BAM — all
  intervals are honored.
- `too_few_sites` (the `decided_by` label in `conversion-matrix.tsv`) now marks
  only templates the applied test genuinely could not evaluate. For count-based
  modes (`count`, and the count fallback of `either`/`adaptive`) that means fewer
  monitored sites than `--max-unconverted-count` — where the count can never be
  reached — rather than merely below the proportion floor `--min-sites`.

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

[Unreleased]: https://github.com/fulcrumgenomics/methylsieve/compare/v0.2.1...HEAD
[0.2.1]: https://github.com/fulcrumgenomics/methylsieve/compare/v0.2.0...v0.2.1
[0.2.0]: https://github.com/fulcrumgenomics/methylsieve/compare/v0.1.0...v0.2.0
[0.1.0]: https://github.com/fulcrumgenomics/methylsieve/releases/tag/v0.1.0
