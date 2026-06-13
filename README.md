[![Build](https://github.com/fulcrumgenomics/methylsieve/actions/workflows/check.yml/badge.svg)](https://github.com/fulcrumgenomics/methylsieve/actions/workflows/check.yml)
[![Version at crates.io](https://img.shields.io/crates/v/methylsieve)](https://crates.io/crates/methylsieve)
[![Bioconda](https://img.shields.io/conda/vn/bioconda/methylsieve.svg?label=bioconda)](https://bioconda.github.io/recipes/methylsieve/README.html)
[![License](http://img.shields.io/badge/license-MIT-blue.svg)](https://github.com/fulcrumgenomics/methylsieve/blob/main/LICENSE)

# methylsieve

Fast, streaming per-template tagging and filtering of **unconverted reads** in
directional bisulfite-sequencing (WGBS) and EM-seq SAM/BAM files.

<p>
<a href="https://fulcrumgenomics.com">
<picture>
  <source media="(prefers-color-scheme: dark)" srcset="https://raw.githubusercontent.com/fulcrumgenomics/methylsieve/main/.github/logos/fulcrumgenomics-dark.svg">
  <source media="(prefers-color-scheme: light)" srcset="https://raw.githubusercontent.com/fulcrumgenomics/methylsieve/main/.github/logos/fulcrumgenomics-light.svg">
  <img alt="Fulcrum Genomics" src="https://raw.githubusercontent.com/fulcrumgenomics/methylsieve/main/.github/logos/fulcrumgenomics-light.svg" height="100">
</picture>
</a>
</p>

[Visit us at Fulcrum Genomics](https://www.fulcrumgenomics.com) to learn more about how we can power your Bioinformatics with methylsieve and beyond.

<a href="mailto:contact@fulcrumgenomics.com?subject=[GitHub inquiry]"><img src="https://img.shields.io/badge/Email_us-%2338b44a.svg?&style=for-the-badge&logo=gmail&logoColor=white"/></a>
<a href="https://www.fulcrumgenomics.com"><img src="https://img.shields.io/badge/Visit_Us-%2326a8e0.svg?&style=for-the-badge&logo=wordpress&logoColor=white"/></a>

In bisulfite sequencing and EM-seq, unmethylated cytosines are converted to
uracil and read as `T`; a read that escaped conversion still reads `C` at
unmethylated positions and silently inflates methylation estimates. methylsieve
reads query-grouped SAM/BAM, makes a single **per-template** decision from all of
a QNAME's primary and supplementary records (using CpH cytosines), and propagates
that decision — an aux tag, the QC-fail flag, or outright removal — to every record
of the template. It also emits a per-context (CpA/CpC/CpT/CpG) and per-spike-in
conversion-rate TSV suitable for MultiQC.

methylsieve runs inline in the alignment pipe with negligible overhead:

```bash
bwa-meth.py --reference genome.fa R1.fq.gz R2.fq.gz \
  | methylsieve --reference genome.fa -o - \
  | dupblaster \
  | samtools sort -o final.bam
```

methylsieve reads SAM or BAM (auto-detected) and writes BAM. Input must be
**query-grouped** — all records for a QNAME adjacent, as produced directly by
the aligner. See `methylsieve --help` for the full option list.

## Motivation

Calling conversion failure *looks* trivial — count the cytosines a read retained,
apply a threshold — and the scripts most pipelines reach for do exactly that. They
are right most of the time. methylsieve is built for the rest of the time: long
reads, single-end and supplementary alignments, high-failure libraries, and
soft-masked references, where a quick per-read count quietly does the wrong thing.
It takes the functionally correct approach, so its calls hold at the edges and not
just in the common case.

- **One decision per template, from all the evidence.** A read pair is one
  molecule. methylsieve pools the monitored cytosines across both mates
  (de-duplicated where they overlap) and every supplementary alignment, decides
  once, and stamps that call on every record — rather than judging each read on
  its own, double-counting the overlap, or calling a pair's two mates differently.
- **Scales to long reads.** A fixed "≥3 unconverted" cutoff over-penalizes longer
  reads, which carry more sites and more sequencing-error noise. The adaptive rate
  test keeps precision high as reads grow — PPV holds on 1×400 reads where an
  absolute count loses ground — and throughput stays flat.
- **No sharp edges.** It evaluates single-end reads and orphans, and stays fast
  even when a large fraction of reads fail conversion — cases where simpler
  per-read approaches tend to break down.
- **Disappears into a pipe.** Input and output each run on a dedicated IO thread
  with a large ring buffer ahead of the worker (`--read-buffer-mb`,
  `--write-buffer-mb`). Those buffers soak up bursty output from the aligner
  upstream and brief stalls from the sorter downstream, so methylsieve keeps the
  pipe flowing instead of becoming the stage everything waits on — and the
  per-record tag it adds stays hidden behind the aligner's latency.

## Benchmarks

methylsieve ships with a reproducible [Snakemake + pixi benchmark pipeline][bench]
that simulates labelled bisulfite/EM-seq libraries — holodeck golden BAMs carrying
a per-read `cf` conversion-failure truth tag — and scores each tool's calls
against that truth, measuring throughput and accuracy against NEB
`mark-nonconverted-reads`, biscuit `bsconv`, and bismark `filter_non_conversion`.
Five datasets span paired- and single-end layouts, short cfDNA and di-nucleosome
inserts, read lengths to 400 bp, and a high (20%) conversion-failure library.

<p align="center">
  <img src="https://raw.githubusercontent.com/fulcrumgenomics/methylsieve/main/docs/img/throughput.png"
       alt="Throughput by tool and dataset — methylsieve is fastest on four of five datasets and stays fast where the others fall off"
       width="700">
</p>

methylsieve is the fastest caller on four of the five datasets and competitive on
the fifth, and it stays fast where the others struggle: on the 20%-failure library
NEB falls to ~34k reads/s — the slowest there — while methylsieve holds ~460k, and
bismark's Perl runs an order of magnitude behind throughout. biscuit is fast, but
only on coordinate-sorted input; the sort it requires is not counted in its time
here.

**Accuracy — sensitivity / PPV (mean of 2 replicates):**

| dataset | layout | fail | methylsieve | biscuit | bismark | NEB |
|---|:---:|:---:|---|---|---|---|
| 350 bp insert | PE | 1% | 0.9996 / **0.994** | 0.9996 / 0.984 | 0.9996 / 0.984 | 0.9996 / 0.973 |
| cfDNA | PE | 1% | 0.9992 / **0.998** | 0.9992 / 0.987 | 0.9992 / 0.987 | 0.9992 / 0.977 |
| 350 bp insert | PE | 20% | 0.9997 / **1.000** | 0.9996 / 0.999 | 0.9996 / 0.999 | 0.9996 / 0.999 |
| 350 bp insert | SE | 1% | 0.9985 / **0.998** | 0.9988 / 0.995 | 0.9988 / 0.995 | **0.000** / — |
| cfDNA di-nucleosome | SE | 1% | 0.9997 / **0.980** | 0.9997 / 0.897 | 0.9997 / 0.897 | **0.000** / — |

Sensitivity is near-identical across tools on paired-end data, so precision (PPV)
is the differentiator — methylsieve holds the highest PPV on every dataset, most
visibly on the noisy 1×400 di-nucleosome reads (0.980 vs 0.897). On single-end
data the difference is categorical: NEB evaluates only proper pairs, so it skips
single-end reads entirely and detects nothing.

See the [benchmark pipeline][bench] for the full methodology, dataset matrix,
fairness notes, and a per-tool "sharp edges" matrix of which mapped reads each
caller actually evaluates.

[bench]: https://github.com/fulcrumgenomics/methylsieve/tree/main/benchmark-pipeline

## How the decision works

For each template, methylsieve tallies the **monitored cytosines** in the
threshold contexts (`--contexts`, default CpH = CpA,CpC,CpT) across all mapped
primary and supplementary records, counting each overlapped reference position
once. From two numbers — `u` unconverted of `n` total monitored sites — it applies
up to two tests:

- **Count test** — flag when `u ≥ --max-unconverted-count` (default 3). Simple
  and absolute; matches NEB's `mark-nonconverted-reads`.
- **Proportion test** — flag when `u / n > --max-unconverted-fraction`
  (default 0.05), *but only* when `n ≥ --min-sites` (default 40). Below that
  floor the proportion is unestimable and the test **abstains**.

The `--mode` option selects how they combine:

| mode | rule | use it when |
|------|------|-------------|
| `count` | count test only | NEB-compatible behavior, or a simple absolute cutoff |
| `proportion` | proportion test only (abstains below `--min-sites`) | not recommended; use `adaptive` instead |
| `either` | flag if **either** test fires | the most aggressive catch |
| `adaptive` *(default)* | proportion at/above `--min-sites`, count below it | read/insert lengths vary — the usual case |

### Why `adaptive` is the default

An absolute count over-penalizes long reads and read pairs: a 100-site read pair
with 4 unconverted C's (4%) is well converted, yet a count threshold of 3 would
flag it. A proportion handles this correctly but cannot be estimated from a read
with only a handful of sites. `adaptive` applies the rate where there is enough
signal (`n ≥ min_sites`) and falls back to the absolute count below that, so
short and long templates are both judged sensibly.

This means `adaptive` can be **more lenient than `count`** on long reads, which is
the intended, statistically correct behavior — but may be unexpected for someone
who expects `--max-unconverted-count` to always fire. The `ch` tag (below) makes
any individual call explainable.

### The `proportion`-mode blind spot

Because the proportion test abstains below `--min-sites`, **`proportion` mode
never flags a template with fewer than `--min-sites` sites**, no matter how badly
converted — those reads pass through. `count`, `either`, and `adaptive` still
evaluate them via the count test. In `proportion` mode methylsieve prints a stderr
warning with the count of templates that escaped, and the
`below_min_sites_templates` column in `--stats` exposes that population in every
mode.

## Defaults and tuning

| flag | default | meaning |
|------|---------|---------|
| `--mode` | `adaptive` | how the two tests combine |
| `--max-unconverted-count` | `3` | count-test threshold |
| `--max-unconverted-fraction` | `0.05` | proportion-test threshold |
| `--min-sites` | `40` | proportion floor / count↔proportion switch point |
| `--min-base-quality` | `20` | skip lower-quality bases when tallying |
| `--contexts` | `CpA,CpC,CpT` | contexts counted toward the decision |

These are coupled, not independent. The crossover where the count and proportion
tests agree is roughly `min_sites ≈ count / fraction`. With count 3 and fraction
0.05 that crossover lands at ~40–60 sites, and **40 is the smallest floor that
keeps the `adaptive` switch continuous** — so a read does not flip its call as it
crosses the threshold by one site. If you lower `--min-sites`, also raise
`--max-unconverted-fraction` or lower `--max-unconverted-count` to keep them
consistent.

Two caveats:

- **Short-insert / cfDNA libraries** have few monitored sites per template (a
  130 bp fragment yields ~25 CpH sites after quality filtering and overlap
  dedup), so most templates fall below `--min-sites = 40` and lean on the count
  fallback. Watch `below_min_sites_templates` in `--stats`; consider a lower
  floor (with the coupling above) if it dominates.
- **Samples with real non-CpG methylation** (plant, neuronal, embryonic-stem
  tissues, where mCH can reach a few percent) will read genuine methylation as
  "unconverted" at fraction 0.05. Use `--mode count` or a higher fraction there.

These are evidence-informed starting points — a survey of peer tools plus the
continuity constraint above — not fixed prescriptions; tune them against your own
controls.

## Outputs

### Marking and removing unconverted reads

A template judged unconverted receives, on *every* one of its records (including
secondaries and supplementaries):

- the aux tag `XX:Z:UC` (configurable via `--tag TAG:Z:VALUE`), and
- the `0x200` QC-fail flag (disable with `--no-qc-fail`).

Use `--remove-unconverted` to drop unconverted templates from the output
entirely instead of marking them.

### Evidence tag (`ch:Z:u/n`)

On by default, every output record is annotated with `u/n` — the unconverted
count and total monitored sites in the `--contexts` subset, i.e. the exact
numerator and denominator behind the decision. It is a **per-template aggregate**
(combining both mates after overlap dedup and trimming) stamped on every read, so
you can see *why* any single read was or was not flagged — particularly useful for
the adaptive-leniency case above. Rename it with `--count-tag <NAME>` or disable
it with `--no-count-tag`. It adds ~8 bytes per record and ~10% CPU; in the
rate-limited alignment pipe that cost is hidden behind the aligner.

### Stats TSV (`--stats`)

One row for the `genome` scope (everything off a control contig) plus one per
`--control-contig`. All four contexts are reported regardless of the decision
subset, so CpG retention on a methylated control reads as a built-in sanity
check. The last four columns are whole-run diagnostics, populated on the
`genome` row only.

| column | description |
|--------|-------------|
| `sample` | Sample name — `--sample` if given, else the unique `@RG SM:` values comma-joined. |
| `methylsieve_version` | Version of methylsieve that wrote the row. |
| `scope` | `genome`, or a `--control-contig` name. |
| `n_templates` | Templates routed to this scope (mapped, unmapped, and zero-site alike). |
| `n_evaluated` | Templates that produced at least one monitored site. |
| `n_unconverted` | Evaluated templates called unconverted (always 0 for control scopes). |
| `n_removed` | Unconverted templates dropped from the output (`--remove-unconverted`). |
| `frac_unconverted` | `n_unconverted / n_evaluated` (blank when none evaluated). |
| `CA_unconv` / `CA_total` | Unconverted and total monitored CpA sites in the scope. |
| `CC_unconv` / `CC_total` | The same, for CpC. |
| `CT_unconv` / `CT_total` | The same, for CpT. |
| `CG_unconv` / `CG_total` | The same, for CpG (reported, but excluded from the default decision). |
| `conv_rate_CpA` | CpA conversion rate, `1 − CA_unconv/CA_total` (blank when no sites). |
| `conv_rate_CpC` | CpC conversion rate. |
| `conv_rate_CpT` | CpT conversion rate. |
| `conv_rate_CpG` | CpG conversion rate (genuine methylation lives here; high retention on a methylated control is the sanity check). |
| `chimeric_to_control_templates` | Genome templates with a supplementary alignment on a control contig. |
| `unmapped_templates` | Templates whose primary R1 was unmapped (never tallied or decided). |
| `zero_site_templates` | Templates with no monitored sites (decided converted by default). |
| `below_min_sites_templates` | Genome templates with evidence but fewer than `--min-sites` sites (the proportion test cannot evaluate them). |

### Conversion matrix (`--conversion-matrix`)

The per-template decision histogram for the `genome` scope — one row per observed
`(checked_sites, unconverted_sites)` cell over the decision contexts. Each cell's
verdict is replayed from the live thresholds, so it never drifts from the calls
written to the BAM, which makes the decision boundary legible: the converted and
conversion-failed populations are visible in the (sites, retained) plane, along
with exactly where the `--mode`/threshold cut falls.

| column | description |
|--------|-------------|
| `sample` | Sample name (as in `--stats`). |
| `checked_sites` | Monitored sites examined per template in this cell — the decision denominator (`n`). |
| `unconverted_sites` | Unconverted monitored sites — the decision numerator (`u`). |
| `conversion_rate` | `1 − unconverted_sites/checked_sites` (blank when `checked_sites` is 0). |
| `n_templates` | Templates with this exact `(checked_sites, unconverted_sites)` pair. |
| `decision` | The cell's verdict: `converted` or `unconverted`. |
| `decided_by` | Which rule decided it: `count`, `proportion`, `min_sites_floor`, `zero_sites`, or `either`. |

## Other behavior

### Overlapping read pairs

For a proper pair whose mates overlap on the reference (common with short inserts
/ cfDNA), each overlapped reference position is counted once. The overlap is split
at its midpoint and each mate keeps the half nearer its own 5' end (where base
quality is highest), so neither read's calls dominate the shared region.

### Template-end trimming

`--ignore-template-ends N` ignores the outermost N bases at each end of the
original fragment — the positions prone to end-repair fill-in and A-tailing
artifacts. For a mapped pair these are the 5' sequenced ends of R1 and R2, trimmed
by *genomic* position, so an overlapped end is dropped in both mates while the
reads' interior ends remain; single-end / orphan reads have both ends trimmed.
Default 0 (off).

### Single-end reads

Single-end reads are evaluated; the monitored strand is chosen per record, so
reverse-mapped reads and supplementaries are handled correctly.

### Spike-in controls

`--control-contig NAME` (repeatable) routes reads whose primary maps to a control
(e.g. unmethylated lambda, methylated pUC19) into their own `--stats` row and
excludes them from tagging, so the conversion rate can be read directly off the
control.

### Reference memory

The genome is preloaded. `--ref-encoding` selects the layout: `twobit` (2-bit,
default) uses ~¼ the memory (~0.8 GB for a human genome), `nibble` (4-bit) ~½, and
`bytes` (1 byte/base) is fastest. `twobit` folds non-ACGT bases (N, IUPAC) to A,
which never changes a conversion call and only relabels the context of a monitored
C/G adjacent to a former N; `nibble` and `bytes` preserve context labeling
exactly.

## Examples

Inline in the alignment pipe, default adaptive mode, uncompressed BAM for speed:

```bash
bwa-meth.py --reference genome.fa R1.fq.gz R2.fq.gz \
  | methylsieve --reference genome.fa -o - \
  | dupblaster | samtools sort -o final.bam
```

NEB-compatible "3 or more unconverted CpH" rule (count only):

```bash
methylsieve -r genome.fa -i in.bam -o out.bam --mode count --max-unconverted-count 3
```

Strict rate-based filtering for uniformly long reads, removing failures outright:

```bash
methylsieve -r genome.fa -i in.bam -o out.bam \
  --mode proportion --max-unconverted-fraction 0.03 --remove-unconverted
```

With spike-in controls and a stats TSV:

```bash
methylsieve -r genome.fa -i in.bam -o out.bam \
  --control-contig lambda --control-contig pUC19 --stats stats.tsv
```

Inspect why a specific read was or was not flagged — the `ch` tag shows the
evidence, the `XX` tag and `0x200` flag show the verdict:

```bash
samtools view out.bam | grep '<read-name>'
# ...  ch:Z:6/52   →   6 unconverted of 52 monitored sites
# flagged reads also carry XX:Z:UC and the QCFAIL (0x200) bit
```

## License

MIT — see [LICENSE](https://github.com/fulcrumgenomics/methylsieve/blob/main/LICENSE).
