[![Build](https://github.com/fulcrumgenomics/methylsieve/actions/workflows/check.yml/badge.svg)](https://github.com/fulcrumgenomics/methylsieve/actions/workflows/check.yml)
[![Version at crates.io](https://img.shields.io/crates/v/methylsieve)](https://crates.io/crates/methylsieve)
[![Bioconda](https://img.shields.io/conda/vn/bioconda/methylsieve.svg?label=bioconda)](https://bioconda.github.io/recipes/methylsieve/README.html)
[![License](http://img.shields.io/badge/license-MIT-blue.svg)](https://github.com/fulcrumgenomics/methylsieve/blob/main/LICENSE)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21445804.svg)](https://doi.org/10.5281/zenodo.21445804)

# methylsieve

Fast, streaming cleanup of directional bisulfite-sequencing (WGBS) and EM-seq
SAM/BAM: per-template **unconverted-read filtering** and **M-bias masking**, done
correctly, in a single pass that disappears into the alignment pipe.

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
uracil and read as `T`. Two artifacts impact methylation estimates: reads
that *escaped conversion* still read `C` at unmethylated positions, and end-repair
fill-in can cause a drop in observed methylation where unmethylated Cs were 
incorporated at recessed 3' ends (**M-bias**).  Both processes affect _templates_
yet most tools evaluate _reads_ independently. Methylsieve addresses both artifacts,
correctly, in one streaming pass over query-grouped SAM/BAM:

- **Unconverted-read filtering.** It makes a single **per-template** decision from
  all of a QNAME's primary and supplementary records (using CpH cytosines) and
  propagates that call — an aux tag, the QC-fail flag, or outright removal — to
  every record of the template, rather than judging each read on its own.
- **M-bias masking** *(opt-in, `--mbias-mask`)*. It learns the per-cycle
  methylation ramp on the fly, freezes the biased 5' (and, for single-end, 3')
  lengths, and sets those bases' qualities to Q2 so base-quality-aware callers
  drop them — no clip, no coordinate/CIGAR rewrite.  Read's 3' ends are also 
  masked if/when they overlap their mate's 5' mask region.

Both run in the same pass, and methylsieve emits per-context (CpA/CpC/CpT/CpG) and
per-spike-in conversion metrics — TSVs plus M-bias and decision-matrix plots —
suitable for QC.

methylsieve runs inline in the alignment pipe with negligible overhead:

```bash
bwa-mem3 --meth genome.fa R1.fq.gz R2.fq.gz \
  | methylsieve --reference genome.fa -o - \
  | dupblaster \
  | samtools sort -o final.bam
```

methylsieve reads SAM or BAM (auto-detected) and writes BAM. Input must be
**query-grouped** — all records for a QNAME adjacent, as produced directly by
the aligner — and should be **adapter-trimmed** first: untrimmed adapter
read-through on short inserts can be force-aligned and read as spurious
unconverted Cytosines.

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
- **M-bias handled automatically, per read class.** The usual fix is to eyeball a
  per-cycle plot and hand a single fixed trim to the downstream caller, applied to
  read 1 and read 2 alike. methylsieve learns the ramp from the data and freezes an
  independent mask for read 1, read 2, and single-end reads (read 1's bias is a
  couple of cycles; read 2's can run past 20), so the right number of bases is
  neutralized — in the same pass, with no measure-then-rerun loop.
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

With `--metrics-prefix`, every template's `(monitored sites, converted sites)` pair
is binned into a density hexbin and the live decision boundary is drawn over it, so
the calls written to the BAM are legible at a glance — the fully-converted mass
hugs the diagonal, and flagged templates fall away below the boundary:

<p align="center">
  <img src="https://raw.githubusercontent.com/fulcrumgenomics/methylsieve/main/docs/img/conversion-matrix.png"
       alt="Conversion-matrix hexbin for an EM-seq library — converted templates cluster on the diagonal, the converted/unconverted boundary is drawn, and flagged templates fall below it"
       width="560">
</p>

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
warning with the count of templates that escaped.

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
  fallback. In `proportion` mode the stderr warning reports how many templates
  escaped; consider a lower floor (with the coupling above) if it dominates.
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

### Metric files (`--metrics-prefix`)

`--metrics-prefix PREFIX` writes a set of files in a single streaming pass (the
output BAM is unchanged):

- **`PREFIX.summary.tsv`** — per-context conversion summary, one row per scope (below).
- **`PREFIX.conversion-matrix.tsv`** — the per-template decision histogram.
- **`PREFIX.mbias.tsv`** — per-read-cycle methylation (the M-bias curve).
- **`PREFIX.mbias.pdf`** — the M-bias curves plotted, with any applied mask shaded.
- **`PREFIX.conversion-matrix.pdf`** — the decision histogram as a density hexbin with the converted/unconverted boundary drawn.

All rates are fractions in `[0, 1]`, never percentages.

**`summary.tsv`** has one row per scope: the `genome` scope first, then one per
`--control-contig`. Every row is *decision-basis* — overlap-deduped, end-trimmed,
and including supplementary evidence — i.e. exactly the sites the unconverted call
acted on. The applied per-read mask lengths ride along as run-level columns
(`r1_mask_5p` / `r2_mask_5p` / `se_mask_5p` / `se_mask_3p`, in sequencing cycles),
blank when masking was not run. All four contexts are reported regardless of the
decision subset, so CpG retention on a methylated control reads as a built-in
sanity check. (Per-read, per-cycle conversion lives in `mbias.tsv`.)

| column | description |
|--------|-------------|
| `sample` | Sample name — `--sample` if given, else the unique `@RG SM:` values comma-joined, else the input file stem. |
| `methylsieve_version` | Version of methylsieve that wrote the row. |
| `scope` | `genome`, or a `--control-contig` name. |
| `r1_mask_5p` / `r2_mask_5p` | Applied 5' mask length (cycles) for read 1 / read 2; blank when masking was not run. |
| `se_mask_5p` / `se_mask_3p` | Applied 5' / 3' mask length for single-end reads; blank when masking was not run. |
| `n_templates` | Templates routed to this scope. |
| `n_mapped` | Templates with at least one mapped primary alignment. |
| `n_evaluated` | Templates that produced at least one monitored site. |
| `n_unconverted` | Evaluated templates called unconverted (always 0 for control scopes). |
| `frac_unconverted` | `n_unconverted / n_evaluated` (blank when none evaluated). |
| `chimeric_to_control_templates` | Genome templates with a supplementary alignment on a control contig (genome row only). |
| `CpA_obs` … `CpG_obs` | Total monitored sites observed per context (`CpH_obs` = CpA + CpC + CpT). |
| `CpA_conv_rate` … `CpG_conv_rate` | Per-context conversion rate, `1 − unconv/total` (blank when no sites). High CpG retention on a methylated control is the sanity check. |
| `CpG_meth_rate` | CpG methylation rate, `unconv/total` (= `1 − CpG_conv_rate`) — the headline biological readout. |

**`mbias.tsv`** has one row per `(read, end, context, cycle)` with coverage:
`sample, read, end, context, cycle, n_methylated, n_total, frac_methylation,
ci_lo, ci_hi` (95% Agresti–Coull interval). `end` is `5p` for paired/orphan
reads; single-end reads also report `3p`. The chosen mask lengths are reported as
the `*_mask_*` columns of `summary.tsv` (and shaded in `mbias.pdf`).

**`conversion-matrix.tsv`** is the per-template decision histogram for the
`genome` scope — one row per observed `(checked_sites, converted_sites)` cell over
the decision contexts. Each cell's verdict is replayed from the live thresholds,
so it never drifts from the calls written to the BAM, making the decision boundary
legible.

| column | description |
|--------|-------------|
| `sample` | Sample name. |
| `checked_sites` | Monitored sites examined per template in this cell — the decision denominator. |
| `converted_sites` | Converted monitored sites (`checked_sites − unconverted`). |
| `conversion_rate` | `converted_sites / checked_sites` (blank when `checked_sites` is 0). |
| `n_templates` | Templates with this exact `(checked_sites, converted_sites)` pair. |
| `decision` | The cell's verdict: `converted` or `unconverted`. |
| `decided_by` | Which rule decided it: `too_few_sites`, `count`, `proportion`, or `either`. |

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
Default 0 (off). This trims the *tally* only; the emitted reads are unchanged.
To alter the reads themselves, see M-bias masking below.

`--mbias-mask` **supersedes** this option: both neutralize the same fragment-end
bias, so when masking is enabled `--ignore-template-ends` is forced to 0 (a
warning is logged if you set a non-zero value explicitly). Use one or the other.

### M-bias-aware masking (`--mbias-mask`)

End-repair fill-in skews the methylation calls at the first sequencing cycles of
a read — especially read 2 — so the per-cycle methylation rate ramps up to a
plateau over the first 10–25 bp. `--mbias-mask` measures that ramp and neutralizes
the biased bases:

<p align="center">
  <img src="https://raw.githubusercontent.com/fulcrumgenomics/methylsieve/main/docs/img/mbias-cfdna.png"
       alt="M-bias curves for a cfDNA EM-seq library — read 2 ramps over ~20 cycles, read 1 needs no 5' mask but declines at the 3' end from adapter read-through"
       width="760">
  <br>
  <img src="https://raw.githubusercontent.com/fulcrumgenomics/methylsieve/main/docs/img/mbias-emseq.png"
       alt="M-bias curves for a standard EM-seq library — read 1 needs no mask, read 2 only three cycles"
       width="760">
</p>

The same option, two libraries, two answers: the short-insert cfDNA library (top)
needs a 20-cycle read-2 mask and also trims read-1 3' bases that read through into
the adapter, while the standard EM-seq library (bottom) is well-behaved enough to
need no read-1 mask and only three read-2 cycles. No single `--ignore-template-ends`
value is right for both — so methylsieve learns the mask per library and per read
class instead of asking you to guess one. The dashed line marks the chosen mask
length; the shaded band is what gets masked.

1. **Learn.** Buffer the first `--mbias-buffer-templates` templates (default
   500,000) and accumulate the per-cycle CpG methylation curve for R1, R2, and
   single-end reads.
2. **Freeze.** For each read class, the mask length is the first 5' cycle whose
   smoothed methylation reaches `--mbias-plateau-fraction` of the plateau
   (default 0.90), minus one — capped at `--mbias-max-mask` (default 30). Single-end
   reads also learn a 3' length (their far template end is unknown).
3. **Mask.** Set the biased bases' qualities to `--mbias-mask-quality` (default 2)
   on every record carrying SEQ and qualities except secondary alignments (so
   primary, supplementary, and unmapped records are all masked) — the first *K*
   cycles from the 5' end, and:
   - **single-end:** also the last *K₃'* cycles;
   - **orphans** (mate unmapped): the 3' end by the mate role's length, in case the
     read ran through the whole template;
   - **proper pairs:** any 3' bases extending past the mate's post-mask 5' end.

Nothing else about the alignment changes — no clip, no coordinate/CIGAR/tag/mate
rewrite. Masked bases fall below `--min-base-quality`, so they also drop out of
methylsieve's own tally.

> **Effective only for base-quality-aware callers.** Masking lowers qualities; a
> downstream methylation caller must honor base quality for it to matter
> (MethylDackel's `-q` does; Bismark's methylation extractor does not filter by
> base quality). The mode is **not idempotent** — don't re-run it on an
> already-masked BAM, since the second pass would learn the curve from masked data.

Pair with `--metrics-prefix` to inspect the learned curve (`PREFIX.mbias.tsv`,
plotted in `PREFIX.mbias.pdf`); the chosen mask lengths are reported as the
`r1_mask_5p` / `r2_mask_5p` / `se_mask_5p` / `se_mask_3p` columns of
`PREFIX.summary.tsv`.

### Single-end reads

Single-end reads are evaluated; the monitored strand is chosen per record, so
reverse-mapped reads and supplementaries are handled correctly.

### Spike-in controls

`--control-contig NAME` (repeatable) routes reads whose primary maps to a control
(e.g. unmethylated lambda, methylated pUC19) into their own metrics-summary row and
excludes them from tagging, so the conversion rate can be read directly off the
control.

### Reference memory

The genome is preloaded once at startup, 2-bit packed (~0.8 GB for a human
genome). Packing folds non-ACGT bases (N, IUPAC ambiguity) to A — this never
changes a conversion call (only genuine C/G positions are monitored, and those
are exact) and only relabels the context of a monitored C/G immediately adjacent
to a former N (assembly-gap edges, below measurement noise).

Loading is index-driven: with a `samtools faidx` `.fai` beside the FASTA, each
contig is read by its byte span in one pass that strips newlines and packs in a
single sweep. Without a `.fai` it falls back to a slower sequential read (and
says so) — indexing is recommended for the fastest startup.

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

With spike-in controls and metrics output:

```bash
methylsieve -r genome.fa -i in.bam -o out.bam \
  --control-contig lambda --control-contig pUC19 --metrics-prefix run
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
