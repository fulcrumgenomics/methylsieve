# methylsieve conversion-failure benchmark

A reproducible Snakemake + pixi pipeline that measures both the **runtime**
(wall/CPU/RSS) and the **detection accuracy** (sensitivity / specificity / PPV)
of [methylsieve](../) against peer bisulfite/EM-seq conversion-failure callers,
across a matrix of dataset geometries.

## How it works

[holodeck](https://github.com/fg-labs/holodeck) simulates EM-seq/bisulfite reads in which a tunable
fraction of whole molecules **fail conversion** (their non-CpG cytosines are
coherently retained). With `--golden-bam` it writes a perfectly-aligned BAM that
carries, on every record, a per-molecule **`cf:i:{0|1}`** truth tag (`1` =
conversion failure) plus Bismark-style `XM`/`XR`/`XG` tags, `MD`, and `NM`.

Because the golden BAM already carries correct alignments, **we skip alignment
entirely** and run every tool directly on it. This isolates conversion-filtering
from alignment quality and gives every tool the same input; `cf` is the
per-template truth label.

```
fetch_reference ─▶ methylate (once) ─▶ simulate (per dataset)
                                            │  results/sim/{id}/{id}.golden.bam  (carries cf:i)
                                            ▼
                       run_<tool> (timed) ─▶ score ─▶ aggregate ─▶ results/benchmark.tsv
```

## Quick start

```bash
cd benchmark-pipeline
./install.sh                          # pixi env (tools) + build methylsieve
./run.sh config/smoke.config.yaml     # fast end-to-end smoke (a few minutes)
./run.sh                              # full run over config/datasets.tsv
./run.sh --dry-run                    # preview the job graph
```

`install.sh` builds **methylsieve** (this repo) from source; everything else —
the **holodeck** simulator (≥0.3) and the peer callers biscuit, NEB's
`mark-nonconverted-reads`, and (optional) bismark — comes from conda via `pixi`.

## Configuration

`config/config.yaml` holds the globals; `config/datasets.tsv` is the matrix, one
row per simulated dataset:

| column | meaning |
|---|---|
| `id` | dataset name (used in paths/table) |
| `layout` | `PE` or `SE` |
| `read_length` | holodeck `-l` |
| `fragment_mean` / `fragment_stddev` | holodeck `-d` / `-s` |
| `coverage` | holodeck `-c` |
| `conversion_rate` | holodeck `--methylation-conversion-rate` |
| `failure_rate` | holodeck `--methylation-failure-rate` (the key axis) |
| `seed` | holodeck `--seed` (per row, for independence) |

Key `config.yaml` knobs: `reference_url` (default UCSC hg38 chr21), `region_bed`
(restrict simulation for fast smoke runs), `methylation_rate` (one shared
methylome), `replicates`, `tools`, and `mixed_policy`. Override any key inline,
e.g. `./run.sh -- --config replicates=1`.

## Tools

| tool | invocation | "called unconverted" signal |
|---|---|---|
| **methylsieve** | `methylsieve -i golden.bam -o out.bam -r ref.fa --mode adaptive` | `XX:Z:UC` tag + `0x200` QC-fail on every record of the template |
| **NEB** | `mark-nonconverted-reads.py --bam golden.bam --reference ref.fa --flag_reads` | `XX:Z:UC` tag (+`0x200`) per read |
| **biscuit** | `biscuit bsconv -m <thr-1> -v ref.fa golden.bam` | membership: the `-v` output *is* the unconverted reads |
| **bismark** *(optional)* | `filter_non_conversion --paired golden.bam` | membership in `*_removed_seqs.bam` |

Each timed `run_<tool>` rule wraps **only** the tool in `env time -v`; any
normalization (e.g. `samtools view -b` on biscuit's stdout) runs untimed. Timed
rules claim the whole `bench` resource pool so no two measurements overlap.

**Tool caveats the table surfaces honestly:**
- **methylsieve** decides per-template and propagates the call to every record,
  so its `mixed` count is structurally 0.
- **NEB** only evaluates proper pairs — it **skips single-end reads entirely**
  (writes them through untagged), so it scores ~0 sensitivity on SE datasets.
- **biscuit**'s threshold is applied by biscuit itself (`-m`), so its call is its
  native conversion-failure logic, not an externally imposed one.
- **bismark** is optional (off by default): `filter_non_conversion` needs
  Bismark-native `XM`/`XR`/`XG` tags, which holodeck's golden BAM carries — it
  runs directly and scored ~1.0 sensitivity in validation. It is opt-in only
  because the conda package pulls in heavy aligners (bowtie2/hisat2/minimap2):
  enable with `./install.sh --with-bismark` (or `pixi install -e bismark`) and
  add `bismark` to `tools`.

## Scoring

`workflow/scripts/score_conversion.py` reads the golden BAM (truth `cf:i`) and a
tool's output, placing each template in a **2×3 contingency table**:

```
                    called converted   called non-converted   called mixed
truth converted (cf=0)     TN                 FP                  cf0_mixed
truth non-conv. (cf=1)     FN                 TP                  cf1_mixed
```

A PE template is `non-converted` when both mates are flagged, `converted` when
neither is, `mixed` when they disagree. SE collapses to 2×2. The raw 6 cells are
emitted so the `mixed_policy` (`converted` | `non-converted` | `separate`,
default `non-converted`) stays a downstream choice; `collate.py` folds them into
`sensitivity` / `specificity` / `ppv` / `f1`. Scoring is memory-bounded: only
flagged QNAMEs are held in a dict; the query-grouped golden BAM is streamed.

## Output: `results/benchmark.tsv`

One row per (dataset × tool × rep) with: dataset params · `tool`, `rep` ·
runtime (`wall_s`, `user_s`, `sys_s`, `cpu_percent`, `max_rss_kb`,
`exit_status`) · throughput (`reads_per_s`, `templates_per_s`, `bases_per_s`,
derived from `wall_s`) · the 6 contingency cells · derived (`n_templates`,
`sensitivity`, `specificity`, `ppv`, `f1`, `mixed_rate`, `mixed_policy`) ·
provenance (`host`, `arch`, `cpu_model`, `tool_version`, `holodeck_version`).

See **Alternative tools & benchmark results** below for headline numbers from a full
run, or `results/benchmark.tsv` for every (dataset × tool × rep) row.

## Alternative tools & benchmark results

methylsieve is benchmarked against the peer conversion-failure callers — NEB
`mark-nonconverted-reads`, biscuit `bsconv`, and bismark `filter_non_conversion`
— each run on holodeck's golden BAM (see **Tools** above for invocations).

<!-- BENCHMARK-RESULTS:START -->
Full run: UCSC hg38 chr21, **32× coverage**, hard-unmasked reference, uncompressed
pipe-realistic input, 2 replicates (median shown). Raw rows in `results/benchmark.tsv`.

### Accuracy — sensitivity / PPV

| dataset | layout | fail | methylsieve | biscuit | bismark | NEB |
|---|---|---|---|---|---|---|
| `ins350_2x150` | PE | 0.01 | 0.9996 / **0.994** | 0.9996 / 0.984 | 0.9996 / 0.984 | 0.9996 / 0.973 |
| `cfdna_2x150` | PE | 0.01 | 0.9992 / **0.998** | 0.9992 / 0.987 | 0.9992 / 0.987 | 0.9992 / 0.977 |
| `ins350_2x150_hifail` | PE | 0.20 | 0.9997 / **1.000** | 0.9996 / 0.999 | 0.9996 / 0.999 | 0.9996 / 0.999 |
| `ins350_1x150` | SE | 0.01 | 0.9985 / 0.998 | 0.9988 / 0.995 | 0.9988 / 0.995 | **0.000** / — |
| `cfdna-2nuc_1x400` | SE | 0.01 | 0.9997 / **0.980** | 0.9997 / 0.897 | 0.9997 / 0.897 | **0.000** / — |

### Performance — wall time · throughput (median of 2 reps)

| dataset | methylsieve | biscuit | bismark | NEB |
|---|---|---|---|---|
| `ins350_2x150` | 22s · 457k/s | 19s · 527k/s | 281s · 36k/s | 84s · 119k/s |
| `cfdna_2x150` | 11s · 874k/s | 18s · 541k/s | 272s · 37k/s | 81s · 123k/s |
| `ins350_2x150_hifail` | 22s · 460k/s | 24s · 417k/s | 242s · 41k/s | 293s · 34k/s |
| `ins350_1x150` | 13s · 786k/s | 19s · 534k/s | 257s · 39k/s | 48s · 208k/s |
| `cfdna-2nuc_1x400` | 9s · 431k/s | 17s · 224k/s | 259s · 14k/s | 37s · 100k/s |

**Read it as:** with a fair (uppercased) reference, every tool detects conversion failures at
~0.999 sensitivity on paired-end data — so the differentiators are:
- **precision on noisy/long reads:** methylsieve holds PPV **0.98** on the Q23→Q17 SBX-style
  `cfdna-2nuc` set where biscuit/bismark drop to **0.90**.
- **single-end support:** NEB scores **0** on both SE rows (structural — it only evaluates proper
  pairs). On a *soft-masked* reference NEB also drops to ~0.69 on PE (see warning); these numbers
  use an uppercased reference.
- **throughput:** methylsieve and biscuit are the fast tier; NEB is mid; bismark (Perl) is ~10–30×
  slower. (Memory, not shown: biscuit streams via faidx at 3–4 MB; methylsieve/NEB load chr21
  at ~150 MB.)
- biscuit and bismark apply the same "≥3 non-CpG retained C" rule and produce *identical* call sets
  on every dataset.
<!-- BENCHMARK-RESULTS:END -->

> [!WARNING]
> **NEB `mark-nonconverted-reads` silently mis-handles soft-masked (mixed-case) reference FASTAs.**
>
> It verifies each retained cytosine with a **case-sensitive** reference check
> (`reference_base == "C"`). Standard references are **soft-masked** — repeats are lowercase
> (`a/c/g/t`), ~44% of UCSC hg38 chr21 — so an unconverted C sitting on a lowercase base is
> **silently not counted**. NEB then misses roughly **one third** of true conversion-failure
> templates whenever their reads fall in repeats. Measured in this harness: NEB sensitivity
> **0.69 on a soft-masked reference vs 1.00 on the same reference uppercased** (specificity
> unchanged at 0.9999).
>
> **If you run NEB, hard-unmask (uppercase) your reference first**, e.g.
> `awk '/^>/{print;next}{print toupper($0)}' ref.fa > ref.upper.fa`. methylsieve and biscuit
> canonicalize case internally and are unaffected. This pipeline uppercases the reference in
> `fetch_reference`, so its numbers reflect conversion-detection skill rather than
> reference-case handling.

### Sharp edges — which mapped reads each tool evaluates

Each caller decides, per read, whether to use it as conversion evidence. They differ
sharply in what they silently skip, drop, or mis-handle (verified against each tool's
source and confirmed empirically):

| mapped-read class | methylsieve | NEB | biscuit | bismark |
|---|---|---|---|---|
| secondary (0x100) | not tallied | skip (untagged) | **counted** | **counted** |
| supplementary (0x800) | counted¹ | skip (untagged) | **counted** | **counted** |
| duplicate (0x400) | **counted** | skip (untagged) | **counted** | **counted** |
| qcfail (0x200) | **counted** | skip (untagged) | dropped | **counted** |
| one-end-mapped singleton | evaluated | skip (untagged) | evaluated | mis-pairs / aborts |
| single-end | evaluated | **can't** | evaluated | evaluated |
| soft-clipped bases | correct | **frame-shift bug** | correct | n/a (XM tag only) |
| low MAPQ | no filter | no filter | no filter | no filter |
| lowercase (soft-masked) ref | case-folded | **silently missed** | case-folded | n/a (no ref) |

¹ methylsieve counts supplementaries by default (`--ignore-supplementary-evidence` to
exclude). **None** of these tools filter PCR duplicates or low-MAPQ reads.

**The spectrum:** NEB skips the *most* — it only evaluates primary, proper-paired, non-dup,
non-qcfail reads; everything else is written through untagged, and unmapped reads are dropped.
biscuit and bismark skip the *least* — both **count PCR duplicates, secondary, and supplementary
alignments toward the conversion call with no dedup safeguard** (biscuit's secondary/dup/supp flag
checks are literally commented out in `bsconv.c`; bismark is entirely flag-blind, working only off
the precomputed `XM` tag). methylsieve sits between: it excludes unmapped and secondary records from
the tally and uniquely handles single-end, one-end-mapped orphans, and soft-masked references — but,
like the others, does not filter PCR duplicates or input-qcfail reads.

Other gotchas the research surfaced:
- **NEB** — besides the soft-masked-reference bug above: its strand dispatch has no `else`, so a
  proper-pair read in a non-FR orientation falls through and is **silently dropped** (not even
  written); soft-clipped reads are **mis-scored** (a trimmed coordinate list is indexed into the
  full read sequence, frame-shifting the base↔reference comparison).
- **biscuit** — the `-m/-f/-a/-c/-t` retention thresholds filter on **non-CpG retention only**;
  CpG retention never triggers a filter. `-v` prints unmapped/qcfail reads with an all-zero
  retention tag (their state is never computed). `-u` (drop unclear bs-strand) is off by default.
- **bismark** — pairs reads **strictly by consecutive line order** with no flag check, so any
  interleaved secondary/supplementary/singleton record desynchronizes and corrupts every later pair
  (or aborts); a read lacking an `XM` tag is passed through in SE but is **fatal** in PE.
- **methylsieve** — the per-template verdict is stamped on *every* record of the QNAME (including
  non-contributing mates); supplementaries are counted but **not deduped** against an overlapping
  primary; templates with zero monitored sites (or, in `proportion` mode, fewer than `--min-sites`)
  pass unflagged. Coordinate-sorted input is rejected with a hard error rather than mis-handled.

## Reuse

Structure, the `bench=100` serialized-timing idiom, `host_info.py`,
`parse_gnu_time.py`, and `install.sh`/`run.sh` follow the sibling
`benchmark-pipeline/` harnesses in `../../chelae` and `../../riker`.
