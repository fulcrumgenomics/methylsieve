#!/usr/bin/env bash
# Run the methylsieve conversion-failure benchmark pipeline.
#
# Defaults to config/config.yaml (which points at config/datasets.tsv via the
# `datasets_tsv` key). Pass an alternate config / datasets sheet positionally:
#
#   ./run.sh                                       # full run (default config)
#   ./run.sh config/smoke.config.yaml              # fast smoke test
#   ./run.sh path/to/cfg.yaml path/to/datasets.tsv # alt config + datasets
#   ./run.sh --dry-run                             # preview the job graph
#   ./run.sh --cores 16                            # cap cores (default: all)
#
# Anything after `--` is forwarded verbatim to snakemake.
set -euo pipefail
cd "$(dirname "$0")"

CONFIG_FILE=""
DATASETS_FILE=""
DRY_RUN=0
CORES="all"
EXTRA_ARGS=()

while [[ $# -gt 0 ]]; do
  case "$1" in
    --dry-run|-n) DRY_RUN=1; shift ;;
    --cores)      CORES="$2"; shift 2 ;;
    --cores=*)    CORES="${1#--cores=}"; shift ;;
    --) shift; EXTRA_ARGS+=("$@"); break ;;
    -h|--help) sed -n '2,14p' "$0"; exit 0 ;;
    -*) echo "unknown flag: $1" >&2; exit 2 ;;
    *)
      if [[ -z "$CONFIG_FILE" ]]; then CONFIG_FILE="$1"
      elif [[ -z "$DATASETS_FILE" ]]; then DATASETS_FILE="$1"
      else echo "too many positional args: $1" >&2; exit 2; fi
      shift ;;
  esac
done

CONFIG_FILE="${CONFIG_FILE:-config/config.yaml}"
[[ -f "$CONFIG_FILE" ]] || { echo "config not found: $CONFIG_FILE" >&2; exit 1; }

# bench=100 sizes the exclusive-lock pool: every timed run rule claims all 100,
# so no two timed measurements overlap; untimed prep/score/aggregate rules
# claim 1 and run freely up to --cores between timed jobs.
SNAKE_ARGS=(
  --configfile "$CONFIG_FILE"
  --cores "$CORES"
  --resources bench=100
  --rerun-incomplete
)
if [[ -n "$DATASETS_FILE" ]]; then
  [[ -f "$DATASETS_FILE" ]] || { echo "datasets not found: $DATASETS_FILE" >&2; exit 1; }
  ABS_DATASETS="$(cd "$(dirname "$DATASETS_FILE")" && pwd)/$(basename "$DATASETS_FILE")"
  SNAKE_ARGS+=(--config "datasets_tsv=$ABS_DATASETS")
fi
[[ "$DRY_RUN" -eq 1 ]] && SNAKE_ARGS+=(-n -p)
# Bash 3.2 (macOS) errors under `set -u` on empty-array expansion; the
# ${name[@]+...} form expands only when the array is set.
SNAKE_ARGS+=(${EXTRA_ARGS[@]+"${EXTRA_ARGS[@]}"})

echo "==> pixi run snakemake ${SNAKE_ARGS[*]}"
exec pixi run snakemake "${SNAKE_ARGS[@]}"
