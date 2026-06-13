#!/usr/bin/env bash
# One-shot setup for the methylsieve conversion-failure benchmark pipeline.
#
# Idempotent: re-running only redoes out-of-date steps. It materializes the pixi
# environment — which provides the simulator (holodeck >=0.3) and the peer
# callers (biscuit, NEB's mark-nonconverted-reads, optionally bismark) plus
# samtools — and builds the one tool under test:
#   - methylsieve  the tool under test, built from this repo (../)
#
#   ./install.sh                  # pixi env + build methylsieve
#   ./install.sh --skip-build     # only materialize the pixi env(s)
#   ./install.sh --with-bismark   # also materialize the optional bismark env
#
# After this succeeds, run benchmarks with ./run.sh.
set -euo pipefail

cd "$(dirname "$0")"
PIPELINE_DIR="$PWD"
METHYLSIEVE_ROOT="$(cd .. && pwd)"

DO_BUILD=1
WITH_BISMARK=0
for arg in "$@"; do
  case "$arg" in
    --skip-build)   DO_BUILD=0 ;;
    --with-bismark) WITH_BISMARK=1 ;;
    -h|--help) sed -n '2,18p' "$0"; exit 0 ;;
    *) echo "unknown flag: $arg" >&2; exit 2 ;;
  esac
done

log() { printf '\033[1;34m==>\033[0m %s\n' "$*"; }

# ---- build toolchain ------------------------------------------------------
# Bail loudly if a C/C++ toolchain or cmake is missing — rustc build scripts
# and some of methylsieve's deps need them. We don't sudo-install.
missing=()
{ command -v cc  >/dev/null 2>&1 || command -v gcc >/dev/null 2>&1; } || missing+=("cc/gcc")
command -v cmake >/dev/null 2>&1 || missing+=("cmake")
if (( ${#missing[@]} > 0 )); then
  echo "ERROR: required build tools missing: ${missing[*]}" >&2
  echo "These are needed to build methylsieve and holodeck." >&2
  case "$(uname -s)" in
    Darwin) echo "  Install with: xcode-select --install && brew install cmake" >&2 ;;
    Linux)  echo "  e.g. sudo apt-get install -y build-essential cmake  (or dnf/yum/pacman equivalent)" >&2 ;;
  esac
  exit 1
fi

# ---- pixi -----------------------------------------------------------------
if ! command -v pixi >/dev/null 2>&1; then
  log "Installing pixi to ~/.pixi"
  curl -fsSL https://pixi.sh/install.sh | bash
  export PATH="$HOME/.pixi/bin:$PATH"
else
  log "pixi found: $(pixi --version)"
fi

log "Materializing pixi environment (default)"
pixi install
if [[ "$WITH_BISMARK" -eq 1 ]]; then
  log "Materializing optional bismark environment"
  pixi install -e bismark
fi

# ---- methylsieve (built from this repo) -----------------------------------
if ! command -v cargo >/dev/null 2>&1; then
  echo "ERROR: cargo not found on PATH. Install rustup: https://rustup.rs" >&2
  exit 1
fi
log "Using cargo: $(cargo --version)"

if [[ "$DO_BUILD" -eq 1 ]]; then
  log "Building methylsieve (release) from $METHYLSIEVE_ROOT"
  cargo build --release --manifest-path "$METHYLSIEVE_ROOT/Cargo.toml"
fi

cat <<EOF

Setup complete.

Next:
  1. (Optional) edit config/config.yaml and config/datasets.tsv.
  2. Smoke test first (fast: one small dataset over a chr21 sub-region):
       ./run.sh config/smoke.config.yaml
  3. Full run:
       ./run.sh
     Preview the job graph any time with:  ./run.sh --dry-run
EOF
