# Contributing to methylsieve

Thanks for your interest in improving methylsieve.

## Development

methylsieve targets the stable Rust toolchain pinned in `rust-toolchain.toml`
(edition 2024). Its IO uses dedicated `ringbuf`-backed threads with
`fgumi-raw-bam` zero-copy records; the bisulfite/EM-seq logic lives in
`src/sieve.rs` and `src/reference.rs`.

### Verification suite

Run all three before sending a change (these mirror CI):

```bash
cargo ci-fmt    # rustfmt --check
cargo ci-lint   # clippy -D warnings
cargo ci-test   # nextest (or `cargo test`)
```

### NEB concordance

`dev/neb_concordance.py` cross-checks methylsieve against NEB's
`mark-nonconverted-reads` on synthetic data. It needs `pysam` and `samtools`:

```bash
python3 -m venv venv && venv/bin/pip install pysam
cargo build --release
venv/bin/python dev/neb_concordance.py
```

## Style & testing

- Idiomatic Rust; meaningful names; doc comments on public items; comments
  explain *why*, not *what*.
- Generate test data programmatically — never commit fixture files. Build SAM
  inputs and tiny indexed FASTAs inline via the helpers in `tests/helpers`.
- Prefer many small focused tests over table-driven ones. Cover expected
  results, error conditions, and boundary cases.
- When matching the behavior of an existing tool, match correctness but feel
  free to improve the interface and output.
