//! Integration-test helpers: build SAM inputs and tiny indexed FASTAs
//! programmatically, run the `methylsieve` binary, and decode its BAM output
//! back to records (via noodles, no `samtools` shell-out).
#![allow(dead_code)] // helpers used by some tests but not all

use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::Command;

use noodles_sam::{
    self as sam,
    alignment::{RecordBuf, record::data::field::Tag},
};
use noodles_util::alignment::io::reader::Builder as AlignmentReaderBuilder;
use tempfile::TempDir;

/// Path to the methylsieve binary built by cargo.
pub fn methylsieve_binary() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_methylsieve"))
}

// ── SAM building ────────────────────────────────────────────────────────────

/// Builder for a SAM file in memory; writes to a temp file on demand.
pub struct SamBuilder {
    pub header: String,
    pub records: Vec<String>,
}

impl SamBuilder {
    pub fn new() -> Self {
        Self { header: String::from("@HD\tVN:1.6\tSO:unsorted\n"), records: Vec::new() }
    }

    /// Add an `@SQ` line.
    pub fn sq(mut self, name: &str, length: usize) -> Self {
        use std::fmt::Write as _;
        let _ = writeln!(self.header, "@SQ\tSN:{name}\tLN:{length}");
        self
    }

    /// Add an `@RG` line with a sample.
    pub fn rg(mut self, id: &str, sample: &str) -> Self {
        use std::fmt::Write as _;
        let _ = writeln!(self.header, "@RG\tID:{id}\tSM:{sample}");
        self
    }

    /// Append a SAM record (11 standard fields; optional trailing aux fields).
    #[allow(clippy::too_many_arguments)]
    pub fn record(
        mut self,
        qname: &str,
        flag: u16,
        rname: &str,
        pos: u32,
        cigar: &str,
        seq: &str,
        qual: &str,
    ) -> Self {
        let seq = if seq.is_empty() { "*" } else { seq };
        let qual = if qual.is_empty() { "*" } else { qual };
        self.records
            .push(format!("{qname}\t{flag}\t{rname}\t{pos}\t60\t{cigar}\t*\t0\t0\t{seq}\t{qual}"));
        self
    }

    pub fn write_to(&self, path: &Path) {
        let mut f = fs::File::create(path).expect("create sam");
        f.write_all(self.header.as_bytes()).unwrap();
        for r in &self.records {
            f.write_all(r.as_bytes()).unwrap();
            f.write_all(b"\n").unwrap();
        }
    }
}

/// Build a Phred-40 (`I`) quality string of the given length.
pub fn q40(len: usize) -> String {
    "I".repeat(len)
}

// ── Reference building ──────────────────────────────────────────────────────

/// Builds a tiny single-line-per-contig FASTA and a matching `.fai` index.
pub struct RefBuilder {
    contigs: Vec<(String, String)>,
}

impl RefBuilder {
    pub fn new() -> Self {
        Self { contigs: Vec::new() }
    }

    /// Add a contig with the given uppercase ACGT sequence.
    pub fn contig(mut self, name: &str, seq: &str) -> Self {
        self.contigs.push((name.to_string(), seq.to_string()));
        self
    }

    /// Write `<path>` and `<path>.fai`. Each contig is a single unwrapped line.
    pub fn write_to(&self, path: &Path) {
        let mut fa = String::new();
        let mut fai = String::new();
        let mut offset: u64 = 0;
        for (name, seq) in &self.contigs {
            let header_line_len = 1 + name.len() + 1; // '>' + name + '\n'
            let seq_offset = offset + header_line_len as u64;
            // name, length, offset, linebases, linewidth(+newline)
            fai.push_str(&format!(
                "{name}\t{}\t{seq_offset}\t{}\t{}\n",
                seq.len(),
                seq.len(),
                seq.len() + 1
            ));
            fa.push('>');
            fa.push_str(name);
            fa.push('\n');
            fa.push_str(seq);
            fa.push('\n');
            offset = seq_offset + seq.len() as u64 + 1; // + sequence + '\n'
        }
        fs::write(path, fa).expect("write fasta");
        let fai_path = format!("{}.fai", path.display());
        fs::write(&fai_path, fai).expect("write fai");
    }
}

// ── Running methylsieve ─────────────────────────────────────────────────────

/// A temp dir with conventional input/output/reference paths.
pub struct TestEnv {
    pub _tmp: TempDir,
    pub input: PathBuf,
    pub reference: PathBuf,
    pub output: PathBuf,
    /// `--metrics-prefix` value passed to the binary.
    pub metrics_prefix: PathBuf,
    /// The `PREFIX.summary.tsv` the metrics prefix produces (the per-context
    /// conversion summary; named `stats` for brevity in assertions).
    pub stats: PathBuf,
}

impl TestEnv {
    pub fn new() -> Self {
        let tmp = TempDir::new().expect("temp dir");
        let p = tmp.path();
        Self {
            input: p.join("in.sam"),
            reference: p.join("ref.fa"),
            output: p.join("out.bam"),
            metrics_prefix: p.join("metrics"),
            stats: p.join("metrics.summary.tsv"),
            _tmp: tmp,
        }
    }

    /// `--metrics-prefix` as a `&str` for arg lists.
    pub fn metrics_prefix_arg(&self) -> String {
        self.metrics_prefix.to_str().unwrap().to_string()
    }
}

/// Outcome of a methylsieve run.
pub struct RunResult {
    pub status_ok: bool,
    pub stderr: String,
    pub output: PathBuf,
}

/// Run methylsieve on `sam`/`reference`, writing BAM to `out`. Returns the run
/// result without asserting success (so error-path tests can inspect stderr).
pub fn run_methylsieve(
    sam: &SamBuilder,
    reference: &RefBuilder,
    env: &TestEnv,
    extra: &[&str],
) -> RunResult {
    sam.write_to(&env.input);
    reference.write_to(&env.reference);
    let out = Command::new(methylsieve_binary())
        .arg("-i")
        .arg(&env.input)
        .arg("-o")
        .arg(&env.output)
        .arg("-r")
        .arg(&env.reference)
        .args(extra)
        .output()
        .expect("methylsieve ran");
    RunResult {
        status_ok: out.status.success(),
        stderr: String::from_utf8_lossy(&out.stderr).into_owned(),
        output: env.output.clone(),
    }
}

/// Run methylsieve, asserting success, and return decoded output records.
pub fn run_ok(
    sam: &SamBuilder,
    reference: &RefBuilder,
    env: &TestEnv,
    extra: &[&str],
) -> Vec<RecordBuf> {
    let r = run_methylsieve(sam, reference, env, extra);
    assert!(r.status_ok, "methylsieve failed: {}", r.stderr);
    read_records(&r.output)
}

// ── Decoding output ─────────────────────────────────────────────────────────

/// Read a SAM/BAM file and return its header and records as owned `RecordBuf`s.
pub fn read_recs_and_header(path: &Path) -> (sam::Header, Vec<RecordBuf>) {
    let mut reader =
        AlignmentReaderBuilder::default().build_from_path(path).expect("open alignment file");
    let header = reader.read_header().expect("read header");
    let records = reader
        .records(&header)
        .map(|r| {
            let r = r.expect("read record");
            RecordBuf::try_from_alignment_record(&header, &*r).expect("record -> RecordBuf")
        })
        .collect();
    (header, records)
}

pub fn read_records(path: &Path) -> Vec<RecordBuf> {
    read_recs_and_header(path).1
}

/// `(qname, flag)` tuples in output order.
pub fn qname_flags(records: &[RecordBuf]) -> Vec<(String, u16)> {
    records
        .iter()
        .map(|r| (r.name().map(|n| n.to_string()).unwrap_or_default(), u16::from(r.flags())))
        .collect()
}

/// The `@PG` lines of a SAM/BAM, as SAM text.
pub fn read_pg_lines(path: &Path) -> Vec<String> {
    let (header, _records) = read_recs_and_header(path);
    let mut buf = Vec::new();
    sam::io::Writer::new(&mut buf).write_header(&header).expect("serialize header");
    String::from_utf8(buf)
        .expect("header is utf-8")
        .lines()
        .filter(|l| l.starts_with("@PG"))
        .map(String::from)
        .collect()
}

/// Value of a string-typed aux tag on a record (e.g. `XX`), or `None`.
pub fn tag_string(rec: &RecordBuf, tag: [u8; 2]) -> Option<String> {
    use noodles_sam::alignment::record_buf::data::field::Value;
    match rec.data().get(&Tag::new(tag[0], tag[1]))? {
        Value::String(s) => Some(s.to_string()),
        _ => None,
    }
}

/// Whether a record carries an aux tag with the given two bytes.
pub fn has_tag(rec: &RecordBuf, tag: [u8; 2]) -> bool {
    rec.data().get(&Tag::new(tag[0], tag[1])).is_some()
}

/// Phred quality scores of a record (raw values, not ASCII-offset).
pub fn quality_scores(rec: &RecordBuf) -> Vec<u8> {
    rec.quality_scores().as_ref().to_vec()
}

/// Count of leading bases at Phred quality `q` (e.g. the Q2 mask window).
pub fn leading_quality_run(rec: &RecordBuf, q: u8) -> usize {
    rec.quality_scores().as_ref().iter().take_while(|&&b| b == q).count()
}

// ── Stats TSV parsing ───────────────────────────────────────────────────────

/// Parse a methylsieve `--stats` TSV into one column→value map per data row
/// (the genome row first, then any control rows).
pub fn read_stats_rows(path: &Path) -> Vec<std::collections::HashMap<String, String>> {
    let text = fs::read_to_string(path).expect("read stats tsv");
    let mut lines = text.lines();
    let header: Vec<&str> = lines.next().expect("header row").split('\t').collect();
    lines
        .map(|line| {
            header
                .iter()
                .zip(line.split('\t'))
                .map(|(h, v)| ((*h).to_string(), v.to_string()))
                .collect()
        })
        .collect()
}

/// Convenience: the genome (first) row of a stats TSV.
pub fn genome_stats(path: &Path) -> std::collections::HashMap<String, String> {
    read_stats_rows(path).into_iter().next().expect("genome row")
}

/// Observed (total monitored) count for a context, from the `<ctx>_obs` column.
pub fn ctx_obs(row: &std::collections::HashMap<String, String>, ctx: &str) -> u64 {
    row[&format!("{ctx}_obs")].parse().unwrap()
}

/// Unconverted count for a context, derived from `<ctx>_obs` and
/// `<ctx>_conv_rate` (the summary reports the rate, not the raw count).
pub fn ctx_unconv(row: &std::collections::HashMap<String, String>, ctx: &str) -> u64 {
    let obs = ctx_obs(row, ctx);
    if obs == 0 {
        return 0;
    }
    let conv: f64 = row[&format!("{ctx}_conv_rate")].parse().unwrap();
    obs - (obs as f64 * conv).round() as u64
}

/// `ctx_obs` against the genome row of a stats TSV.
pub fn genome_ctx_obs(path: &Path, ctx: &str) -> u64 {
    ctx_obs(&genome_stats(path), ctx)
}

/// `ctx_unconv` against the genome row of a stats TSV.
pub fn genome_ctx_unconv(path: &Path, ctx: &str) -> u64 {
    ctx_unconv(&genome_stats(path), ctx)
}

/// SAM FLAG bit constants for fixtures and assertions.
pub const FLAG_PAIRED: u16 = 0x1;
pub const FLAG_PROPER_PAIR: u16 = 0x2;
pub const FLAG_UNMAPPED: u16 = 0x4;
pub const FLAG_MATE_UNMAPPED: u16 = 0x8;
pub const FLAG_REVERSE: u16 = 0x10;
pub const FLAG_MATE_REVERSE: u16 = 0x20;
pub const FLAG_FIRST_SEGMENT: u16 = 0x40;
pub const FLAG_LAST_SEGMENT: u16 = 0x80;
pub const FLAG_SECONDARY: u16 = 0x100;
pub const FLAG_QC_FAIL: u16 = 0x200;
pub const FLAG_SUPPLEMENTARY: u16 = 0x800;
