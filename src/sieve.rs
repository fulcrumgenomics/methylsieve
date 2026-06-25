//! Per-template record processing — the heart of methylsieve.
//!
//! Each "block" is the run of consecutive records sharing a QNAME (methylsieve
//! requires query-grouped input). For each block we:
//!
//! 1. Classify the template by its primary R1's contig: **control** (a
//!    `--control-contig`) or **main** (the genome).
//! 2. Tally per-context converted/unconverted cytosines across the evaluated
//!    records (primary R1, primary R2, and — unless suppressed —
//!    supplementaries; secondaries never contribute). The reference base
//!    monitored is decided **per record** (`monitor_C = (R1 or unpaired) XOR
//!    reverse`), so reverse-mapped supplementaries flip correctly. For an
//!    overlapping proper pair, reference positions covered by both mates are
//!    counted once: the overlap is split at its midpoint and each mate keeps the
//!    half nearer its own 5' end (higher base quality) — see
//!    [`RecordProcessor::overlap_skip`]. With `--ignore-template-ends`, the
//!    outermost bases of each fragment terminus are skipped by genomic position
//!    in every record that covers them — see [`RecordProcessor::template_termini`].
//! 3. For main templates, make one unconverted/converted decision from the
//!    aggregated counts and propagate it — tag and/or QC-fail flag — to *every*
//!    record of the template (including secondaries/supplementaries), or drop
//!    them all with `--remove-unconverted`. Control templates are passed through
//!    untouched but still tallied into their own stats scope.

use std::collections::BTreeMap;

use anyhow::Result;
use fgumi_raw_bam::RawRecord;

use crate::mbias::{MbiasAccumulator, ReadEnd, ReadRole};
use crate::record::{
    FLAG_FIRST_SEGMENT, FLAG_LAST_SEGMENT, FLAG_PAIRED, FLAG_QC_FAIL, FLAG_REVERSE, FLAG_SECONDARY,
    FLAG_SUPPLEMENTARY, FLAG_UNMAPPED, five_prime_ref_span, has, is_primary_mapped, monitor_c_of,
    read_role, three_prime_ref_span,
};
use crate::reference::{BASE_A, BASE_C, BASE_G, BASE_T, Context, RefCodes, RefEncoding, Reference};
use crate::{ContextMask, TagSpec};

// ── Processor ───────────────────────────────────────────────────────────────

/// Per-block driver holding the reference and resolved options.
pub(crate) struct RecordProcessor {
    reference: Reference,
    opts: ProcessorOptions,
    /// Whether any `--control-contig` is configured. When false, the
    /// per-template chimeric-to-control scan is skipped entirely.
    has_controls: bool,
}

impl RecordProcessor {
    /// Build from a loaded reference and resolved options.
    #[must_use]
    pub(crate) fn new(reference: Reference, opts: ProcessorOptions) -> Self {
        let has_controls = opts.scope_of_tid.iter().any(Option::is_some);
        Self { reference, opts, has_controls }
    }

    /// Process one QNAME block: classify, tally, decide, propagate, emit.
    ///
    /// # Errors
    /// Returns an error if the block is not query-grouped (no primary present)
    /// or the output write fails.
    /// Classify, tally, accumulate M-bias, decide, and **stamp** every record of
    /// the template with the conversion tag / flag (no write). Returns whether
    /// the template should be emitted ([`Disposition::Keep`]) or dropped
    /// (`--remove-unconverted`). The caller owns emission, so it can buffer the
    /// stamped records (M-bias learn phase), mask them, or write them straight
    /// through.
    pub(crate) fn process_block(
        &self,
        block: &mut [RawRecord],
        stats: &mut Stats,
        mut mbias: Option<&mut MbiasAccumulator>,
    ) -> Result<Disposition> {
        stats.total_templates += 1;

        let class_idx = self.classification_index(block)?;
        let class_rec = &block[class_idx];

        // Unmapped primary R1 → pass through, no tally, no decision.
        if has(class_rec.flags(), FLAG_UNMAPPED) {
            stats.genome.n_templates += 1;
            stats.unmapped_templates += 1;
            self.stamp(block, Action::PassThrough, (0, 0));
            return Ok(Disposition::Keep);
        }

        let tid = class_rec.ref_id();
        let scope_idx = self.opts.scope_of_tid.get(tid as usize).copied().flatten();

        // Tally per-context counts for this template across its evidence
        // records. For overlapping proper pairs, the reference positions covered
        // by both mates are counted once — the overlap is split between the mates.
        // With `--ignore-template-ends`, the outermost bases of each fragment terminus
        // are skipped in every record that covers them (see `template_termini`).
        let overlap = self.overlap_skip(block);
        let trimming = self.opts.ignore_template_ends > 0;
        let termini =
            if trimming { self.template_termini(block) } else { TemplateTermini::default() };
        let mut counters = PerContextCounters::default();
        for (i, rec) in block.iter().enumerate() {
            if self.is_evidence_record(rec) {
                // Per-record trim: own 5' (or both ends for SE/orphan) via the read
                // window, plus one interior genomic skip. The overlap dedup region
                // takes precedence; an intruding mate terminus is used only when
                // there is no overlap (it is otherwise contained within it). The
                // terminus machinery is bypassed entirely when not trimming.
                let (trim_lo, trim_hi, mate_terminus) = if trimming {
                    let (lo, hi) = termini.own_trim_for(rec, self.opts.ignore_template_ends);
                    (lo, hi, termini.mate_skip_for(rec))
                } else {
                    (0, 0, None)
                };
                // Overlap handling: split the overlap at its midpoint and let each
                // mate keep the half nearer its own 5' end (where its base quality
                // is higher), so neither read's calls dominate the whole overlap.
                let overlap_iv = overlap.and_then(|(i1, i2, (os, oe))| {
                    if i != i1 && i != i2 {
                        return None;
                    }
                    if trimming {
                        // Keep legacy R1-takes-all so a 5' template-end trim still
                        // occupies the single skip slot (split and trimming aren't
                        // combined).
                        return (i == i2).then_some((os, oe));
                    }
                    let mid = os + (oe - os) / 2;
                    // A forward read's 5' end is at the low genomic coord, so it
                    // keeps [os, mid) and skips [mid, oe); a reverse read mirrors
                    // it. Proper FR pairs have opposite strands, so the two halves
                    // partition the overlap with no double count.
                    Some(if has(rec.flags(), FLAG_REVERSE) { (os, mid) } else { (mid, oe) })
                });
                let trim = RecordTrim { trim_lo, trim_hi, skip: overlap_iv.or(mate_terminus) };
                // When M-bias is being collected, a primary mapped record's
                // decision tally and M-bias accumulation share a single reference
                // scan: the M-bias sites (every called cytosine, no trim/dedup) are
                // a superset of the decision's, so one walk feeds both. Supplementary
                // evidence — and every record when M-bias is off — takes the plain
                // decision walk, leaving the masking-off hot path untouched.
                match mbias.as_deref_mut() {
                    Some(acc) if is_primary_mapped(rec.flags()) => {
                        self.tally_and_accumulate(rec, trim, &mut counters, acc);
                    }
                    _ => self.tally_record(rec, trim, &mut counters),
                }
            }
        }

        // Diagnostic (main templates only): did a *supplementary* land on a
        // control contig? Computed over all records independent of whether
        // supplementary evidence is suppressed, since it reflects mapping, not
        // tallying.
        let saw_control_supp = self.has_controls
            && scope_idx.is_none()
            && block.iter().any(|rec| {
                let f = rec.flags();
                has(f, FLAG_SUPPLEMENTARY) && !has(f, FLAG_UNMAPPED) && {
                    let rtid = rec.ref_id();
                    rtid >= 0
                        && self.opts.scope_of_tid.get(rtid as usize).copied().flatten().is_some()
                }
            });

        let monitored = counters.monitored_total();
        // Decision numerator/denominator over the threshold contexts; also the
        // `--count-tag` u/n, so compute once and reuse.
        let counts =
            (counters.unconv_in(self.opts.contexts), counters.total_in(self.opts.contexts));

        match scope_idx {
            Some(ci) => {
                // Control template: tally, never decide, never tag.
                let scope = &mut stats.controls[ci];
                scope.n_templates += 1;
                if monitored > 0 {
                    scope.n_evaluated += 1;
                }
                scope.counters.add(&counters);
                self.stamp(block, Action::PassThrough, counts);
                Ok(Disposition::Keep)
            }
            None => {
                // Main (genome) template.
                stats.genome.n_templates += 1;
                stats.genome.counters.add(&counters);
                if self.opts.record_matrix {
                    // Histogram cell keyed by (checked, unconverted) over the
                    // decision contexts; the decision/decided_by are replayed
                    // per cell at output time via `classify`.
                    *stats.conversion_matrix.entry((counts.1, counts.0)).or_insert(0) += 1;
                }
                if saw_control_supp {
                    stats.chimeric_to_control_templates += 1;
                }
                if monitored == 0 {
                    stats.zero_site_templates += 1;
                } else {
                    stats.genome.n_evaluated += 1;
                    // Flag the proportion-test blind spot: some threshold-context
                    // evidence, but below the site floor (uses the decision's own
                    // subset denominator `counts.1`, not the all-context `monitored`).
                    if counts.1 > 0 && counts.1 < u64::from(self.opts.min_sites) {
                        stats.below_min_sites_templates += 1;
                    }
                }

                let unconverted = self.decide(&counters);
                let action = if unconverted {
                    stats.genome.n_unconverted += 1;
                    if self.opts.remove_unconverted {
                        stats.genome.n_removed += 1;
                        Action::Remove
                    } else {
                        Action::Mark
                    }
                } else {
                    Action::PassThrough
                };
                if action == Action::Remove {
                    return Ok(Disposition::Drop);
                }
                self.stamp(block, action, counts);
                Ok(Disposition::Keep)
            }
        }
    }

    /// Index of the record whose contig classifies the template: primary R1 if
    /// present, else the sole unpaired primary, else any primary. Bails when no
    /// primary exists (the signature of non-query-grouped input).
    fn classification_index(&self, block: &[RawRecord]) -> Result<usize> {
        let mut r1_primary = None;
        let mut unpaired_primary = None;
        let mut any_primary = None;
        for (i, rec) in block.iter().enumerate() {
            let f = rec.flags();
            if has(f, FLAG_SECONDARY | FLAG_SUPPLEMENTARY) {
                continue;
            }
            any_primary.get_or_insert(i);
            if !has(f, FLAG_PAIRED) {
                unpaired_primary.get_or_insert(i);
            } else if has(f, FLAG_FIRST_SEGMENT) {
                r1_primary.get_or_insert(i);
            }
        }
        r1_primary.or(unpaired_primary).or(any_primary).ok_or_else(|| {
            let qname = block.first().map(|r| r.read_name().to_vec()).unwrap_or_default();
            anyhow::anyhow!(
                "QNAME {} appeared with {} record(s) but no primary alignment — the primary must \
                 be elsewhere in the stream. This almost always means the input is not \
                 query-grouped (e.g. coordinate-sorted). Re-sort with `samtools sort -n` and \
                 re-run.",
                String::from_utf8_lossy(&qname),
                block.len(),
            )
        })
    }

    /// Whether a record contributes conversion evidence: mapped, not secondary,
    /// and (unless suppressed) supplementaries are allowed.
    #[inline]
    fn is_evidence_record(&self, rec: &RawRecord) -> bool {
        let f = rec.flags();
        if has(f, FLAG_UNMAPPED) || has(f, FLAG_SECONDARY) {
            return false;
        }
        if has(f, FLAG_SUPPLEMENTARY) && self.opts.ignore_supplementary_evidence {
            return false;
        }
        true
    }

    /// For a proper paired-end template whose two primary mates overlap on the
    /// reference, return `(r1_index, r2_index, (start, end))`: the two primary
    /// records and the reference interval where they overlap. The caller splits
    /// that interval at its midpoint and assigns each half to the mate whose 5'
    /// end is nearer it, so each overlapped reference position is counted once.
    ///
    /// Returns `None` when the template isn't a qualifying pair or the mates
    /// don't overlap. Dedup is applied only when both mates monitor the **same**
    /// strand (`monitor_c(R1) == monitor_c(R2)`), which holds for proper FR
    /// pairs — the only case where an overlapped reference position is the same
    /// monitored base for both mates. In any other orientation the mates tally
    /// distinct positions, so skipping by interval would wrongly drop evidence;
    /// there we leave both mates intact.
    fn overlap_skip(&self, block: &[RawRecord]) -> Option<(usize, usize, (usize, usize))> {
        let mut r1 = None;
        let mut r2 = None;
        for (i, rec) in block.iter().enumerate() {
            let f = rec.flags();
            if has(f, FLAG_SECONDARY | FLAG_SUPPLEMENTARY | FLAG_UNMAPPED) || !has(f, FLAG_PAIRED) {
                continue;
            }
            if has(f, FLAG_FIRST_SEGMENT) {
                r1.get_or_insert(i);
            } else if has(f, FLAG_LAST_SEGMENT) {
                r2.get_or_insert(i);
            }
        }
        let (i1, i2) = (r1?, r2?);
        let (a, b) = (&block[i1], &block[i2]);
        if a.ref_id() < 0 || a.ref_id() != b.ref_id() {
            return None;
        }
        if monitor_c_of(a.flags()) != monitor_c_of(b.flags()) {
            return None;
        }
        // Reference spans [pos, pos + reference_length).
        let a1 = a.pos() as usize;
        let b1 = a1 + a.reference_length().max(0) as usize;
        let a2 = b.pos() as usize;
        let b2 = a2 + b.reference_length().max(0) as usize;
        let start = a1.max(a2);
        let end = b1.min(b2);
        if start < end { Some((i1, i2, (start, end))) } else { None }
    }

    /// Resolve the template's two fragment termini (in reference coordinates)
    /// for `--ignore-template-ends`. The end-repair fill-in and A-tailing
    /// artifacts sit at the physical ends of the original fragment, so we trim by
    /// genomic position, not read position:
    ///
    /// - **Mapped pair:** the two termini are the 5' sequenced ends of primary R1
    ///   and R2 — for a standard FR pair these are exactly the fragment's left and
    ///   right ends. Each is skipped in *both* mates wherever they overlap, so an
    ///   end seen by both reads (short insert / cfDNA) is trimmed in each.
    /// - **Single-end or orphan** (one mapped mate): the far end can't be located,
    ///   so both ends of the lone read are trimmed instead.
    ///
    /// Returns an empty (all-`None`) value when `--ignore-template-ends` is 0.
    /// Chimeric/supplementary segments are handled approximately: termini come
    /// from the primaries, and a supplementary contributes only via the genomic
    /// skips where it covers a primary-derived terminus.
    ///
    /// Kept out of line so it never bloats the hot per-block path: it is only
    /// reached when `--ignore-template-ends` is set.
    #[inline(never)]
    fn template_termini(&self, block: &[RawRecord]) -> TemplateTermini {
        let n = self.opts.ignore_template_ends as usize;
        if n == 0 {
            return TemplateTermini::default();
        }
        let mut r1: Option<&RawRecord> = None;
        let mut r2: Option<&RawRecord> = None;
        let mut unpaired: Option<&RawRecord> = None;
        for rec in block {
            let f = rec.flags();
            if has(f, FLAG_SECONDARY | FLAG_SUPPLEMENTARY | FLAG_UNMAPPED) {
                continue;
            }
            if !has(f, FLAG_PAIRED) {
                unpaired.get_or_insert(rec);
            } else if has(f, FLAG_FIRST_SEGMENT) {
                r1.get_or_insert(rec);
            } else if has(f, FLAG_LAST_SEGMENT) {
                r2.get_or_insert(rec);
            }
        }

        let tag =
            |rec: &RawRecord, span: Option<(usize, usize)>| span.map(|(s, e)| (rec.ref_id(), s, e));
        let five = |rec: &RawRecord| tag(rec, five_prime_ref_span(rec, n));
        let three = |rec: &RawRecord| tag(rec, three_prime_ref_span(rec, n));

        // Single-end primary: trim both ends of the one read.
        if let Some(rec) = unpaired {
            return TemplateTermini { r1: five(rec), r2: three(rec), single: true };
        }
        match (r1, r2) {
            // Mapped pair: a terminus from each mate's 5' sequenced end.
            (Some(a), Some(b)) => TemplateTermini { r1: five(a), r2: five(b), single: false },
            // Orphan (mate unmapped/absent): trim both ends of the lone read.
            (Some(m), None) | (None, Some(m)) => {
                TemplateTermini { r1: five(m), r2: three(m), single: true }
            }
            (None, None) => TemplateTermini::default(),
        }
    }

    /// Walk one record's aligned positions and add its monitored cytosines to
    /// `counters`, dispatching once on the reference encoding so the inner walk
    /// is monomorphized per [`RefCodes`] layout.
    fn tally_record(&self, rec: &RawRecord, trim: RecordTrim, counters: &mut PerContextCounters) {
        let tid = rec.ref_id();
        match self.reference.encoding() {
            RefEncoding::Bytes => {
                if let Some(c) = self.reference.byte_codes(tid) {
                    self.tally_aligned(rec, c, trim, counters);
                }
            }
            RefEncoding::Nibble => {
                if let Some(c) = self.reference.nibble_codes(tid) {
                    self.tally_aligned(rec, c, trim, counters);
                }
            }
            RefEncoding::TwoBit => {
                if let Some(c) = self.reference.twobit_codes(tid) {
                    self.tally_aligned(rec, c, trim, counters);
                }
            }
        }
    }

    /// Walk one record's aligned positions over a concrete reference encoding.
    ///
    /// The inner per-aligned-base work (ref-base check → context → BQ →
    /// read-base compare → counter bump) is the hot path, run by [`tally_span`].
    fn tally_aligned<R: RefCodes + Copy>(
        &self,
        rec: &RawRecord,
        refc: R,
        trim: RecordTrim,
        counters: &mut PerContextCounters,
    ) {
        let seq_len = rec.l_seq() as usize;
        if seq_len == 0 {
            return;
        }
        let ref_len = refc.len();

        // Per-record monitored strand (MethylDackel getStrand): treat single-end
        // and R1 the same; XOR with the record's own reverse bit.
        let f = rec.flags();
        let monitor_c = monitor_c_of(f);

        // Template-end trim window over *stored* read positions [keep_lo, keep_hi).
        // `trim.trim_lo`/`trim.trim_hi` are already resolved (by the caller) for
        // this record's strand and role: a paired read trims only its own 5' end,
        // a single-end/orphan read trims both ends. Stored soft-clips occupy the
        // read-position extremes, so they naturally count toward the budget. The
        // interior genomic skip (`trim.skip`) carries the mate's terminus or the
        // PE-overlap dedup region, which fall in this read's interior.
        let keep_lo = trim.trim_lo;
        let keep_hi = seq_len.saturating_sub(trim.trim_hi);

        let min_bq = self.opts.min_base_quality;
        // The monitored reference base is fixed for the whole record by strand.
        let monitored_base = if monitor_c { BASE_C } else { BASE_G };
        let mut read_pos: usize = 0;
        let pos = rec.pos();
        if pos < 0 {
            return;
        }
        let mut ref_pos: usize = pos as usize;

        let params = SpanParams {
            monitor_c,
            monitored_base,
            min_bq,
            keep_lo,
            keep_hi,
            seq_len,
            ref_len,
            skip: trim.skip,
        };
        for op in rec.cigar_ops_iter() {
            let len = (op >> 4) as usize;
            let code = op & 0xf;
            match code {
                // M, =, X — aligned; both read and reference advance.
                0 | 7 | 8 => {
                    tally_span(rec, refc, read_pos, ref_pos, len, &params, counters);
                    read_pos += len;
                    ref_pos += len;
                }
                // I — insertion: consumes read only.
                1 => read_pos += len,
                // S — soft clip: present in SEQ, consumes read only.
                4 => read_pos += len,
                // D, N — deletion / skip: consume reference only.
                2 | 3 => ref_pos += len,
                // H (5) hard clip, P (6) padding: consume neither stored read nor
                // reference.
                _ => {}
            }
        }
    }

    /// Tally one primary mapped record's decision counts **and** its per-cycle
    /// M-bias in a single reference scan, dispatching once on the encoding.
    ///
    /// Used in place of a separate [`Self::tally_aligned`] + M-bias pass whenever
    /// M-bias is being collected (`--metrics-prefix` or, pre-freeze, masking):
    /// the M-bias sites — every called cytosine at its true cycle, no end-trim and
    /// no overlap dedup — are a superset of the decision's, so one walk feeds both
    /// and `classify_site` runs once per site instead of twice. The masking-off
    /// decision path never reaches here; it keeps using [`tally_aligned`]
    /// untouched, so that codegen-sensitive hot loop is unaffected.
    fn tally_and_accumulate(
        &self,
        rec: &RawRecord,
        trim: RecordTrim,
        counters: &mut PerContextCounters,
        acc: &mut MbiasAccumulator,
    ) {
        let tid = rec.ref_id();
        match self.reference.encoding() {
            RefEncoding::Bytes => {
                if let Some(c) = self.reference.byte_codes(tid) {
                    self.fused_walk(rec, c, trim, counters, acc);
                }
            }
            RefEncoding::Nibble => {
                if let Some(c) = self.reference.nibble_codes(tid) {
                    self.fused_walk(rec, c, trim, counters, acc);
                }
            }
            RefEncoding::TwoBit => {
                if let Some(c) = self.reference.twobit_codes(tid) {
                    self.fused_walk(rec, c, trim, counters, acc);
                }
            }
        }
    }

    /// Single-scan decision tally + M-bias accumulation for one primary mapped
    /// record over a concrete reference encoding (monomorphized per [`RefCodes`]).
    ///
    /// The scan visits the full aligned range — every called cytosine — and
    /// records each into the M-bias accumulator at its 5' cycle (reverse records
    /// store SEQ forward-genomic, so their 5' end is the high stored position;
    /// single-end reads also record from the 3' end). It additionally bumps the
    /// decision counters for the subset of sites the decision walk would have
    /// counted: stored read position inside the trim window `[keep_lo, keep_hi)`
    /// and genomic position outside the interior skip — exactly the range
    /// [`tally_span`] carves out, evaluated per site here (cheap, since the site
    /// is already classified). `classify_site` and the per-base reference probe
    /// thus run once instead of once per walk.
    fn fused_walk<R: RefCodes + Copy>(
        &self,
        rec: &RawRecord,
        refc: R,
        trim: RecordTrim,
        counters: &mut PerContextCounters,
        acc: &mut MbiasAccumulator,
    ) {
        let seq_len = rec.l_seq() as usize;
        let pos = rec.pos();
        if seq_len == 0 || pos < 0 {
            return;
        }
        let ref_len = refc.len();
        let f = rec.flags();
        let monitor_c = monitor_c_of(f);
        let monitored_base = if monitor_c { BASE_C } else { BASE_G };
        let role = read_role(f);
        let is_se = role == ReadRole::Se;
        let reverse = has(f, FLAG_REVERSE);
        let min_bq = self.opts.min_base_quality;

        // Decision-counter gate, mirroring `tally_span` exactly: a site counts
        // toward the decision iff its stored read position is in `[keep_lo,
        // keep_hi)` and its genomic position is outside the interior skip. M-bias
        // records every site regardless. An absent skip is an empty interval.
        let keep_lo = trim.trim_lo;
        let keep_hi = seq_len.saturating_sub(trim.trim_hi);
        let (skip_s, skip_e) = trim.skip.unwrap_or((usize::MAX, usize::MAX));

        let mut read_pos = 0usize;
        let mut ref_pos = pos as usize;
        for op in rec.cigar_ops_iter() {
            let len = (op >> 4) as usize;
            match op & 0xf {
                // M, =, X — aligned; both read and reference advance.
                0 | 7 | 8 => {
                    // Clip the op to positions valid in both SEQ and the contig,
                    // so the inner loop needs no per-base bounds test.
                    let k1 = len
                        .min(seq_len.saturating_sub(read_pos))
                        .min(ref_len.saturating_sub(ref_pos));
                    for k in 0..k1 {
                        let gp = ref_pos + k;
                        if !refc.monitors(gp, monitored_base) {
                            continue;
                        }
                        let rp = read_pos + k;
                        if let Some((ctx, unconverted)) =
                            classify_site(rec, refc, rp, gp, monitor_c, min_bq, ref_len)
                        {
                            let cycle_5p = if reverse { seq_len - 1 - rp } else { rp };
                            acc.record(role, ReadEnd::FivePrime, ctx, cycle_5p, unconverted);
                            if is_se {
                                acc.record(
                                    role,
                                    ReadEnd::ThreePrime,
                                    ctx,
                                    seq_len - 1 - cycle_5p,
                                    unconverted,
                                );
                            }
                            if rp >= keep_lo && rp < keep_hi && !(gp >= skip_s && gp < skip_e) {
                                counters.record(ctx, unconverted);
                            }
                        }
                    }
                    read_pos += len;
                    ref_pos += len;
                }
                1 | 4 => read_pos += len,
                2 | 3 => ref_pos += len,
                _ => {}
            }
        }
    }

    /// Decide whether a template is unconverted from its aggregated counts.
    /// Thin wrapper over [`Self::classify`] on the threshold-context totals.
    fn decide(&self, counters: &PerContextCounters) -> bool {
        let unconv = counters.unconv_in(self.opts.contexts);
        let monitored = counters.total_in(self.opts.contexts);
        self.classify(unconv, monitored).0
    }

    /// Classify a template from its `(unconverted, monitored)` site counts over
    /// the threshold contexts: returns whether it is unconverted and which arm
    /// of the decision logic applied. Pure in `(unconv, monitored)` given the
    /// configured mode/thresholds, so `--conversion-matrix` can replay it per
    /// cell without re-tallying; `decide` is this `.0`.
    ///
    /// `min_sites` is a **floor**: the proportion test is unestimable below it
    /// and abstains. A template with no monitored sites is never unconverted.
    pub(crate) fn classify(&self, unconv: u64, monitored: u64) -> (bool, DecidedBy) {
        if monitored == 0 {
            return (false, DecidedBy::ZeroSites);
        }
        let count_hit = unconv >= u64::from(self.opts.max_unconverted_count);
        // The proportion is only estimable at/above the site floor; below it the
        // proportion test abstains.
        let at_floor = monitored >= u64::from(self.opts.min_sites);
        let frac_hit =
            at_floor && (unconv as f64) / (monitored as f64) > self.opts.max_unconverted_fraction;
        match self.opts.mode {
            DecisionMode::Count => (count_hit, DecidedBy::Count),
            DecisionMode::Proportion => {
                if at_floor {
                    (frac_hit, DecidedBy::Proportion)
                } else {
                    (false, DecidedBy::MinSitesFloor)
                }
            }
            DecisionMode::Either => (count_hit || frac_hit, DecidedBy::Either),
            // Trust the rate at/above the floor (an absolute count over-penalizes
            // long reads); fall back to the count below it, where the rate can't
            // be estimated.
            DecisionMode::Adaptive => {
                if at_floor {
                    (frac_hit, DecidedBy::Proportion)
                } else {
                    (count_hit, DecidedBy::Count)
                }
            }
        }
    }

    /// Stamp every record of the template with the conversion tag / QC-fail flag
    /// (for [`Action::Mark`]) and the optional per-record count tag. Mutates in
    /// place; does not write. `counts` is the template's `(unconverted, total)`
    /// over the decision contexts (for `--count-tag`). Not called for
    /// [`Action::Remove`] (the caller drops the template instead).
    fn stamp(&self, block: &mut [RawRecord], action: Action, counts: (u64, u64)) {
        // The count tag is a per-template aggregate (u/n over the decision
        // contexts): build the value once, stamp every record. Applied on every
        // template, flagged or not, so a user can inspect surprises either way.
        let count_value = self.opts.count_tag.map(|_| format!("{}/{}", counts.0, counts.1));
        for rec in block.iter_mut() {
            if action == Action::Mark {
                if self.opts.qc_fail {
                    rec.set_flags(rec.flags() | FLAG_QC_FAIL);
                }
                // Idempotent: don't append a second copy on a re-run.
                if rec.tags().find_string(&self.opts.tag.tag).is_none() {
                    rec.tags_editor().append_string(&self.opts.tag.tag, &self.opts.tag.value);
                }
            }
            if let (Some(tag), Some(value)) = (&self.opts.count_tag, &count_value)
                && rec.tags().find_string(tag).is_none()
            {
                rec.tags_editor().append_string(tag, value.as_bytes());
            }
        }
    }
}

/// Whether a processed template should be emitted or dropped.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Disposition {
    /// Emit every record (already stamped).
    Keep,
    /// Drop the whole template (`--remove-unconverted`).
    Drop,
}

/// How the per-template unconverted decision combines the count and proportion
/// tests. `min_sites` is a floor that disables the proportion test below it.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub(crate) enum DecisionMode {
    /// Count test only: flag when unconverted ≥ `max_unconverted_count`.
    Count,
    /// Proportion test only: flag when sites ≥ `min_sites` AND fraction >
    /// `max_unconverted_fraction`. Templates with fewer than `min_sites` sites
    /// are NEVER flagged (the proportion is unestimable) — they pass through.
    Proportion,
    /// Flag when EITHER the count or the proportion test fires.
    Either,
    /// Proportion test at/above `min_sites`, count test below it. The count
    /// fallback still evaluates low-site templates (unlike `Proportion`), while
    /// high-site templates are judged on rate rather than an absolute count that
    /// over-penalizes long reads / read pairs.
    Adaptive,
}

/// Which arm of [`RecordProcessor::classify`] produced a template's verdict —
/// surfaced per cell by `--conversion-matrix`.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub(crate) enum DecidedBy {
    /// No monitored sites in the threshold contexts → converted by default.
    ZeroSites,
    /// The count test was the operative arm (count mode, or the sub-floor
    /// fallback in adaptive mode).
    Count,
    /// The proportion test was the operative arm (at/above `min_sites`).
    Proportion,
    /// Proportion mode below the `min_sites` floor: the test abstains and the
    /// template passes through converted.
    MinSitesFloor,
    /// Either mode: the count and proportion tests OR'd together.
    Either,
}

impl DecidedBy {
    pub(crate) fn as_str(self) -> &'static str {
        match self {
            DecidedBy::ZeroSites => "zero_sites",
            DecidedBy::Count => "count",
            DecidedBy::Proportion => "proportion",
            DecidedBy::MinSitesFloor => "min_sites_floor",
            DecidedBy::Either => "either",
        }
    }
}

/// Runtime options the processor uses, resolved from [`crate::Args`].
pub(crate) struct ProcessorOptions {
    /// Contexts counted toward the unconverted threshold.
    pub(crate) contexts: ContextMask,
    /// How the count and proportion tests combine (`--mode`).
    pub(crate) mode: DecisionMode,
    /// Count threshold (`--max-unconverted-count`).
    pub(crate) max_unconverted_count: u32,
    /// Fraction threshold (`--max-unconverted-fraction`).
    pub(crate) max_unconverted_fraction: f64,
    /// Minimum monitored sites for the proportion test to apply (its floor, and
    /// the count↔proportion switch point in `adaptive`).
    pub(crate) min_sites: u32,
    /// Skip read bases below this base quality.
    pub(crate) min_base_quality: u8,
    /// Ignore the outermost N bases at each end of the template (fragment) when
    /// tallying — the end-repair / A-tailing–prone positions. See
    /// [`RecordProcessor::template_termini`].
    pub(crate) ignore_template_ends: u32,
    /// Exclude supplementaries from tallying (still tagged/flagged).
    pub(crate) ignore_supplementary_evidence: bool,
    /// The aux tag to set on unconverted templates.
    pub(crate) tag: TagSpec,
    /// If set, stamp every record with this aux tag carrying the template's
    /// `unconverted/total` site counts (over the decision contexts).
    pub(crate) count_tag: Option<[u8; 2]>,
    /// Whether to OR the QC-fail flag into unconverted records.
    pub(crate) qc_fail: bool,
    /// Drop unconverted templates from the output entirely.
    pub(crate) remove_unconverted: bool,
    /// Scope per BAM tid: `None` → genome, `Some(i)` → `controls[i]`.
    pub(crate) scope_of_tid: Vec<Option<usize>>,
    /// Accumulate the per-`(checked, unconverted)` decision histogram for
    /// `--conversion-matrix`. Off by default to avoid per-template map cost.
    pub(crate) record_matrix: bool,
}

/// What to do with a template's records on emit.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Action {
    /// Emit every record unchanged.
    PassThrough,
    /// Tag and/or QC-fail every record (unconverted, kept).
    Mark,
    /// Drop every record (unconverted, `--remove-unconverted`).
    Remove,
}

/// The two template (fragment) termini in reference coordinates, derived from
/// the 5' sequenced ends of the primary mates. Each is `(tid, start, end)` for
/// the half-open reference span of the outermost trimmed bases; `None` when that
/// end is absent or fully soft-clipped. Empty when `--ignore-template-ends` is 0.
#[derive(Debug, Clone, Copy, Default)]
struct TemplateTermini {
    /// 5' terminus of the primary first segment (R1), or the lone read.
    r1: Option<(i32, usize, usize)>,
    /// 5' terminus of the primary last segment (R2). For a single-end/orphan
    /// read this holds that read's *3'* end (its other, unknown-far terminus).
    r2: Option<(i32, usize, usize)>,
    /// True when only one mate is mapped (single-end or orphan): both ends of the
    /// single read are trimmed, since the far template end can't be located.
    single: bool,
}

impl TemplateTermini {
    /// Stored-read bases to trim from this record's (low, high) ends. A paired
    /// read trims only its own 5' end (low for forward, high for reverse); a
    /// single-end/orphan read trims both ends. `n` is `--ignore-template-ends`.
    #[inline(never)]
    fn own_trim_for(&self, rec: &RawRecord, n: u32) -> (usize, usize) {
        let n = n as usize;
        if n == 0 {
            return (0, 0);
        }
        if self.single {
            (n, n)
        } else if has(rec.flags(), FLAG_REVERSE) {
            (0, n)
        } else {
            (n, 0)
        }
    }

    /// The genomic terminus belonging to the *other* mate, if it intrudes into
    /// this record's reference span (so non-overlapping mates contribute no
    /// per-base skip). `None` for single-end/orphan reads (both ends already
    /// handled by [`own_trim_for`]).
    #[inline(never)]
    fn mate_skip_for(&self, rec: &RawRecord) -> Option<(usize, usize)> {
        if self.single {
            return None;
        }
        let f = rec.flags();
        let mate = if has(f, FLAG_FIRST_SEGMENT) {
            self.r2
        } else if has(f, FLAG_LAST_SEGMENT) {
            self.r1
        } else {
            return None;
        };
        let (ttid, s, e) = mate?;
        if ttid != rec.ref_id() {
            return None;
        }
        let rstart = rec.pos().max(0) as usize;
        let rend = rstart + rec.reference_length().max(0) as usize;
        (s < rend && rstart < e).then_some((s, e))
    }
}

/// Per-record trimming resolved at the template level and handed to the tally.
#[derive(Debug, Clone, Copy, Default)]
struct RecordTrim {
    /// Stored read bases to drop from the low-position end.
    trim_lo: usize,
    /// Stored read bases to drop from the high-position end.
    trim_hi: usize,
    /// Single genomic interval to skip in this record's interior: the PE-overlap
    /// dedup region, or the mate's intruding template terminus (the overlap
    /// subsumes the terminus when both apply).
    skip: Option<(usize, usize)>,
}

// ── Stats ───────────────────────────────────────────────────────────────────

/// Run-wide statistics, at the **template** level (one per QNAME block).
#[derive(Debug, Clone)]
pub(crate) struct Stats {
    /// The main genome scope (everything not on a control contig).
    pub(crate) genome: ScopeStats,
    /// One scope per `--control-contig`, in the order they were passed.
    pub(crate) controls: Vec<ScopeStats>,
    /// Main templates that had a supplementary on a control contig (diagnostic).
    pub(crate) chimeric_to_control_templates: u64,
    /// Templates whose primary R1 was unmapped (never tallied or decided).
    pub(crate) unmapped_templates: u64,
    /// Templates that produced zero monitored sites (decided converted).
    pub(crate) zero_site_templates: u64,
    /// Genome templates with some threshold-context evidence but fewer than
    /// `min_sites` sites — the proportion test can't evaluate them. In
    /// `proportion` mode these pass through unflagged (the blind spot); in
    /// `either`/`adaptive` the count test still covers them.
    pub(crate) below_min_sites_templates: u64,
    /// Total templates processed.
    pub(crate) total_templates: u64,
    /// Genome decision histogram keyed by `(checked_sites, unconverted_sites)`
    /// over the threshold contexts → template count. Populated only when
    /// `record_matrix` is set (drives the `--conversion-matrix` output).
    pub(crate) conversion_matrix: BTreeMap<(u64, u64), u64>,
}

impl Stats {
    /// Build empty stats with one genome scope plus a scope per control name.
    #[must_use]
    pub(crate) fn new(control_names: &[String]) -> Self {
        Self {
            genome: ScopeStats::new("genome".to_string()),
            controls: control_names.iter().map(|n| ScopeStats::new(n.clone())).collect(),
            chimeric_to_control_templates: 0,
            unmapped_templates: 0,
            zero_site_templates: 0,
            below_min_sites_templates: 0,
            total_templates: 0,
            conversion_matrix: BTreeMap::new(),
        }
    }
}

/// Stats for one reporting scope (the genome, or one control contig).
#[derive(Debug, Clone)]
pub(crate) struct ScopeStats {
    /// Scope name: `"genome"` or the control contig's name.
    pub(crate) name: String,
    /// Templates routed to this scope (mapped, unmapped, and zero-site alike).
    pub(crate) n_templates: u64,
    /// Templates that produced at least one monitored site (in any context).
    pub(crate) n_evaluated: u64,
    /// Templates decided unconverted (always 0 for control scopes).
    pub(crate) n_unconverted: u64,
    /// Unconverted templates dropped under `--remove-unconverted`.
    pub(crate) n_removed: u64,
    /// Aggregated per-context tallies across the scope.
    pub(crate) counters: PerContextCounters,
}

impl ScopeStats {
    fn new(name: String) -> Self {
        Self {
            name,
            n_templates: 0,
            n_evaluated: 0,
            n_unconverted: 0,
            n_removed: 0,
            counters: PerContextCounters::default(),
        }
    }
}

// ── Per-context counters ────────────────────────────────────────────────────

/// Converted/unconverted tallies per methylation context, indexed by
/// [`Context::index`] (CpA, CpC, CpG, CpT).
#[derive(Debug, Default, Clone, Copy)]
pub(crate) struct PerContextCounters {
    /// Unconverted cytosines (ref C read as C / ref G read as G) per context.
    unconv: [u64; 4],
    /// Total monitored & called cytosines (converted + unconverted) per context.
    total: [u64; 4],
}

impl PerContextCounters {
    #[inline]
    pub(crate) fn record(&mut self, ctx: Context, unconverted: bool) {
        let i = ctx.index();
        self.total[i] += 1;
        if unconverted {
            self.unconv[i] += 1;
        }
    }

    /// Add `unconv` unconverted of `total` monitored sites to `ctx` (bulk form of
    /// [`Self::record`], for aggregating from already-counted sources).
    pub(crate) fn add_counts(&mut self, ctx: Context, unconv: u64, total: u64) {
        self.unconv[ctx.index()] += unconv;
        self.total[ctx.index()] += total;
    }

    /// Unconverted count for one context.
    #[must_use]
    pub(crate) fn unconv_for(&self, ctx: Context) -> u64 {
        self.unconv[ctx.index()]
    }

    /// Total monitored count for one context.
    #[must_use]
    pub(crate) fn total_for(&self, ctx: Context) -> u64 {
        self.total[ctx.index()]
    }

    /// Accumulate `other` into `self`.
    pub(crate) fn add(&mut self, other: &PerContextCounters) {
        for i in 0..4 {
            self.unconv[i] += other.unconv[i];
            self.total[i] += other.total[i];
        }
    }

    /// Sum of unconverted counts over the contexts in `mask`.
    #[must_use]
    pub(crate) fn unconv_in(&self, mask: ContextMask) -> u64 {
        Context::ALL.iter().filter(|c| mask.contains(**c)).map(|c| self.unconv[c.index()]).sum()
    }

    /// Sum of total monitored counts over the contexts in `mask`.
    #[must_use]
    pub(crate) fn total_in(&self, mask: ContextMask) -> u64 {
        Context::ALL.iter().filter(|c| mask.contains(**c)).map(|c| self.total[c.index()]).sum()
    }

    /// Total monitored counts across all four contexts.
    #[must_use]
    pub(crate) fn monitored_total(&self) -> u64 {
        self.total.iter().sum()
    }
}

// ── Tally kernel ────────────────────────────────────────────────────────────
//
// The per-record hot path walks each aligned (M/=/X) span and, for every
// reference position equal to the strand's monitored base (C for top, G for
// bottom), classifies the read base as converted/unconverted in its reference
// context. The scan over the reference span ([`tally_span`]) is the dominant
// work; the per-site work after a match is shared via [`tally_site`].

/// Fixed-per-record parameters threaded into the span kernels.
#[derive(Debug, Clone, Copy)]
struct SpanParams {
    /// Whether this record monitors reference C (top) or G (bottom).
    pub(crate) monitor_c: bool,
    /// The monitored 4-bit reference base code (C=2 or G=4).
    pub(crate) monitored_base: u8,
    /// Minimum base quality to tally a site.
    pub(crate) min_bq: u8,
    /// Inclusive lower bound of the kept read-position window (the record's own
    /// 5'-end template trim folds into this; see [`RecordProcessor::tally_aligned`]).
    pub(crate) keep_lo: usize,
    /// Exclusive upper bound of the kept read-position window.
    pub(crate) keep_hi: usize,
    /// Stored SEQ length (read positions are valid in `0..seq_len`).
    pub(crate) seq_len: usize,
    /// Contig length (reference positions are valid in `0..ref_len`).
    pub(crate) ref_len: usize,
    /// A single reference half-open interval `[start, end)` to skip even when it
    /// falls inside the kept window: the PE-overlap dedup region, or the mate's
    /// template terminus that intrudes into this read. One interval suffices
    /// because an intruding mate terminus is always contained in the overlap
    /// region when both apply — see [`RecordProcessor::process_block`]. `None` in
    /// the common case (one cheap discriminant check on the per-base hot path).
    pub(crate) skip: Option<(usize, usize)>,
}

/// Classify one already-matched monitored cytosine at read position `rp` /
/// reference position `gp`: returns its `(context, unconverted)`, or `None` to
/// drop it (base quality below `min_bq`, a chromosome edge with no neighbor, or
/// a read base that is neither the unconverted nor converted call — SNP / N).
/// `refc[gp]` is assumed to equal the monitored base. Shared by the decision
/// tally ([`tally_site`]) and the fused decision+M-bias walk
/// ([`RecordProcessor::fused_walk`]) so the two agree on context/edge/drop rules;
/// `#[inline]` so it folds into the decision hot loop with no call overhead.
#[inline]
fn classify_site<R: RefCodes>(
    rec: &RawRecord,
    refc: R,
    rp: usize,
    gp: usize,
    monitor_c: bool,
    min_bq: u8,
    ref_len: usize,
) -> Option<(Context, bool)> {
    if rec.get_qual(rp) < min_bq {
        return None;
    }
    // Context from the reference neighbor (chrom-end safe), decoded in the
    // encoding's native space (no per-neighbor 4-bit decode for packed layouts).
    let ctx = if monitor_c {
        if gp + 1 >= ref_len {
            return None;
        }
        refc.ctx_top(gp + 1)
    } else {
        if gp == 0 {
            return None;
        }
        refc.ctx_bottom(gp - 1)
    }?;
    let unconverted = match (monitor_c, rec.get_base(rp)) {
        (true, BASE_C) | (false, BASE_G) => true,  // unconverted
        (true, BASE_T) | (false, BASE_A) => false, // converted
        _ => return None,                          // SNP / N — drop
    };
    Some((ctx, unconverted))
}

/// Classify a single already-matched monitored cytosine and record it into the
/// decision counters. `refc[gp]` is assumed to equal the monitored base.
#[inline]
fn tally_site<R: RefCodes>(
    rec: &RawRecord,
    refc: R,
    rp: usize,
    gp: usize,
    p: &SpanParams,
    counters: &mut PerContextCounters,
) {
    if let Some((ctx, unconverted)) =
        classify_site(rec, refc, rp, gp, p.monitor_c, p.min_bq, p.ref_len)
    {
        counters.record(ctx, unconverted);
    }
}

/// Tally a contiguous run of aligned positions `[lo, hi)` (in local `k`
/// coordinates) with no per-base skip test. Pulled out so the hot caller stays
/// small; `#[inline]` lets it fuse into [`tally_span`]'s common path. Args are
/// threaded as scalars (not a `Range`) so they stay in registers on the hot path.
#[inline]
#[allow(clippy::too_many_arguments)]
fn tally_run<R: RefCodes + Copy>(
    rec: &RawRecord,
    refc: R,
    rp_start: usize,
    gp_start: usize,
    lo: usize,
    hi: usize,
    p: &SpanParams,
    counters: &mut PerContextCounters,
) {
    for k in lo..hi {
        let gp = gp_start + k;
        // Reference check first: rejects ~79% of bases (only ~21% are the
        // monitored C/G) before the gather work. `monitors` compares in the
        // encoding's native space (no per-base decode for packed layouts).
        if !refc.monitors(gp, p.monitored_base) {
            continue;
        }
        tally_site(rec, refc, rp_start + k, gp, p, counters);
    }
}

/// Scalar reference-span tally generic over the reference encoding: walk `len`
/// aligned positions from read position `rp_start` / reference position
/// `gp_start`. Monomorphized per [`RefCodes`] impl, so each encoding gets its
/// own branch-free-on-layout loop.
///
/// Any genomic skip (mate terminus / PE-overlap dedup) is carved out of the scan
/// *range* rather than tested per base: the common no-skip case is a single tight
/// loop with zero per-base skip overhead — strictly less work than a per-base
/// skip test — and the rare skip case scans the two sub-ranges around the
/// excluded middle. (Keeping the skip handling inline matters: any out-of-line
/// helper here stops `tally_span` from fusing into the caller and regresses the
/// hot path.)
fn tally_span<R: RefCodes + Copy>(
    rec: &RawRecord,
    refc: R,
    rp_start: usize,
    gp_start: usize,
    len: usize,
    p: &SpanParams,
    counters: &mut PerContextCounters,
) {
    // Hoist the read-position window and contig/SEQ bounds out of the per-base
    // loop: across an M-span `rp = rp_start + k` and `gp = gp_start + k` move
    // together, so the four per-base bounds checks collapse to a single `k` range.
    let k0 = p.keep_lo.saturating_sub(rp_start);
    let k1 = len
        .min(p.keep_hi.saturating_sub(rp_start))
        .min(p.seq_len.saturating_sub(rp_start))
        .min(p.ref_len.saturating_sub(gp_start));
    match p.skip {
        None => tally_run(rec, refc, rp_start, gp_start, k0, k1, p, counters),
        Some((s, e)) => {
            // Skip applies where `gp ∈ [s, e)`, i.e. `k ∈ [s-gp_start, e-gp_start)`,
            // clamped into `[k0, k1)`. A skip outside the span leaves one side empty.
            let sk0 = s.saturating_sub(gp_start).clamp(k0, k1);
            let sk1 = e.saturating_sub(gp_start).clamp(k0, k1);
            tally_run(rec, refc, rp_start, gp_start, k0, sk0, p, counters);
            tally_run(rec, refc, rp_start, gp_start, sk1, k1, p, counters);
        }
    }
}

#[cfg(test)]
mod tests {
    use std::io::{BufRead, Cursor};

    use super::*;
    use crate::reference::{Reference, encode_ref_base};
    use crate::sam_reader::SamReader;
    use crate::{TagSpec, parse_contexts};

    /// Encode an ASCII reference string into 4-bit codes for a tiny contig.
    fn enc(seq: &str) -> Vec<u8> {
        seq.bytes().map(encode_ref_base).collect()
    }

    fn cph_mask() -> ContextMask {
        parse_contexts("CpA,CpC,CpT").unwrap()
    }

    fn opts(contexts: ContextMask, count: u32) -> ProcessorOptions {
        ProcessorOptions {
            contexts,
            mode: DecisionMode::Either,
            max_unconverted_count: count,
            max_unconverted_fraction: 1.0,
            min_sites: 5,
            min_base_quality: 0,
            ignore_template_ends: 0,
            ignore_supplementary_evidence: false,
            tag: TagSpec { tag: [b'X', b'X'], value: b"UC".to_vec() },
            count_tag: None,
            qc_fail: true,
            remove_unconverted: false,
            scope_of_tid: vec![None],
            record_matrix: false,
        }
    }

    /// CpA counters with `unconv` unconverted out of `total` monitored sites.
    fn cph_counters(unconv: u64, total: u64) -> PerContextCounters {
        let mut c = PerContextCounters::default();
        for i in 0..total {
            c.record(Context::CpA, i < unconv);
        }
        c
    }

    /// A processor whose `decide` uses the given mode/thresholds (the reference
    /// is irrelevant to the decision, so a 1-base stub suffices).
    fn decider(mode: DecisionMode, count: u32, fraction: f64, min_sites: u32) -> RecordProcessor {
        let mut o = opts(cph_mask(), count);
        o.mode = mode;
        o.max_unconverted_fraction = fraction;
        o.min_sites = min_sites;
        RecordProcessor::new(Reference::from_encoded_contigs(vec![enc("C")]), o)
    }

    #[test]
    fn adaptive_below_floor_uses_count() {
        // 3 unconverted of 4 sites, floor 40 → proportion abstains, count (≥3)
        // decides → flagged.
        let d = decider(DecisionMode::Adaptive, 3, 0.05, 40);
        assert!(d.decide(&cph_counters(3, 4)));
        assert!(!d.decide(&cph_counters(2, 4)), "2 < count threshold 3 → not flagged");
    }

    #[test]
    fn adaptive_above_floor_uses_proportion_not_count() {
        // 4 unconverted of 100 sites: count (≥3) would fire, but at/above the
        // floor adaptive trusts the rate (4% < 5%) → NOT flagged. This is the
        // leniency on long reads that distinguishes adaptive from either/count.
        let d = decider(DecisionMode::Adaptive, 3, 0.05, 40);
        assert!(!d.decide(&cph_counters(4, 100)));
        assert!(d.decide(&cph_counters(10, 100)), "10% > 5% → flagged");
    }

    #[test]
    fn proportion_below_floor_never_flags() {
        // The blind spot: 4 of 4 unconverted (100%) but below the 40-site floor →
        // proportion can't evaluate → passes through.
        let d = decider(DecisionMode::Proportion, 3, 0.05, 40);
        assert!(!d.decide(&cph_counters(4, 4)));
        assert!(d.decide(&cph_counters(3, 40)), "3/40 = 7.5% > 5% at the floor → flagged");
    }

    #[test]
    fn count_mode_ignores_fraction_and_floor() {
        // 4 of 100 (4% < 5%): count mode flags purely on the absolute count.
        let d = decider(DecisionMode::Count, 3, 0.05, 40);
        assert!(d.decide(&cph_counters(4, 100)));
        // And it fires below the floor too, where proportion would abstain.
        assert!(d.decide(&cph_counters(3, 4)));
    }

    #[test]
    fn either_fires_when_only_one_test_hits() {
        // count threshold 10 so the count test misses; proportion (8/100 = 8% >
        // 5%) carries it → either flags, but count-only would not.
        let either = decider(DecisionMode::Either, 10, 0.05, 40);
        assert!(either.decide(&cph_counters(8, 100)));
        assert!(!decider(DecisionMode::Count, 10, 0.05, 40).decide(&cph_counters(8, 100)));
        // Below the floor, only the count arm can fire; it does at 10/12.
        assert!(either.decide(&cph_counters(10, 12)));
    }

    #[test]
    fn adaptive_switch_is_continuous_at_the_default_floor() {
        // Documents why min_sites=40 pairs with count=3, fraction=0.05: at the
        // switch point the two tests agree at 3 unconverted, so crossing the
        // floor doesn't flip the call. Just below (39 sites) count decides; just
        // above (40 sites) proportion does — both flag at 3, neither at 2.
        let d = decider(DecisionMode::Adaptive, 3, 0.05, 40);
        assert!(d.decide(&cph_counters(3, 39)) && d.decide(&cph_counters(3, 40)));
        assert!(!d.decide(&cph_counters(2, 39)) && !d.decide(&cph_counters(2, 40)));
    }

    #[test]
    fn classify_reports_decided_by_and_matches_decide() {
        // Zero monitored sites → converted, regardless of mode.
        let d = decider(DecisionMode::Adaptive, 3, 0.05, 40);
        assert_eq!(d.classify(0, 0), (false, DecidedBy::ZeroSites));

        // Count mode always reports the count arm.
        let c = decider(DecisionMode::Count, 3, 0.05, 40);
        assert_eq!(c.classify(3, 10), (true, DecidedBy::Count));
        assert_eq!(c.classify(2, 10), (false, DecidedBy::Count));

        // Proportion mode abstains below the floor, applies the rate at/above it.
        let p = decider(DecisionMode::Proportion, 3, 0.05, 40);
        assert_eq!(p.classify(5, 10), (false, DecidedBy::MinSitesFloor));
        assert_eq!(p.classify(5, 50), (true, DecidedBy::Proportion));

        // Adaptive: count arm below the floor, proportion arm at/above it.
        assert_eq!(d.classify(3, 10), (true, DecidedBy::Count));
        assert_eq!(d.classify(1, 50), (false, DecidedBy::Proportion));

        // Either OR's the two tests under one label.
        let e = decider(DecisionMode::Either, 3, 0.05, 40);
        assert_eq!(e.classify(3, 10), (true, DecidedBy::Either));
        assert_eq!(e.classify(2, 10), (false, DecidedBy::Either));

        // classify(.0) can never drift from decide() over the same counts.
        for &(u, t) in &[(0, 0), (3, 10), (2, 10), (5, 50), (1, 50), (40, 100)] {
            assert_eq!(d.classify(u, t).0, d.decide(&cph_counters(u, t)), "u={u} t={t}");
        }
    }

    /// Parse a one-contig SAM document into `RawRecord`s via methylsieve's own
    /// SAM reader (so the records carry a correct BAM layout: packed 4-bit SEQ,
    /// quals, CIGAR, flags, tid). `contig_len` sizes the single `@SQ` line.
    fn parse_sam_records(records: &[&str], contig_len: usize) -> Vec<RawRecord> {
        let mut sam = format!("@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:chr1\tLN:{contig_len}\n");
        for r in records {
            sam.push_str(r);
            sam.push('\n');
        }
        let boxed: Box<dyn BufRead> = Box::new(Cursor::new(sam.into_bytes()));
        let mut reader = SamReader::new(boxed);
        reader.read_header().expect("read header");
        let mut out = Vec::new();
        loop {
            let mut rec = RawRecord::new();
            if !reader.read_record(&mut rec).expect("read record") {
                break;
            }
            out.push(rec);
        }
        out
    }

    /// A single SAM record line with the standard 11 fields, no aux.
    fn sam_line(flag: u16, pos: u32, cigar: &str, seq: &str, qual: &str) -> String {
        format!("q\t{flag}\tchr1\t{pos}\t60\t{cigar}\t*\t0\t0\t{seq}\t{qual}")
    }

    #[test]
    fn forward_read_top_strand_counts_unconverted_cph() {
        // Reference C A C A ... ; monitored ref C (top). ref[i+1]=A → CpA.
        // Read identical to ref → all C's unconverted.
        let reference = Reference::from_encoded_contigs(vec![enc("CACACACACA")]);
        let proc = RecordProcessor::new(reference, opts(cph_mask(), 3));
        let recs = parse_sam_records(&[&sam_line(0, 1, "10M", "CACACACACA", "IIIIIIIIII")], 10);
        let mut counters = PerContextCounters::default();
        proc.tally_record(&recs[0], RecordTrim::default(), &mut counters);
        assert_eq!(counters.unconv[Context::CpA.index()], 5);
        assert_eq!(counters.total[Context::CpA.index()], 5);
        assert!(proc.decide(&counters), "5 unconverted CpA ≥ 3 → unconverted");
    }

    #[test]
    fn converted_read_is_not_flagged() {
        let reference = Reference::from_encoded_contigs(vec![enc("CACACACACA")]);
        let proc = RecordProcessor::new(reference, opts(cph_mask(), 3));
        let recs = parse_sam_records(&[&sam_line(0, 1, "10M", "TATATATATA", "IIIIIIIIII")], 10);
        let mut counters = PerContextCounters::default();
        proc.tally_record(&recs[0], RecordTrim::default(), &mut counters);
        assert_eq!(counters.unconv[Context::CpA.index()], 0);
        assert_eq!(counters.total[Context::CpA.index()], 5);
        assert!(!proc.decide(&counters));
    }

    #[test]
    fn reverse_single_end_read_monitors_ref_g() {
        // Reference G T G T ...; a reverse SE read (flag 0x10) monitors ref G
        // (monitor_C = treat_as_read1(true) XOR reverse(true) = false → G).
        // ref[i-1]=T → bottom context CpA. Read G at those positions → unconv.
        let reference = Reference::from_encoded_contigs(vec![enc("TGTGTGTGTG")]);
        let proc = RecordProcessor::new(reference, opts(cph_mask(), 3));
        // SEQ as stored (already reverse-complemented for a reverse alignment):
        // matches ref so the monitored G's read as G (unconverted).
        let recs =
            parse_sam_records(&[&sam_line(FLAG_REVERSE, 1, "10M", "TGTGTGTGTG", "IIIIIIIIII")], 10);
        let mut counters = PerContextCounters::default();
        proc.tally_record(&recs[0], RecordTrim::default(), &mut counters);
        // ref G at positions 1,3,5,7,9 (0-based); each has ref[i-1]=T → CpA.
        assert_eq!(counters.unconv[Context::CpA.index()], 5);
        assert_eq!(counters.total[Context::CpA.index()], 5);
    }
}
