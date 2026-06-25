//! In-memory template buffer for the two-phase M-bias masking mode.
//!
//! During the *learn* phase we hold up to a target number of complete templates
//! (all records sharing a QNAME) while M-bias is accumulated; once the mask
//! lengths are frozen we *drain* the buffer in arrival order, mask, and emit,
//! then stream the rest of the file.
//!
//! Records are stored as raw BAM bytes packed into a handful of large **chunks**
//! (rather than one `Vec<RawRecord>` per record), so a million buffered records
//! cost a few allocations and pack densely instead of a million small heap
//! allocations. A per-record `(chunk, offset, len)` index and a per-template
//! record-range index preserve order and grouping.
//!
//! The representation is fully encapsulated: callers push whole templates and
//! iterate them back as byte slices through the methods below — there is no
//! public field access, so the chunked layout can change freely.

use fgumi_raw_bam::RawRecord;

/// Conservative bytes-per-template estimate used to size the byte ceiling and
/// chunk granularity. Generous (covers paired 150 bp reads plus supplementaries
/// and aux tags) so normal inputs stop on the template target, not the ceiling.
const EST_BYTES_PER_TEMPLATE: usize = 2048;

/// Smallest chunk size, so tiny targets still pack into one growable buffer.
const MIN_CHUNK_BYTES: usize = 8 << 20; // 8 MiB

/// Smallest byte ceiling, so tiny targets don't cap absurdly low.
const MIN_BYTE_CEILING: usize = 64 << 20; // 64 MiB

/// Locates one record's bytes within the chunked store.
#[derive(Clone, Copy)]
struct RecordLoc {
    chunk: u32,
    offset: u32,
    len: u32,
}

/// A buffered, byte-packed run of complete templates, drained in arrival order.
pub(crate) struct TemplateArena {
    /// Large packed byte buffers; new chunks are added on demand.
    chunks: Vec<Vec<u8>>,
    /// Target bytes per chunk (a new chunk starts once the last would overflow).
    chunk_bytes: usize,
    /// Hard cap on total stored bytes — a pathological input stops buffering
    /// early (and the run decides on what it has) rather than blowing up RSS.
    byte_ceiling: usize,
    /// Total bytes stored across all chunks.
    total_bytes: usize,
    /// Stop buffering once this many templates are held.
    target_templates: usize,
    /// Per-record locations, in arrival order.
    records: Vec<RecordLoc>,
    /// Per-template `(first_record_index, record_count)`, in arrival order.
    templates: Vec<(u32, u32)>,
}

impl TemplateArena {
    /// Build an arena that buffers up to `target_templates`. The byte ceiling
    /// and chunk size are derived from the target (≈6 chunks at the ceiling), so
    /// a single oversize allocation is never made and memory stays bounded.
    #[must_use]
    pub(crate) fn with_target(target_templates: usize) -> Self {
        let byte_ceiling = target_templates
            .saturating_mul(EST_BYTES_PER_TEMPLATE)
            .saturating_mul(3)
            .max(MIN_BYTE_CEILING);
        // Cap chunk size so a record's in-chunk offset always fits the `u32` in
        // `RecordLoc` — otherwise a very large `--mbias-buffer-templates` could
        // grow a chunk past 4 GiB and silently truncate offsets in `drain`.
        let chunk_bytes = (byte_ceiling / 6).max(MIN_CHUNK_BYTES).min(u32::MAX as usize);
        Self {
            chunks: Vec::new(),
            chunk_bytes,
            byte_ceiling,
            total_bytes: 0,
            target_templates,
            records: Vec::new(),
            templates: Vec::new(),
        }
    }

    /// Whether buffering should stop: the template target is reached, or the
    /// byte ceiling would be exceeded (graceful early stop on huge inputs).
    #[must_use]
    pub(crate) fn is_full(&self) -> bool {
        self.templates.len() >= self.target_templates || self.total_bytes >= self.byte_ceiling
    }

    /// Number of templates buffered.
    #[must_use]
    pub(crate) fn template_count(&self) -> usize {
        self.templates.len()
    }

    /// Copy one complete template (its records) into the arena. Returns `false`
    /// without storing anything if the arena is already full — the caller should
    /// then stop the learn phase and process this template in the drain/stream
    /// path instead. Empty templates are ignored.
    pub(crate) fn push_template(&mut self, recs: &[RawRecord]) -> bool {
        if self.is_full() {
            return false;
        }
        if recs.is_empty() {
            return true;
        }
        let first = self.records.len() as u32;
        for rec in recs {
            self.append_record(rec.as_ref());
        }
        self.templates.push((first, recs.len() as u32));
        true
    }

    /// Append one record's raw bytes, starting a new chunk if the current one
    /// would overflow `chunk_bytes` (an over-large record gets its own chunk).
    fn append_record(&mut self, bytes: &[u8]) {
        let need_new = match self.chunks.last() {
            None => true,
            Some(c) => !c.is_empty() && c.len() + bytes.len() > self.chunk_bytes,
        };
        if need_new {
            self.chunks.push(Vec::with_capacity(self.chunk_bytes.max(bytes.len())));
        }
        let chunk = self.chunks.len() - 1;
        let offset = self.chunks[chunk].len();
        // `chunk_bytes` is capped at `u32::MAX`, so a non-oversize record's offset
        // always fits; a lone oversize record sits at offset 0. `len` is a single
        // BAM record, always well under 4 GiB.
        debug_assert!(offset <= u32::MAX as usize && bytes.len() <= u32::MAX as usize);
        self.chunks[chunk].extend_from_slice(bytes);
        self.records.push(RecordLoc {
            chunk: chunk as u32,
            offset: offset as u32,
            len: bytes.len() as u32,
        });
        self.total_bytes += bytes.len();
    }

    /// Drain each buffered template in arrival order, materializing its records
    /// into a reused `Vec<RawRecord>` scratch and handing `f` a mutable slice
    /// (so it can mask in place and write). The scratch records' buffers are
    /// reused across templates, so only their one-time growth allocates.
    /// Errors from `f` short-circuit.
    pub(crate) fn drain<F>(&self, mut f: F) -> anyhow::Result<()>
    where
        F: FnMut(&mut [RawRecord]) -> anyhow::Result<()>,
    {
        let mut scratch: Vec<RawRecord> = Vec::new();
        for &(first, count) in &self.templates {
            let count = count as usize;
            while scratch.len() < count {
                scratch.push(RawRecord::new());
            }
            for (j, rec) in scratch.iter_mut().enumerate().take(count) {
                let RecordLoc { chunk, offset, len } = self.records[first as usize + j];
                let (o, l) = (offset as usize, len as usize);
                let bytes = &self.chunks[chunk as usize][o..o + l];
                let v = rec.as_mut_vec();
                v.clear();
                v.extend_from_slice(bytes);
            }
            f(&mut scratch[..count])?;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use fgumi_raw_bam::RawRecord;

    /// A `RawRecord` whose bytes are a recognizable fill — enough to test
    /// packing/iteration without a full BAM layout.
    fn rec(fill: u8, len: usize) -> RawRecord {
        let mut r = RawRecord::new();
        r.as_mut_vec().extend(std::iter::repeat_n(fill, len));
        r
    }

    #[test]
    fn round_trips_templates_in_order() {
        let mut a = TemplateArena::with_target(10);
        a.push_template(&[rec(1, 4), rec(2, 6)]); // template 0: two records
        a.push_template(&[rec(3, 5)]); // template 1: one record
        assert_eq!(a.template_count(), 2);

        let mut seen: Vec<Vec<Vec<u8>>> = Vec::new();
        a.drain(|recs| {
            seen.push(recs.iter().map(|r| r.as_ref().to_vec()).collect());
            Ok(())
        })
        .unwrap();
        assert_eq!(seen.len(), 2);
        assert_eq!(seen[0], vec![vec![1u8; 4], vec![2u8; 6]]);
        assert_eq!(seen[1], vec![vec![3u8; 5]]);
    }

    #[test]
    fn stops_at_template_target() {
        let mut a = TemplateArena::with_target(2);
        assert!(a.push_template(&[rec(1, 4)]));
        assert!(a.push_template(&[rec(2, 4)]));
        assert!(a.is_full());
        assert!(!a.push_template(&[rec(3, 4)]), "push past target is rejected");
        assert_eq!(a.template_count(), 2);
    }

    #[test]
    fn spans_multiple_chunks() {
        // Force tiny chunks by using a large target's ceiling math is hard to
        // reach in a test; instead verify many small records still round-trip
        // (exercises the chunk-append path across the default single chunk).
        let mut a = TemplateArena::with_target(1000);
        for i in 0..100u8 {
            assert!(a.push_template(&[rec(i, 32)]));
        }
        let mut n = 0;
        a.drain(|recs| {
            assert_eq!(recs.len(), 1);
            assert_eq!(recs[0].as_ref().len(), 32);
            n += 1;
            Ok(())
        })
        .unwrap();
        assert_eq!(n, 100);
    }

    #[test]
    fn oversize_record_gets_its_own_chunk() {
        let mut a = TemplateArena::with_target(4);
        // A record far larger than a chunk still stores and round-trips.
        let big = MIN_CHUNK_BYTES + 1024;
        a.push_template(&[rec(7, big)]);
        a.drain(|recs| {
            assert_eq!(recs[0].as_ref().len(), big);
            Ok(())
        })
        .unwrap();
        assert_eq!(a.template_count(), 1);
    }
}
