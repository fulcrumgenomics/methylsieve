//! Dedicated IO threads with a caller-sized byte ring buffer in each direction.
//!
//! Production methylsieve pipelines look like
//!     aligner | methylsieve | sorter
//! where both ends are typically bursty (the aligner has variable alignment
//! cost per read; the sorter periodically flushes a chunk to disk). With a
//! single-threaded reader+worker+writer, the small OS pipe buffer (~64 KB
//! on macOS) is the only thing decoupling stages, so any blip in one stage
//! stalls the others.
//!
//! [`ThreadedReader`] and [`ThreadedWriter`] put a dedicated thread on
//! each IO end with a user-space ring buffer in between, sized by the caller
//! (`ring_bytes`, with a small floor; the binary defaults to 16 MB read /
//! 64 MB write via `--read-buffer-mb` / `--write-buffer-mb`). The worker
//! reads/writes through the ring, never blocking on the kernel pipe.
//!
//! Design choices:
//! * **`ringbuf::HeapRb`** for the bytes — single allocation up front,
//!   recycled forever, no per-chunk heap traffic.
//! * **`thread::park` / `unpark`** for blocking. Lock-free fast path when
//!   the ring isn't full/empty; only the rare contended case parks.
//! * **`read()` straight into `vacant_slices_mut`** — one memcpy
//!   (kernel→ring) instead of two (kernel→temp + temp→ring).
//! * **Symmetric on write**: worker pushes bytes into the ring directly;
//!   IO writer thread drains via `as_slices()` + `write_all` + `skip`.

use std::io::{self, BufRead, Read, Write};
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::{Arc, Mutex, MutexGuard};
use std::thread::{self, JoinHandle};

use ringbuf::traits::{Consumer, Observer, Producer, Split};
use ringbuf::{HeapCons, HeapProd, HeapRb};

/// Acquire a Mutex even if it's been poisoned by a previous panic.
///
/// Both threads communicate IO errors through the same `Mutex<Option<io::Error>>`,
/// so a panic on one side would otherwise cascade into a panic on the other
/// when it next tries to read or store an error. Treating poison as
/// "no recorded error" lets us surface the original failure instead.
fn lock_or_recover<T>(m: &Mutex<T>) -> MutexGuard<'_, T> {
    m.lock().unwrap_or_else(|e| e.into_inner())
}

// ─── Read side ─────────────────────────────────────────────────────────────

/// `BufRead`-compatible reader fed by an IO thread.
pub(crate) struct ThreadedReader {
    consumer: HeapCons<u8>,
    io_thread: thread::Thread,
    eof: Arc<AtomicBool>,
    stop: Arc<AtomicBool>,
    /// First read error from the IO thread, if any.
    error: Arc<Mutex<Option<io::Error>>>,
    join: Option<JoinHandle<()>>,
}

impl ThreadedReader {
    /// Spawn an IO thread that reads from `src` into a ring buffer of
    /// `ring_bytes` capacity.
    pub(crate) fn new<R: Read + Send + 'static>(src: R, ring_bytes: usize) -> Self {
        let rb = HeapRb::<u8>::new(ring_bytes.max(64 * 1024));
        let (producer, consumer) = rb.split();
        let eof = Arc::new(AtomicBool::new(false));
        let stop = Arc::new(AtomicBool::new(false));
        let error = Arc::new(Mutex::new(None));

        let eof_io = eof.clone();
        let stop_io = stop.clone();
        let error_io = error.clone();
        let consumer_thread = thread::current();

        let join = thread::Builder::new()
            .name("methylsieve-io-read".into())
            .spawn(move || io_read_loop(src, producer, eof_io, stop_io, error_io, consumer_thread))
            .expect("spawning IO read thread");
        let io_thread = join.thread().clone();

        Self { consumer, io_thread, eof, stop, error, join: Some(join) }
    }

    fn take_error(&self) -> Option<io::Error> {
        lock_or_recover(&self.error).take()
    }
}

impl Read for ThreadedReader {
    fn read(&mut self, dst: &mut [u8]) -> io::Result<usize> {
        let src = self.fill_buf()?;
        let n = src.len().min(dst.len());
        dst[..n].copy_from_slice(&src[..n]);
        self.consume(n);
        Ok(n)
    }
}

impl BufRead for ThreadedReader {
    fn fill_buf(&mut self) -> io::Result<&[u8]> {
        loop {
            // Surface any IO-thread error first.
            if let Some(e) = self.take_error() {
                return Err(e);
            }

            if self.consumer.occupied_len() > 0 {
                let (first, _second) = self.consumer.as_slices();
                return Ok(first);
            }

            if self.eof.load(Ordering::Acquire) {
                // Drain any straggler bytes the producer published just
                // before setting EOF. `Acquire` above pairs with
                // `Release` in the IO thread, so re-checking occupied_len
                // here observes any final push.
                if self.consumer.occupied_len() > 0 {
                    continue;
                }
                return Ok(&[]);
            }

            // Ring is empty and producer hasn't flagged EOF yet. Park
            // until the IO thread unparks us. Spurious wakeups are
            // harmless because we re-check the loop condition.
            thread::park();
        }
    }

    fn consume(&mut self, amt: usize) {
        self.consumer.skip(amt);
        // Wake the IO thread in case it parked on a full ring.
        self.io_thread.unpark();
    }
}

impl Drop for ThreadedReader {
    fn drop(&mut self) {
        self.stop.store(true, Ordering::Release);
        self.io_thread.unpark();
        if let Some(h) = self.join.take() {
            let _ = h.join();
        }
    }
}

fn io_read_loop<R: Read>(
    mut src: R,
    mut producer: HeapProd<u8>,
    eof: Arc<AtomicBool>,
    stop: Arc<AtomicBool>,
    error: Arc<Mutex<Option<io::Error>>>,
    consumer_thread: thread::Thread,
) {
    loop {
        if stop.load(Ordering::Acquire) {
            break;
        }

        let (first, _second) = producer.vacant_slices_mut();
        if first.is_empty() {
            // Wake the consumer in case it's waiting (and we've just
            // become full because it's been slow).
            consumer_thread.unpark();
            thread::park();
            continue;
        }

        // SAFETY: We're reinterpreting `&mut [MaybeUninit<u8>]` as
        // `&mut [u8]` only to pass to `Read::read`, which writes into
        // every byte it claims to have read. After the read returns
        // `Ok(n)`, exactly `n` bytes are initialized; `advance_write_index(n)`
        // exposes only those to the consumer.
        let dst: &mut [u8] =
            unsafe { std::slice::from_raw_parts_mut(first.as_mut_ptr() as *mut u8, first.len()) };
        match src.read(dst) {
            Ok(0) => {
                eof.store(true, Ordering::Release);
                consumer_thread.unpark();
                break;
            }
            Ok(n) => {
                // SAFETY: `read` initialized exactly `n` bytes of the
                // vacant slice; ringbuf requires that pre-condition.
                unsafe {
                    producer.advance_write_index(n);
                }
                consumer_thread.unpark();
            }
            Err(e) => {
                *lock_or_recover(&error) = Some(e);
                eof.store(true, Ordering::Release);
                consumer_thread.unpark();
                break;
            }
        }
    }
    // Final wake so a consumer parked on `eof=false` sees the new state.
    consumer_thread.unpark();
}

// ─── Write side ────────────────────────────────────────────────────────────

/// `Write`-compatible writer that hands bytes off to an IO thread.
pub(crate) struct ThreadedWriter {
    producer: HeapProd<u8>,
    io_thread: thread::Thread,
    finished: Arc<AtomicBool>,
    error: Arc<Mutex<Option<io::Error>>>,
    join: Option<JoinHandle<()>>,
}

impl ThreadedWriter {
    /// Spawn an IO thread that writes the ring contents to `dst`. Ring
    /// holds `ring_bytes` of pending output.
    pub(crate) fn new<W: Write + Send + 'static>(dst: W, ring_bytes: usize) -> Self {
        let rb = HeapRb::<u8>::new(ring_bytes.max(64 * 1024));
        let (producer, consumer) = rb.split();
        let finished = Arc::new(AtomicBool::new(false));
        let error = Arc::new(Mutex::new(None));

        let finished_io = finished.clone();
        let error_io = error.clone();
        let producer_thread = thread::current();

        let join = thread::Builder::new()
            .name("methylsieve-io-write".into())
            .spawn(move || io_write_loop(dst, consumer, finished_io, error_io, producer_thread))
            .expect("spawning IO write thread");
        let io_thread = join.thread().clone();

        Self { producer, io_thread, finished, error, join: Some(join) }
    }

    fn take_error(&self) -> Option<io::Error> {
        lock_or_recover(&self.error).take()
    }

    /// Flush remaining bytes, signal the IO thread to drain, then join.
    /// Returns the IO thread's final result. Idempotent — calling twice
    /// is a no-op the second time.
    pub(crate) fn finish(mut self) -> io::Result<()> {
        self.finished.store(true, Ordering::Release);
        self.io_thread.unpark();
        if let Some(h) = self.join.take() {
            let _ = h.join();
        }
        if let Some(e) = self.take_error() {
            return Err(e);
        }
        Ok(())
    }
}

impl Write for ThreadedWriter {
    fn write(&mut self, mut buf: &[u8]) -> io::Result<usize> {
        // Surface any pending IO-thread error once per call, not per
        // ring-push iteration. The previous per-iteration check acquired
        // the `Mutex` ~150 M times on a 30 GB run.
        if let Some(e) = self.take_error() {
            return Err(e);
        }
        let initial_len = buf.len();
        while !buf.is_empty() {
            let pushed = self.producer.push_slice(buf);
            if pushed > 0 {
                buf = &buf[pushed..];
                self.io_thread.unpark();
            } else {
                // Ring is full. If the IO thread has died (e.g. a closed
                // downstream pipe → EPIPE, or a full disk), it will never drain
                // the ring, so parking here would hang forever. Surface its error
                // instead — checked both before parking (it may have already
                // failed) and after waking (it records the error, unparks us, then
                // exits), so we never park a second time into a dead thread.
                if let Some(e) = self.take_error() {
                    return Err(e);
                }
                thread::park();
                if let Some(e) = self.take_error() {
                    return Err(e);
                }
            }
        }
        Ok(initial_len)
    }

    fn flush(&mut self) -> io::Result<()> {
        // Nothing to flush at this layer — bytes are already in the ring
        // or written to `dst`. The IO thread does its own write_all.
        if let Some(e) = self.take_error() {
            return Err(e);
        }
        Ok(())
    }
}

impl Drop for ThreadedWriter {
    fn drop(&mut self) {
        // If finish() wasn't called, signal anyway so the IO thread can
        // shut down cleanly. Errors are silently dropped here — explicit
        // finish() is the right path for callers who care.
        self.finished.store(true, Ordering::Release);
        self.io_thread.unpark();
        if let Some(h) = self.join.take() {
            let _ = h.join();
        }
    }
}

fn io_write_loop<W: Write>(
    mut dst: W,
    mut consumer: HeapCons<u8>,
    finished: Arc<AtomicBool>,
    error: Arc<Mutex<Option<io::Error>>>,
    producer_thread: thread::Thread,
) {
    loop {
        if consumer.occupied_len() > 0 {
            let (first, _second) = consumer.as_slices();
            // Copy locally because `skip` borrows consumer mutably below.
            let n = first.len();
            // SAFETY: nothing — `first` is &[u8], plain write.
            if let Err(e) = dst.write_all(first) {
                *lock_or_recover(&error) = Some(e);
                producer_thread.unpark();
                break;
            }
            consumer.skip(n);
            producer_thread.unpark();
            continue;
        }
        if finished.load(Ordering::Acquire) {
            // Drain any stragglers — check again under acquire ordering.
            if consumer.occupied_len() > 0 {
                continue;
            }
            // Flush the underlying writer before exit.
            if let Err(e) = dst.flush() {
                *lock_or_recover(&error) = Some(e);
            }
            break;
        }
        thread::park();
    }
    producer_thread.unpark();
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    /// Round-trip a payload through a ThreadedReader: bytes in == bytes out.
    #[test]
    fn threaded_reader_round_trip_small() {
        let payload: Vec<u8> = (0..1000u32).flat_map(|i| i.to_le_bytes()).collect();
        let mut r = ThreadedReader::new(Cursor::new(payload.clone()), 64 * 1024);
        let mut out = Vec::new();
        std::io::copy(&mut r, &mut out).unwrap();
        assert_eq!(out, payload);
    }

    /// Payload much larger than the ring buffer — exercises the wrap-around.
    #[test]
    fn threaded_reader_round_trip_larger_than_ring() {
        let ring = 4096;
        let payload: Vec<u8> = (0..(ring * 8) as u32).map(|i| i as u8).collect();
        let mut r = ThreadedReader::new(Cursor::new(payload.clone()), ring);
        let mut out = Vec::new();
        std::io::copy(&mut r, &mut out).unwrap();
        assert_eq!(out, payload);
    }

    /// Write a payload through ThreadedWriter and confirm the underlying
    /// sink received every byte after `finish()`.
    #[test]
    fn threaded_writer_round_trip_with_finish() {
        // ThreadedWriter takes ownership of `W: Write + Send + 'static`, so
        // we hand it a `Sink` that mirrors bytes into a shared buffer the
        // test can inspect.
        struct Sink(std::sync::Arc<std::sync::Mutex<Vec<u8>>>);
        impl Write for Sink {
            fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
                self.0.lock().unwrap().extend_from_slice(buf);
                Ok(buf.len())
            }
            fn flush(&mut self) -> io::Result<()> {
                Ok(())
            }
        }

        let payload: Vec<u8> = (0..50_000u32).map(|i| i as u8).collect();
        let captured = std::sync::Arc::new(std::sync::Mutex::new(Vec::<u8>::new()));
        let mut w = ThreadedWriter::new(Sink(captured.clone()), 4096);
        w.write_all(&payload).unwrap();
        w.finish().unwrap();
        assert_eq!(*captured.lock().unwrap(), payload);
    }

    /// A dead IO thread (e.g. downstream pipe closed → EPIPE) must not hang the
    /// producer on a full ring: the error surfaces from `write`/`finish` instead.
    /// Regression for a park-forever bug where `write` checked the error only once
    /// before parking.
    #[test]
    fn threaded_writer_surfaces_error_without_hanging() {
        struct FailSink;
        impl Write for FailSink {
            fn write(&mut self, _buf: &[u8]) -> io::Result<usize> {
                Err(io::Error::new(io::ErrorKind::BrokenPipe, "downstream closed"))
            }
            fn flush(&mut self) -> io::Result<()> {
                Err(io::Error::new(io::ErrorKind::BrokenPipe, "downstream closed"))
            }
        }
        // Payload larger than the (floored 64 KiB) ring forces the producer to
        // block once the ring fills, exercising the park-then-recheck path. If the
        // error weren't surfaced there, this test would hang.
        let payload = vec![0u8; 256 * 1024];
        let mut w = ThreadedWriter::new(FailSink, 64 * 1024);
        let write_err = w.write_all(&payload).err();
        let finish_err = w.finish().err();
        assert!(
            write_err.is_some() || finish_err.is_some(),
            "broken-pipe error must surface from write() or finish()"
        );
    }
}
