/*************************************************************************************
                           The MIT License

   BWA-MEM2  (Sequence alignment using Burrows-Wheeler Transform),
   Copyright (C) 2019  Intel Corporation, Heng Li.

   Rust port of kseq.h — streaming FASTA/FASTQ reader.

   Original kseq.h Copyright (c) 2008, 2009, 2011 Attractive Chaos
   <attractor@live.co.uk>; The MIT License.

Authors: Rust port based on original C by Heng Li <hli@jimmy.harvard.edu>
*****************************************************************************************/

use std::fs::File;
use std::io::{self, BufRead, BufReader, Read};

use anyhow::{bail, Result};
use flate2::read::GzDecoder;

// ---------------------------------------------------------------------------
// Public record type
// ---------------------------------------------------------------------------

/// A single FASTA or FASTQ record.
#[derive(Debug, Clone, Default)]
pub struct KSeq {
    /// Sequence name — everything between the `>` / `@` sigil and the first
    /// whitespace on the header line.
    pub name: String,
    /// Optional comment: the part of the header line after the name, with the
    /// leading whitespace character stripped.
    pub comment: String,
    /// Nucleotide sequence bytes, whitespace stripped and concatenated across
    /// all FASTA continuation lines.
    pub seq: Vec<u8>,
    /// Base-quality scores (Phred+33 encoded); empty for FASTA records.
    pub qual: Vec<u8>,
}

// ---------------------------------------------------------------------------
// Streaming parser
// ---------------------------------------------------------------------------

/// Buffered FASTA / FASTQ parser.  Yields one [`KSeq`] record per [`next`] call.
///
/// `R` must implement [`BufRead`].  Obtain a suitable reader for plain or
/// gzip-compressed files with [`open_seq_file`].
pub struct KSeqReader<R: BufRead> {
    inner: R,
    /// The most recently read `>` or `@` byte that begins the *current* or
    /// *next* record's header.  Set to `0` before the first record has been
    /// found or after successful FASTQ quality parsing, when the next header
    /// char has not yet been read.
    last_char: u8,
    /// Set once the underlying stream reports EOF.
    eof: bool,
}

impl<R: BufRead> KSeqReader<R> {
    /// Wrap an existing buffered reader.
    pub fn new(reader: R) -> Self {
        KSeqReader {
            inner: reader,
            last_char: 0,
            eof: false,
        }
    }

    // -----------------------------------------------------------------------
    // Low-level byte helpers
    // -----------------------------------------------------------------------

    /// Consume and return the next byte, or `None` at EOF.
    #[inline]
    fn read_byte(&mut self) -> io::Result<Option<u8>> {
        let buf = self.inner.fill_buf()?;
        if buf.is_empty() {
            self.eof = true;
            return Ok(None);
        }
        let b = buf[0];
        self.inner.consume(1);
        Ok(Some(b))
    }

    /// Peek at the next byte without consuming it, or `None` at EOF.
    #[inline]
    fn peek_byte(&mut self) -> io::Result<Option<u8>> {
        let buf = self.inner.fill_buf()?;
        if buf.is_empty() {
            self.eof = true;
            Ok(None)
        } else {
            Ok(Some(buf[0]))
        }
    }

    /// Read bytes into `dst` until `delimiter` (the delimiter itself is consumed
    /// but not appended) or until EOF.
    ///
    /// Returns `true` when the delimiter was found; `false` on EOF with no data.
    fn read_until_byte(&mut self, delimiter: u8, dst: &mut Vec<u8>) -> io::Result<bool> {
        let mut found_any = false;
        loop {
            let available = self.inner.fill_buf()?;
            if available.is_empty() {
                self.eof = true;
                return Ok(found_any);
            }
            match available.iter().position(|&b| b == delimiter) {
                Some(pos) => {
                    dst.extend_from_slice(&available[..pos]);
                    self.inner.consume(pos + 1); // consume through the delimiter
                    return Ok(true);
                }
                None => {
                    let len = available.len();
                    dst.extend_from_slice(available);
                    self.inner.consume(len);
                    found_any = true;
                }
            }
        }
    }

    /// Read one line (up to and including `\n`) into `dst`.
    ///
    /// Trailing `\r` (Windows CRLF) is stripped.  Returns `true` if at least
    /// one byte was available before EOF; `false` on immediate EOF.
    fn read_line_into(&mut self, dst: &mut Vec<u8>) -> io::Result<bool> {
        let found = self.read_until_byte(b'\n', dst)?;
        // Strip Windows-style \r if present.
        if dst.last() == Some(&b'\r') {
            dst.pop();
        }
        Ok(found)
    }

    // -----------------------------------------------------------------------
    // Record parser
    // -----------------------------------------------------------------------

    /// Return the next FASTA or FASTQ record, or `None` at EOF.
    ///
    /// # Errors
    /// Propagates I/O errors and returns an error for malformed FASTQ input
    /// (quality line length != sequence length).
    pub fn next(&mut self) -> Result<Option<KSeq>> {
        if self.eof {
            return Ok(None);
        }

        // ---- Find the opening `>` / `@` of the next record ----------------
        //
        // After a FASTA record we leave the `>` of the next record *in the
        // underlying buffer* (we peeked but did not consume it), so the scan
        // below will find it immediately.  After FASTQ quality parsing we also
        // check for `@` before returning, restoring it via `last_char` only
        // when we had to consume it.
        let header_char: u8 = if self.last_char != 0 {
            let hc = self.last_char;
            self.last_char = 0;
            hc
        } else {
            // Scan forward until we find `>` or `@`.
            loop {
                match self.read_byte()? {
                    None => return Ok(None),
                    Some(b'>') => break b'>',
                    Some(b'@') => break b'@',
                    Some(_) => continue,
                }
            }
        };

        // ---- Header line: name + optional comment -------------------------
        // Read the rest of the header line (everything after `>` / `@`).
        let mut header_line: Vec<u8> = Vec::with_capacity(256);
        self.read_line_into(&mut header_line)?;

        // Split at the first whitespace.
        let ws_pos = header_line
            .iter()
            .position(|&b| b == b' ' || b == b'\t');

        let (name_bytes, comment_bytes) = match ws_pos {
            None => (header_line.as_slice(), &b""[..]),
            Some(p) => {
                let (n, rest) = header_line.split_at(p);
                let c = rest.get(1..).unwrap_or(&b""[..]);
                (n, c)
            }
        };

        let name = String::from_utf8_lossy(name_bytes).into_owned();
        let comment = String::from_utf8_lossy(comment_bytes).into_owned();

        // ---- Sequence -----------------------------------------------------
        let is_fastq = header_char == b'@';
        let mut seq: Vec<u8> = Vec::with_capacity(512);

        if is_fastq {
            // FASTQ: read sequence lines until the `+` separator line.
            loop {
                match self.peek_byte()? {
                    None => break, // truncated file — return what we have
                    Some(b'+') => {
                        // Consume and discard the entire `+…` separator line.
                        let mut discard: Vec<u8> = Vec::new();
                        self.read_line_into(&mut discard)?;
                        break;
                    }
                    Some(b'\n') | Some(b'\r') => {
                        self.read_byte()?; // skip blank lines
                    }
                    Some(_) => {
                        let mut line: Vec<u8> = Vec::new();
                        self.read_line_into(&mut line)?;
                        for &b in &line {
                            if b != b' ' && b != b'\t' {
                                seq.push(b);
                            }
                        }
                    }
                }
            }
        } else {
            // FASTA: read sequence lines until a new header sigil or EOF.
            //
            // When we peek and see `>` or `@`, we do NOT consume it — we leave
            // it in the underlying buffer.  The next call to `next()` will
            // consume it during the header-scan phase above.
            loop {
                match self.peek_byte()? {
                    None => break, // EOF
                    Some(b'>') | Some(b'@') => {
                        // Leave the sigil byte in the buffer for the next record.
                        break;
                    }
                    Some(b'\n') | Some(b'\r') => {
                        self.read_byte()?; // skip blank/empty lines
                    }
                    Some(_) => {
                        let mut line: Vec<u8> = Vec::new();
                        self.read_line_into(&mut line)?;
                        for &b in &line {
                            if b != b' ' && b != b'\t' {
                                seq.push(b);
                            }
                        }
                    }
                }
            }
            // last_char stays 0 — the sigil is still in the buffer.
        }

        // ---- Quality string (FASTQ only) ----------------------------------
        let mut qual: Vec<u8> = Vec::new();
        if is_fastq {
            // Accumulate quality lines until we reach seq.len() bytes.
            while qual.len() < seq.len() {
                match self.peek_byte()? {
                    None => break,
                    Some(b'\n') | Some(b'\r') => {
                        self.read_byte()?; // skip blank lines
                    }
                    Some(_) => {
                        let mut line: Vec<u8> = Vec::new();
                        self.read_line_into(&mut line)?;
                        for &b in &line {
                            if b != b'\r' {
                                qual.push(b);
                            }
                        }
                    }
                }
            }
            // Guard against over-reading in degenerate multi-line quality.
            qual.truncate(seq.len());

            // Validate.
            if !seq.is_empty() && qual.len() != seq.len() {
                bail!(
                    "FASTQ record '{}': quality length {} != sequence length {}",
                    name,
                    qual.len(),
                    seq.len()
                );
            }

            // After the quality block the next record's `@` (or `>`) may
            // follow immediately.  We peek without consuming so the next
            // call to `next()` finds it during the normal header scan.
            // (No special action needed — we simply don't consume the byte.)
        }

        Ok(Some(KSeq {
            name,
            comment,
            seq,
            qual,
        }))
    }
}

// ---------------------------------------------------------------------------
// File opening helper
// ---------------------------------------------------------------------------

/// Open a FASTA or FASTQ file, transparently decompressing gzip data.
///
/// Gzip format is detected by the `.gz` / `.gzip` extension *or* by the two
/// magic bytes `0x1f 0x8b` at the beginning of the file.  A boxed [`BufRead`]
/// trait object is returned so callers do not need to know the concrete type.
pub fn open_seq_file(path: &str) -> Result<Box<dyn BufRead>> {
    let file = File::open(path)
        .map_err(|e| anyhow::anyhow!("cannot open '{}': {}", path, e))?;

    // Decide by extension first; if ambiguous, peek at the magic bytes.
    let gz_by_ext = path.ends_with(".gz") || path.ends_with(".gzip");

    let is_gz = if gz_by_ext {
        true
    } else {
        // Peek via a BufReader wrapping a shared reference — no bytes consumed.
        let mut probe = BufReader::new(&file);
        let header = probe.fill_buf()?;
        header.len() >= 2 && header[0] == 0x1f && header[1] == 0x8b
    };

    if is_gz {
        // Re-open a fresh file handle for GzDecoder (the probe consumed nothing).
        let file2 = File::open(path)
            .map_err(|e| anyhow::anyhow!("cannot re-open '{}': {}", path, e))?;
        Ok(Box::new(BufReader::new(GzDecoder::new(file2))))
    } else {
        Ok(Box::new(BufReader::new(file)))
    }
}

// ---------------------------------------------------------------------------
// Batch reader
// ---------------------------------------------------------------------------

/// Read up to `batch_size` records from `reader`.
///
/// Returns fewer records (possibly zero) only when EOF is reached.  May be
/// called repeatedly on the same `reader` to stream through a file in batches.
///
/// Note: the `reader` must be a persistent, stateful [`BufRead`] across calls.
/// Each call creates an internal [`KSeqReader`] that only borrows the underlying
/// reader transiently; because FASTA sigil bytes (`>` / `@`) are left in the
/// buffer rather than being consumed into a separate field, the parser state is
/// fully captured in the underlying byte stream.
pub fn read_batch(reader: &mut dyn BufRead, batch_size: usize) -> Result<Vec<KSeq>> {
    // Thin adapter so `&mut dyn BufRead` satisfies `R: BufRead`.
    struct DynAdapter<'a>(&'a mut dyn BufRead);

    impl Read for DynAdapter<'_> {
        fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
            self.0.read(buf)
        }
    }

    impl BufRead for DynAdapter<'_> {
        fn fill_buf(&mut self) -> io::Result<&[u8]> {
            self.0.fill_buf()
        }
        fn consume(&mut self, amt: usize) {
            self.0.consume(amt);
        }
    }

    let mut parser = KSeqReader::new(DynAdapter(reader));
    let mut batch = Vec::with_capacity(batch_size.min(1024));

    while batch.len() < batch_size {
        match parser.next()? {
            None => break,
            Some(rec) => batch.push(rec),
        }
    }

    Ok(batch)
}

// ---------------------------------------------------------------------------
// Unit tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    fn make_reader(s: &str) -> KSeqReader<Cursor<Vec<u8>>> {
        KSeqReader::new(Cursor::new(s.as_bytes().to_vec()))
    }

    // -----------------------------------------------------------------------
    // FASTA
    // -----------------------------------------------------------------------

    #[test]
    fn single_fasta_record() {
        let input = ">seq1 a comment\nACGTACGT\n";
        let mut r = make_reader(input);
        let rec = r.next().unwrap().expect("expected a record");
        assert_eq!(rec.name, "seq1");
        assert_eq!(rec.comment, "a comment");
        assert_eq!(&rec.seq, b"ACGTACGT");
        assert!(rec.qual.is_empty());
        assert!(r.next().unwrap().is_none());
    }

    #[test]
    fn fasta_no_comment() {
        let input = ">simple\nNNNN\n";
        let mut r = make_reader(input);
        let rec = r.next().unwrap().unwrap();
        assert_eq!(rec.name, "simple");
        assert_eq!(rec.comment, "");
        assert_eq!(&rec.seq, b"NNNN");
    }

    #[test]
    fn multiline_fasta() {
        // Three lines: "ACGT" + "TGCA" + "AAAA" = 12 bases.
        let input = ">chr1\nACGT\nTGCA\nAAAA\n";
        let mut r = make_reader(input);
        let rec = r.next().unwrap().expect("expected a record");
        assert_eq!(&rec.seq, b"ACGTTGCAAAAA");
    }

    #[test]
    fn two_fasta_records() {
        let input = ">s1\nAAAA\n>s2\nCCCC\n";
        let mut r = make_reader(input);
        let r1 = r.next().unwrap().unwrap();
        let r2 = r.next().unwrap().unwrap();
        assert_eq!(r1.name, "s1");
        assert_eq!(&r1.seq, b"AAAA");
        assert_eq!(r2.name, "s2");
        assert_eq!(&r2.seq, b"CCCC");
        assert!(r.next().unwrap().is_none());
    }

    #[test]
    fn fasta_no_trailing_newline() {
        let input = ">s1\nACGT";
        let mut r = make_reader(input);
        let rec = r.next().unwrap().unwrap();
        assert_eq!(&rec.seq, b"ACGT");
    }

    // -----------------------------------------------------------------------
    // FASTQ
    // -----------------------------------------------------------------------

    #[test]
    fn single_fastq_record() {
        let input = "@read1 comment\nACGT\n+\nIIII\n";
        let mut r = make_reader(input);
        let rec = r.next().unwrap().unwrap();
        assert_eq!(rec.name, "read1");
        assert_eq!(rec.comment, "comment");
        assert_eq!(&rec.seq, b"ACGT");
        assert_eq!(&rec.qual, b"IIII");
        assert!(r.next().unwrap().is_none());
    }

    #[test]
    fn two_fastq_records() {
        let input = "@r1\nACGT\n+\nIIII\n@r2\nTTTT\n+\nJJJJ\n";
        let mut r = make_reader(input);
        let r1 = r.next().unwrap().unwrap();
        let r2 = r.next().unwrap().unwrap();
        assert_eq!(r1.name, "r1");
        assert_eq!(&r1.seq, b"ACGT");
        assert_eq!(&r1.qual, b"IIII");
        assert_eq!(r2.name, "r2");
        assert_eq!(&r2.seq, b"TTTT");
        assert_eq!(&r2.qual, b"JJJJ");
        assert!(r.next().unwrap().is_none());
    }

    #[test]
    fn fastq_quality_length_mismatch_is_error() {
        // Quality is shorter than the sequence — must produce an error.
        let input = "@r1\nACGTACGT\n+\nIII\n";
        let mut r = make_reader(input);
        assert!(r.next().is_err());
    }

    // -----------------------------------------------------------------------
    // Batch reading (exercises inter-batch state preservation)
    // -----------------------------------------------------------------------

    #[test]
    fn read_batch_collects_and_continues() {
        let input = ">a\nAAAA\n>b\nCCCC\n>c\nGGGG\n";
        let mut reader: Box<dyn BufRead> = Box::new(Cursor::new(input.as_bytes().to_vec()));

        let batch1 = read_batch(reader.as_mut(), 2).unwrap();
        assert_eq!(batch1.len(), 2);
        assert_eq!(batch1[0].name, "a");
        assert_eq!(batch1[1].name, "b");

        let batch2 = read_batch(reader.as_mut(), 2).unwrap();
        assert_eq!(batch2.len(), 1);
        assert_eq!(batch2[0].name, "c");

        let batch3 = read_batch(reader.as_mut(), 2).unwrap();
        assert!(batch3.is_empty());
    }

    #[test]
    fn read_batch_fastq_continues_across_calls() {
        let input = "@r1\nACGT\n+\nIIII\n@r2\nTTTT\n+\nJJJJ\n@r3\nGGGG\n+\nKKKK\n";
        let mut reader: Box<dyn BufRead> = Box::new(Cursor::new(input.as_bytes().to_vec()));

        let batch1 = read_batch(reader.as_mut(), 2).unwrap();
        assert_eq!(batch1.len(), 2);
        assert_eq!(batch1[0].name, "r1");
        assert_eq!(batch1[1].name, "r2");

        let batch2 = read_batch(reader.as_mut(), 2).unwrap();
        assert_eq!(batch2.len(), 1);
        assert_eq!(batch2[0].name, "r3");
    }

    // -----------------------------------------------------------------------
    // Edge cases
    // -----------------------------------------------------------------------

    #[test]
    fn empty_input_returns_none() {
        let mut r = make_reader("");
        assert!(r.next().unwrap().is_none());
    }

    #[test]
    fn fasta_with_windows_crlf() {
        let input = ">seq\r\nACGT\r\n";
        let mut r = make_reader(input);
        let rec = r.next().unwrap().unwrap();
        assert_eq!(&rec.seq, b"ACGT");
    }

    #[test]
    fn fasta_header_with_tab_separator() {
        let input = ">seq\tsome comment\nACGT\n";
        let mut r = make_reader(input);
        let rec = r.next().unwrap().unwrap();
        assert_eq!(rec.name, "seq");
        assert_eq!(rec.comment, "some comment");
        assert_eq!(&rec.seq, b"ACGT");
    }

    #[test]
    fn ten_consecutive_fasta_records() {
        let mut input = String::new();
        for i in 0..10usize {
            input.push_str(&format!(">seq{}\nACGT\n", i));
        }
        let mut r = make_reader(&input);
        let mut count = 0usize;
        while let Some(_rec) = r.next().unwrap() {
            count += 1;
        }
        assert_eq!(count, 10);
    }

    #[test]
    fn fasta_with_blank_lines_between_records() {
        let input = ">s1\nAAAA\n\n>s2\nCCCC\n";
        let mut r = make_reader(input);
        let r1 = r.next().unwrap().unwrap();
        let r2 = r.next().unwrap().unwrap();
        assert_eq!(&r1.seq, b"AAAA");
        assert_eq!(&r2.seq, b"CCCC");
    }
}
