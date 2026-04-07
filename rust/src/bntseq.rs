/*************************************************************************************
                           The MIT License

   BWA-MEM2  (Sequence alignment using Burrows-Wheeler Transform),
   Copyright (C) 2019  Intel Corporation, Heng Li.

   Rust port of reference sequence management from bntseq.cpp / bntseq.h.

Authors: Rust port based on original C/C++ by Heng Li <hli@jimmy.harvard.edu>
*****************************************************************************************/

use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Read, Write};

use anyhow::{bail, Context, Result};

// ---------------------------------------------------------------------------
// Structures
// ---------------------------------------------------------------------------

/// Annotation record for one reference sequence (one line pair in .ann file).
#[derive(Debug, Clone, Default)]
pub struct BntAnn {
    /// Offset (0-based) of this sequence in the packed reference.
    pub offset: i64,
    /// Length of this sequence in bases.
    pub len: i32,
    /// Number of ambiguous bases in this sequence.
    pub n_ambs: i32,
    /// GenInfo identifier (gi number), often 0.
    pub gi: u32,
    /// Non-zero if this is an ALT/alternative contig.
    pub is_alt: i32,
    /// Sequence name (RNAME in SAM / FASTA header up to first whitespace).
    pub name: String,
    /// Full FASTA header annotation (everything after the name).
    pub anno: String,
}

/// One contiguous run of ambiguous (non-ACGT) bases in the packed reference.
#[derive(Debug, Clone, Default)]
pub struct BntAmb {
    /// Start offset (0-based) in the packed reference.
    pub offset: i64,
    /// Length of this ambiguous run.
    pub len: i32,
    /// The IUPAC character that was replaced (e.g. `b'N'`).
    pub amb: u8,
}

/// Full reference metadata: sequence annotations and ambiguity table.
#[derive(Debug, Clone)]
pub struct BntSeq {
    /// Total number of packed bases (forward strand only).
    pub l_pac: i64,
    /// Number of reference sequences.
    pub n_seqs: i32,
    /// Random seed used when replacing ambiguous bases.
    pub seed: u32,
    /// Per-sequence annotations.
    pub anns: Vec<BntAnn>,
    /// Number of ambiguous regions.
    pub n_holes: i32,
    /// Ambiguous regions.
    pub ambs: Vec<BntAmb>,
}

// ---------------------------------------------------------------------------
// I/O
// ---------------------------------------------------------------------------

impl BntSeq {
    /// Load metadata from `{prefix}.ann` and `{prefix}.amb`.
    pub fn restore(prefix: &str) -> Result<Self> {
        // --- .ann file ---
        let ann_path = format!("{}.ann", prefix);
        let ann_file = File::open(&ann_path)
            .with_context(|| format!("cannot open '{}'", ann_path))?;
        let mut ann_reader = BufReader::new(ann_file);

        // First line: l_pac n_seqs seed
        let mut line = String::new();
        ann_reader.read_line(&mut line).context("reading .ann header")?;
        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() < 3 {
            bail!("malformed .ann header: '{}'", line.trim());
        }
        let l_pac: i64 = parts[0].parse().context("parsing l_pac")?;
        let n_seqs: i32 = parts[1].parse().context("parsing n_seqs")?;
        let seed: u32 = parts[2].parse().context("parsing seed")?;

        let mut anns: Vec<BntAnn> = Vec::with_capacity(n_seqs as usize);

        for i in 0..n_seqs {
            // First line: "gi name [anno]"  (C format: gi and name first, optional anno)
            line.clear();
            ann_reader.read_line(&mut line)
                .with_context(|| format!("reading .ann record {} line 1", i))?;
            let p: Vec<&str> = line.splitn(3, char::is_whitespace).collect();
            let gi: u32 = p.first().and_then(|s| s.trim().parse().ok()).unwrap_or(0);
            let name = p.get(1).map(|s| s.trim().to_string()).unwrap_or_default();
            let raw_anno = p.get(2).map(|s| s.trim().to_string()).unwrap_or_default();
            // C writes "(null)" when anno is empty
            let anno = if raw_anno == "(null)" { String::new() } else { raw_anno };

            // Second line: "offset len n_ambs"
            line.clear();
            ann_reader.read_line(&mut line)
                .with_context(|| format!("reading .ann record {} line 2", i))?;
            let p: Vec<&str> = line.split_whitespace().collect();
            if p.len() < 3 {
                bail!("malformed .ann record {}: '{}'", i, line.trim());
            }
            let offset: i64 = p[0].parse().context("parsing offset")?;
            let len: i32 = p[1].parse().context("parsing len")?;
            let n_ambs: i32 = p[2].parse().context("parsing n_ambs")?;

            anns.push(BntAnn { offset, len, n_ambs, gi, is_alt: 0, name, anno });
        }

        // --- .amb file ---
        let amb_path = format!("{}.amb", prefix);
        let amb_file = File::open(&amb_path)
            .with_context(|| format!("cannot open '{}'", amb_path))?;
        let mut amb_reader = BufReader::new(amb_file);

        line.clear();
        amb_reader.read_line(&mut line).context("reading .amb header")?;
        let p: Vec<&str> = line.split_whitespace().collect();
        let n_holes: i32 = if p.len() >= 3 {
            p[2].parse().unwrap_or(0)
        } else if p.len() >= 2 {
            p[1].parse().unwrap_or(0)
        } else {
            0
        };

        let mut ambs: Vec<BntAmb> = Vec::with_capacity(n_holes as usize);
        for _ in 0..n_holes {
            line.clear();
            if amb_reader.read_line(&mut line)? == 0 { break; }
            let p: Vec<&str> = line.split_whitespace().collect();
            if p.len() < 3 { continue; }
            let offset: i64 = p[0].parse().unwrap_or(0);
            let len: i32 = p[1].parse().unwrap_or(0);
            let amb: u8 = p[2].bytes().next().unwrap_or(b'N');
            ambs.push(BntAmb { offset, len, amb });
        }

        Ok(BntSeq { l_pac, n_seqs, seed, anns, n_holes, ambs })
    }

    /// Write metadata to `{prefix}.ann` and `{prefix}.amb`.
    pub fn dump(&self, prefix: &str) -> Result<()> {
        let ann_path = format!("{}.ann", prefix);
        let ann_file = File::create(&ann_path)
            .with_context(|| format!("cannot create '{}'", ann_path))?;
        let mut w = BufWriter::new(ann_file);
        writeln!(w, "{} {} {}", self.l_pac, self.n_seqs, self.seed)?;
        for ann in &self.anns {
            // Match C format: "gi name [anno]" then "offset len n_ambs"
            if ann.anno.is_empty() {
                writeln!(w, "{} {}", ann.gi, ann.name)?;
            } else {
                writeln!(w, "{} {} {}", ann.gi, ann.name, ann.anno)?;
            }
            writeln!(w, "{} {} {}", ann.offset, ann.len, ann.n_ambs)?;
        }
        w.flush()?;

        let amb_path = format!("{}.amb", prefix);
        let amb_file = File::create(&amb_path)
            .with_context(|| format!("cannot create '{}'", amb_path))?;
        let mut w = BufWriter::new(amb_file);
        writeln!(w, "{} {} {}", self.l_pac, self.n_seqs, self.n_holes)?;
        for amb in &self.ambs {
            writeln!(w, "{} {} {}", amb.offset, amb.len, amb.amb as char)?;
        }
        w.flush()?;
        Ok(())
    }

    /// Binary search: return the index of the reference sequence containing
    /// packed position `pos` (0-based, forward strand). Returns -1 if not found.
    pub fn pos2rid(&self, pos: i64) -> i32 {
        if pos < 0 || self.anns.is_empty() { return -1; }
        let mut lo: i32 = 0;
        let mut hi: i32 = self.n_seqs - 1;
        while lo < hi {
            let mid = (lo + hi + 1) / 2;
            if self.anns[mid as usize].offset <= pos { lo = mid; } else { hi = mid - 1; }
        }
        let ann = &self.anns[lo as usize];
        if pos < ann.offset + ann.len as i64 { lo } else { -1 }
    }

    /// Build a `BntSeq` and packed reference from a FASTA file, writing
    /// `.pac`, `.ann`, `.amb` files with the given prefix.
    /// Returns `(bntseq, packed_forward_sequence)`.
    pub fn fasta_to_bntseq(fasta_path: &str, prefix: &str) -> Result<(Self, Vec<u8>)> {
        let file = File::open(fasta_path)
            .with_context(|| format!("cannot open FASTA '{}'", fasta_path))?;
        let reader = BufReader::new(file);

        let mut anns: Vec<BntAnn> = Vec::new();
        let mut ambs: Vec<BntAmb> = Vec::new();
        let mut pac: Vec<u8> = Vec::new();
        let mut l_pac: i64 = 0;

        // Simple LCG for N-substitution (matches bwa seed=11)
        let mut rng: u32 = 11;
        let lcg = |s: u32| -> u32 { s.wrapping_mul(1103515245).wrapping_add(12345) };

        let mut cur_name = String::new();
        let mut cur_anno = String::new();
        let mut cur_offset: i64 = 0;
        let mut cur_len: i32 = 0;
        let mut cur_n_ambs: i32 = 0;
        let mut in_amb = false;
        let mut amb_start: i64 = 0;
        let mut amb_char: u8 = b'N';
        let mut in_seq = false;

        for raw in reader.lines() {
            let line = raw.context("reading FASTA")?;
            let line = line.trim_end();
            if line.starts_with('>') {
                if in_seq && cur_len > 0 {
                    if in_amb {
                        let run = (l_pac - amb_start) as i32;
                        if run > 0 { ambs.push(BntAmb { offset: amb_start, len: run, amb: amb_char }); }
                        in_amb = false;
                    }
                    anns.push(BntAnn {
                        offset: cur_offset, len: cur_len, n_ambs: cur_n_ambs,
                        gi: 0, is_alt: 0, name: cur_name.clone(), anno: cur_anno.clone(),
                    });
                }
                let hdr = &line[1..];
                let mut it = hdr.splitn(2, char::is_whitespace);
                cur_name  = it.next().unwrap_or("").to_string();
                cur_anno  = it.next().unwrap_or("").to_string();
                cur_offset = l_pac;
                cur_len = 0;
                cur_n_ambs = 0;
                in_seq = true;
            } else if in_seq {
                for &b in line.as_bytes() {
                    let two_bit = nuc_to_2bit(b);
                    let base: u8 = if two_bit < 4 {
                        if in_amb {
                            let run = (l_pac - amb_start) as i32;
                            if run > 0 { ambs.push(BntAmb { offset: amb_start, len: run, amb: amb_char }); }
                            in_amb = false;
                        }
                        two_bit
                    } else {
                        if !in_amb {
                            amb_start = l_pac;
                            amb_char = b.to_ascii_uppercase();
                            in_amb = true;
                            cur_n_ambs += 1;
                        }
                        rng = lcg(rng);
                        (rng >> 24) as u8 & 3
                    };
                    let byte_pos = (l_pac >> 2) as usize;
                    let shift = ((3 - (l_pac & 3)) * 2) as u32;
                    if byte_pos >= pac.len() { pac.push(0); }
                    pac[byte_pos] |= base << shift;
                    l_pac += 1;
                    cur_len += 1;
                }
            }
        }
        // Flush last sequence
        if in_seq && cur_len > 0 {
            if in_amb {
                let run = (l_pac - amb_start) as i32;
                if run > 0 { ambs.push(BntAmb { offset: amb_start, len: run, amb: amb_char }); }
            }
            anns.push(BntAnn {
                offset: cur_offset, len: cur_len, n_ambs: cur_n_ambs,
                gi: 0, is_alt: 0, name: cur_name, anno: cur_anno,
            });
        }

        let n_seqs = anns.len() as i32;
        let n_holes = ambs.len() as i32;
        let bnt = BntSeq { l_pac, n_seqs, seed: 11, anns, n_holes, ambs };

        // Write .pac
        let pac_path = format!("{}.pac", prefix);
        let mut pf = File::create(&pac_path)
            .with_context(|| format!("cannot create '{}'", pac_path))?;
        pf.write_all(&pac)?;
        pf.flush()?;

        bnt.dump(prefix)?;
        Ok((bnt, pac))
    }
}

// ---------------------------------------------------------------------------
// Standalone helpers
// ---------------------------------------------------------------------------

/// Map an ASCII nucleotide to its 2-bit encoding (A=0, C=1, G=2, T=3, other=4).
#[inline]
pub fn nuc_to_2bit(c: u8) -> u8 {
    match c {
        b'A' | b'a' => 0,
        b'C' | b'c' => 1,
        b'G' | b'g' => 2,
        b'T' | b't' => 3,
        _ => 4,
    }
}

/// Load the packed reference from `{prefix}.pac`.
pub fn restore_pac(prefix: &str) -> Result<Vec<u8>> {
    let pac_path = format!("{}.pac", prefix);
    let mut file = File::open(&pac_path)
        .with_context(|| format!("cannot open '{}'", pac_path))?;
    let mut buf = Vec::new();
    file.read_to_end(&mut buf)?;
    Ok(buf)
}

/// Extract a single 2-bit base from a packed array.
#[inline]
pub fn pac2nt(pac: &[u8], pos: i64) -> u8 {
    let byte_idx = (pos >> 2) as usize;
    let shift = ((3 - (pos & 3)) * 2) as u32;
    if byte_idx < pac.len() { (pac[byte_idx] >> shift) & 3 } else { 0 }
}

/// Set a single 2-bit base in a packed array.
#[inline]
fn pack_base(pac: &mut [u8], pos: i64, base: u8) {
    let byte_idx = (pos >> 2) as usize;
    let shift = ((3 - (pos & 3)) * 2) as u32;
    pac[byte_idx] &= !(3u8 << shift);
    pac[byte_idx] |= (base & 3) << shift;
}

/// Build the forward+reverse-complement packed sequence (2*l_pac bases).
pub fn restore_pac_rev(fwd_pac: &[u8], l_pac: i64) -> Vec<u8> {
    let total = 2 * l_pac;
    let nbytes = ((total + 3) / 4) as usize;
    let mut pac = vec![0u8; nbytes];
    for i in 0..l_pac {
        pack_base(&mut pac, i, pac2nt(fwd_pac, i));
    }
    for i in l_pac..total {
        let fwd_base = pac2nt(fwd_pac, total - 1 - i);
        pack_base(&mut pac, i, 3 - fwd_base);
    }
    pac
}

// ---------------------------------------------------------------------------
// Unit tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn nuc_to_2bit_basic() {
        assert_eq!(nuc_to_2bit(b'A'), 0);
        assert_eq!(nuc_to_2bit(b'C'), 1);
        assert_eq!(nuc_to_2bit(b'G'), 2);
        assert_eq!(nuc_to_2bit(b'T'), 3);
        assert_eq!(nuc_to_2bit(b'N'), 4);
    }

    #[test]
    fn pac2nt_roundtrip() {
        let mut pac = vec![0u8; 2];
        pack_base(&mut pac, 0, 0);
        pack_base(&mut pac, 1, 1);
        pack_base(&mut pac, 2, 2);
        pack_base(&mut pac, 3, 3);
        assert_eq!(pac2nt(&pac, 0), 0);
        assert_eq!(pac2nt(&pac, 1), 1);
        assert_eq!(pac2nt(&pac, 2), 2);
        assert_eq!(pac2nt(&pac, 3), 3);
    }

    #[test]
    fn pos2rid_basic() {
        let bnt = BntSeq {
            l_pac: 20, n_seqs: 2, seed: 0, n_holes: 0,
            anns: vec![
                BntAnn { offset: 0,  len: 10, ..Default::default() },
                BntAnn { offset: 10, len: 10, ..Default::default() },
            ],
            ambs: vec![],
        };
        assert_eq!(bnt.pos2rid(0),  0);
        assert_eq!(bnt.pos2rid(9),  0);
        assert_eq!(bnt.pos2rid(10), 1);
        assert_eq!(bnt.pos2rid(19), 1);
        assert_eq!(bnt.pos2rid(20), -1);
    }
}
