/*************************************************************************************
                           The MIT License

   BWA-MEM2  (Sequence alignment using Burrows-Wheeler Transform),
   Copyright (C) 2019  Intel Corporation, Heng Li.

   Rust port of paired-end processing from bwamem_pair.cpp.

Authors: Rust port based on original C/C++ by
         Vasimuddin Md <vasimuddin.md@intel.com>; Sanchit Misra <sanchit.misra@intel.com>;
         Heng Li <hli@jimmy.harvard.edu>
*****************************************************************************************/

//! Paired-end insert-size estimation and pair-scoring.

use crate::bntseq::BntSeq;
use crate::sam::{mem_aln2sam, mem_reg2aln, write_unmapped_sam};
use crate::types::{MemAln, MemAlnReg, MemOpt, MemPeStat, BSeq, MEM_F_PE};

// ---------------------------------------------------------------------------
// Insert-size statistics
// ---------------------------------------------------------------------------

/// Orientation modes for paired-end reads.
///
/// | Mode | Read1 | Read2 | Typical use         |
/// |------|-------|-------|---------------------|
/// | 0 FR | fwd   | rev   | Illumina paired-end |
/// | 1 RF | rev   | fwd   | mate-pair           |
/// | 2 FF | fwd   | fwd   | (rare)              |
/// | 3 RR | rev   | rev   | (rare)              |
const N_ORI: usize = 4;

/// Estimate insert-size distributions from a batch of alignment pairs.
///
/// For each of the four orientation modes, collects insert sizes from pairs
/// where both ends have a single unambiguous primary alignment, then computes
/// mean and standard deviation.
pub fn mem_pestat(
    opt: &MemOpt,
    bns: &BntSeq,
    pairs: &[(Vec<MemAlnReg>, Vec<MemAlnReg>)],
) -> [MemPeStat; 4] {
    let mut sizes: [Vec<i64>; N_ORI] = Default::default();

    for (regs1, regs2) in pairs {
        // Both ends must have exactly one primary alignment
        let primaries1: Vec<&MemAlnReg> = regs1.iter().filter(|r| r.secondary < 0).collect();
        let primaries2: Vec<&MemAlnReg> = regs2.iter().filter(|r| r.secondary < 0).collect();
        if primaries1.len() != 1 || primaries2.len() != 1 { continue; }
        let r1 = primaries1[0];
        let r2 = primaries2[0];
        if r1.rid != r2.rid { continue; }

        let rev1 = r1.rb >= bns.l_pac;
        let rev2 = r2.rb >= bns.l_pac;

        // Map positions to forward-strand coordinates
        let p1 = if rev1 { 2 * bns.l_pac - r1.re } else { r1.rb };
        let p2 = if rev2 { 2 * bns.l_pac - r2.re } else { r2.rb };

        let isize = (p1 - p2).abs();
        if isize > opt.max_ins as i64 { continue; }

        let ori = match (rev1, rev2) {
            (false, true)  if p1 <= p2 => 0, // FR
            (true,  false) if p1 >  p2 => 1, // RF
            (false, false)              => 2, // FF
            (true,  true)               => 3, // RR
            _ => continue,
        };
        sizes[ori].push(isize);
    }

    let mut stats: [MemPeStat; 4] = Default::default();
    for i in 0..N_ORI {
        let v = &mut sizes[i];
        if v.len() < 20 {
            stats[i].failed = true;
            continue;
        }
        v.sort_unstable();
        // Trim 2.5% from each tail to get robust mean/std
        let lo = v.len() / 40;
        let hi = v.len() - lo;
        let trimmed = &v[lo..hi];
        let sum: i64 = trimmed.iter().sum();
        let avg = sum as f64 / trimmed.len() as f64;
        let var = trimmed.iter().map(|&x| (x as f64 - avg).powi(2)).sum::<f64>()
            / trimmed.len() as f64;
        let std = var.sqrt();
        let low  = (avg - 4.0 * std).max(1.0) as i32;
        let high = (avg + 4.0 * std).min(opt.max_ins as f64) as i32;
        stats[i] = MemPeStat { low, high, failed: false, avg, std };
    }
    stats
}

// ---------------------------------------------------------------------------
// Pair scoring
// ---------------------------------------------------------------------------

/// Check whether two alignment regions form a proper pair under orientation
/// `ori` with insert size in `[pes.low, pes.high]`.
fn is_proper_pair(r1: &MemAlnReg, r2: &MemAlnReg, pes: &MemPeStat, l_pac: i64) -> bool {
    if pes.failed { return false; }
    if r1.rid != r2.rid { return false; }

    let rev1 = r1.rb >= l_pac;
    let rev2 = r2.rb >= l_pac;
    let p1 = if rev1 { 2 * l_pac - r1.re } else { r1.rb };
    let p2 = if rev2 { 2 * l_pac - r2.re } else { r2.rb };
    let isize = (p1 - p2).abs() as i32;

    if isize < pes.low || isize > pes.high { return false; }

    // Orientation check: FR = (fwd, rev, p1 <= p2)
    matches!((rev1, rev2, p1 <= p2), (false, true, true))
        || matches!((rev1, rev2, p1 > p2), (true, false, true))
}

// ---------------------------------------------------------------------------
// Pair processing
// ---------------------------------------------------------------------------

/// Process one read pair: score pairing, set SAM flags, and generate SAM lines.
///
/// Returns `(sam_read1, sam_read2)`.
pub fn mem_pair(
    opt: &MemOpt,
    bns: &BntSeq,
    pac: &[u8],
    pes: &[MemPeStat; 4],
    seq1: &BSeq,
    seq2: &BSeq,
    regs1: &mut Vec<MemAlnReg>,
    regs2: &mut Vec<MemAlnReg>,
    _n_processed: i64,
) -> (String, String) {
    // Find the best-scoring proper pair
    let mut best_pair: Option<(usize, usize, i32)> = None;
    let mut best_pair_score = i32::MIN;

    for (i, r1) in regs1.iter().enumerate() {
        for (j, r2) in regs2.iter().enumerate() {
            for pe in pes.iter() {
                if is_proper_pair(r1, r2, pe, bns.l_pac) {
                    let pair_score = r1.score + r2.score;
                    if pair_score > best_pair_score {
                        best_pair_score = pair_score;
                        best_pair = Some((i, j, pair_score));
                    }
                }
            }
        }
    }

    // Set SAM flags for paired-end
    let pe_flag: u16 = 0x1; // paired

    let (sam1, sam2) = if let Some((i1, i2, _score)) = best_pair {
        // Proper pair found
        let r1 = &regs1[i1];
        let r2 = &regs2[i2];

        let mut aln1 = mem_reg2aln(opt, bns, pac, seq1.seq.len() as i32, &seq1.seq, r1);
        let mut aln2 = mem_reg2aln(opt, bns, pac, seq2.seq.len() as i32, &seq2.seq, r2);

        aln1.flag |= pe_flag | 0x2 | 0x40; // paired, proper, read1
        aln2.flag |= pe_flag | 0x2 | 0x80; // paired, proper, read2
        if aln2.is_rev { aln1.flag |= 0x20; }
        if aln1.is_rev { aln2.flag |= 0x20; }

        let s1 = mem_aln2sam(opt, bns, seq1, &aln1, &[], Some(&aln2));
        let s2 = mem_aln2sam(opt, bns, seq2, &aln2, &[], Some(&aln1));
        (s1, s2)
    } else {
        // No proper pair: output best individual alignments
        let aln1_opt = regs1.iter().find(|r| r.secondary < 0)
            .map(|r| mem_reg2aln(opt, bns, pac, seq1.seq.len() as i32, &seq1.seq, r));
        let aln2_opt = regs2.iter().find(|r| r.secondary < 0)
            .map(|r| mem_reg2aln(opt, bns, pac, seq2.seq.len() as i32, &seq2.seq, r));

        let s1 = match &aln1_opt {
            Some(aln) => {
                let mut a = aln.clone();
                a.flag |= pe_flag | 0x40;
                if let Some(m) = &aln2_opt { if m.is_rev { a.flag |= 0x20; } }
                mem_aln2sam(opt, bns, seq1, &a, &[], aln2_opt.as_ref())
            }
            None => write_unmapped_sam(seq1, pe_flag | 0x4 | 0x40, aln2_opt.as_ref(), bns),
        };
        let s2 = match &aln2_opt {
            Some(aln) => {
                let mut a = aln.clone();
                a.flag |= pe_flag | 0x80;
                if let Some(m) = &aln1_opt { if m.is_rev { a.flag |= 0x20; } }
                mem_aln2sam(opt, bns, seq2, &a, &[], aln1_opt.as_ref())
            }
            None => write_unmapped_sam(seq2, pe_flag | 0x4 | 0x80, aln1_opt.as_ref(), bns),
        };
        (s1, s2)
    };

    (sam1, sam2)
}

// ---------------------------------------------------------------------------
// Unit tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn pe_stat_default_failed() {
        let s: [MemPeStat; 4] = Default::default();
        for i in 0..4 {
            assert!(s[i].failed);
        }
    }

    #[test]
    fn is_proper_pair_failed_stat() {
        let pes = MemPeStat { failed: true, low: 100, high: 500, avg: 300.0, std: 50.0 };
        let r1 = MemAlnReg { rb: 100, re: 250, rid: 0, ..Default::default() };
        let r2 = MemAlnReg { rb: 300, re: 450, rid: 0, ..Default::default() };
        assert!(!is_proper_pair(&r1, &r2, &pes, 1000));
    }
}
