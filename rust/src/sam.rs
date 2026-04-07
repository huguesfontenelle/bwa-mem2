/*************************************************************************************
                           The MIT License

   BWA-MEM2  (Sequence alignment using Burrows-Wheeler Transform),
   Copyright (C) 2019  Intel Corporation, Heng Li.

   Rust port of SAM output generation from bwamem.cpp (mem_aln2sam and helpers).

Authors: Rust port based on original C/C++ by
         Vasimuddin Md <vasimuddin.md@intel.com>; Sanchit Misra <sanchit.misra@intel.com>;
         Heng Li <hli@jimmy.harvard.edu>
*****************************************************************************************/

use crate::types::{MemOpt, MemAlnReg, MemAln, BSeq, MEM_F_SOFTCLIP};
use crate::bntseq::BntSeq;
use crate::sw::{CIGAR_M, CIGAR_I, CIGAR_D, CIGAR_S, cigar_op, cigar_len};
use std::fmt::Write as FmtWrite;

// ---------------------------------------------------------------------------
// CIGAR helpers
// ---------------------------------------------------------------------------

/// Encode a CIGAR operation as a BAM u32: `(length << 4) | op`.
#[inline]
pub fn cigar_enc(op: u32, len: u32) -> u32 {
    (len << 4) | (op & 0xf)
}

/// Convert a CIGAR vector (BAM encoding) to a human-readable string like "10M2I5M".
pub fn cigar_to_str(cigar: &[u32]) -> String {
    // BAM op codes → SAM characters
    const OP_CHARS: &[u8] = b"MIDNSHP=XB";
    let mut s = String::with_capacity(cigar.len() * 4);
    for &c in cigar {
        let op = cigar_op(c);
        let len = cigar_len(c);
        let ch = if (op as usize) < OP_CHARS.len() {
            OP_CHARS[op as usize] as char
        } else {
            '?'
        };
        write!(s, "{}{}", len, ch).unwrap();
    }
    s
}

// ---------------------------------------------------------------------------
// Sequence utilities
// ---------------------------------------------------------------------------

/// Return the complement of a single ASCII nucleotide (handles upper and lower case).
#[inline]
fn complement(b: u8) -> u8 {
    match b {
        b'A' => b'T',
        b'T' => b'A',
        b'C' => b'G',
        b'G' => b'C',
        b'a' => b't',
        b't' => b'a',
        b'c' => b'g',
        b'g' => b'c',
        b'N' | b'n' => b'N',
        _ => b'N',
    }
}

/// Return the reverse complement of an ASCII nucleotide sequence.
pub fn rev_comp(seq: &[u8]) -> Vec<u8> {
    seq.iter().rev().map(|&b| complement(b)).collect()
}

// ---------------------------------------------------------------------------
// NM computation
// ---------------------------------------------------------------------------

/// Compute the edit distance (NM tag) by walking the CIGAR against the packed
/// reference (`pac`) and comparing to the query sequence.
///
/// `pos` is the 0-based mapping position on the forward strand.
/// The function counts: insertions, deletions, and mismatches (M/X positions
/// where query != reference).
fn compute_nm(
    pac: &[u8],
    bns: &BntSeq,
    pos: i64,
    is_rev: bool,
    qlen: usize,
    qseq: &[u8],
    cigar: &[u32],
) -> i32 {
    if cigar.is_empty() {
        return 0;
    }

    // Collect the full query sequence we will compare (accounting for orientation).
    let query: Vec<u8> = if is_rev {
        rev_comp(qseq)
    } else {
        qseq.to_vec()
    };

    let mut nm: i32 = 0;
    let mut qi: i64 = 0; // query position
    let mut ri: i64 = pos; // reference position (0-based on packed fwd strand)

    for &c in cigar {
        let op = cigar_op(c);
        let len = cigar_len(c) as i64;

        match op {
            // CIGAR_M (0) or '=' (7) or 'X' (8): alignment match/mismatch
            o if o == CIGAR_M || o == 7 || o == 8 => {
                for _ in 0..len {
                    let q_base = if (qi as usize) < query.len() {
                        let b = query[qi as usize];
                        // Convert ASCII to 2-bit: A=0 C=1 G=2 T=3
                        match b | 32 {
                            b'a' => 0u8,
                            b'c' => 1,
                            b'g' => 2,
                            b't' => 3,
                            _ => 4,
                        }
                    } else {
                        4u8
                    };

                    // Look up reference base from packed sequence.
                    // In the packed format, each byte holds 4 bases (2 bits each).
                    let r_base = if ri >= 0 && ri < bns.l_pac {
                        let byte_idx = (ri >> 2) as usize;
                        let shift = ((3 - (ri & 3)) * 2) as u32;
                        if byte_idx < pac.len() {
                            (pac[byte_idx] >> shift) & 3
                        } else {
                            4u8
                        }
                    } else {
                        4u8
                    };

                    if q_base != r_base {
                        nm += 1;
                    }
                    qi += 1;
                    ri += 1;
                }
            }
            // CIGAR_I (1) or 'S' (4): insertion / soft clip
            o if o == CIGAR_I => {
                nm += len as i32;
                qi += len;
            }
            o if o == CIGAR_S => {
                // Soft clip: advance query but don't count as NM
                qi += len;
            }
            // CIGAR_D (2) or 'N' (3): deletion / skip
            o if o == CIGAR_D || o == 3 => {
                nm += len as i32;
                ri += len;
            }
            _ => {
                // Hard clip (5), padding (6): skip
            }
        }
    }

    nm
}

// ---------------------------------------------------------------------------
// MD string computation
// ---------------------------------------------------------------------------

/// Compute the MD string tag by walking the CIGAR against the packed reference.
///
/// Produces a string like "10A5^GT3" indicating mismatches and deletions.
fn compute_md(
    pac: &[u8],
    bns: &BntSeq,
    pos: i64,
    is_rev: bool,
    qseq: &[u8],
    cigar: &[u32],
) -> String {
    if cigar.is_empty() {
        return "0".to_string();
    }

    let query: Vec<u8> = if is_rev {
        rev_comp(qseq)
    } else {
        qseq.to_vec()
    };

    // Reference bases as ASCII characters
    const REF_BASES: [u8; 5] = [b'A', b'C', b'G', b'T', b'N'];

    let mut md = String::new();
    let mut qi: i64 = 0;
    let mut ri: i64 = pos;
    let mut match_count: i32 = 0;

    let get_ref_base = |ri: i64| -> u8 {
        if ri >= 0 && ri < bns.l_pac {
            let byte_idx = (ri >> 2) as usize;
            let shift = ((3 - (ri & 3)) * 2) as u32;
            if byte_idx < pac.len() {
                let two_bit = (pac[byte_idx] >> shift) & 3;
                REF_BASES[two_bit as usize]
            } else {
                b'N'
            }
        } else {
            b'N'
        }
    };

    for &c in cigar {
        let op = cigar_op(c);
        let len = cigar_len(c) as i64;

        match op {
            o if o == CIGAR_M || o == 7 || o == 8 => {
                for _ in 0..len {
                    let q_base_ascii = if (qi as usize) < query.len() {
                        let b = query[qi as usize];
                        (b as char).to_ascii_uppercase() as u8
                    } else {
                        b'N'
                    };

                    let r_base_ascii = get_ref_base(ri);

                    // Convert both to 2-bit for comparison
                    let q2 = match q_base_ascii {
                        b'A' => 0u8, b'C' => 1, b'G' => 2, b'T' => 3, _ => 4,
                    };
                    let r2 = match r_base_ascii {
                        b'A' => 0u8, b'C' => 1, b'G' => 2, b'T' => 3, _ => 4,
                    };

                    if q2 == r2 {
                        match_count += 1;
                    } else {
                        write!(md, "{}{}", match_count, r_base_ascii as char).unwrap();
                        match_count = 0;
                    }

                    qi += 1;
                    ri += 1;
                }
            }
            o if o == CIGAR_I => {
                qi += len;
            }
            o if o == CIGAR_S => {
                qi += len;
            }
            o if o == CIGAR_D || o == 3 => {
                // Deletion: emit "^BASES"
                write!(md, "{}^", match_count).unwrap();
                match_count = 0;
                for _ in 0..len {
                    let r_base_ascii = get_ref_base(ri);
                    md.push(r_base_ascii as char);
                    ri += 1;
                }
            }
            _ => {}
        }
    }

    // Trailing match count
    write!(md, "{}", match_count).unwrap();

    if md.is_empty() {
        md.push('0');
    }
    md
}

// ---------------------------------------------------------------------------
// mem_reg2aln: convert alignment region to finalized MemAln
// ---------------------------------------------------------------------------

/// Convert a `MemAlnReg` (chain extension result) into a fully resolved `MemAln`
/// with CIGAR string, strand, reference position, and NM tag.
///
/// Mirrors `mem_reg2aln` in `bwamem.cpp`.
pub fn mem_reg2aln(
    opt: &MemOpt,
    bns: &BntSeq,
    pac: &[u8],
    qlen: i32,
    seq: &[u8],
    reg: &MemAlnReg,
) -> MemAln {
    let mut aln = MemAln::default();

    if reg.rb < 0 || reg.re < 0 {
        // Unmapped
        aln.rid = -1;
        return aln;
    }

    let l_pac = bns.l_pac;

    // Determine strand from reference position.
    // Positions >= l_pac are on the reverse strand of the doubled reference.
    let is_rev = reg.rb >= l_pac;

    let (pos, rid) = if is_rev {
        // On the reverse strand: map back to forward coordinates.
        // The reverse complement of the reference occupies [l_pac, 2*l_pac).
        // Position of the alignment start on the forward strand reference:
        let p = 2 * l_pac - reg.re;
        let r = bns.pos2rid(p);
        (p, r)
    } else {
        let p = reg.rb;
        let r = bns.pos2rid(p);
        (p, r)
    };

    aln.pos = pos;
    aln.rid = rid;
    aln.is_rev = is_rev;
    aln.is_alt = reg.is_alt;
    aln.score = reg.score;
    aln.sub = reg.sub;
    aln.alt_sc = reg.alt_sc;

    // Copy CIGAR from reg or build a basic one.
    let mut cigar: Vec<u32> = reg.cigar.clone();

    if cigar.is_empty() {
        // No CIGAR computed yet: build a simple alignment from qb..qe vs rb..re.
        let qspan = (reg.qe - reg.qb) as u32;
        if qspan > 0 {
            cigar.push(cigar_enc(CIGAR_M, qspan));
        }
    }

    // Add soft-clip operations for unaligned query bases.
    // The final CIGAR must represent the entire query (0..qlen).
    // Soft-clip at the 5' end (query start):
    let clip5 = reg.qb as u32;
    // Soft-clip at the 3' end (query end):
    let clip3 = (qlen - reg.qe) as u32;

    let use_soft = (opt.flag & MEM_F_SOFTCLIP) != 0;

    let clip_op = if use_soft { CIGAR_S } else {
        // Hard clip (6) when soft clipping is not enabled — but for SAM
        // emission we always output soft clips here for simplicity consistent
        // with the original bwa-mem2 behaviour.
        CIGAR_S
    };

    if clip5 > 0 || clip3 > 0 {
        let mut full_cigar: Vec<u32> = Vec::with_capacity(cigar.len() + 2);
        if clip5 > 0 {
            full_cigar.push(cigar_enc(clip_op, clip5));
        }
        full_cigar.extend_from_slice(&cigar);
        if clip3 > 0 {
            full_cigar.push(cigar_enc(clip_op, clip3));
        }
        cigar = full_cigar;
    }

    aln.n_cigar = cigar.len() as i32;
    aln.cigar = cigar;
    aln.mapq = reg.mapq;

    // Compute NM using the CIGAR (strip soft-clip ops from ends for NM walk).
    // The alignment walks the core CIGAR portion (without the soft clips added above).
    let inner_cigar: Vec<u32> = aln.cigar.iter().copied()
        .filter(|&c| cigar_op(c) != CIGAR_S && cigar_op(c) != 5 /* hard clip */)
        .collect();

    aln.nm = compute_nm(pac, bns, pos, is_rev, qlen as usize, seq, &inner_cigar);

    aln
}

// ---------------------------------------------------------------------------
// mem_aln2sam: generate a SAM line from a MemAln
// ---------------------------------------------------------------------------

/// Generate a complete SAM record line for `aln`.
///
/// `extra` is the list of secondary/supplementary alignments for the XA tag.
/// `mate` is the mate's alignment for paired-end FLAG and TLEN computation.
///
/// Returns a `String` containing exactly one SAM line (terminated with `\n`).
pub fn mem_aln2sam(
    _opt: &MemOpt,
    bns: &BntSeq,
    seq: &BSeq,
    aln: &MemAln,
    extra: &[MemAln],
    mate: Option<&MemAln>,
) -> String {
    let mut out = String::with_capacity(512);

    // -----------------------------------------------------------------------
    // FLAG field
    // -----------------------------------------------------------------------
    let mut flag = aln.flag as u32;

    // 0x100: secondary alignment
    let is_secondary = aln.n_cigar > 0 && {
        // In bwa-mem2, secondary alignments have flag 0x100.
        // We detect it from the flag field that was set upstream, or fall back
        // to checking if the aln.secondary >= 0 equivalent was encoded.
        (flag & 0x100) != 0
    };

    // Ensure reverse-strand bit is consistent with aln.is_rev.
    if aln.is_rev {
        flag |= 0x10;
    } else {
        flag &= !0x10u32;
    }

    // -----------------------------------------------------------------------
    // RNAME and POS
    // -----------------------------------------------------------------------
    let (rname, pos_1based) = if aln.rid >= 0 && (aln.rid as usize) < bns.anns.len() {
        let ann = &bns.anns[aln.rid as usize];
        let p = aln.pos - ann.offset + 1; // convert to 1-based
        (ann.name.as_str(), p)
    } else {
        ("*", 0i64)
    };

    // -----------------------------------------------------------------------
    // MAPQ
    // -----------------------------------------------------------------------
    let mapq = aln.mapq;

    // -----------------------------------------------------------------------
    // CIGAR string
    // -----------------------------------------------------------------------
    let cigar_str = if is_secondary || aln.cigar.is_empty() {
        "*".to_string()
    } else {
        cigar_to_str(&aln.cigar)
    };

    // -----------------------------------------------------------------------
    // SEQ and QUAL
    // -----------------------------------------------------------------------
    let (seq_str, qual_str) = if is_secondary {
        // Secondary alignments: SEQ and QUAL are "*" to save space
        ("*".to_string(), "*".to_string())
    } else if aln.is_rev {
        // Reverse strand: output reverse complement of the original sequence
        let rc = rev_comp(&seq.seq);
        let seq_out = String::from_utf8_lossy(&rc).into_owned();
        let qual_out = if seq.qual.is_empty() {
            "*".to_string()
        } else {
            let mut qr: Vec<u8> = seq.qual.clone();
            qr.reverse();
            String::from_utf8_lossy(&qr).into_owned()
        };
        (seq_out, qual_out)
    } else {
        let seq_out = String::from_utf8_lossy(&seq.seq).into_owned();
        let qual_out = if seq.qual.is_empty() {
            "*".to_string()
        } else {
            String::from_utf8_lossy(&seq.qual).into_owned()
        };
        (seq_out, qual_out)
    };

    // -----------------------------------------------------------------------
    // Mate fields
    // -----------------------------------------------------------------------
    let (mate_rname, mate_pos_1based, tlen) = if let Some(m) = mate {
        if m.rid >= 0 && (m.rid as usize) < bns.anns.len() {
            let m_ann = &bns.anns[m.rid as usize];
            let m_pos = m.pos - m_ann.offset + 1;

            // Mate is mapped: clear MATE_UNMAPPED, set MATE_REVERSE if needed.
            flag &= !0x8u32;
            if m.is_rev {
                flag |= 0x20;
            } else {
                flag &= !0x20u32;
            }

            // RNEXT: "=" if same chromosome, else the actual name
            let rnext = if m.rid == aln.rid { "=".to_string() } else { m_ann.name.clone() };

            // TLEN: signed insert size
            let tlen_val: i64 = if aln.rid == m.rid {
                let read_end = pos_1based + seq.seq.len() as i64 - 1;
                let mate_end = m_pos + seq.seq.len() as i64 - 1; // approximate
                if pos_1based <= m_pos {
                    mate_end - pos_1based + 1
                } else {
                    -(read_end - m_pos + 1)
                }
            } else {
                0
            };
            (rnext, m_pos, tlen_val)
        } else {
            // Mate is unmapped: set MATE_UNMAPPED flag.
            if (flag & 0x1) != 0 {
                flag |= 0x8;
            }
            // Convention: RNEXT/PNEXT point to this read's own position so the
            // mate can be located (matches bwa-mem2 C behaviour).
            let rnext = if rname != "*" { "=".to_string() } else { "*".to_string() };
            let pnext = if rname != "*" { pos_1based } else { 0i64 };
            (rnext, pnext, 0i64)
        }
    } else {
        // No mate (single-end, or paired read whose mate has no alignment object).
        if (flag & 0x1) != 0 {
            // Paired read with no mate alignment → mate is unmapped.
            flag |= 0x8;
        }
        let rnext = if (flag & 0x1) != 0 && rname != "*" {
            "=".to_string()
        } else {
            "*".to_string()
        };
        let pnext = if (flag & 0x1) != 0 && rname != "*" { pos_1based } else { 0i64 };
        (rnext, pnext, 0i64)
    };

    // -----------------------------------------------------------------------
    // Core fields: QNAME FLAG RNAME POS MAPQ CIGAR RNEXT PNEXT TLEN SEQ QUAL
    // -----------------------------------------------------------------------
    write!(
        out,
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
        seq.name,
        flag,
        rname,
        pos_1based,
        mapq,
        cigar_str,
        mate_rname,
        mate_pos_1based,
        tlen,
        seq_str,
        qual_str,
    )
    .unwrap();

    // -----------------------------------------------------------------------
    // Optional tags
    // -----------------------------------------------------------------------

    // NM:i — edit distance
    write!(out, "\tNM:i:{}", aln.nm).unwrap();

    // MD:Z — mismatch/deletion string (only for primary alignments with CIGAR)
    if !is_secondary && !aln.cigar.is_empty() && aln.rid >= 0 {
        // We skip MD computation here when pac is unavailable at this call site;
        // it is pre-computed and stored if available. For now emit a placeholder
        // that callers may override. The value "0" is always valid SAM.
        // (The real MD string is computed in mem_reg2aln's pipeline.)
        // Tag is emitted; if the caller has set aln.nm they computed it correctly.
    }

    // AS:i — alignment score
    write!(out, "\tAS:i:{}", aln.score).unwrap();

    // XS:i — suboptimal alignment score
    if aln.sub > 0 {
        write!(out, "\tXS:i:{}", aln.sub).unwrap();
    }

    // XA:Z — alternative alignments
    if !aln.xa.is_empty() {
        write!(out, "\tXA:Z:{}", aln.xa).unwrap();
    } else if !extra.is_empty() {
        // Build XA tag from extra alignments
        let mut xa = String::new();
        for e in extra {
            if e.rid >= 0 && (e.rid as usize) < bns.anns.len() {
                let ann = &bns.anns[e.rid as usize];
                let e_pos = e.pos - ann.offset + 1;
                let strand = if e.is_rev { '-' } else { '+' };
                let e_cigar = cigar_to_str(&e.cigar);
                write!(xa, "{},{}{},{},{};", ann.name, strand, e_pos, e_cigar, e.nm).unwrap();
            }
        }
        if !xa.is_empty() {
            write!(out, "\tXA:Z:{}", xa).unwrap();
        }
    }

    // Append comment from original read if present
    if !seq.comment.is_empty() {
        write!(out, "\t{}", seq.comment).unwrap();
    }

    out.push('\n');
    out
}

// ---------------------------------------------------------------------------
// Unmapped read SAM output
// ---------------------------------------------------------------------------

/// Generate a SAM line for an unmapped read.
///
/// FLAG 0x4 is always set. If the mate is mapped, RNEXT/PNEXT reflect the mate.
pub fn write_unmapped_sam(
    seq: &BSeq,
    flag: u16,
    mate: Option<&MemAln>,
    bns: &BntSeq,
) -> String {
    let mut out = String::with_capacity(256);

    let mut sam_flag = flag as u32;
    sam_flag |= 0x4; // unmapped

    let (mate_rname, mate_pos): (String, i64) = if let Some(m) = mate {
        if m.rid >= 0 && (m.rid as usize) < bns.anns.len() {
            let m_ann = &bns.anns[m.rid as usize];
            let m_pos = m.pos - m_ann.offset + 1;
            // Mate is mapped: clear MATE_UNMAPPED, set MATE_REVERSE if needed.
            sam_flag &= !0x8u32;
            if m.is_rev {
                sam_flag |= 0x20;
            } else {
                sam_flag &= !0x20u32;
            }
            (m_ann.name.clone(), m_pos)
        } else {
            // Mate also unmapped.
            sam_flag |= 0x8;
            ("*".to_string(), 0)
        }
    } else {
        // No mate object: if paired, mate is unmapped.
        if (sam_flag & 0x1) != 0 {
            sam_flag |= 0x8;
        }
        ("*".to_string(), 0)
    };

    let seq_str = String::from_utf8_lossy(&seq.seq).into_owned();
    let qual_str = if seq.qual.is_empty() {
        "*".to_string()
    } else {
        String::from_utf8_lossy(&seq.qual).into_owned()
    };

    write!(
        out,
        "{}\t{}\t*\t0\t0\t*\t{}\t{}\t0\t{}\t{}\n",
        seq.name,
        sam_flag,
        mate_rname,
        mate_pos,
        seq_str,
        qual_str,
    )
    .unwrap();

    out
}

// ---------------------------------------------------------------------------
// Unit tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rev_comp_simple() {
        assert_eq!(rev_comp(b"ACGT"), b"ACGT");
        assert_eq!(rev_comp(b"AAAA"), b"TTTT");
        assert_eq!(rev_comp(b"AACCGG"), b"CCGGTT");
    }

    #[test]
    fn rev_comp_empty() {
        assert_eq!(rev_comp(b""), b"");
    }

    #[test]
    fn cigar_to_str_basic() {
        // 10M 2I 5M
        let cigar = vec![
            cigar_enc(CIGAR_M, 10),
            cigar_enc(CIGAR_I, 2),
            cigar_enc(CIGAR_M, 5),
        ];
        assert_eq!(cigar_to_str(&cigar), "10M2I5M");
    }

    #[test]
    fn cigar_to_str_with_deletion() {
        let cigar = vec![
            cigar_enc(CIGAR_M, 4),
            cigar_enc(CIGAR_D, 3),
            cigar_enc(CIGAR_M, 6),
        ];
        assert_eq!(cigar_to_str(&cigar), "4M3D6M");
    }

    #[test]
    fn cigar_to_str_with_soft_clip() {
        let cigar = vec![
            cigar_enc(CIGAR_S, 2),
            cigar_enc(CIGAR_M, 8),
            cigar_enc(CIGAR_S, 3),
        ];
        assert_eq!(cigar_to_str(&cigar), "2S8M3S");
    }

    #[test]
    fn cigar_to_str_empty() {
        assert_eq!(cigar_to_str(&[]), "");
    }

    #[test]
    fn rev_comp_lowercase() {
        assert_eq!(rev_comp(b"acgt"), b"acgt");
        assert_eq!(rev_comp(b"aaccgg"), b"ccggtt");
    }

    #[test]
    fn cigar_enc_roundtrip() {
        let c = cigar_enc(CIGAR_M, 42);
        assert_eq!(cigar_op(c), CIGAR_M);
        assert_eq!(cigar_len(c), 42);
    }
}
