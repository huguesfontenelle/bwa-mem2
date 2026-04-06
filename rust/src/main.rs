/*************************************************************************************
                           The MIT License

   BWA-MEM2  (Sequence alignment using Burrows-Wheeler Transform),
   Copyright (C) 2019  Intel Corporation, Heng Li.

   Rust port — CLI entry point.

Authors: Rust port based on original C/C++ by
         Vasimuddin Md <vasimuddin.md@intel.com>; Sanchit Misra <sanchit.misra@intel.com>;
         Heng Li <hli@jimmy.harvard.edu>
*****************************************************************************************/

mod bntseq;
mod bwt;
mod fmi_search;
mod index;
mod kseq;
mod mem;
mod pe;
mod sam;
mod sw;
mod types;

use std::io::{self, BufWriter, Write};
use std::sync::{Arc, Mutex};

use anyhow::{Context, Result};
use clap::{Args, Parser, Subcommand};
use log::info;
use rayon::prelude::*;

use bntseq::{restore_pac, restore_pac_rev, BntSeq};
use bwt::FmIndex;
use index::build_index;
use kseq::{open_seq_file, KSeqReader};
use mem::mem_align1;
use pe::{mem_pair, mem_pestat};
use sam::{mem_aln2sam, mem_reg2aln, write_unmapped_sam};
use types::{fill_scmat, BSeq, MemOpt, MemPeStat, MEM_F_PE, MEM_F_SOFTCLIP};

const VERSION: &str = "2.2.1";

// ---------------------------------------------------------------------------
// CLI definition
// ---------------------------------------------------------------------------

/// BWA-MEM2: fast and accurate short-read aligner (Rust port)
#[derive(Parser)]
#[command(name = "bwa-mem2", version = VERSION, about = "BWA-MEM2 Rust port")]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand)]
enum Command {
    /// Build FM-index for a reference genome
    Index(IndexArgs),
    /// Align reads to a reference genome
    Mem(MemArgs),
    /// Print version
    Version,
}

#[derive(Args)]
struct IndexArgs {
    /// Reference FASTA file
    reference: String,
    /// Index prefix [default: same as reference]
    #[arg(short = 'p', long)]
    prefix: Option<String>,
}

#[derive(Args)]
struct MemArgs {
    /// Reference genome (index prefix)
    reference: String,
    /// Read file (FASTA or FASTQ, optionally gzip-compressed)
    reads: String,
    /// Mate read file for paired-end alignment
    mates: Option<String>,

    // Alignment options
    /// Number of threads [1]
    #[arg(short = 't', default_value_t = 1)]
    threads: i32,
    /// Minimum seed length [19]
    #[arg(short = 'k', default_value_t = 19)]
    min_seed_len: i32,
    /// Band width for banded alignment [100]
    #[arg(short = 'w', default_value_t = 100)]
    band_width: i32,
    /// Off-diagonal X-dropoff (Z-drop) [100]
    #[arg(short = 'd', default_value_t = 100)]
    zdrop: i32,
    /// Skip seeds with more than INT occurrences [500]
    #[arg(short = 'c', default_value_t = 500)]
    max_occ: i32,
    /// Match score [1]
    #[arg(short = 'A', default_value_t = 1)]
    match_score: i32,
    /// Mismatch penalty [4]
    #[arg(short = 'B', default_value_t = 4)]
    mismatch_pen: i32,
    /// Gap open penalty [6]
    #[arg(short = 'O', default_value_t = 6)]
    gap_open: i32,
    /// Gap extension penalty [1]
    #[arg(short = 'E', default_value_t = 1)]
    gap_ext: i32,
    /// Clipping penalty [5]
    #[arg(short = 'L', default_value_t = 5)]
    clip_pen: i32,
    /// Penalty for unpaired reads [17]
    #[arg(short = 'U', default_value_t = 17)]
    unpaired_pen: i32,
    /// Treat ALT contigs as part of the primary assembly
    #[arg(short = 'j', default_value_t = false)]
    no_alt: bool,
    /// Mark shorter split hits as secondary (for compatibility with BWA-SW/Bowtie2)
    #[arg(short = 'M', default_value_t = false)]
    mark_secondary: bool,
    /// Use soft clipping CIGAR operation for supplementary alignments
    #[arg(short = 'Y', default_value_t = false)]
    soft_clip_sup: bool,
    /// Output all found alignments for single-end or unpaired paired-end reads
    #[arg(short = 'a', default_value_t = false)]
    output_all: bool,
    /// Read group header line (e.g. '@RG\tID:foo\tSM:bar')
    #[arg(short = 'R')]
    rg: Option<String>,
    /// Insert SAM header lines
    #[arg(short = 'H')]
    header: Option<String>,
    /// Output file [stdout]
    #[arg(short = 'o')]
    output: Option<String>,
    /// Chunk size for processing (bases per batch) [10000000]
    #[arg(long, default_value_t = 10_000_000)]
    chunk_size: i64,
}

// ---------------------------------------------------------------------------
// SAM header
// ---------------------------------------------------------------------------

fn write_sam_header(
    bns: &BntSeq,
    rg_line: Option<&str>,
    extra_headers: Option<&str>,
    pg_cmd: &str,
    out: &mut dyn Write,
) -> Result<()> {
    writeln!(out, "@HD\tVN:1.6\tSO:unsorted")?;
    for ann in &bns.anns {
        writeln!(out, "@SQ\tSN:{}\tLN:{}", ann.name, ann.len)?;
    }
    if let Some(rg) = rg_line {
        // Normalise: replace literal \t with tab
        let rg = rg.replace("\\t", "\t");
        if rg.starts_with("@RG") {
            writeln!(out, "{}", rg)?;
        } else {
            writeln!(out, "@RG\t{}", rg)?;
        }
    }
    if let Some(h) = extra_headers {
        for line in h.lines() {
            writeln!(out, "{}", line)?;
        }
    }
    writeln!(out, "@PG\tID:bwa-mem2\tPN:bwa-mem2\tVN:{}\tCL:{}", VERSION, pg_cmd)?;
    Ok(())
}

// ---------------------------------------------------------------------------
// Index command
// ---------------------------------------------------------------------------

fn cmd_index(args: &IndexArgs) -> Result<()> {
    let prefix = args.prefix.as_deref().unwrap_or(&args.reference);
    eprintln!("[bwa-mem2] Building index for '{}' (prefix '{}')", args.reference, prefix);
    build_index(args.reference.as_str(), prefix)?;
    eprintln!("[bwa-mem2] Index built successfully.");
    Ok(())
}

// ---------------------------------------------------------------------------
// Mem command
// ---------------------------------------------------------------------------

fn cmd_mem(args: &MemArgs) -> Result<()> {
    // Build options
    let mut opt = MemOpt::default();
    opt.n_threads      = args.threads;
    opt.min_seed_len   = args.min_seed_len;
    opt.w              = args.band_width;
    opt.zdrop          = args.zdrop;
    opt.max_occ        = args.max_occ;
    opt.a              = args.match_score;
    opt.b              = args.mismatch_pen;
    opt.o_del          = args.gap_open;
    opt.o_ins          = args.gap_open;
    opt.e_del          = args.gap_ext;
    opt.e_ins          = args.gap_ext;
    opt.pen_clip5      = args.clip_pen;
    opt.pen_clip3      = args.clip_pen;
    opt.pen_unpaired   = args.unpaired_pen;
    opt.chunk_size     = args.chunk_size;
    if args.soft_clip_sup { opt.flag |= MEM_F_SOFTCLIP; }
    if args.mates.is_some() { opt.flag |= MEM_F_PE; }

    // Recompute scoring matrix
    types::fill_scmat(opt.a, opt.b, &mut opt.mat);

    // Load index
    let prefix = &args.reference;
    eprintln!("[bwa-mem2] Loading index '{}'...", prefix);
    let fmi = FmIndex::load(prefix).with_context(|| format!("loading FM-index from '{}'", prefix))?;
    let bns = BntSeq::restore(prefix).with_context(|| format!("loading BNT from '{}'", prefix))?;
    let fwd_pac = restore_pac(prefix).with_context(|| format!("loading PAC from '{}'", prefix))?;
    let pac = restore_pac_rev(&fwd_pac, bns.l_pac);
    eprintln!("[bwa-mem2] Index loaded. Reference: {} sequences, {} Mbp total.",
        bns.n_seqs, bns.l_pac / 1_000_000);

    // Set up thread pool
    rayon::ThreadPoolBuilder::new()
        .num_threads(opt.n_threads.max(1) as usize)
        .build_global()
        .ok();

    // Set up output
    let stdout = io::stdout();
    let mut out: Box<dyn Write> = match &args.output {
        Some(path) => Box::new(BufWriter::new(
            std::fs::File::create(path).with_context(|| format!("creating output '{}'", path))?
        )),
        None => Box::new(BufWriter::new(stdout.lock())),
    };

    // SAM header
    let pg_cmd = format!("bwa-mem2 mem {}", prefix);
    write_sam_header(&bns, args.rg.as_deref(), args.header.as_deref(), &pg_cmd, &mut out)?;

    // Arc-wrap immutable shared data
    let opt   = Arc::new(opt);
    let fmi   = Arc::new(fmi);
    let bns   = Arc::new(bns);
    let pac   = Arc::new(pac);
    let out   = Arc::new(Mutex::new(out));

    let is_pe = args.mates.is_some();
    let mut n_processed: i64 = 0;

    // Open read file(s)
    let mut reader1_box = open_seq_file(&args.reads)?;
    let reader1 = reader1_box.as_mut();
    let mut kseq1 = KSeqReader::new(reader1);

    let (mut reader2_opt, mut kseq2_opt) = if let Some(mates) = &args.mates {
        let mut r = open_seq_file(mates)?;
        let kseq = KSeqReader::new(unsafe {
            // SAFETY: we're extending the lifetime of the reader's reference to match
            // the outer scope. This is safe because reader2_opt owns the Box.
            &mut *(r.as_mut() as *mut dyn io::BufRead)
        });
        (Some(r), Some(kseq))
    } else {
        (None, None)
    };

    // Insert-size stats (for PE mode)
    let mut pes: [MemPeStat; 4] = Default::default();
    let mut pe_stats_computed = false;

    loop {
        // Read a batch of sequences
        let batch_bases = opt.chunk_size;
        let mut seqs1: Vec<BSeq> = Vec::new();
        let mut seqs2: Vec<BSeq> = Vec::new();
        let mut total_bases: i64 = 0;

        while total_bases < batch_bases {
            match kseq1.next()? {
                None => break,
                Some(ks) => {
                    total_bases += ks.seq.len() as i64;
                    seqs1.push(BSeq {
                        id: seqs1.len() as i32,
                        name: ks.name,
                        comment: ks.comment,
                        seq: ks.seq,
                        qual: ks.qual,
                        sam: String::new(),
                    });
                    if let Some(ref mut kseq2) = kseq2_opt {
                        match kseq2.next()? {
                            None => {}
                            Some(ks2) => seqs2.push(BSeq {
                                id: seqs2.len() as i32,
                                name: ks2.name,
                                comment: ks2.comment,
                                seq: ks2.seq,
                                qual: ks2.qual,
                                sam: String::new(),
                            }),
                        }
                    }
                }
            }
        }

        if seqs1.is_empty() { break; }

        let batch_size = seqs1.len();
        eprintln!("[bwa-mem2] Aligning batch of {} reads...", batch_size);

        // Align reads in parallel
        let opt_ref  = Arc::clone(&opt);
        let fmi_ref  = Arc::clone(&fmi);
        let bns_ref  = Arc::clone(&bns);
        let pac_ref  = Arc::clone(&pac);

        let results: Vec<Vec<_>> = seqs1.par_iter().map(|seq| {
            mem_align1(&opt_ref, &fmi_ref, &bns_ref, &pac_ref, seq)
        }).collect();

        // Output SAM
        let mut out_guard = out.lock().unwrap();

        if is_pe && !seqs2.is_empty() {
            // Estimate insert-size stats from first batch
            if !pe_stats_computed {
                let pairs: Vec<_> = results.iter().zip(seqs2.iter().map(|s| {
                    mem_align1(&opt_ref, &fmi_ref, &bns_ref, &pac_ref, s)
                })).map(|(r1, r2)| (r1.clone(), r2)).collect();
                pes = mem_pestat(&opt_ref, &bns_ref, &pairs.iter().map(|(a,b)| (a.clone(), b.clone())).collect::<Vec<_>>());
                pe_stats_computed = true;
                for (i, pe) in pes.iter().enumerate() {
                    if !pe.failed {
                        eprintln!("[bwa-mem2] Insert size distribution (mode {}): avg={:.1} std={:.1} [{}, {}]",
                            i, pe.avg, pe.std, pe.low, pe.high);
                    }
                }
            }

            // For each pair, compute alignments and output
            for (seq1, regs1) in seqs1.iter().zip(results.iter()) {
                let idx = seq1.id as usize;
                if idx >= seqs2.len() { continue; }
                let seq2 = &seqs2[idx];
                let mut regs1_local = regs1.clone();
                let mut regs2_local = mem_align1(&opt_ref, &fmi_ref, &bns_ref, &pac_ref, seq2);

                let (s1, s2) = mem_pair(
                    &opt_ref, &bns_ref, &pac_ref, &pes,
                    seq1, seq2,
                    &mut regs1_local, &mut regs2_local,
                    n_processed,
                );
                out_guard.write_all(s1.as_bytes())?;
                out_guard.write_all(s2.as_bytes())?;
            }
        } else {
            // Single-end output
            for (seq, regs) in seqs1.iter().zip(results.iter()) {
                if regs.is_empty() {
                    let s = write_unmapped_sam(seq, 0x4, None, &bns_ref);
                    out_guard.write_all(s.as_bytes())?;
                } else {
                    for reg in regs {
                        if reg.secondary >= 0 && !args.output_all { continue; }
                        let aln = mem_reg2aln(&opt_ref, &bns_ref, &pac_ref, seq.seq.len() as i32, &seq.seq, reg);
                        let s = mem_aln2sam(&opt_ref, &bns_ref, seq, &aln, &[], None);
                        out_guard.write_all(s.as_bytes())?;
                    }
                }
            }
        }

        n_processed += batch_size as i64;
        drop(out_guard);
    }

    eprintln!("[bwa-mem2] Done. Processed {} reads.", n_processed);
    Ok(())
}

// ---------------------------------------------------------------------------
// Entry point
// ---------------------------------------------------------------------------

fn main() {
    env_logger::init();
    let cli = Cli::parse();
    let result = match &cli.command {
        Command::Index(args) => cmd_index(args),
        Command::Mem(args)   => cmd_mem(args),
        Command::Version     => {
            println!("{}", VERSION);
            Ok(())
        }
    };
    if let Err(e) = result {
        eprintln!("[bwa-mem2] Error: {:#}", e);
        std::process::exit(1);
    }
}
