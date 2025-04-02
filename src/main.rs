/*
ml Clang/18.1.8-GCCcore-13.3.0
samtools view data/HLA-A_reads_ds.bam | head -n 1
samtools view data/HLA-A_reads.snc.bam | head -n 2
cargo build --release && cp ~/develop/splitncigar/target/release/splitncigar ~/.local/bin/splitncigar && cp ~/develop/splitncigar/target/release/slowview ~/.local/bin/slowview && cp ~/develop/splitncigar/target/release/bamsummary ~/.local/bin/bamsummary
hg38=/Users/sfurlan/refs/GRCh38/GRCh38.p13.genome.fa
splitncigar --input data/HLA-A_reads_ds.bam \
            --output HLA-A_reads.snc.bam \
            --reference /fh/fast/furlan_s/grp/refs/GRCh38/GRCh38.p13.genome.fa
splitncigar --debug-splits \
            --input data/HLA-A_reads_ds.bam \
            --output HLA-A_reads.snc.bam \
            --reference /fh/fast/furlan_s/grp/refs/GRCh38/GRCh38.p13.genome.fa
samtools view HLA-A_reads.snc.bam | head -n 2
samtools sort HLA-A_reads.snc.bam -o HLA-A_reads.snc.sorted.bam
samtools index HLA-A_reads.snc.sorted.bam
slowview HLA-A_reads.snc.sorted.bam
bamsummary HLA-A_reads.snc.sorted.bam
bamsummary data/HLA-A_reads.snc.snc_fc.bam
rm HLA-A_reads*
*/


use clap::{Arg, Command};
use rayon::prelude::*;
use rust_htslib::bam;
use rust_htslib::bam::{Read, Writer};
use rust_htslib::bam::record::Aux;
use std::error::Error;
use std::io::{stdin, stdout, Write};

/// Control how we handle an existing SA tag when adding new sub-alignments.
#[derive(Debug, Clone, Copy)]
enum SaTagMode {
    Overwrite,
    Merge,
    Preserve,
}

fn main() -> Result<(), Box<dyn Error>> {
    let matches = Command::new("split_ncigar_reads")
        .version("0.1")
        .author("Your Name")
        .about("Splits reads with N in the CIGAR string (spliced reads) by chromosome.")
        .arg(
            Arg::new("input")
                .short('i')
                .long("input")
                .help("Input BAM file (sorted and indexed)")
                .num_args(1)
                .required(true),
        )
        .arg(
            Arg::new("output")
                .short('o')
                .long("output")
                .help("Output (unsorted) BAM file")
                .num_args(1)
                .required(true),
        )
        .arg(
            Arg::new("reference")
                .short('r')
                .long("reference")
                .help("Reference FASTA file")
                .num_args(1)
                .required(true),
        )
        .arg(
            Arg::new("refactor_cigar_string")
                .long("refactor-cigar-string")
                .help("Refactor CIGAR strings with NDN elements")
                .action(clap::ArgAction::SetTrue),
        )
        .arg(
            Arg::new("skip_mq_transform")
                .long("skip-mq-transform")
                .help("Skip the mapping quality 255 -> 60 transformation")
                .action(clap::ArgAction::SetTrue),
        )
        .arg(
            Arg::new("debug_splits")
                .long("debug-splits")
                .help("Print debug info and process sequentially")
                .action(clap::ArgAction::SetTrue),
        )
        .arg(
            Arg::new("no_flag_correction")
                .long("no-flag-correction")
                .help("Disable flag correction (default: enabled)")
                .action(clap::ArgAction::SetTrue),
        )
        .arg(
            Arg::new("sa_mode")
                .long("sa-mode")
                .help("How to handle SA tags: overwrite, merge, or preserve (default: merge)")
                .num_args(1)
                .default_value("merge"),
        )
        .arg(
            Arg::new("threads")
                .long("threads")
                .help("Number of cores to use (default: 1)")
                .num_args(1)
                .default_value("1"),
        )
        .get_matches();

    let input_path = matches.get_one::<String>("input").unwrap();
    let output_path = matches.get_one::<String>("output").unwrap();
    let _reference = matches.get_one::<String>("reference").unwrap();

    let refactor_cigar_string = matches.get_flag("refactor_cigar_string");
    let skip_mq_transform = matches.get_flag("skip_mq_transform");
    let debug_splits = matches.get_flag("debug_splits");
    let flag_correction = !matches.get_flag("no_flag_correction");

    // Parse SA mode.
    let sa_mode_str = matches.get_one::<String>("sa_mode").unwrap().to_lowercase();
    let sa_mode = match sa_mode_str.as_str() {
        "overwrite" => SaTagMode::Overwrite,
        "merge" => SaTagMode::Merge,
        "preserve" => SaTagMode::Preserve,
        other => {
            eprintln!("Unknown SA mode '{}'; using default 'merge'", other);
            SaTagMode::Merge
        }
    };

    // Parse number of threads.
    let threads: usize = matches.get_one::<String>("threads").unwrap().parse()?;

    // Open a BAM reader to retrieve header and target names.
    let mut bam_reader = bam::Reader::from_path(input_path)?;
    let header = bam::Header::from_template(bam_reader.header());
    let headerview = bam::HeaderView::from_header(&header);
    let target_names: Vec<String> = headerview
        .target_names()
        .iter()
        .map(|n| String::from_utf8_lossy(n).into_owned())
        .collect();

    // Process each chromosome. If debug mode is enabled or only one thread is used, process sequentially.
    let processed: Vec<Vec<bam::Record>> = if debug_splits || threads == 1 {
        target_names
            .iter()
            .map(|chrom| {
                process_chromosome(
                    input_path,
                    chrom,
                    &headerview,
                    flag_correction,
                    refactor_cigar_string,
                    skip_mq_transform,
                    debug_splits,
                    sa_mode,
                )
            })
            .collect::<Result<Vec<_>, Box<dyn Error>>>()?
    } else {
        // Build a thread pool with the requested number of threads.
        let pool = rayon::ThreadPoolBuilder::new().num_threads(threads).build()?;
        pool.install(|| {
            target_names
                .par_iter()
                .map(|chrom| {
                    process_chromosome(
                        input_path,
                        chrom,
                        &headerview,
                        flag_correction,
                        refactor_cigar_string,
                        skip_mq_transform,
                        debug_splits,
                        sa_mode,
                    )
                })
                .collect::<Result<Vec<_>, Box<dyn Error>>>()
        })?
    };

    // Open a writer for the single output (unsorted) BAM file.
    let mut bam_writer = Writer::from_path(output_path, &header, bam::Format::Bam)?;
    // Write all processed records (order is not guaranteed).
    for record_vec in processed {
        for rec in record_vec {
            bam_writer.write(&rec)?;
        }
    }

    println!("Output written to {} (unsorted)", output_path);
    println!("You can sort the BAM file afterwards (e.g., using samtools sort).");
    Ok(())
}

/// Process all records for a given chromosome using an IndexedReader.
/// Returns a vector of processed (split) records.
fn process_chromosome(
    input_path: &str,
    chrom: &str,
    headerview: &bam::HeaderView,
    flag_correction: bool,
    refactor_cigar_string: bool,
    skip_mq_transform: bool,
    debug_splits: bool,
    sa_mode: SaTagMode,
) -> Result<Vec<bam::Record>, Box<dyn Error>> {
    // Open an IndexedReader (the BAM must be sorted and indexed).
    let mut idx_reader = bam::IndexedReader::from_path(input_path)?;
    // Fetch all records for the given chromosome.
    idx_reader.fetch(chrom, 0, None)?;

    let mut processed_records = Vec::new();
    for result in idx_reader.records() {
        let mut record = result?;
        let processed = process_record(
            &mut record,
            flag_correction,
            refactor_cigar_string,
            skip_mq_transform,
            debug_splits,
            headerview,
            sa_mode,
        )?;
        processed_records.extend(processed);
    }
    Ok(processed_records)
}

/// Process a single record: transform MAPQ, optionally refactor CIGAR,
/// split the read, correct flags, and annotate SA tags.
fn process_record(
    record: &mut bam::Record,
    flag_correction: bool,
    _refactor_cigar_string: bool, // placeholder for future CIGAR refactoring
    skip_mq_transform: bool,
    debug_splits: bool,
    headerview: &bam::HeaderView,
    sa_mode: SaTagMode,
) -> Result<Vec<bam::Record>, Box<dyn Error>> {
    if !skip_mq_transform && record.mapq() == 255 {
        record.set_mapq(60);
    }
    if _refactor_cigar_string {
        // Insert any desired CIGAR refactoring logic here.
    }
    if debug_splits {
        eprintln!(
            "[DEBUG] Original read:\n  QNAME={}\n  pos={}\n  seq_len={}\n  CIGAR={}\n  seq={}",
            String::from_utf8_lossy(record.qname()),
            record.pos(),
            record.seq().len(),
            record.cigar().to_string(),
            decode_seq(record)
        );
    }
    let orig_flag = record.flags();
    let mut split_records = split_ncigar_read(record)?;
    if flag_correction {
        for rec in split_records.iter_mut() {
            rec.set_flags(orig_flag);
        }
    }
    if split_records.len() > 1 {
        repair_supplementary_tags(&mut split_records);
        let rname = tid_to_refname(headerview, record.tid());
        annotate_sa_tag(&mut split_records, &rname, sa_mode)?;
    }
    if debug_splits && split_records.len() > 1 {
        eprintln!("    -> split into {} records", split_records.len());
    }
    Ok(split_records)
}

/// Return the reference name for a given tid.
fn tid_to_refname(header_view: &bam::HeaderView, tid: i32) -> String {
    if tid < 0 {
        return "*".to_string();
    }
    let idx = tid as u32;
    if idx >= header_view.target_count() {
        return "*".to_string();
    }
    let bytes = header_view.tid2name(idx);
    String::from_utf8_lossy(bytes).to_string()
}

/// Convert a read's sequence into a String.
fn decode_seq(record: &bam::Record) -> String {
    String::from_utf8_lossy(&record.seq().as_bytes()).to_string()
}

/// Splits a read with 'N' operators in its CIGAR string into subrecords.
fn split_ncigar_read(record: &bam::Record) -> Result<Vec<bam::Record>, Box<dyn Error>> {
    let cigar = record.cigar();
    let cigar_elems: Vec<_> = cigar.iter().collect();
    let mut has_n = false;
    for c in &cigar_elems {
        if c.char() == 'N' {
            has_n = true;
            break;
        }
    }
    if !has_n {
        return Ok(vec![record.clone()]);
    }
    let mut splits = Vec::new();
    let mut section_has_match = false;
    let mut start_idx = 0;
    for (i, c) in cigar_elems.iter().enumerate() {
        let op = c.char();
        if matches!(op, 'M' | '=' | 'X' | 'I' | 'D') {
            section_has_match = true;
        }
        if op == 'N' {
            if section_has_match {
                let split = split_read_softclip(record, start_idx, i)?;
                splits.push(split);
            }
            start_idx = i + 1;
            section_has_match = false;
        }
    }
    if section_has_match && start_idx < cigar_elems.len() {
        let split = split_read_softclip(record, start_idx, cigar_elems.len())?;
        splits.push(split);
    }
    Ok(splits)
}

/// Creates a subrecord by soft-clipping outside the specified CIGAR subrange.
fn split_read_softclip(
    record: &bam::Record,
    cigar_start: usize,
    cigar_end: usize,
) -> Result<bam::Record, Box<dyn Error>> {
    let binding = record.cigar();
    let cigar_elems: Vec<_> = binding.into_iter().collect();
    let mut start = cigar_start;
    let mut end = cigar_end;
    while start < end && cigar_elems[start].char() == 'D' {
        start += 1;
    }
    while end > start && cigar_elems[end - 1].char() == 'D' {
        end -= 1;
    }
    if start >= end {
        return Err(format!(
            "Empty subrange between Ns: {}",
            record.cigar().to_string()
        )
        .into());
    }
    let mut new_start = record.pos();
    for c in &cigar_elems[0..start] {
        if consumes_reference(c.char()) {
            new_start += c.len() as i64;
        }
    }
    let sub_ops: Vec<bam::record::Cigar> = cigar_elems[start..end]
        .iter()
        .map(|&c| c.clone())
        .collect();
    let prefix_qbases: usize = cigar_elems[0..start]
        .iter()
        .filter(|c| consumes_query(c.char()))
        .map(|c| c.len() as usize)
        .sum();
    let middle_qbases: usize = cigar_elems[start..end]
        .iter()
        .filter(|c| consumes_query(c.char()))
        .map(|c| c.len() as usize)
        .sum();
    let total_len = record.seq().len();
    let leading_s = prefix_qbases;
    let trailing_s = total_len.saturating_sub(leading_s + middle_qbases);
    let mut new_ops = Vec::new();
    if leading_s > 0 {
        new_ops.push(bam::record::Cigar::SoftClip(leading_s as u32));
    }
    new_ops.extend(sub_ops);
    if trailing_s > 0 {
        new_ops.push(bam::record::Cigar::SoftClip(trailing_s as u32));
    }
    let mut new_rec = record.clone();
    new_rec.set_pos(new_start);
    new_rec.set(
        record.qname(),
        Some(&bam::record::CigarString::from(new_ops)),
        &record.seq().as_bytes(),
        record.qual(),
    );
    Ok(new_rec)
}

/// Returns true if the CIGAR operator consumes reference.
fn consumes_reference(op: char) -> bool {
    matches!(op, 'M' | 'D' | 'N' | '=' | 'X')
}

/// Returns true if the CIGAR operator consumes query.
fn consumes_query(op: char) -> bool {
    matches!(op, 'M' | 'I' | 'S' | '=' | 'X')
}

/// Remove NM, MD, NH tags and mark non-primary (supplementary) records.
fn repair_supplementary_tags(records: &mut [bam::Record]) {
    let remove_tags = [b"NM", b"MD", b"NH"];
    for rec in records.iter_mut() {
        for tag in &remove_tags {
            let _ = rec.remove_aux(*tag);
        }
    }
    if records.len() > 1 {
        for rec in records.iter_mut().skip(1) {
            let f = rec.flags();
            rec.set_flags(f | 0x800);
        }
    }
}

/// Annotate records with an SA tag containing sub-alignments.
fn annotate_sa_tag(
    records: &mut [bam::Record],
    rname: &str,
    sa_mode: SaTagMode,
) -> Result<(), Box<dyn Error>> {
    if records.len() < 2 {
        return Ok(());
    }
    let sub_aligns: Vec<String> = records
        .iter()
        .map(|rec| {
            let pos1 = rec.pos() + 1;
            let strand = if (rec.flags() & 0x10) != 0 { '-' } else { '+' };
            let cigar_str = rec.cigar().to_string();
            let mapq = rec.mapq();
            let nm = 0;
            format!("{},{},{},{},{},{}", rname, pos1, strand, cigar_str, mapq, nm)
        })
        .collect();
    for (i, rec) in records.iter_mut().enumerate() {
        let new_entries: Vec<String> = sub_aligns
            .iter()
            .enumerate()
            .filter(|(j, _)| *j != i)
            .map(|(_, sa)| sa.clone())
            .collect();
        handle_sa_tag(rec, &new_entries, sa_mode)?;
    }
    Ok(())
}

/// Handle the SA tag according to the specified mode.
fn handle_sa_tag(
    rec: &mut bam::Record,
    new_subaligns: &[String],
    mode: SaTagMode,
) -> Result<(), Box<dyn Error>> {
    let old_sa = match rec.aux(b"SA") {
        Ok(Aux::String(s)) => s.to_string(),
        _ => String::new(),
    };
    match mode {
        SaTagMode::Overwrite => {
            let _ = rec.remove_aux(b"SA");
            if !new_subaligns.is_empty() {
                let mut merged = new_subaligns.join(";");
                merged.push(';');
                rec.push_aux(b"SA", Aux::String(&merged))?;
            }
        }
        SaTagMode::Preserve => {
            if old_sa.is_empty() && !new_subaligns.is_empty() {
                let mut new_sa = new_subaligns.join(";");
                new_sa.push(';');
                rec.push_aux(b"SA", Aux::String(&new_sa))?;
            }
        }
        SaTagMode::Merge => {
            let mut old_list: Vec<String> = old_sa
                .split(';')
                .filter(|s| !s.is_empty())
                .map(|s| s.to_string())
                .collect();
            for ent in new_subaligns {
                if !old_list.contains(ent) {
                    old_list.push(ent.clone());
                }
            }
            let _ = rec.remove_aux(b"SA");
            if !old_list.is_empty() {
                let mut merged = old_list.join(";");
                merged.push(';');
                rec.push_aux(b"SA", Aux::String(&merged))?;
            }
        }
    }
    Ok(())
}
