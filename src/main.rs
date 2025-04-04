/*
ml Clang/18.1.8-GCCcore-13.3.0
samtools view data/HLA-A_reads_ds.bam | head -n 1
samtools view data/HLA-A_reads.snc.bam | head -n 2
cargo build --release \
    && cp ~/develop/splitncigar/target/release/splitncigar ~/.local/bin/splitncigar \
    && cp ~/develop/splitncigar/target/release/slowview ~/.local/bin/slowview \
    && cp ~/develop/splitncigar/target/release/bamsummary ~/.local/bin/bamsummary \
    && cp ~/develop/splitncigar/target/release/cosmictovcf ~/.local/bin/cosmictovcf
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

hg38=/Users/sfurlan/refs/GRCh38/GRCh38.p13.genome.fa
cosmic=/Users/sfurlan/Downloads/Cosmic_MutantCensus_Tsv_v101_GRCh38/Cosmic_MutantCensus_v101_GRCh38.tsv
cosmictovcf $cosmic cosmic.vcf $hg38
head $cosmic
*/

use clap::{Arg, Command};
use rust_htslib::bam;
use rust_htslib::bam::{Read, Writer};
use rust_htslib::bam::record::Aux;
use std::error::Error;
use std::io::{stdin, stdout, Write}; // for reading user input & flushing

/// Control how we handle an existing SA tag when adding new sub-alignments.
#[derive(Debug, Clone, Copy)]
enum SaTagMode {
    /// Remove any old SA, then write only your new subalignments.
    Overwrite,
    /// Parse + merge with old SA (avoid duplicates, etc.).
    Merge,
    /// If the record already has SA, keep it. If not, add your new ones.
    Preserve,
}

fn main() -> Result<(), Box<dyn Error>> {
    let matches = Command::new("split_ncigar_reads")
        .version("0.1")
        .author("Your Name")
        .about("Splits reads with N in the CIGAR string (spliced reads).")
        .arg(
            Arg::new("input")
                .short('i')
                .long("input")
                .help("Input BAM file")
                .num_args(1)
                .required(true),
        )
        .arg(
            Arg::new("output")
                .short('o')
                .long("output")
                .help("Output BAM file")
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
                .help("Print debug info for each read and prompt for user input.")
                .action(clap::ArgAction::SetTrue),
        )
        // Option to disable flag correction.
        .arg(
            Arg::new("no_flag_correction")
                .long("no-flag-correction")
                .help("Disable flag correction (default: enabled)")
                .action(clap::ArgAction::SetTrue),
        )
        // New: Option to choose how to handle SA tags.
        .arg(
            Arg::new("sa_mode")
                .long("sa-mode")
                .help("Specify how to handle SA tags: overwrite, merge, or preserve (default: merge)")
                .num_args(1)
                .default_value("merge"),
        )
        .get_matches();

    let input_path = matches.get_one::<String>("input").unwrap();
    let output_path = matches.get_one::<String>("output").unwrap();
    // let reference_path = matches.get_one::<String>("reference").unwrap(); // Not used in code below.

    let refactor_cigar_string = matches.get_flag("refactor_cigar_string");
    let skip_mq_transform = matches.get_flag("skip_mq_transform");
    let debug_splits = matches.get_flag("debug_splits");
    // Determine whether to perform flag correction (enabled by default).
    let flag_correction = !matches.get_flag("no_flag_correction");

    // Parse the SA mode from command line.
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

    // Open the input BAM.
    let mut bam_reader = bam::Reader::from_path(input_path)?;
    let header = bam::Header::from_template(bam_reader.header());
    let headerview = bam::HeaderView::from_header(&header);
    let mut bam_writer = Writer::from_path(output_path, &header, bam::Format::Bam)?;

    // Process each read record.
    for result in bam_reader.records() {
        let mut record = result?;

        // Optionally apply the mapping quality transform (255 -> 60).
        if !skip_mq_transform && record.mapq() == 255 {
            record.set_mapq(60);
        }

        // Optionally refactor the CIGAR (placeholder).
        if refactor_cigar_string {
            // e.g. merging adjacent N-D-N segments if needed.
        }

        // Debug: print the original read.
        if debug_splits {
            eprintln!(
                "[DEBUG] Original read:\n  QNAME={}\n  pos={}\n  seq_len={}\n  CIGAR={}\n  seq={}",
                String::from_utf8_lossy(record.qname()),
                record.pos(),
                record.seq().len(),
                record.cigar().to_string(),
                decode_seq(&record)
            );
        }

        // Save the original flag from the unsplit record.
        let orig_flag = record.flags();

        // Split the read if it has N operators.
        let mut split_records = split_ncigar_read(&record)?;

        // If flag correction is enabled (default), then override the flags
        // of each split record with the original flag.
        if flag_correction {
            for rec in split_records.iter_mut() {
                rec.set_flags(orig_flag);
            }
        }

        // If multiple subreads, mark them supplementary and set SA:Z.
        if split_records.len() > 1 {
            repair_supplementary_tags(&mut split_records);

            // Obtain reference name from tid.
            let rname = tid_to_refname(&headerview, record.tid());

            // Build and add sub-alignments using the chosen SA mode.
            annotate_sa_tag(&mut split_records, &rname, sa_mode)?;
        }

        // Debug: print splitted records if more than one.
        if debug_splits && split_records.len() > 1 {
            eprintln!("    -> splitted into {} records:", split_records.len());
            for (i, rec) in split_records.iter().enumerate() {
                eprintln!(
                    "       split#{}:\n         QNAME={}\n         pos={}\n         seq_len={}\n         CIGAR={}\n         seq={}",
                    i,
                    String::from_utf8_lossy(rec.qname()),
                    rec.pos(),
                    rec.seq().len(),
                    rec.cigar().to_string(),
                    decode_seq(rec)
                );
            }
        }

        // Write out each subrecord.
        for rec in split_records {
            bam_writer.write(&rec)?;
        }

        // If debugging, prompt user for 'q' or Enter.
        if debug_splits {
            let should_quit = pause_for_debug()?;
            if should_quit {
                eprintln!("User quit debug mode. Stopping early...");
                break;
            }
        }
    }

    Ok(())
}

/// Prompt user to press Enter to continue, or 'q' to quit. 
fn pause_for_debug() -> Result<bool, Box<dyn Error>> {
    eprint!("Press Enter to continue, or 'q' then Enter to quit: ");
    stdout().flush()?; 
    let mut input = String::new();
    stdin().read_line(&mut input)?;
    Ok(input.trim().eq_ignore_ascii_case("q"))
}

/// If tid < 0 or out-of-range, return "*", else return the actual ref name.
fn tid_to_refname(header_view: &bam::HeaderView, tid: i32) -> String {
    if tid < 0 {
        return "*".to_string();
    }
    let idx = tid as u32;
    if idx >= header_view.target_count() {
        return "*".to_string();
    }
    let bytes: &[u8] = header_view.tid2name(idx);
    String::from_utf8_lossy(bytes).to_string()
}

/// Convert a read's 4-bit or ASCII-coded bases into a String.
fn decode_seq(record: &bam::Record) -> String {
    String::from_utf8_lossy(&record.seq().as_bytes()).to_string()
}

/// If the read's CIGAR has 'N', create multiple subrecords. 
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
    // last chunk
    if section_has_match && start_idx < cigar_elems.len() {
        let split = split_read_softclip(record, start_idx, cigar_elems.len())?;
        splits.push(split);
    }
    Ok(splits)
}

/// Soft-clip outside the subrange [cigar_start_index..cigar_end_index].
fn split_read_softclip(
    record: &bam::Record,
    cigar_start: usize,
    cigar_end: usize,
) -> Result<bam::Record, Box<dyn Error>> {
    let binding = record.cigar();
    let cigar_elems: Vec<_> = binding.into_iter().collect();

    // skip leading/trailing 'D' in subrange
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

    // alignment start
    let mut new_start = record.pos();
    for c in &cigar_elems[0..start] {
        if consumes_reference(c.char()) {
            new_start += c.len() as i64;
        }
    }

    // collect middle portion
    let sub_ops: Vec<bam::record::Cigar> = cigar_elems[start..end]
        .iter()
        .map(|&c| c.clone())
        .collect();

    // figure out how many bases in prefix, subrange
    let mut prefix_qbases = 0;
    for c in &cigar_elems[0..start] {
        if consumes_query(c.char()) {
            prefix_qbases += c.len() as usize;
        }
    }
    let mut middle_qbases = 0;
    for c in &cigar_elems[start..end] {
        if consumes_query(c.char()) {
            middle_qbases += c.len() as usize;
        }
    }
    let total_len = record.seq().len();
    let leading_s = prefix_qbases;
    let trailing_s = total_len.saturating_sub(leading_s + middle_qbases);

    // new CIGAR
    let mut new_ops = Vec::new();
    if leading_s > 0 {
        new_ops.push(bam::record::Cigar::SoftClip(leading_s as u32));
    }
    new_ops.extend(sub_ops);
    if trailing_s > 0 {
        new_ops.push(bam::record::Cigar::SoftClip(trailing_s as u32));
    }

    // build a new record
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

/// returns true if the operator consumes reference
fn consumes_reference(op: char) -> bool {
    matches!(op, 'M' | 'D' | 'N' | '=' | 'X')
}

/// returns true if the operator consumes query
fn consumes_query(op: char) -> bool {
    matches!(op, 'M' | 'I' | 'S' | '=' | 'X')
}

/// remove NM, MD, NH and set 0x800 on non-primary
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
            rec.set_flags(f | 0x800); // supplementary flag
        }
    }
}

/// Decide how to handle the existing SA tag.
fn annotate_sa_tag(
    records: &mut [bam::Record],
    rname: &str,
    sa_mode: SaTagMode,
) -> Result<(), Box<dyn Error>> {
    if records.len() < 2 {
        return Ok(());
    }

    // Build new sub-align lines for each record.
    let sub_aligns: Vec<String> = records
        .iter()
        .map(|rec| {
            let pos1 = rec.pos() + 1;
            let strand = if (rec.flags() & 0x10) != 0 { '-' } else { '+' };
            let cigar_str = rec.cigar().to_string();
            let mapq = rec.mapq();
            let nm = 0; // or parse rec.aux(b"NM")
            format!("{},{},{},{},{},{}", rname, pos1, strand, cigar_str, mapq, nm)
        })
        .collect();

    // For each record, gather sub-align lines of the *other* split records.
    for (i, rec) in records.iter_mut().enumerate() {
        let mut new_entries = Vec::new();
        for (j, sa) in sub_aligns.iter().enumerate() {
            if j != i {
                new_entries.push(sa.clone());
            }
        }
        handle_sa_tag(rec, &new_entries, sa_mode)?;
    }

    Ok(())
}

/// Actually handle Overwrite, Merge, or Preserve existing SA:Z.
fn handle_sa_tag(
    rec: &mut bam::Record,
    new_subaligns: &[String],
    mode: SaTagMode,
) -> Result<(), Box<dyn Error>> {
    // fetch existing SA if any.
    let old_sa = match rec.aux(b"SA") {
        Ok(Aux::String(s)) => s.to_string(),
        _ => String::new(),
    };

    match mode {
        SaTagMode::Overwrite => {
            // Remove old, then write new.
            let _ = rec.remove_aux(b"SA");
            if !new_subaligns.is_empty() {
                let mut merged = new_subaligns.join(";");
                merged.push(';');
                rec.push_aux(b"SA", Aux::String(&merged))?;
            }
        }
        SaTagMode::Preserve => {
            // If old SA is present, do nothing; else add new.
            if old_sa.is_empty() && !new_subaligns.is_empty() {
                let mut new_sa = new_subaligns.join(";");
                new_sa.push(';');
                rec.push_aux(b"SA", Aux::String(&new_sa))?;
            }
        }
        SaTagMode::Merge => {
            // Parse old; add new if not present; rewrite.
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
