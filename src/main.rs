/*
ml Clang/18.1.8-GCCcore-13.3.0
samtools view data/HLA-A_reads_ds.bam | head -n 1
samtools view data/HLA-A_reads.snc.bam | head -n 2
cargo build --release && cp ~/develop/splitncigar/target/release/splitncigar ~/.local/bin/splitncigar
hg38=/Users/sfurlan/refs/GRCh38/GRCh38.p13.genome.fa
splitncigar --input data/HLA-A_reads_ds.bam --output HLA-A_reads.snc.bam --reference /fh/fast/furlan_s/grp/refs/GRCh38/GRCh38.p13.genome.fa
samtools view HLA-A_reads.snc.bam | head -n 1
rm HLA-A_reads.snc.bam
*/


use clap::{Arg, Command};
use rust_htslib::bam;
use rust_htslib::bam::{Read, Record, Writer};
use std::error::Error;
use std::io::{stdin, stdout, Write}; // for reading user input & flushing

fn main() -> Result<(), Box<dyn Error>> {
    let matches = Command::new("split_ncigar_reads")
        .version("0.1")
        .author("Your Name")
        .about("Splits reads with N in the CIGAR string (e.g. spanning splicing events).")
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
                .short('R')
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
                .help("Print debug info for each read and then prompt user to press Enter or 'q'.")
                .action(clap::ArgAction::SetTrue),
        )
        .get_matches();

    let input_path = matches.get_one::<String>("input").unwrap();
    let output_path = matches.get_one::<String>("output").unwrap();
    // let reference_path = matches.get_one::<String>("reference").unwrap(); // not used below

    let refactor_cigar_string = matches.get_flag("refactor_cigar_string");
    let skip_mq_transform = matches.get_flag("skip_mq_transform");
    let debug_splits = matches.get_flag("debug_splits");

    // Open the input BAM and create the writer using its header.
    let mut bam_reader = bam::Reader::from_path(input_path)?;
    let header = bam::Header::from_template(bam_reader.header());
    let mut bam_writer = Writer::from_path(output_path, &header, bam::Format::Bam)?;

    // Process each read record.
    for result in bam_reader.records() {
        let mut record = result?;

        // Optionally apply a mapping quality transformation:
        if !skip_mq_transform && record.mapq() == 255 {
            record.set_mapq(60);
        }

        // Optionally refactor the CIGAR string (stub).
        if refactor_cigar_string {
            // In a complete implementation, you might merge adjacent N-D-N segments, etc.
        }

        // Debug: print the original read before splitting.
        if debug_splits {
            eprintln!(
                "[DEBUG] Original read: QNAME={}, pos={}, seq_len={}, CIGAR={}",
                String::from_utf8_lossy(record.qname()),
                record.pos(),
                record.seq().len(),
                record.cigar().to_string()
            );
        }

        // Split the read if it contains N operators.
        let mut split_records = split_ncigar_read(&record)?;

        // If the read was split into multiple records, update auxiliary tags.
        if split_records.len() > 1 {
            repair_supplementary_tags(&mut split_records);
        }

        // Debug: print info about the splitted records if more than one.
        if debug_splits && split_records.len() > 1 {
            eprintln!("    -> splitted into {} records:", split_records.len());
            for (i, rec) in split_records.iter().enumerate() {
                eprintln!(
                    "       split#{}: QNAME={}, pos={}, seq_len={}, CIGAR={}",
                    i,
                    String::from_utf8_lossy(rec.qname()),
                    rec.pos(),
                    rec.seq().len(),
                    rec.cigar().to_string()
                );
            }
        }

        // Write out each resulting record.
        for rec in split_records {
            bam_writer.write(&rec)?;
        }

        // If debugging, prompt user to continue or quit after this read’s debug info.
        if debug_splits {
            let should_quit = pause_for_debug()?;
            if should_quit {
                eprintln!("User quit debug mode. Stopping early...");
                break; // stops reading/writing the rest of the file
            }
        }
    }

    Ok(())
}

/// Prompt user to press Enter to continue or 'q' to quit.
/// Returns `true` if user typed 'q', otherwise `false`.
fn pause_for_debug() -> Result<bool, Box<dyn Error>> {
    eprint!("Press Enter to continue, or 'q' then Enter to quit: ");
    stdout().flush()?; // ensure the prompt is displayed immediately

    let mut input = String::new();
    stdin().read_line(&mut input)?;
    Ok(input.trim().eq_ignore_ascii_case("q"))
}

/// Splits a read into one or more records if its CIGAR string contains N operators.
/// If no N is present, returns a vector with a single clone of the original record.
fn split_ncigar_read(record: &Record) -> Result<Vec<Record>, Box<dyn Error>> {
    let cigar = record.cigar();
    let cigar_elements: Vec<_> = cigar.iter().collect();
    let mut has_n = false;
    for cig in &cigar_elements {
        if cig.char() == 'N' {
            has_n = true;
            break;
        }
    }
    if !has_n {
        return Ok(vec![record.clone()]);
    }

    let mut split_records = Vec::new();
    let mut section_has_match = false;
    let mut first_cigar_index = 0;
    let num_cigar_elements = cigar_elements.len();

    for (i, cig) in cigar_elements.iter().enumerate() {
        let op = cig.char();
        // Consider these operators as “real” bases.
        if op == 'M' || op == '=' || op == 'X' || op == 'I' || op == 'D' {
            section_has_match = true;
        }
        if op == 'N' {
            if section_has_match {
                let split = split_read_based_on_cigar(record, first_cigar_index, i)?;
                split_records.push(split);
            }
            first_cigar_index = i + 1;
            section_has_match = false;
        }
    }
    // Process the last section (if any).
    if first_cigar_index < num_cigar_elements && section_has_match {
        let split = split_read_based_on_cigar(record, first_cigar_index, num_cigar_elements)?;
        split_records.push(split);
    }
    Ok(split_records)
}

/// Creates a new record by “clipping” the input record to the region defined by
/// the CIGAR elements between `cigar_start_index` and `cigar_end_index`.
/// This function adjusts the CIGAR, the alignment start position,
/// and also the read bases + quality arrays.
fn split_read_based_on_cigar(
    record: &Record,
    cigar_start_index: usize,
    cigar_end_index: usize,
) -> Result<Record, Box<dyn Error>> {
    let original_cigar = record.cigar();
    let cigar_elements: Vec<_> = original_cigar.iter().collect();

    // 1) Adjust start/end so that we skip leading or trailing D in this subset.
    let mut start = cigar_start_index;
    let mut end = cigar_end_index;
    while start < cigar_elements.len() && cigar_elements[start].char() == 'D' {
        start += 1;
    }
    while end > start && cigar_elements[end - 1].char() == 'D' {
        end -= 1;
    }
    if start >= end {
        return Err(format!(
            "Empty section between Ns in CIGAR: {}",
            record.cigar().to_string()
        )
        .into());
    }

    // 2) Compute the new alignment start.
    // In BAM, `pos()` is 0-based. For each reference-consuming operator
    // before `start`, add up the operator lengths to shift the alignment.
    let mut new_start = record.pos();
    let prefix = &cigar_elements[0..start];
    for cig in prefix {
        if consumes_reference(cig.char()) {
            new_start += cig.len() as i64;
        }
    }

    // 3) Build the sub-CIGAR by taking the slice [start..end].
    let new_cigar: Vec<bam::record::Cigar> = cigar_elements[start..end]
        .iter()
        .map(|&c| c.clone())
        .collect();

    // 4) Figure out how many query bases we skip in the prefix (so we can trim the sequence).
    //    We also need to know how many query bases are in this subrange,
    //    so we can slice out exactly that part of the read.
    let mut skipped_query_bases = 0;
    for cig in &cigar_elements[0..start] {
        if consumes_query(cig.char()) {
            skipped_query_bases += cig.len() as usize;
        }
    }

    let mut subrange_query_bases = 0;
    for cig in &cigar_elements[start..end] {
        if consumes_query(cig.char()) {
            subrange_query_bases += cig.len() as usize;
        }
    }

    // 5) Slice the read’s sequence and quality accordingly.
    let full_seq = record.seq().as_bytes();
    let full_qual = record.qual();
    if skipped_query_bases + subrange_query_bases > full_seq.len() {
        return Err(format!(
            "Mismatch splitting read: skipping {} + subrange {} > seq len {}",
            skipped_query_bases,
            subrange_query_bases,
            full_seq.len()
        )
        .into());
    }

    let new_seq = &full_seq[skipped_query_bases..skipped_query_bases + subrange_query_bases];
    let new_qual = &full_qual[skipped_query_bases..skipped_query_bases + subrange_query_bases];

    // 6) Create the new record from the old record. Then set the new CIGAR, new pos, and new sequence.
    let mut new_record = record.clone();
    new_record.set(
        record.qname(),
        Some(&bam::record::CigarString::from(new_cigar)),
        new_seq,
        new_qual,
    );
    new_record.set_pos(new_start);

    Ok(new_record)
}

/// Returns true if the operator consumes reference bases.
fn consumes_reference(op: char) -> bool {
    matches!(op, 'M' | 'D' | 'N' | '=' | 'X')
}

/// Returns true if the operator consumes query bases.
/// Include `S` so soft-clipped bases are counted in the read length.
fn consumes_query(op: char) -> bool {
    matches!(op, 'M' | 'I' | 'S' | '=' | 'X')
}

/// Clears certain auxiliary tags and marks non-primary split reads as supplementary.
fn repair_supplementary_tags(records: &mut Vec<Record>) {
    // List of tags to remove.
    let tags_to_remove = [b"NM", b"MD", b"NH"];
    for record in records.iter_mut() {
        for tag in &tags_to_remove {
            let _ = record.remove_aux(*tag);
        }
    }
    // If there are multiple records, designate the first as primary and mark the rest as supplementary.
    if records.len() > 1 {
        for record in records.iter_mut().skip(1) {
            let flags = record.flags();
            record.set_flags(flags | 0x800); // 0x800 is the supplementary flag
        }
    }
}
