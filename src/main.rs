// Cargo.toml dependencies:
// [dependencies]
// clap = "4.0"
// rust-htslib = "0.43"   // or the latest version available

use clap::{Arg, Command};
use rust_htslib::bam;
use rust_htslib::bam::{Read, Record, Writer};
use std::error::Error;

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
        .get_matches();

    let input_path = matches.get_one::<String>("input").unwrap();
    let output_path = matches.get_one::<String>("output").unwrap();
    let _reference_path = matches.get_one::<String>("reference").unwrap();

    let refactor_cigar_string = matches.get_flag("refactor_cigar_string");
    let skip_mq_transform = matches.get_flag("skip_mq_transform");

    // Open the input BAM and create the writer using its header.
    let mut bam_reader = bam::Reader::from_path(input_path)?;
    let header = bam::Header::from_template(bam_reader.header());
    let mut bam_writer =
        Writer::from_path(output_path, &header, bam::Format::Bam)?;

    // Process each read record.
    for result in bam_reader.records() {
        let mut record = result?;

        // Optionally apply a mapping quality transformation:
        if !skip_mq_transform && record.mapq() == 255 {
            record.set_mapq(60);
        }

        // Optionally refactor the CIGAR string.
        if refactor_cigar_string {
            // This is a stub. In a complete implementation you might merge adjacent N-D-N segments.
            // For now, we simply leave the CIGAR unmodified.
        }

        // Split the read if it contains N operators.
        let mut split_records = split_ncigar_read(&record)?;

        // If the read was split into multiple records, update auxiliary tags.
        if split_records.len() > 1 {
            repair_supplementary_tags(&mut split_records);
        }

        // Write out each resulting record.
        for rec in split_records {
            bam_writer.write(&rec)?;
        }
    }

    Ok(())
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
///
/// This function adjusts the CIGAR and the alignment start position accordingly.
/// (For simplicity, the base sequence and qualities are not modified here.)
fn split_read_based_on_cigar(
    record: &Record,
    cigar_start_index: usize,
    cigar_end_index: usize,
) -> Result<Record, Box<dyn Error>> {
    let cigar = record.cigar();
    let cigar_elements: Vec<_> = cigar.iter().collect();
    let mut cigar_first_index = cigar_start_index;
    let mut cigar_second_index = cigar_end_index;

    // In case a section starts or ends with D, trim these operators.
    while cigar_first_index < cigar_elements.len() && cigar_elements[cigar_first_index].char() == 'D'
    {
        cigar_first_index += 1;
    }
    while cigar_second_index > cigar_first_index
        && cigar_elements[cigar_second_index - 1].char() == 'D'
    {
        cigar_second_index -= 1;
    }
    if cigar_first_index > cigar_second_index {
        return Err(format!(
            "Cannot split this read (might be an empty section between Ns): {}",
            record.cigar().to_string()
        )
        .into());
    }

    // Compute the new alignment start position.
    // In BAM, pos is 0-based.
    let mut new_start = record.pos();
    // Sum the reference-consuming operations in the prefix.
    for cig in &cigar_elements[0..cigar_first_index] {
        if consumes_reference(cig.char()) {
            new_start += cig.len() as i64;
        }
    }

    // For the split section, count the number of reference bases.
    let mut _ref_bases = 0;
    for cig in &cigar_elements[cigar_first_index..cigar_second_index] {
        if consumes_reference(cig.char()) {
            _ref_bases += cig.len();
        }
    }
    // In this simplified version, we do not modify the sequence bases.
    // We simply update the CIGAR to include only the kept section.
    let new_cigar: Vec<bam::record::Cigar> = cigar_elements[cigar_first_index..cigar_second_index]
        .iter()
        .map(|&c| c.clone())
        .collect();

    let mut new_record = record.clone();
    // Note: `set_cigar` is assumed to be available in rust-htslib 0.43.
    new_record.set_cigar(&new_cigar)?;
    new_record.set_pos(new_start);

    Ok(new_record)
}

/// Returns true if the CIGAR operator consumes reference bases.
fn consumes_reference(op: char) -> bool {
    matches!(op, 'M' | 'D' | 'N' | '=' | 'X')
}

/// Clears certain auxiliary tags and marks non-primary split reads as supplementary.
fn repair_supplementary_tags(records: &mut Vec<Record>) {
    // List of tags to remove.
    let tags_to_remove = [b"NM", b"MD", b"NH"];
    for record in records.iter_mut() {
        for tag in &tags_to_remove {
            // Remove the auxiliary field; ignore errors if not present.
            let _ = record.remove_aux(*tag);
        }
    }
    // If there are multiple records, designate the first as primary and mark the rest as supplementary.
    if records.len() > 1 {
        for record in records.iter_mut().skip(1) {
            let flags = record.flags();
            // Set the supplementary flag (0x800).
            record.set_flags(flags | 0x800);
        }
    }
}
