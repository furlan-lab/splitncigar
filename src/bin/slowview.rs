use rust_htslib::bam::{self, Read, Record};
use rust_htslib::bam::record::Aux;
use std::error::Error;
use std::io::{stdin, stdout, Write};
// use std::path::PathBuf;

pub fn main() -> Result<(), Box<dyn Error>> {
    // Parse command-line argument for the BAM file path
    let bam_path = match std::env::args().nth(1) {
        Some(p) => p,
        None => {
            eprintln!("Usage: slow_bam_viewer <input.bam>");
            std::process::exit(1);
        }
    };

    // Open the BAM file for reading
    let mut reader = bam::Reader::from_path(&bam_path)?;
    let header = bam::Header::from_template(reader.header());

    // We'll convert the header into a `HeaderView` for fetching reference names
    let header_view = bam::HeaderView::from_header(&header);

    // Iterate over each read in the BAM
    for (i, result) in reader.records().enumerate() {
        let record = result?;

        // Print an index for reference
        println!("================ READ {} ================", i + 1);

        // 1) QNAME
        let qname = String::from_utf8_lossy(record.qname());
        println!("QNAME: {}", qname);

        // 2) Flags (integer + binary format)
        let flags = record.flags();
        println!("FLAGS: {} (binary: {:b})", flags, flags);

        // 3) Reference name
        let tid = record.tid();
        if tid < 0 {
            println!("RNAME: * (unmapped)");
        } else {
            // Convert TID to the actual reference name
            let refname_bytes = header_view.tid2name(tid as u32);
            let refname = String::from_utf8_lossy(refname_bytes);
            println!("RNAME: {}", refname);
        }

        // 4) Position (0-based in rust-htslib)
        println!("POS: {}", record.pos());

        // 5) MAPQ
        println!("MAPQ: {}", record.mapq());

        // 6) CIGAR
        let cigar_string = record.cigar().to_string();
        println!("CIGAR: {}", cigar_string);

        // 7) Sequence
        //    In many short-read BAMs, rust-htslib stores sequences in nibble-coded form,
        //    but your data may have ASCII. We'll decode it anyway:
        let seq_str = decode_seq(&record);
        println!("SEQ: {}", seq_str);

        // 8) Qualities (ASCII phred)
        // let qual = record.qual();
        // Build a string of ASCII characters offset by 33
        // or just print them as numeric for clarity
        // let qual_string: Vec<String> = qual.iter().map(|&q| q.to_string()).collect();
        let qual_string = decode_qual(&record);
        println!("QUAL: [{}]", qual_string);

        // 9) All auxiliary tags
        //    We'll list each "TAG:VALUE" in alphabetical order by tag
        let mut aux_data = Vec::new();
        for aux in record.aux_iter() {
            match aux {
                Ok((tag, value)) => {
                    aux_data.push((tag, aux_value_to_string(value)));
                }
                Err(e) => {
                    eprintln!("Error parsing aux: {}", e);
                }
            }
        }
        // Sort by the two-byte tag, e.g. b"NM", b"MD"
        aux_data.sort_by_key(|(tag, _)| (*tag));
        println!("TAGS:");
        for (tag, val_str) in aux_data {
            let t_str = String::from_utf8_lossy(&tag);
            println!("  {}: {}", t_str, val_str);
        }

        // Pause for user input
        if !prompt_for_next_read()? {
            println!("Quitting...");
            break;
        }
    }

    Ok(())
}

/// Convert a record's nibble-coded or ASCII-coded bases to a string
fn decode_seq(record: &Record) -> String {
    record.seq().as_bytes().iter().map(|&b| {
        match b {
            1 => 'A',
            2 => 'C',
            4 => 'G',
            8 => 'T',
            15 => 'N',
            _ => b as char, // If it's ASCII, or an IUPAC code, etc.
        }
    }).collect()
}

/// Convert a record's nibble-coded or ASCII-coded bases to a string
fn decode_qual(record: &Record) -> String {
    record.qual().iter().map(|f| *f as char).collect()
}

/// Convert an Aux value into a readable String
fn aux_value_to_string(value: Aux) -> String {
    match value {
        Aux::Char(c) => format!("Char({})", c as char),
        Aux::I8(v) => format!("I8({})", v),
        Aux::U8(v) => format!("U8({})", v),
        Aux::I16(v) => format!("I16({})", v),
        Aux::U16(v) => format!("U16({})", v),
        Aux::I32(v) => format!("I32({})", v),
        Aux::U32(v) => format!("U32({})", v),
        Aux::Float(f) => format!("Float({:.3})", f),
        Aux::Double(d) => format!("Double({:.3})", d),
        Aux::String(s) => format!("String({})", s),
        Aux::HexByteArray(h) => format!("Hex({})", String::from_utf8_lossy(h.as_ref())),
        Aux::ArrayFloat(floats) => format!("ArrayFloat({:?})", floats),
        // Aux::ArrayDouble(a) => format!("ArrayDouble({:?})", a),
        Aux::ArrayI8(a) => format!("ArrayI8({:?})", a),
        Aux::ArrayU8(a) => format!("ArrayU8({:?})", a),
        Aux::ArrayI16(a) => format!("ArrayI16({:?})", a),
        Aux::ArrayU16(a) => format!("ArrayU16({:?})", a),
        Aux::ArrayI32(a) => format!("ArrayI32({:?})", a),
        Aux::ArrayU32(a) => format!("ArrayU32({:?})", a),
    }
}

/// Prompt user to press Enter for next read, or 'q' to quit.
fn prompt_for_next_read() -> Result<bool, Box<dyn Error>> {
    eprint!("Press Enter for next read, or 'q' then Enter to quit: ");
    stdout().flush()?;

    let mut input = String::new();
    stdin().read_line(&mut input)?;
    if input.trim().eq_ignore_ascii_case("q") {
        return Ok(false);
    }
    Ok(true)
}
