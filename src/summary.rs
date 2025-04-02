use rust_htslib::bam::{self, Read};
use rust_htslib::bam::record::Aux;
use std::collections::HashMap;
use std::env;

fn main() {
    let args: Vec<String> = env::args().collect();
    let bam_path = &args[1];

    let mut reader = bam::Reader::from_path(bam_path).expect("Failed to open BAM file");

    let mut total_reads = 0;
    let mut total_seq_len = 0;
    let mut total_mapq = 0;
    let mut tag_counts: HashMap<String, usize> = HashMap::new();
    let mut tag_bytes: HashMap<String, usize> = HashMap::new();

    for result in reader.records() {
        let record = result.expect("Failed to read record");
        total_reads += 1;
        total_seq_len += record.seq().len();
        total_mapq += record.mapq() as usize;

        let aux_iter = record.aux_iter();
        for aux in aux_iter {
            if let Ok((tag, value)) = aux {
                let tag_str = String::from_utf8_lossy(&tag).to_string();
                *tag_counts.entry(tag_str.clone()).or_insert(0) += 1;
                *tag_bytes.entry(tag_str.clone()).or_insert(0) += aux_size(&value);
            }
        }
    }

    println!("Summary for BAM file: {}", bam_path);
    println!("Total reads: {}", total_reads);
    println!("Average read length: {:.2}", total_seq_len as f64 / total_reads as f64);
    println!("Average MAPQ: {:.2}", total_mapq as f64 / total_reads as f64);
    println!("Tag sizes (approx. byte contribution):");

    let mut tags_sorted: Vec<_> = tag_bytes.iter().collect();
    tags_sorted.sort_by(|a, b| b.1.cmp(a.1));
    for (tag, bytes) in tags_sorted {
        let count = tag_counts.get(tag).unwrap_or(&0);
        println!("  {:<2}  {:>6} times, ~{:>8} bytes", tag, count, bytes);
    }
}

fn aux_size(aux: &Aux) -> usize {
    match aux {
        Aux::Char(_) => 1,
        Aux::I8(_) | Aux::U8(_) => 1,
        Aux::I16(_) | Aux::U16(_) => 2,
        Aux::I32(_) | Aux::U32(_) => 4,
        Aux::Float(_) => 4,
        Aux::Double(_) => 8,
        Aux::String(s) => s.len() + 1, // null-terminated
        Aux::HexByteArray(v) => v.len(),
        Aux::ArrayI8(a) => a.len(),
        Aux::ArrayU8(a) => a.len(),
        Aux::ArrayI16(a) => a.len() * 2,
        Aux::ArrayU16(a) => a.len() * 2,
        Aux::ArrayI32(a) => a.len() * 4,
        Aux::ArrayU32(a) => a.len() * 4,
        Aux::ArrayFloat(a) => a.len() * 4,
    }
}
