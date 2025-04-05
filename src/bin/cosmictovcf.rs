#![allow(non_snake_case)]
#![allow(dead_code)]
use std::collections::HashSet;
use std::env;
use std::error::Error;
use std::io::{BufRead, BufReader, Read};

use csv::ReaderBuilder;
use rust_htslib::{bcf, bcf::header::Header, bcf::Writer, faidx};
use serde::Deserialize;
use flate2::read::GzDecoder;

/// A record representing one line of the COSMIC TSV.
#[derive(Debug, Deserialize)]
struct CosmicRecord {
    GENE_SYMBOL: String,
    COSMIC_GENE_ID: String,
    TRANSCRIPT_ACCESSION: String,
    COSMIC_SAMPLE_ID: String,
    SAMPLE_NAME: String,
    COSMIC_PHENOTYPE_ID: String,
    GENOMIC_MUTATION_ID: String,
    LEGACY_MUTATION_ID: String,
    MUTATION_ID: String,
    MUTATION_CDS: String,
    MUTATION_AA: String,
    MUTATION_DESCRIPTION: String,
    MUTATION_ZYGOSITY: String,
    LOH: String,
    CHROMOSOME: String,
    GENOME_START: String,
    GENOME_STOP: String,
    STRAND: String,
    PUBMED_PMID: String,
    COSMIC_STUDY_ID: String,
    HGVSP: String,
    HGVSC: String,
    HGVSG: String,
    GENOMIC_WT_ALLELE: String,
    GENOMIC_MUT_ALLELE: String,
    MUTATION_SOMATIC_STATUS: String,
}

/// Left-normalize an indel variant.
fn left_normalize_variant(
    chrom: &str,
    mut pos: usize,
    mut ref_allele: String,
    mut alt_allele: String,
    fasta: &faidx::Reader,
) -> Result<(usize, String, String), Box<dyn Error>> {
    if ref_allele.is_empty() || alt_allele.is_empty() {
        return Err("Reference or alternate allele is empty during normalization".into());
    }
    loop {
        if pos <= 1 {
            break;
        }
        let last_ref = ref_allele.chars().last().unwrap();
        let last_alt = alt_allele.chars().last().unwrap();
        if last_ref != last_alt {
            break;
        }
        let prev_base = fasta.fetch_seq(chrom, pos - 2, pos - 1)?.to_ascii_uppercase();
        if prev_base.len() != 1 {
            break;
        }
        let prev_base_str = std::str::from_utf8(&prev_base)?.to_uppercase();
        let prev_char = prev_base_str.chars().next().unwrap();
        if last_ref != prev_char {
            break;
        }
        pos -= 1;
        ref_allele = format!("{}{}", prev_char, &ref_allele[..ref_allele.len() - 1]);
        alt_allele = format!("{}{}", prev_char, &alt_allele[..alt_allele.len() - 1]);
    }
    Ok((pos, ref_allele, alt_allele))
}

/// Normalize a variant (prepend left base if needed and then left-normalize).
fn normalize_variant(
    chrom: &str,
    pos: usize,
    ref_allele: &str,
    alt_allele: &str,
    fasta: &faidx::Reader,
) -> Result<(usize, String, String), Box<dyn Error>> {
    let mut pos = pos;
    let mut norm_ref = ref_allele.to_string();
    let mut norm_alt = alt_allele.to_string();

    if norm_ref == "." || norm_ref.is_empty() || norm_ref.len() != norm_alt.len() {
        if pos > 1 {
            let left_base = fasta
                .fetch_seq(chrom, pos - 2, pos - 1)
                .map(|seq| seq.to_ascii_uppercase())
                .map_err(|_| format!("Failed to fetch sequence for {}:{}-{}", chrom, pos - 2, pos - 1))?;
            let left_base_str = String::from_utf8(left_base.clone())?;
            if norm_ref == "." || norm_ref.is_empty() {
                norm_ref = left_base_str.clone();
                norm_alt = format!("{}{}", left_base_str, norm_alt);
            } else {
                norm_ref = format!("{}{}", left_base_str.clone(), norm_ref);
                norm_alt = format!("{}{}", left_base_str, norm_alt);
            }
            pos -= 1;
        }
    }
    left_normalize_variant(chrom, pos, norm_ref, norm_alt, fasta)
}

fn main() -> Result<(), Box<dyn Error>> {
    // Expect 3 or 4 arguments: input TSV, output VCF, reference FASTA, [compress]
    let args: Vec<String> = env::args().collect();
    if args.len() < 4 {
        eprintln!(
            "Usage: {} <input_cosmic_tsv> <output_vcf> <reference_fasta> [compress]",
            args[0]
        );
        std::process::exit(1);
    }
    let input_tsv = &args[1];
    let output_vcf = &args[2];
    let ref_fasta_path = &args[3];

    // Parse compression option (default true).
    let compress_flag: bool = if args.len() > 4 {
        args[4].parse().unwrap_or(true)
    } else {
        true
    };

    // Open the FASTA.
    let fasta = faidx::Reader::from_path(ref_fasta_path)?;

    // Build contig set and VCF header.
    let mut fasta_contigs = HashSet::new();
    let mut header = Header::new();
    header.push_record(b"##fileformat=VCFv4.2");
    header.push_record(format!("##reference={}", ref_fasta_path).as_bytes());
    header.push_record(b"##FILTER=<ID=PASS,Description=\"All filters passed\">");
    header.push_record(b"##INFO=<ID=GeneSymbol,Number=1,Type=String,Description=\"Gene symbol\">");
    header.push_record(b"##INFO=<ID=Transcript,Number=1,Type=String,Description=\"Transcript Accession\">");
    header.push_record(b"##INFO=<ID=SampleName,Number=1,Type=String,Description=\"Sample Name\">");
    header.push_record(b"##INFO=<ID=PubMed,Number=1,Type=String,Description=\"PubMed ID\">");
    header.push_record(b"##INFO=<ID=StudyID,Number=1,Type=String,Description=\"COSMIC Study ID\">");
    header.push_record(b"##INFO=<ID=CDS,Number=1,Type=String,Description=\"CDS Mutation\">");
    header.push_record(b"##INFO=<ID=AA,Number=1,Type=String,Description=\"Protein Mutation\">");
    header.push_record(b"##INFO=<ID=Description,Number=1,Type=String,Description=\"Mutation Description\">");
    header.push_record(b"##INFO=<ID=Zygosity,Number=1,Type=String,Description=\"Zygosity\">");
    header.push_record(b"##INFO=<ID=LOH,Number=1,Type=String,Description=\"Loss of Heterozygosity\">");
    header.push_record(b"##INFO=<ID=SomaticStatus,Number=1,Type=String,Description=\"Somatic Status\">");

    let fai_path = format!("{}.fai", ref_fasta_path);
    let fai_file = std::fs::File::open(&fai_path)?;
    let reader = BufReader::new(fai_file);
    for line in reader.lines() {
        let line = line?;
        let mut parts = line.split('\t');
        if let Some(contig) = parts.next() {
            fasta_contigs.insert(contig.to_string());
            if let Some(length_str) = parts.next() {
                let contig_length: u64 = length_str.parse()?;
                let vcf_contig = if contig.starts_with("chr") {
                    contig.to_string()
                } else {
                    format!("chr{}", contig)
                };
                header.push_record(format!("##contig=<ID={},length={}>", vcf_contig, contig_length).as_bytes());
            }
        }
    }
    header.push_record(b"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO");

    // Create VCF writer (invert compress_flag per rust-htslib behavior).
    let mut writer = Writer::from_path(output_vcf, &header, !compress_flag, bcf::Format::Vcf)?;

    // Open the TSV file (handle gzipped input if needed).
    let file = std::fs::File::open(input_tsv)?;
    let reader: Box<dyn Read> = if input_tsv.ends_with(".gz") {
        Box::new(GzDecoder::new(file))
    } else {
        Box::new(file)
    };
    let mut rdr = ReaderBuilder::new().delimiter(b'\t').from_reader(reader);

    // Counters for summary.
    let mut total_records = 0;
    let mut processed_records = 0;
    let mut skipped_records = 0;
    let mut ref_mismatch_errors = 0;
    let mut corrected_records = 0;

    for result in rdr.deserialize() {
        total_records += 1;
        let record: CosmicRecord = match result {
            Ok(r) => r,
            Err(e) => {
                eprintln!("Skipping invalid record: {}", e);
                skipped_records += 1;
                continue;
            }
        };

        if record.CHROMOSOME.trim().is_empty()
            || record.CHROMOSOME.starts_with('#')
            || record.GENOME_START.trim().is_empty()
        {
            skipped_records += 1;
            continue;
        }

        let original_chrom = record.CHROMOSOME.trim();
        let fasta_contig = if fasta_contigs.contains(original_chrom) {
            original_chrom.to_string()
        } else if fasta_contigs.contains(&format!("chr{}", original_chrom)) {
            format!("chr{}", original_chrom)
        } else {
            eprintln!("Sequence {} not found in FASTA", original_chrom);
            skipped_records += 1;
            continue;
        };
        let vcf_chrom = if fasta_contig.starts_with("chr") {
            fasta_contig.clone()
        } else {
            format!("chr{}", fasta_contig)
        };

        // Parse GENOME_START and GENOME_STOP.
        let mut pos: usize = match record.GENOME_START.trim().parse() {
            Ok(p) => p,
            Err(_) => {
                eprintln!("Skipping record with invalid start position: {}", record.GENOME_START);
                skipped_records += 1;
                continue;
            }
        };
        let stop: usize = match record.GENOME_STOP.trim().parse() {
            Ok(s) => s,
            Err(_) => {
                eprintln!("Skipping record with invalid stop position: {}", record.GENOME_STOP);
                skipped_records += 1;
                continue;
            }
        };

        let var_id = record.GENOMIC_MUTATION_ID.as_str();
        let raw_ref_orig = record.GENOMIC_WT_ALLELE.trim();
        let raw_alt_orig = record.GENOMIC_MUT_ALLELE.trim();

        // Declare mutable raw_ref and raw_alt so that we can update them if needed.
        let mut raw_ref = if raw_ref_orig.is_empty() {
            ".".to_string()
        } else {
            raw_ref_orig.to_string()
        };
        let mut raw_alt = raw_alt_orig.to_string();

        // Full reference check:
        // If TSV provides a nonempty REF (not "."), fetch the reference sequence for [GENOME_START, GENOME_STOP].
        // For one-base substitutions, if REF length is 1 and stop == pos+1, subtract 1 from stop.
        if !raw_ref_orig.is_empty() && raw_ref_orig != "." {
            // if raw_ref_orig.len() == 1 && stop == pos + 1 {
            //     stop -= 1;
            // }
            let expected_full = fasta.fetch_seq(&fasta_contig, pos - 1, stop -1 )?;
            let expected_full_str = String::from_utf8(expected_full)?.to_uppercase();
            if expected_full_str != raw_ref_orig.to_uppercase() {
                // Attempt correction for one-base REF.
                if raw_ref_orig.len() == 1 && expected_full_str.len() == 2 && pos > 1 {
                    let expected_extra = fasta.fetch_seq(&fasta_contig, pos - 2, pos - 2 + raw_ref_orig.len() + 1)?;
                    let expected_extra_str = String::from_utf8(expected_extra)?.to_uppercase();
                    if expected_extra_str.ends_with(&raw_ref_orig.to_uppercase()) {
                        let extra = expected_extra_str.len() - raw_ref_orig.len();
                        pos = pos.saturating_sub(extra);
                        raw_ref = expected_extra_str.clone();
                        if !raw_alt_orig.is_empty() {
                            raw_alt = format!("{}{}", expected_extra_str.chars().next().unwrap(), raw_alt_orig);
                        }
                        corrected_records += 1;
                    } else {
                        eprintln!(
                            "Full reference mismatch for record {}: TSV REF '{}' does not match reference '{}'",
                            var_id, raw_ref_orig, expected_full_str
                        );
                        ref_mismatch_errors += 1;
                        continue;
                    }
                } else {
                    eprintln!(
                        "Full reference mismatch for record {}: TSV REF '{}' does not match reference '{}'",
                        var_id, raw_ref_orig, expected_full_str
                    );
                    ref_mismatch_errors += 1;
                    continue;
                }
            }
        }

        // Handle missing ALT allele as deletion.
        if raw_alt.is_empty() {
            if pos > 1 {
                let left_base = fasta.fetch_seq(&fasta_contig, pos - 2, pos - 1)?.to_ascii_uppercase();
                let left_base_str = String::from_utf8(left_base)?;
                let new_pos = pos - 1;
                let new_ref = format!("{}{}", left_base_str, raw_ref);
                let new_alt = left_base_str.clone();
                pos = new_pos;
                raw_ref = new_ref;
                raw_alt = new_alt;
            } else {
                eprintln!("Skipping deletion variant at beginning: REF='{}', ALT='{}'", raw_ref, raw_alt);
                skipped_records += 1;
                continue;
            }
        }

        let (norm_pos, norm_ref, norm_alt) =
            match normalize_variant(&fasta_contig, pos, &raw_ref, &raw_alt, &fasta) {
                Ok(v) => v,
                Err(e) => {
                    eprintln!("Skipping record due to normalization error: {}", e);
                    skipped_records += 1;
                    continue;
                }
            };

        let mut vcf_record = writer.empty_record();
        let rid = writer.header().name2rid(vcf_chrom.as_bytes())?;
        vcf_record.set_rid(Some(rid));
        vcf_record.set_pos((norm_pos - 1) as i64);
        let _ = vcf_record.set_id(var_id.as_bytes());
        if let Err(e) = vcf_record.set_alleles(&[norm_ref.as_bytes(), norm_alt.as_bytes()]) {
            eprintln!("Failed to set alleles for record {}: {}", var_id, e);
            skipped_records += 1;
            continue;
        }
        vcf_record.set_qual(-1.0);
        vcf_record.set_filters(&[&b"PASS"[..]])?;
        vcf_record.push_info_string(b"GeneSymbol", &[record.GENE_SYMBOL.as_bytes()])?;
        vcf_record.push_info_string(b"Transcript", &[record.TRANSCRIPT_ACCESSION.as_bytes()])?;
        vcf_record.push_info_string(b"SampleName", &[record.SAMPLE_NAME.as_bytes()])?;
        vcf_record.push_info_string(b"PubMed", &[record.PUBMED_PMID.as_bytes()])?;
        vcf_record.push_info_string(b"StudyID", &[record.COSMIC_STUDY_ID.as_bytes()])?;
        vcf_record.push_info_string(b"CDS", &[record.MUTATION_CDS.as_bytes()])?;
        vcf_record.push_info_string(b"AA", &[record.MUTATION_AA.as_bytes()])?;
        vcf_record.push_info_string(b"Description", &[record.MUTATION_DESCRIPTION.as_bytes()])?;
        vcf_record.push_info_string(b"Zygosity", &[record.MUTATION_ZYGOSITY.as_bytes()])?;
        vcf_record.push_info_string(b"LOH", &[record.LOH.as_bytes()])?;
        vcf_record.push_info_string(b"SomaticStatus", &[record.MUTATION_SOMATIC_STATUS.as_bytes()])?;

        writer.write(&vcf_record)?;
        processed_records += 1;
    }

    println!(
        "Summary: Total records read: {}. Processed: {}. Skipped: {}. Reference mismatches: {}. Corrected records: {}.",
        total_records, processed_records, skipped_records, ref_mismatch_errors, corrected_records
    );

    Ok(())
}
