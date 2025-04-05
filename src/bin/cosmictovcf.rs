/*
Below is a high-level overview of how the code processes the COSMIC TSV and creates the VCF, with special attention to how missing sequences (REF and ALT alleles) are handled:

1. **Reading Inputs and Setting Up:**  
   - **Command-Line Arguments:**  
     The program expects at least three arguments: the input COSMIC TSV file, the output VCF filename, and the reference FASTA file (with its index present). An optional fourth argument ("true" or "false") specifies whether to compress the output.
     
   - **Loading the Reference FASTA:**  
     The code opens the FASTA (using rust-htslib) and reads the corresponding FASTA index (.fai) to build a set of available contig names. This set is later used to map the TSV chromosome to the correct FASTA contig (adding a "chr" prefix if needed).

   - **Building the VCF Header:**  
     Using the contig names and other metadata, the program builds a VCF header (all contig names are ensured to start with "chr").

2. **Processing Each TSV Record:**  
   The TSV file is read record-by-record. For each record:
   
   - **Filtering Incomplete Records:**  
     Records missing essential genomic information (e.g. chromosome or start position) are skipped.
     
   - **Mapping Chromosome to FASTA Contig:**  
     The TSV’s chromosome value is mapped to the appropriate FASTA contig name.
     
   - **Full Reference Check:**  
     If the TSV provides a nonempty REF allele (i.e. not missing or “.”), the code fetches the full reference sequence from the FASTA using the TSV’s GENOME_START and GENOME_STOP positions. It then confirms that this fetched sequence matches the TSV’s REF allele (after converting both to uppercase). If they do not match, the record is skipped and counted as a reference mismatch.
     
   - **Handling Alleles:**  
     *Missing REF Allele:* If the TSV REF is empty, it’s replaced with “.” (triggering later left‑flanking base addition).  
     *Missing ALT Allele (Deletions):* If the ALT allele is empty, the record is assumed to represent a deletion; the left‑flanking base is fetched, the coordinate is adjusted (shifted one base to the left), and the REF and ALT are updated accordingly.
     
   - **Normalization:**  
     The variant is left‑normalized (shifting left while the last bases of REF and ALT match the preceding reference base).
     
   - **Writing the VCF Record:**  
     A VCF record is created (with 0‑based positions) and annotated with INFO fields from the TSV. It’s then written to the output VCF (which is compressed or not according to the optional flag).

3. **Summary:**  
   At the end, the code prints a summary with counts for total records read, processed, skipped, reference mismatches, and any records that were “corrected” (if applicable).

**In summary:**  
The code confirms that each TSV record’s provided reference allele (over the region defined by GENOME_START–GENOME_STOP) matches the FASTA’s sequence before proceeding. Missing REF alleles trigger left‑flanking base addition for insertions, and missing ALT alleles (deletions) are rescued similarly. After normalization, a valid VCF record is written for each record that passes the checks.

*/

#![allow(non_snake_case)]
#![allow(dead_code)]
use std::collections::HashSet;
use std::env;
use std::error::Error;
use std::io::{BufRead, BufReader};

use csv::ReaderBuilder;
use rust_htslib::{bcf, bcf::header::Header, bcf::Writer, faidx};
use serde::Deserialize;

/// A record representing one line of the COSMIC TSV.
/// Field names must match the TSV header exactly.
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

/// Naively left‑normalize an indel variant.
/// Shifts the variant to the left as long as the last base of the REF and ALT alleles is identical to
/// the base immediately preceding the variant.
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
        // Cannot shift if at the beginning.
        if pos <= 1 {
            break;
        }
        let last_ref = ref_allele.chars().last().unwrap();
        let last_alt = alt_allele.chars().last().unwrap();
        if last_ref != last_alt {
            break;
        }
        // FASTA coordinates are 0‑based; fetch region (pos‑2, pos‑1)
        let prev_base = fasta.fetch_seq(chrom, pos - 2, pos - 1)?.to_ascii_uppercase();
        if prev_base.len() != 1 {
            break;
        }
        let prev_base_str = std::str::from_utf8(&prev_base)?.to_uppercase();
        let prev_char = prev_base_str.chars().next().unwrap();
        if last_ref != prev_char {
            break;
        }
        // Shift left: drop the last base from both alleles and prepend the left base.
        pos -= 1;
        ref_allele = format!("{}{}", prev_char, &ref_allele[..ref_allele.len() - 1]);
        alt_allele = format!("{}{}", prev_char, &alt_allele[..alt_allele.len() - 1]);
    }
    Ok((pos, ref_allele, alt_allele))
}

/// Normalize a variant:
/// - If the REF allele is missing (i.e. "." or empty) or if the alleles differ in length (an indel),
///   fetch the left‑flanking base from the reference and prepend it to both alleles.
/// - Then apply left‑normalization.
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
    // Expect three or four command-line arguments:
    // 1. Input COSMIC TSV file
    // 2. Output VCF file
    // 3. Reference FASTA (indexed)
    // 4. (Optional) Compression option ("true" for compressed, "false" for uncompressed)
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

    // Parse compression option. Default: true (compress)
    // NOTE: Due to rust-htslib behavior, passing true means we want compression,
    // so we invert the flag when calling Writer::from_path.
    let compress_flag: bool = if args.len() > 4 {
        args[4].parse().unwrap_or(true)
    } else {
        true
    };

    // Open the FASTA.
    let fasta = faidx::Reader::from_path(ref_fasta_path)?;

    // Build a set of contig names from the FASTA index (.fai)
    let mut fasta_contigs = HashSet::new();
    // Build the VCF header.
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
                // Ensure contig names in header start with "chr".
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

    // Create a VCF writer.
    // Due to rust-htslib behavior, we pass !compress_flag to force_bgzf.
    let mut writer = Writer::from_path(output_vcf, &header, !compress_flag, bcf::Format::Vcf)?;

    let mut rdr = ReaderBuilder::new().delimiter(b'\t').from_path(input_tsv)?;

    // Counters for summary.
    let mut total_records = 0;
    let mut processed_records = 0;
    let mut skipped_records = 0;
    let mut ref_mismatch_errors = 0;
    let corrected_records = 0;

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

        // Check for missing essential genomic info.
        if record.CHROMOSOME.trim().is_empty()
            || record.CHROMOSOME.starts_with('#')
            || record.GENOME_START.trim().is_empty()
        {
            skipped_records += 1;
            continue;
        }

        let original_chrom = record.CHROMOSOME.trim();
        // Map the TSV chromosome to a FASTA contig.
        let fasta_contig = if fasta_contigs.contains(original_chrom) {
            original_chrom.to_string()
        } else if fasta_contigs.contains(&format!("chr{}", original_chrom)) {
            format!("chr{}", original_chrom)
        } else {
            eprintln!("Sequence {} not found in FASTA", original_chrom);
            skipped_records += 1;
            continue;
        };
        // For VCF output, ensure contig name starts with "chr".
        let vcf_chrom = if fasta_contig.starts_with("chr") {
            fasta_contig.clone()
        } else {
            format!("chr{}", fasta_contig)
        };

        // Parse genomic start position.
        let mut pos: usize = match record.GENOME_START.trim().parse() {
            Ok(p) => p,
            Err(_) => {
                eprintln!("Skipping record with invalid start position: {}", record.GENOME_START);
                skipped_records += 1;
                continue;
            }
        };

        let var_id = record.GENOMIC_MUTATION_ID.as_str();

        let raw_ref_orig = record.GENOMIC_WT_ALLELE.trim();
        let raw_alt_orig = record.GENOMIC_MUT_ALLELE.trim();

        // Full reference check: if TSV provides a nonempty REF allele (not "."), confirm that the
        // reference sequence from the FASTA for the region [GENOME_START, GENOME_STOP] matches it.
        if !raw_ref_orig.is_empty() && raw_ref_orig != "." {
            let stop: usize = match record.GENOME_STOP.trim().parse() {
                Ok(s) => s,
                Err(_) => {
                    eprintln!("Skipping record with invalid stop position: {}", record.GENOME_STOP);
                    skipped_records += 1;
                    continue;
                }
            };
            let expected_full = fasta.fetch_seq(&fasta_contig, pos - 1, stop - 1)?;
            let expected_full_str = String::from_utf8(expected_full)?.to_uppercase();
            if expected_full_str != raw_ref_orig.to_uppercase() {
                eprintln!(
                    "Full reference mismatch for record {}: TSV REF '{}' does not match reference '{}'",
                    var_id, raw_ref_orig, expected_full_str
                );
                ref_mismatch_errors += 1;
                continue;
            }
        }

        // If REF allele is empty, treat as insertion by setting it to "."
        let mut raw_ref = if raw_ref_orig.is_empty() {
            ".".to_string()
        } else {
            raw_ref_orig.to_string()
        };
        let mut raw_alt = raw_alt_orig.to_string();

        // If ALT allele is empty, treat as deletion.
        // For a deletion, fetch the left‑flanking base and adjust coordinates and alleles.
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

        // Normalize the variant (left‑normalization).
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
        // VCF positions are 0‑based.
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
