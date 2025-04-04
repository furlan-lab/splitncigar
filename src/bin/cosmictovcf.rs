/*
Below is a high-level overview of how the code processes the COSMIC TSV and creates the VCF, with special attention to how missing sequences (REF and ALT alleles) are handled:

1. **Reading Inputs and Setting Up:**  
   - **Command-Line Arguments:**  
     The program expects three arguments: the input COSMIC TSV file, the output VCF filename (which will be gzipped), and the reference FASTA file (with its index present).  
     
   - **Loading the Reference FASTA:**  
     The code opens the FASTA using a FASTA reader from rust-htslib. It also reads the corresponding FASTA index (.fai) to build a set of available contig names. This set is later used to map the chromosome information from the TSV to the proper FASTA contig (for example, by adding the "chr" prefix if needed).

   - **Building the VCF Header:**  
     Using the contig names (from the FASTA index) and other metadata (such as reference, filter, and INFO fields), the program builds a VCF header. All contig names in the header are formatted to start with "chr".

2. **Processing Each TSV Record:**  
   The TSV file is read record-by-record using a CSV reader configured for tab-delimited input. For each record:
   
   - **Filtering Incomplete Records:**  
     The program first checks if essential genomic information is missing (e.g., if the chromosome field is empty, starts with a "#", or if the genome start position is missing). Such records are counted as skipped.
     
   - **Mapping Chromosome to FASTA Contig:**  
     The code takes the chromosome value from the record and verifies whether it exists in the FASTA contigs set. If the raw value isn’t found, it attempts to add a "chr" prefix. If a matching contig is found, that name is used for FASTA lookups and then re-formatted (if needed) for the VCF header.

3. **Handling Alleles (REF and ALT):**  
   The code carefully handles cases where allele information may be incomplete:
   
   - **Missing REF Allele:**  
     If the REF allele from the TSV is empty, it is substituted with a dot `"."`. This indicates a missing allele in the input and triggers special handling in the normalization step.  
     
   - **Missing ALT Allele (Deletions):**  
     If the ALT allele is empty, the program interprets the event as a deletion. In VCF, deletions are represented by including the left‐flanking base:
     - The code fetches the left base from the reference (using the FASTA) at the position immediately preceding the mutation.
     - It then adjusts the position (shifting it one base to the left) and constructs a new REF allele by concatenating the left base with the original REF allele.
     - The new ALT allele becomes just the left base.  
     
   These steps ensure that even if the TSV does not supply complete allele information, the code can “rescue” the record by reconstructing a valid VCF representation.

4. **Normalization of Variants:**  
   After handling missing alleles:
   - The code calls a normalization function that may further adjust the variant’s position and allele representation (left-normalization).  
   - This normalization repeatedly checks if the last base of both alleles matches the base immediately preceding the variant in the reference. If so, the variant is shifted one base to the left. This standardizes the representation for indels.

5. **Creating and Writing the VCF Record:**  
   With the normalized information:
   - A VCF record is created where the contig is set (using the “chr”-prefixed name), the position is converted to 0-based indexing, and the alleles (REF and ALT) are assigned.
   - The record is annotated with additional INFO fields (such as gene symbol, transcript accession, etc.) taken from the TSV.
   - The record is then written to the VCF output, which is produced in BGZF (gzipped) format.

6. **Record Counters and Final Summary:**  
   Throughout processing, the code maintains counters for:
   - Total records read.
   - Records successfully processed (written to the VCF).
   - Records skipped (due to missing/incomplete genomic information or normalization issues).  
     
   At the end, a summary message is printed displaying these counts.

---

**In summary:**  
The code reads each COSMIC TSV record, checks and maps genomic information to the reference FASTA, and then handles incomplete allele data:
- Missing REF alleles are replaced with `"."` to trigger left-base fetching.
- Missing ALT alleles (interpreted as deletions) are rescued by prepending the left-flanking base from the reference, adjusting the variant accordingly.
After normalizing indels by shifting them left as necessary, the code writes a complete VCF record for each valid entry, all while tracking and reporting on the number of records processed versus skipped.

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

/// Naively left‐normalize an indel variant.
/// Shifts the variant to the left as long as the last base of the reference
/// allele and alternate allele is identical to the base immediately preceding the variant.
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
        // FASTA coordinates are 0-based; fetch region (pos-2, pos-1)
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
///   fetch the left flanking base from the reference and prepend it to both alleles.
/// - Then apply left normalization.
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
    // Expect three command-line arguments:
    // 1. Input COSMIC TSV file
    // 2. Output VCF file (gzipped)
    // 3. Reference FASTA (indexed)
    let args: Vec<String> = env::args().collect();
    if args.len() < 4 {
        eprintln!("Usage: {} <input_cosmic_tsv> <output_vcf.gz> <reference_fasta>", args[0]);
        std::process::exit(1);
    }
    let input_tsv = &args[1];
    let output_vcf = &args[2];
    let ref_fasta_path = &args[3];

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
                // For the VCF header, ensure contig names start with "chr".
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

    // Create a VCF writer with BGZF (gzipped) output.
    let mut writer = Writer::from_path(output_vcf, &header, true, bcf::Format::Vcf)?;

    let mut rdr = ReaderBuilder::new().delimiter(b'\t').from_path(input_tsv)?;

    // Counters for summary.
    let mut total_records = 0;
    let mut processed_records = 0;
    let mut skipped_records = 0;

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

        // Check if genomic information is incomplete.
        if record.CHROMOSOME.trim().is_empty() ||
           record.CHROMOSOME.starts_with('#') ||
           record.GENOME_START.trim().is_empty() {
            skipped_records += 1;
            continue;
        }

        let original_chrom = record.CHROMOSOME.trim();
        // Determine the correct FASTA contig name.
        let fasta_contig = if fasta_contigs.contains(original_chrom) {
            original_chrom.to_string()
        } else if fasta_contigs.contains(&format!("chr{}", original_chrom)) {
            format!("chr{}", original_chrom)
        } else {
            eprintln!("Sequence {} not found in FASTA", original_chrom);
            skipped_records += 1;
            continue;
        };
        // For VCF output, ensure the contig name starts with "chr".
        let vcf_chrom = if fasta_contig.starts_with("chr") {
            fasta_contig.clone()
        } else {
            format!("chr{}", fasta_contig)
        };

        // Parse genomic start position.
        let mut pos: usize = match record.GENOME_START.parse() {
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

        // If REF allele is empty, treat it as insertion (set to ".")
        let mut raw_ref = if raw_ref_orig.is_empty() { ".".to_string() } else { raw_ref_orig.to_string() };
        let mut raw_alt = raw_alt_orig.to_string();

        // If ALT allele is empty, interpret this as a deletion.
        // In VCF a deletion is represented by prepending the left flanking base.
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
        // VCF positions are 0-based.
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
        "Summary: Total records read: {}. Processed: {}. Skipped: {}.",
        total_records, processed_records, skipped_records
    );

    Ok(())
}
