![splitncigar Logo](logo.png)

# splitncigar

version 0.2.0


**splitncigar** is a Rust-based tool for processing BAM files that contain spliced reads. It is modelled after the popular [Java implementation](https://gatk.broadinstitute.org/hc/en-us/articles/360036858811-SplitNCigarReads) with some additional tweaks.  It performs the following functions:

- **Splitting Reads**: Splits reads with 'N' operators in the CIGAR string (indicative of spliced reads) into subreads.
- **Flag Correction**: Optionally corrects the FLAG field in each split subread by copying the original FLAG from the unsplit read. This effort was motivated by this [publication](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-023-02923-y).
- **Direct Inspection**: In our experience, long-read RNA sequencing reads processed with the Java SplitNCigar followed by flag correction are not fully IGV-compatible, as the read details do not show when "clicking" on a read. In this Rust implementation, the output is not subject to this problem.
- **Supplementary Alignment (SA) Tag Handling**: Offers three modes for handling SA tags:
  - **overwrite**: Remove any old SA tag and write only the new subalignments.
  - **merge**: Merge new subalignments with any existing SA tag, avoiding duplicates.
  - **preserve**: Keep any existing SA tag if present; otherwise, add the new subalignments.
- **Additional Utilities**:
  - **slowview**: A lightweight BAM viewer.
  - **bamsummary**: Generates a summary report of BAM file statistics (e.g., read counts, average lengths, mapping qualities, and tag contributions).

## Features

- **Read Splitting**: Identify and split spliced reads (with intronic 'N' operators) into distinct subreads.
- **FLAG Correction**: By default, copies the FLAG from the original unsplit read to all split records. This can be disabled via a command-line flag.
- **Customizable SA Handling**: Choose from three different SA tag modes (overwrite, merge, or preserve) via the `--sa-mode` option.
- **Debug Mode**: Run with `--debug-splits` to print detailed processing information and step through splits interactively.
- **Utility Binaries**: In addition to the main `splitncigar` binary, use `slowview` for viewing BAM files and `bamsummary` for summarizing BAM statistics.  `cosmictovcf` takes a COSMIC tsv file and creates a vcf for easier annotation using VEP of SNPEff

## Installation

Build the project in release mode:
```sh
cargo build --release
```

Install the binaries (optional):
```sh
cp target/release/splitncigar ~/.local/bin/splitncigar
cp target/release/slowview   ~/.local/bin/slowview
cp target/release/bamsummary ~/.local/bin/bamsummary
```

### Prerequisites

- [Rust](https://www.rust-lang.org/) (latest stable release recommended)
- (Optional) [Samtools](http://www.htslib.org/) installed and accessible in your PATH
- (Optional) A reference FASTA file for the `--reference` option

### Building

Clone the repository:

```sh
git clone https://github.com/yourusername/splitncigar.git
cd splitncigar
```

## Usage

#### splitncigar
The primary binary splits spliced reads and optionally corrects FLAG and SA tag fields.

##### Basic Command
```sh
splitncigar --input <input.bam> --output <output.bam> --reference <reference.fa>
```

##### Options
--refactor-cigar-string
Refactor CIGAR strings containing NDN elements.
--skip-mq-transform
Skip the mapping quality transformation (255 → 60).
--debug-splits
Print debug information for each read and interactively prompt for user input.
--no-flag-correction
Disable flag correction (by default, flag correction is enabled).
--sa-mode <mode>

Specify how to handle SA tags. Valid options are:
overwrite
merge (default)
preserve


#### slowview

slowview provides a quick look at the contents of a BAM file.
```sh
slowview <bam_file>
```

#### bamsummary
bamsummary generates a summary report of a BAM file, including read counts, average read length, average mapping quality, and tag statistics.

```sh
bamsummary <bam_file>
```


## Examples

Below is an example workflow:
```sh
# Check input headers
samtools view data/HLA-A_reads_ds.bam | head -n 1
samtools view data/HLA-A_reads.snc.bam | head -n 2

# Build and install the binaries
cargo build --release && \
  cp target/release/splitncigar ~/.local/bin/splitncigar && \
  cp target/release/slowview   ~/.local/bin/slowview && \
  cp target/release/bamsummary ~/.local/bin/bamsummary

# Define reference FASTA
hg38=/path/to/GRCh38.p13.genome.fa

# Run splitncigar with default settings (flag correction enabled, SA mode merge)
splitncigar --input data/HLA-A_reads_ds.bam \
            --output HLA-A_reads.snc.bam \
            --reference $hg38

# Run splitncigar in debug mode
splitncigar --debug-splits \
            --input data/HLA-A_reads_ds.bam \
            --output HLA-A_reads.snc.bam \
            --reference $hg38

# View the resulting BAM file
samtools view HLA-A_reads.snc.bam | head -n 2
samtools sort HLA-A_reads.snc.bam -o HLA-A_reads.snc.sorted.bam
samtools index HLA-A_reads.snc.sorted.bam

# Use slowview to examine the sorted BAM file
slowview HLA-A_reads.snc.sorted.bam

# Generate BAM summary reports comparing output to data generated using the Java implementation 
# of SplitNCigar followed by an Rscript to fix the flags as mentioned in the publication above
bamsummary HLA-A_reads.snc.sorted.bam
bamsummary data/HLA-A_reads.snc.snc_fc.bam

```

## License

This project is licensed under the [MIT License](LICENSE). See the LICENSE file for details.

## Contributing

Contributions are welcome! Please open an issue or submit a pull request on GitHub.

## Contact

For questions or support, please contact us.

## Version History

<ul>
  <li>0.1.0 First Version</li>
  <li>0.2.0 Added cosmictovcf</li>
</ul>
