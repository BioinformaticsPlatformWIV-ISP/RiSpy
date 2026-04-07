# RiSpy

### Fingerprint Detection Workflow

A workflow for **genetic fingerprint detection** using Whole-Genome Sequencing (WGS) data from **Illumina** and/or **Oxford Nanopore Technology (ONT)** platforms.  

- A **fingerprint** is defined as a collection of barcodes.  
- Each **barcode** consists of **two SNVs (single nucleotide variants)**.  

```text
Fingerprint
│
├── Barcode 1 ──> SNV1 + SNV2
├── Barcode 2 ──> SNV3 + SNV4
├── Barcode 3 ──> SNV5 + SNV6
└── Barcode N ──> SNV(2N-1) + SNV(2N)
```
----

## Installation

### Step 1 — Install the workflow
You can install the workflow directly from the repository:

```bash
pip install git+https://to_be_added.git
````

or, if you have cloned the repo locally:

```bash
pip install .
````

### Step 2 — Install additional dependencies with Conda
We provide an environment.yml file with the required tools

```bash
conda env create -f environment.yaml
conda activate fingerprint_environment
````
---
## Dependencies
This workflow was developed and tested using **Python 3.10**.   
In addition to the Python dependencies specified in `pyproject.toml`,  
the **Fingerprint Detection Workflow** requires several bioinformatics tools.  

These are the versions we have **tested and validated** the workflow on:  

- fastqc = 0.11.7  
- trimmomatic = 0.39  
- bwa = 0.7.17  
- samtools = 1.21  
- qualimap = 2.2.2  
- bcftools = 1.21  
- ucsc-liftover = 469  
- bedtools = 2.31.0  
- minimap2 = 2.28  
- transanno = 0.4.5  
- fastp = 1.0.1  
- nanoplot = 1.46.1  
- seqkit = 2.3.1  
- chopper = 0.8.0  
- deepvariant = 1.8.0

Other versions of these tools **may also work**, but they have not been tested by us.

### Manual installation
If you install these tools manually rather than through Conda,  
please ensure that the corresponding **binaries are available in your `PATH`**,  
otherwise the workflow will not be able to run them.

**Note on DeepVariant:** Installing DeepVariant through Conda is **not maintained nor recommended** by the DeepVariant team.  
> The recommended way to install DeepVariant is via **Docker**.  
> Please follow the official installation instructions here: [DeepVariant Quick Start](https://github.com/google/deepvariant/blob/r1.9/docs/deepvariant-quick-start.md)

---

## Command-Line Parameters

The `fingerprint` workflow provides multiple subcommands.  
Each subcommand has its own parameter set.

---

## Illumina-complete-run

### Description:

```markdown
Complete analysis of Illumina sequencing data, including read QC, trimming, mapping, variant calling, and fingerprint detection.
```
### Expected outputs:

```markdown
OUTPUT_DIR/Fingerprint_Detection – fingerprint detection results

OUTPUT_DIR/sample/ – final html report for each sample

OUTPUT_DIR/sample/fastqc_before_trimming/ – fastqc results for raw reads

OUTPUT_DIR/sample/fastqc_after_trimming/ – fastqc results for trimmed reads

OUTPUT_DIR/sample/depth/ – mean and median depth of coverage results

OUTPUT_DIR/sample/breadth/ – breadth of coverage results

OUTPUT_DIR/sample/count/ – number of primary reads mapped

OUTPUT_DIR/sample/fold80/ – fold80 results

OUTPUT_DIR/sample/mapping/ – bam file

OUTPUT_DIR/sample/qualimap/ – bam qc results

OUTPUT_DIR/sample/trimming/ – trimmed reads

OUTPUT_DIR/sample/deep_variant/ – raw and filtered variants results with stats
```

## usage: 

```markdown
fingerprint Illumina-complete-run --reference REF --barcodes BARCODES --ref-vcf REF_VCF --fastq-dir FASTQ_DIR --output-dir OUTPUT_DIR
                                         [--dry-run] [--chain-file CHAIN_FILE] [--chromosome-mapping CHROMOSOME_MAPPING] [--centromere CENTROMERE]
                                         [--threads THREADS] [--mem-mb MEM_MB] [--fastp-window-size FASTP_WINDOW_SIZE] [--fastp-quality FASTP_QUALITY]
                                         [--fastp-min-len FASTP_MIN_LEN] [--breadth-of-coverage-depth-cutoff DEPTH_CUTOFF]
                                         [--deep-variant-bin-version VERSION] [--variant-filtering-min-QUAL QUAL] [--variant-filtering-min-DP DP]
                                         [--fastqc-raw-threads N] [--fastqc-raw-memory MB] [--trim-reads-threads N] [--trim-reads-memory MB]
                                         [--fastqc-trimmed-threads N] [--fastqc-trimmed-memory MB] [--consensus-fasta-index-threads N]
                                         [--consensus-fasta-index-memory MB] [--alignment-threads N] [--alignment-memory MB]
                                         [--qualimap-threads N] [--qualimap-memory MB] [--samtools-depth-threads N] [--samtools-depth-memory MB]
                                         [--calculate-mean-and-median-depth-memory MB] [--breadth-of-coverage-memory MB] [--count-bam-threads N]
                                         [--count-bam-memory MB] [--fold-80-memory MB] [--deep-variant-threads N] [--deep-variant-memory MB]
                                         [--variant-filtering-memory MB] [--unique-SNPs-memory MB] [--bed-maker-memory MB] [--liftover-memory MB]
                                         [--chromosome-name-conversion-memory MB] [--filter-fingerprints-memory MB]
                                         [--validation-fingerprints-alt-memory MB] [--validation-fingerprints-ref-memory MB]
                                         [--barcodes-positions-new-ref-memory MB] [--chromosome-arms-memory MB]

options:
  --reference REF, -r REF
                        Reference genome FASTA [required]
  --barcodes BARCODES, -b BARCODES
                        TSV file containing 2-SNV barcodes [required] (default: data/Nipponbare_Barcodes.tsv)
  --ref-vcf REF_VCF, -v REF_VCF
                        VCF used to generate the barcodes [required]
  --fastq-dir FASTQ_DIR, -f FASTQ_DIR
                        Directory with FASTQ files (FASTQ files must follow the naming convention: SampleName_R1.fastq.gz and SampleName_R2.fastq.gz) [required]
  --output-dir OUTPUT_DIR, -o OUTPUT_DIR
                        Output directory [required]

  --dry-run             Generate the config and perform a dry run without executing
  --chain-file CHAIN_FILE
                        Chain file for liftover (default: data/CHAINFILE.chain)
  --chromosome-mapping CHROMOSOME_MAPPING
                        Chromosome name mapping (default: data/chromosome_mapping.txt)
  --centromere CENTROMERE
                        Centromere TSV file (default: data/Centromere.tsv)
  --threads THREADS, -t THREADS
                        Number of CPU cores (default: 8)
  --mem-mb MEM_MB, -m MEM_MB
                        Total memory in MB (default: 16000)
  --deep-variant-bin-version VERSION
                        DeepVariant binary version (default: 1.8.0)

  --fastp-window-size FASTP_WINDOW_SIZE
                        Window size for Fastp (default: 4)
  --fastp-quality FASTP_QUALITY
                        Quality cutoff for Fastp (default: 20)
  --fastp-min-len FASTP_MIN_LEN
                        Minimum read length for Fastp (default: 40)

  --breadth-of-coverage-depth-cutoff DEPTH_CUTOFF
                        Depth cutoff for breadth of coverage (default: 20)

  --variant-filtering-min-QUAL QUAL
                        Minimum QUAL for variant filtering (default: 20)
  --variant-filtering-min-DP DP
                        Minimum depth for variant filtering (default: 5)

  --fastqc-raw-threads N
  --fastqc-raw-memory MB
  --trim-reads-threads N
  --trim-reads-memory MB
  --fastqc-trimmed-threads N
  --fastqc-trimmed-memory MB
  --consensus-fasta-index-threads N
  --consensus-fasta-index-memory MB
  --alignment-threads N
  --alignment-memory MB
  --qualimap-threads N
  --qualimap-memory MB
  --samtools-depth-threads N
  --samtools-depth-memory MB
  --calculate-mean-and-median-depth-memory MB
  --breadth-of-coverage-memory MB
  --count-bam-threads N
  --count-bam-memory MB
  --fold-80-memory MB
  --deep-variant-threads N
  --deep-variant-memory MB
  --variant-filtering-memory MB
  --unique-SNPs-memory MB
  --bed-maker-memory MB
  --liftover-memory MB
  --chromosome-name-conversion-memory MB
  --filter-fingerprints-memory MB
  --validation-fingerprints-alt-memory MB
  --validation-fingerprints-ref-memory MB
  --barcodes-positions-new-ref-memory MB
  --chromosome-arms-memory MB
```

---

## ONT-complete-run

### Description:

```markdown
Complete analysis of ONT sequencing data, including read QC, trimming, mapping, merging bam files from different flow cells, variant calling, and fingerprint detection.
```
### Expected outputs:

```markdown
OUTPUT_DIR/merged_flowcells – final html report for merged flow cells

OUTPUT_DIR/Fingerprint_Detection – fingerprint detection results

OUTPUT_DIR/merged_flowcells/deep_variant/ – raw and filtered variants results with stats

OUTPUT_DIR/sample/nanoplot_raw/ – nanoplot results for raw reads

OUTPUT_DIR/sample/nanoplot_after_filtering/ – nanoplot results for trimmed reads

OUTPUT_DIR/sample/depth/ – mean and median depth of coverage results

OUTPUT_DIR/sample/breadth/ – breadth of coverage results

OUTPUT_DIR/sample/count/ – number of primary reads mapped

OUTPUT_DIR/sample/fold80/ – fold80 results

OUTPUT_DIR/sample/mapping/ – bam file

OUTPUT_DIR/sample/qualimap/ – bam qc results

OUTPUT_DIR/sample/chopper/ – trimmed reads
```

## usage: 

```markdown
fingerprint ONT-complete-run --reference REF --barcodes BARCODES --ref-vcf REF_VCF --fastq-dir FASTQ_DIR --output-dir OUTPUT_DIR
                                    [--chopper-quality Q] [--chopper-min-len L] [--nanoplot-raw-threads N] [--nanoplot-after-filtering-threads N]
                                    [--threads THREADS] [--mem-mb MEM_MB] [...]

options:
  --reference REF, -r REF
                        Reference genome FASTA [required]
  --barcodes BARCODES, -b BARCODES
                        TSV file containing 2-SNV barcodes [required] (default: data/Nipponbare_Barcodes.tsv)
  --ref-vcf REF_VCF, -v REF_VCF
                        VCF used to generate the barcodes [required]
  --fastq-dir FASTQ_DIR, -f FASTQ_DIR
                        Directory with FASTQ files (FASTQ files must follow the naming convention: SampleName.fastq.gz) [required]
  --output-dir OUTPUT_DIR, -o OUTPUT_DIR
                        Output directory [required]

  --chopper-quality Q   Quality cutoff for Chopper (default: 10)
  --chopper-min-len L   Minimum read length for Chopper (default: 150)
  --nanoplot-raw-threads N
                        Threads for NanoPlot (raw reads)
  --nanoplot-after-filtering-threads N
                        Threads for NanoPlot (after filtering)

  --threads THREADS     Number of CPU cores (default: 8)
  --mem-mb MEM_MB       Total memory in MB (default: 16000)
  [plus step-specific memory/thread options as in Illumina]
```
---

## scan-barcodes

### Description:

```markdown
Performing only fingerprint detection.
```
### Expected outputs:

```markdown
OUTPUT_DIR/sample – final html report for each sample

OUTPUT_DIR/Fingerprint_Detection – fingerprint detection results
```

## usage: 

```markdown
fingerprint scan-barcodes --barcodes BARCODES --ref-vcf REF_VCF --bam BAM --query-vcf QUERY_VCF --output-dir OUTPUT_DIR
                                 [--threads THREADS] [--mem-mb MEM_MB] [...]

options:
  --barcodes BARCODES, -b BARCODES
                        TSV file of 2-SNV barcodes [required] (default: data/Nipponbare_Barcodes.tsv)
  --ref-vcf REF_VCF, -v REF_VCF
                        Reference VCF used to generate the barcodes [required]
  --bam BAM, -B BAM     Input BAM file (indexed) [required]
  --query-vcf QUERY_VCF, -q QUERY_VCF
                        Input VCF file (indexed) [required]
  --output-dir OUTPUT_DIR, -o OUTPUT_DIR
                        Output directory [required]

  --threads THREADS     Number of CPU cores (default: 8)
  --mem-mb MEM_MB       Total memory in MB (default: 16000)
  [plus step-specific memory/thread options]
```

---

## generate-chain

### Description:

```markdown
Generates a liftover chain file between an old and a new reference genome. Required for coordinate conversion when comparing variants across different references.
```
### Expected outputs:

```markdown
OUTPUT_DIR/chain_file/CHAINFILE.chain – ready to use chain file for liftover
```

## usage: 

```markdown
fingerprint generate-chain --old-ref OLD_REF --new-ref NEW_REF --output-dir OUTPUT_DIR
                                  [--threads THREADS] [--mem-mb MEM_MB] [--paf2chain-memory MB] [--chain-sorting-memory MB]

options:
  --old-ref OLD_REF, -O OLD_REF
                        Reference genome FASTA to convert from [required]
  --new-ref NEW_REF, -N NEW_REF
                        Reference genome FASTA to convert to [required]
  --output-dir OUTPUT_DIR, -o OUTPUT_DIR
                        Output directory [required]

  --threads THREADS     Number of CPU cores (default: 8)
  --mem-mb MEM_MB       Total memory in MB (default: 16000)
  --paf2chain-memory MB Memory for PAF-to-chain conversion (default: 8000)
  --chain-sorting-memory MB
                        Memory for chain sorting (default: 8000)
  [plus step-specific memory/thread options]
```

---

> **Note:** The default number of threads is set conservatively and may not be optimal for large FASTQ datasets.  
> For faster analysis, we recommend increasing the thread count for compute-intensive steps   
> particularly **alignment** and **variant calling** up to the maximum supported by your system.

---

## Help
All commands support:  
```bash
fingerprint <command> --help
```

