import typer
from pathlib import Path
import os
import fingerprint_cli
from fingerprint_cli.cli_utils import (
    get_default_centromere_path,
    get_default_chain_file_path,
    get_default_chromosome_mapping_path,
    get_default_Nip_barcodes_path,
    get_default_Nip_ref_vcf_path,
    load_and_write_config,
    run_snakemake,
    set_paths,
    set_resources
)
from fingerprint_cli.validators.illumina import validate_illumina_inputs
from fingerprint_cli.validators.ont import validate_ont_inputs
from fingerprint_cli.validators.detection import validate_detection_inputs
from fingerprint_cli.validators.chain import validate_chain_inputs

app = typer.Typer(help="Detection of 2-SNV barcodes (fingerprint) from WGS data")

""""
-----------------------------
Illumina FULL PIPELINE SUBCOMMAND
-----------------------------
"""

@app.command("Illumina-complete-run")
def full_pipeline_Illumina(
    reference: Path = typer.Option(..., "--reference", "-r", help="Reference genome FASTA"),
    barcodes: Path = typer.Option(..., "--barcodes", "-b", help="tsv file of the fingerprint (2‑SNV barcodes)"),
    ref_vcf: Path = typer.Option(..., "--ref-vcf", "-v", help="VCF used to generate the fingerprint"),
    output_dir: Path = typer.Option(..., "--output-dir", "-o", help="Output directory"),
    fastq_dir: Path = typer.Option(..., "--fastq-dir", "-f", help="Directory with FASTQ files (FASTQ files must follow the naming convention: SampleName_R1.fastq.gz and SampleName_R2.fastq.gz)"),
    dry_run: bool = typer.Option(False, "--dry-run", help="Generate the config file and perform a dry-run without executing rules"),
    chain_file: Path = typer.Option(get_default_chain_file_path(), "--chain-file", "-C", help="Chain file for liftover"),
    chromosome_mapping: Path = typer.Option(get_default_chromosome_mapping_path(), "--chromosome-mapping", "-M", help="Chromosome name mapping"),
    centromere: Path = typer.Option(get_default_centromere_path(), "--centromere", "-c", help="Centromere TSV file"),
    cores: int = typer.Option(8, "--threads", "-t", help="Number of CPU cores"),
    mem_mb: int = typer.Option(16000, "--mem-mb", "-m", help="Memory in MB"),

    fastp_window_size: int = typer.Option(4, "--fastp-window-size", help="window size for fastp"),
    fastp_quality: int = typer.Option(20, "--fastp-quality", help="quality cutoff for fastp"),
    fastp_min_len: int = typer.Option(40, "--fastp-min-len", help="minimum length for fastp"),
    breadth_of_coverage_depth_cutoff: int = typer.Option(20, "--breadth-of-coverage-depth-cutoff", help="depth cutoff for calculating the breadth of coverage"),
    deep_variant_BIN_VERSION: str = typer.Option("1.8.0", "--deep-variant-bin-version", help="BIN VERSION for DeepVariant"),
    variant_filtering_min_QUAL: int = typer.Option(20, "--variant-filtering-min-QUAL", help="minimum QUAL for variant filtering"),
    variant_filtering_min_DP: int = typer.Option(5, "--variant-filtering-min-DP", help="minimum DP for variant filtering"),

    fastqc_raw_threads: int = typer.Option(4, "--fastqc-raw-threads", help="threads for fastqc_raw"),
    fastqc_raw_memory: int = typer.Option(8000, "--fastqc-raw-memory", help="memory (MB) for fastqc_raw"),
    trim_reads_threads: int = typer.Option(8, "--trim-reads-threads", help="threads for trim_reads"),
    trim_reads_memory: int = typer.Option(16000, "--trim-reads-memory", help="memory (MB) for trim_reads"),
    fastqc_trimmed_threads: int = typer.Option(4, "--fastqc-trimmed-threads", help="threads for fastqc_trimmed"),
    fastqc_trimmed_memory: int = typer.Option(8000, "--fastqc-trimmed-memory", help="memory (MB) for fastqc_trimmed"),
    consensus_fasta_index_threads: int = typer.Option(2, "--consensus-fasta-index-threads", help="threads for consensus_fasta_index"),
    consensus_fasta_index_memory: int = typer.Option(4000, "--consensus-fasta-index-memory", help="memory (MB) for consensus_fasta_index"),
    alignment_threads: int = typer.Option(8, "--alignment-threads", help="threads for alignment"),
    alignment_memory: int = typer.Option(16000, "--alignment-memory", help="memory (MB) for alignment"),
    qualimap_threads: int = typer.Option(8, "--qualimap-threads", help="threads for qualimap"),
    qualimap_memory: int = typer.Option(16000, "--qualimap-memory", help="memory (MB) for qualimap"),
    samtools_depth_threads: int = typer.Option(4, "--samtools-depth-threads", help="threads for samtools_depth"),
    samtools_depth_memory: int = typer.Option(2000, "--samtools-depth-memory", help="memory (MB) for samtools_depth"),
    calculate_mean_and_median_depth_memory: int = typer.Option(16000, "--calculate-mean-and-median-depth-memory", help="memory (MB) for calculate_mean_and_median_depth"),
    breadth_of_coverage_memory: int = typer.Option(2000, "--breadth-of-coverage-memory", help="memory (MB) for breadth_of_coverage"),
    count_bam_threads: int = typer.Option(8, "--count-bam-threads", help="threads for count_bam"),
    count_bam_memory: int = typer.Option(2000, "--count-bam-memory", help="memory (MB) for count_bam"),
    fold_80_memory: int = typer.Option(16000, "--fold-80-memory", help="memory (MB) for fold_80"),
    deep_variant_threads: int = typer.Option(8, "--deep-variant-threads", help="threads for deep_variant"),
    deep_variant_memory: int = typer.Option(4000, "--deep-variant-memory", help="memory (MB) for deep_variant"),
    variant_filtering_memory: int = typer.Option(2000, "--variant-filtering-memory", help="memory (MB) for variant_filtering"),
    unique_SNPs_memory: int = typer.Option(4000, "--unique-SNPs-memory", help="memory (MB) for unique_SNPs"),
    bed_maker_memory: int = typer.Option(2000, "--bed-maker-memory", help="memory (MB) for bed_maker"),
    liftover_memory: int = typer.Option(4000, "--liftover-memory", help="memory (MB) for liftover"),
    chromosome_name_conversion_memory: int = typer.Option(4000, "--chromosome-name-conversion-memory", help="memory (MB) for chromosome_name_conversion"),
    filter_fingerprints_memory: int = typer.Option(4000, "--filter-fingerprints-memory", help="memory (MB) for filter_fingerprints"),
    validation_fingerprints_alt_memory: int = typer.Option(4000, "--validation-fingerprints-alt-memory", help="memory (MB) for validation_fingerprints_alt"),
    validation_fingerprints_ref_memory: int = typer.Option(4000, "--validation-fingerprints-ref-memory", help="memory (MB) for validation_fingerprints_ref"),
    barcodes_positions_new_ref_memory: int = typer.Option(4000, "--barcodes-positions-new-ref-memory", help="memory (MB) for barcodes_positions_new_ref"),
    chromosome_arms_memory: int = typer.Option(4000, "--chromosome-arms-memory", help="memory (MB) for chromosome_arms"),

):
    """
    Running the complete Illumina pipeline from FASTQ files to fingerprint detection.
    """
    validate_illumina_inputs(
        reference=reference,
        barcodes=barcodes,
        ref_vcf=ref_vcf,
        fastq_dir=fastq_dir,
        chain_file=chain_file,
        chromosome_mapping=chromosome_mapping,
        centromere=centromere,
    )
    overrides = {
        "paths": set_paths(
            reference,
            barcodes,
            ref_vcf,
            output_dir,
            fastq_dir,
            chain_file,
            chromosome_mapping,
            centromere,
        ),
        "parameters": {
            "fastp": {
                "window_size": int(fastp_window_size),
                "quality": int(fastp_quality),
                "min_len": int(fastp_min_len),
            },
            "breadth_of_coverage": {"depth_cutoff": int(breadth_of_coverage_depth_cutoff)},
            "deep_variant": {"BIN_VERSION": str(deep_variant_BIN_VERSION)},
            "variant_filtering": {
                "min_QUAL": int(variant_filtering_min_QUAL),
                "min_DP": int(variant_filtering_min_DP),
            },
        },
        "resources": set_resources({
            "fastqc_raw": {"threads": fastqc_raw_threads, "memory": fastqc_raw_memory},
            "trim_reads": {"threads": trim_reads_threads, "memory": trim_reads_memory},
            "fastqc_trimmed": {"threads": fastqc_trimmed_threads, "memory": fastqc_trimmed_memory},
            "consensus_fasta_index": {"threads": consensus_fasta_index_threads, "memory": consensus_fasta_index_memory},
            "alignment": {"threads": alignment_threads, "memory": alignment_memory},
            "qualimap": {"threads": qualimap_threads, "memory": qualimap_memory},
            "samtools_depth": {"threads": samtools_depth_threads, "memory": samtools_depth_memory},
            "calculate_mean_and_median_depth": {"memory": calculate_mean_and_median_depth_memory},
            "breadth_of_coverage": {"memory": breadth_of_coverage_memory},
            "count_bam": {"threads": count_bam_threads, "memory": count_bam_memory},
            "fold_80": {"memory": fold_80_memory},
            "deep_variant": {"threads": deep_variant_threads, "memory": deep_variant_memory},
            "variant_filtering": {"memory": variant_filtering_memory},
            "unique_SNPs": {"memory": unique_SNPs_memory},
            "bed_maker": {"memory": bed_maker_memory},
            "liftover": {"memory": liftover_memory},
            "chromosome_name_conversion": {"memory": chromosome_name_conversion_memory},
            "filter_fingerprints": {"memory": filter_fingerprints_memory},
            "validation_fingerprints_alt": {"memory": validation_fingerprints_alt_memory},
            "validation_fingerprints_ref": {"memory": validation_fingerprints_ref_memory},
            "barcodes_positions_new_ref": {"memory": barcodes_positions_new_ref_memory},
            "chromosome_arms": {"memory": chromosome_arms_memory},
        }),
    }
    config_path = load_and_write_config("Illumina_full_pipeline_config_template.yaml", overrides, output_dir)
    os.environ["TMPDIR"] = str(Path(config_path).parent / "tmp")
    snakefiles = [
        "illumina_qc_mapping_variant_calling.smk",
        "mapq.smk",
        "fingerprint_detection.smk"
    ]
    base = Path(__file__).parent / "snakefiles"
    for sf in snakefiles:
        run_snakemake(base / sf, config_path, cores, mem_mb, dry_run)

""" 
-----------------------------
ONT FULL PIPELINE SUBCOMMAND
-----------------------------
"""

@app.command("ONT-complete-run")
def full_pipeline_ONT(
    reference: Path = typer.Option(..., "--reference", "-r", help="Reference genome FASTA"),
    barcodes: Path = typer.Option(..., "--barcodes", "-b", help="tsv file of the fingerprint (2‑SNV barcodes)"),
    ref_vcf: Path = typer.Option(..., "--ref-vcf", "-v", help="VCF used to generate the fingerprint"),
    output_dir: Path = typer.Option(..., "--output-dir", "-o", help="Output directory"),
    fastq_dir: Path = typer.Option(..., "--fastq-dir", "-f", help="Directory with FASTQ files (FASTQ files must follow the naming convention: SampleName.fastq.gz)"),
    dry_run: bool = typer.Option(False, "--dry-run", help="Generate the config file and perform a dry-run without executing rules"),
    chain_file: Path = typer.Option(get_default_chain_file_path(), "--chain-file", "-C", help="Chain file for liftover"),
    chromosome_mapping: Path = typer.Option(get_default_chromosome_mapping_path(), "--chromosome-mapping", "-M", help="Chromosome name mapping"),
    centromere: Path = typer.Option(get_default_centromere_path(), "--centromere", "-c", help="Centromere TSV file"),
    cores: int = typer.Option(8, "--threads", "-t", help="Number of CPU cores"),
    mem_mb: int = typer.Option(16000, "--mem-mb", "-m", help="Memory in MB"),

    chopper_quality: int = typer.Option(10, "--chopper-quality", help="quality cutoff for chopper"),
    chopper_min_len: int = typer.Option(150, "--chopper-min-len", help="minimum length for chopper"),
    breadth_of_coverage_depth_cutoff: int = typer.Option(20, "--breadth-of-coverage-depth-cutoff", help="depth cutoff for calculating the breadth of coverage"),
    breadth_of_coverage_merged_depth_cutoff: int = typer.Option(20, "--breadth-of-coverage-merged-depth-cutoff", help="depth cutoff for calculating the breadth of coverage for merged BAM"),
    deep_variant_BIN_VERSION: str = typer.Option("1.8.0", "--deep-variant-bin-version", help="BIN VERSION for DeepVariant"),
    variant_filtering_min_QUAL: int = typer.Option(10, "--variant-filtering-min-QUAL", help="minimum QUAL for variant filtering"),
    variant_filtering_min_DP: int = typer.Option(5, "--variant-filtering-min-DP", help="minimum DP for variant filtering"),

    nanoplot_raw_threads: int = typer.Option(4, "--nanoplot-raw-threads", help="threads for nanoplot_raw"),
    nanoplot_raw_memory: int = typer.Option(2000, "--nanoplot-raw-memory", help="memory (MB) for nanoplot_raw"),
    extract_headers_threads: int = typer.Option(2, "--extract-headers-threads", help="threads for extract_headers"),
    extract_headers_memory: int = typer.Option(2000, "--extract-headers-memory", help="memory (MB) for extract_headers"),
    filter_reads_threads: int = typer.Option(4, "--filter-reads-threads", help="threads for filter_reads"),
    filter_reads_memory: int = typer.Option(2000, "--filter-reads-memory", help="memory (MB) for filter_reads"),
    trim_reads_threads: int = typer.Option(4, "--trim-reads-threads", help="threads for trim_reads"),
    trim_reads_memory: int = typer.Option(2000, "--trim-reads-memory", help="memory (MB) for trim_reads"),
    nanoplot_after_filtering_threads: int = typer.Option(4, "--nanoplot-after-filtering-threads", help="threads for nanoplot_after_filtering"),
    nanoplot_after_filtering_memory: int = typer.Option(2000, "--nanoplot-after-filtering-memory", help="memory (MB) for nanoplot_after_filtering"),
    index_reference_threads: int = typer.Option(4, "--index-reference-threads", help="threads for index_reference"),
    index_reference_memory: int = typer.Option(5000, "--index-reference-memory", help="memory (MB) for index_reference"),
    alignment_threads: int = typer.Option(8, "--alignment-threads", help="threads for alignment"),
    alignment_memory: int = typer.Option(16000, "--alignment-memory", help="memory (MB) for alignment"),
    qualimap_threads: int = typer.Option(8, "--qualimap-threads", help="threads for qualimap"),
    qualimap_memory: int = typer.Option(16000, "--qualimap-memory", help="memory (MB) for qualimap"),
    samtools_depth_threads: int = typer.Option(2, "--samtools-depth-threads", help="threads for samtools_depth"),
    samtools_depth_memory: int = typer.Option(2000, "--samtools-depth-memory", help="memory (MB) for samtools_depth"),
    calculate_mean_and_median_depth_memory: int = typer.Option(16000, "--calculate-mean-and-median-depth-memory", help="memory (MB) for calculate_mean_and_median_depth"),
    breadth_of_coverage_memory: int = typer.Option(2000, "--breadth-of-coverage-memory", help="memory (MB) for breadth_of_coverage"),
    count_bam_threads: int = typer.Option(8, "--count-bam-threads", help="threads for count_bam"),
    count_bam_memory: int = typer.Option(2000, "--count-bam-memory-memory", help="memory (MB) for count_bam"),
    fold_80_memory: int = typer.Option(16000, "--fold-80-memory", help="memory (MB) for fold_80"),
    merge_bams_threads: int = typer.Option(4, "--merge-bams-threads", help="threads for merge_bams"),
    merge_bams_memory: int = typer.Option(2000, "--merge-bams-memory", help="memory (MB) for merge_bams"),
    qualimap_merged_threads: int = typer.Option(8, "--qualimap-merged-threads", help="threads for qualimap_merged"),
    qualimap_merged_memory: int = typer.Option(16000, "--qualimap-merged-memory", help="memory (MB) for qualimap_merged"),
    samtools_depth_merged_threads: int = typer.Option(4, "--samtools-depth-merged-threads", help="threads for samtools_depth_merged"),
    samtools_depth_merged_memory: int = typer.Option(2000, "--samtools-depth-merged-memory", help="memory (MB) for samtools_depth_merged"),
    calculate_mean_and_median_depth_merged_memory: int = typer.Option(16000, "--calculate-mean-and-median-depth-merged-memory", help="memory (MB) for calculate_mean_and_median_depth_merged"),
    breadth_of_coverage_merged_memory: int = typer.Option(2000, "--breadth-of-coverage-merged-memory", help="memory (MB) for breadth_of_coverage_merged"),
    count_bam_merged_threads: int = typer.Option(8, "--count-bam-merged-threads", help="threads for count_bam_merged"),
    count_bam_merged_memory: int = typer.Option(2000, "--count-bam-merged-memory", help="memory (MB) for count_bam_merged"),
    fold_80_merged_memory: int = typer.Option(16000, "--fold-80-merged-memory", help="memory (MB) for fold_80_merged"),
    deep_variant_threads: int = typer.Option(8, "--deep-variant-threads", help="threads for deep_variant"),
    deep_variant_memory: int = typer.Option(16000, "--deep-variant-memory", help="memory (MB) for deep_variant"),
    variant_filtering_memory: int = typer.Option(2000, "--variant-filtering-memory", help="memory (MB) for variant_filtering"),
    unique_SNPs_memory: int = typer.Option(4000, "--unique-SNPs-memory", help="memory (MB) for unique_SNPs"),
    bed_maker_memory: int = typer.Option(2000, "--bed-maker-memory", help="memory (MB) for bed_maker"),
    liftover_memory: int = typer.Option(8000, "--liftover-memory", help="memory (MB) for liftover"),
    chromosome_name_conversion_memory: int = typer.Option(4000, "--chromosome-name-conversion-memory", help="memory (MB) for chromosome_name_conversion"),
    filter_fingerprints_memory: int = typer.Option(8000, "--filter-fingerprints-memory", help="memory (MB) for filter_fingerprints"),
    validation_fingerprints_alt_memory: int = typer.Option(4000, "--validation-fingerprints-alt-memory", help="memory (MB) for validation_fingerprints_alt"),
    validation_fingerprints_ref_memory: int = typer.Option(4000, "--validation-fingerprints-ref-memory", help="memory (MB) for validation_fingerprints_ref"),
    barcodes_positions_new_ref_memory: int = typer.Option(4000, "--barcodes-positions-new-ref-memory", help="memory (MB) for barcodes_positions_new_ref"),
    chromosome_arms_memory: int = typer.Option(4000, "--chromosome-arms-memory", help="memory (MB) for chromosome_arms"),

):

    """
    Running the complete ONT pipeline from FASTQ files to fingerprint detection.
    """
    validate_ont_inputs(
        reference=reference,
        barcodes=barcodes,
        ref_vcf=ref_vcf,
        fastq_dir=fastq_dir,
        chain_file=chain_file,
        chromosome_mapping=chromosome_mapping,
        centromere=centromere,
    )
    overrides = {
        "paths": set_paths(
            reference,
            barcodes,
            ref_vcf,
            output_dir,
            fastq_dir,
            chain_file,
            chromosome_mapping,
            centromere,
        ),
        "parameters": {
            "chopper": {
                "quality": int(chopper_quality),
                "min_len": int(chopper_min_len),
            },
            "breadth_of_coverage": {
                "depth_cutoff": int(breadth_of_coverage_depth_cutoff),
            },
            "breadth_of_coverage_merged": {
                "depth_cutoff": int(breadth_of_coverage_merged_depth_cutoff),
            },
            "deep_variant": {
                "BIN_VERSION": str(deep_variant_BIN_VERSION),
            },
            "variant_filtering": {
                "min_QUAL": int(variant_filtering_min_QUAL),
                "min_DP": int(variant_filtering_min_DP),
            },
        },
        "resources": set_resources({
            "nanoplot_raw": {"threads": nanoplot_raw_threads, "memory": nanoplot_raw_memory},
            "extract_headers": {"threads": extract_headers_threads, "memory": extract_headers_memory},
            "filter_reads": {"threads": filter_reads_threads, "memory": filter_reads_memory},
            "trim_reads": {"threads": trim_reads_threads, "memory": trim_reads_memory},
            "nanoplot_after_filtering": {"threads": nanoplot_after_filtering_threads,"memory": nanoplot_after_filtering_memory},
            "index_reference": {"threads": index_reference_threads, "memory": index_reference_memory},
            "alignment": {"threads": alignment_threads, "memory": alignment_memory},
            "qualimap": {"threads": qualimap_threads, "memory": qualimap_memory},
            "samtools_depth": {"threads": samtools_depth_threads, "memory": samtools_depth_memory},
            "calculate_mean_and_median_depth": {"memory": calculate_mean_and_median_depth_memory},
            "breadth_of_coverage": {"memory": breadth_of_coverage_memory},
            "count_bam": {"threads": count_bam_threads, "memory": count_bam_memory},
            "fold_80": {"memory": fold_80_memory},
            "merge_bams": {"threads": merge_bams_threads, "memory": merge_bams_memory},
            "qualimap_merged": {"threads": qualimap_merged_threads, "memory": qualimap_merged_memory},
            "samtools_depth_merged": {"threads": samtools_depth_merged_threads, "memory": samtools_depth_merged_memory},
            "calculate_mean_and_median_depth_merged": {"memory": calculate_mean_and_median_depth_merged_memory},
            "breadth_of_coverage_merged": {"memory": breadth_of_coverage_merged_memory},
            "count_bam_merged": {"threads": count_bam_merged_threads, "memory": count_bam_merged_memory},
            "fold_80_merged": {"memory": fold_80_merged_memory},
            "deep_variant": {"threads": deep_variant_threads, "memory": deep_variant_memory},
            "variant_filtering": {"memory": variant_filtering_memory},
            "unique_SNPs": {"memory": unique_SNPs_memory},
            "bed_maker": {"memory": bed_maker_memory},
            "liftover": {"memory": liftover_memory},
            "chromosome_name_conversion": {"memory": chromosome_name_conversion_memory},
            "filter_fingerprints": {"memory": filter_fingerprints_memory},
            "validation_fingerprints_alt": {"memory": validation_fingerprints_alt_memory},
            "validation_fingerprints_ref": {"memory": validation_fingerprints_ref_memory},
            "barcodes_positions_new_ref": {"memory": barcodes_positions_new_ref_memory},
            "chromosome_arms": {"memory": chromosome_arms_memory},
        }),
    }

    config_path = load_and_write_config("ONT_full_pipeline_config_template.yaml", overrides, output_dir)
    os.environ["TMPDIR"] = str(Path(config_path).parent / "tmp")
    snakefiles = [
        "ont_qc_mapping_variant_calling.smk",
        "mapq.smk",
        "fingerprint_detection.smk"
    ]
    base = Path(__file__).parent / "snakefiles"
    for sf in snakefiles:
        run_snakemake(base / sf, config_path, cores, mem_mb, dry_run)

"""
-----------------------------
ONLY DETECTION SUBCOMMAND
-----------------------------
"""

@app.command("scan-barcodes")
def only_detection(
    barcodes: Path = typer.Option(get_default_Nip_barcodes_path(), "--barcodes", "-b", help="tsv file of the fingerprint (2‑SNV barcodes)"),
    ref_vcf: Path = typer.Option(get_default_Nip_ref_vcf_path(), "--ref-vcf", "-v", help="VCF used to generate the fingerprint"),
    output_dir: Path = typer.Option(..., "--output-dir", "-o", help="Output directory"),
    bam: Path = typer.Option(..., "--bam", "-B", help="BAM file for sample to check the fingerprint (must be indexed)"),
    query_vcf: Path = typer.Option(..., "--query-vcf", "-q", help="VCF file for sample to check the fingerprint (must be indexed)"),
    dry_run: bool = typer.Option(False, "--dry-run", help="Generate the config file and perform a dry-run without executing rules"),
    chain_file: Path = typer.Option(get_default_chain_file_path(), "--chain-file", "-C", help="Chain file for liftover"),
    chromosome_mapping: Path = typer.Option(get_default_chromosome_mapping_path(), "--chromosome-mapping", "-M", help="Chromosome name mapping"),
    centromere: Path = typer.Option(get_default_centromere_path(), "--centromere", "-c", help="Centromere TSV file"),
    cores: int = typer.Option(8, "--threads", "-t", help="Number of CPU cores"),
    mem_mb: int = typer.Option(16000, "--mem-mb", "-m", help="Memory in MB"),

    unique_SNPs_memory: int = typer.Option(4000, "--unique-SNPs-memory", help="memory (MB) for unique_SNPs"),
    bed_maker_memory: int = typer.Option(2000, "--bed-maker-memory", help="memory (MB) for bed_maker"),
    liftover_memory: int = typer.Option(4000, "--liftover-memory", help="memory (MB) for liftover"),
    chromosome_name_conversion_memory: int = typer.Option(4000, "--chromosome-name-conversion-memory", help="memory (MB) for chromosome_name_conversion"),
    filter_fingerprints_memory: int = typer.Option(4000, "--filter-fingerprints-memory", help="memory (MB) for filter_fingerprints"),
    validation_fingerprints_alt_memory: int = typer.Option(4000, "--validation-fingerprints-alt-memory", help="memory (MB) for validation_fingerprints_alt"),
    validation_fingerprints_ref_memory: int = typer.Option(4000, "--validation-fingerprints-ref-memory", help="memory (MB) for validation_fingerprints_ref"),
    barcodes_positions_new_ref_memory: int = typer.Option(4000, "--barcodes-positions-new-ref-memory", help="memory (MB) for barcodes_positions_new_ref"),
    chromosome_arms_memory: int = typer.Option(4000, "--chromosome-arms-memory", help="memory (MB) for chromosome_arms"),
):
    """
    Running only the fingerprint detection from provided BAM and VCF files.
    """
    validate_detection_inputs(
        barcodes=barcodes,
        ref_vcf=ref_vcf,
        bam=bam,
        query_vcf=query_vcf,
        chain_file=chain_file,
        chromosome_mapping=chromosome_mapping,
        centromere=centromere,
    )
    overrides = {
        "paths": {
            "barcodes": str(barcodes),
            "ref_vcf": str(ref_vcf),
            "output_dir": str(output_dir),
            "bam": str(bam),
            "query_vcf": str(query_vcf),
            "chain_file": str(chain_file),
            "chromosome_mapping": str(chromosome_mapping),
            "centromere": str(centromere),
        },
        "resources": set_resources({
            "unique_SNPs": {"memory": unique_SNPs_memory},
            "bed_maker": {"memory": bed_maker_memory},
            "liftover": {"memory": liftover_memory},
            "chromosome_name_conversion": {"memory": chromosome_name_conversion_memory},
            "filter_fingerprints": {"memory": filter_fingerprints_memory},
            "validation_fingerprints_alt": {"memory": validation_fingerprints_alt_memory},
            "validation_fingerprints_ref": {"memory": validation_fingerprints_ref_memory},
            "barcodes_positions_new_ref": {"memory": barcodes_positions_new_ref_memory},
            "chromosome_arms": {"memory": chromosome_arms_memory},
        }),
    }
    config_path = load_and_write_config("only_detection_config_template.yaml", overrides, output_dir)
    os.environ["TMPDIR"] = str(Path(config_path).parent / "tmp")
    snakefile = Path(__file__).parent / "snakefiles/fingerprint_detection.smk"
    run_snakemake(snakefile, config_path, cores, mem_mb, dry_run)

"""
-----------------------------
Chain file generation SUBCOMMAND
-----------------------------
"""

@app.command("generate-chain")
def chain_file(
    output_dir: Path = typer.Option(..., "--output-dir", "-o", help="Output directory"),
    old_ref: Path = typer.Option(..., "--old-ref", "-O", help="The reference genome FASTA to convert from"),
    new_ref: Path = typer.Option(..., "--new-ref", "-N", help="The reference genome FASTA to convert to"),
    dry_run: bool = typer.Option(False, "--dry-run", help="Generate the config file and perform a dry-run without executing rules"),
    cores: int = typer.Option(8, "--threads", "-t", help="Number of CPU cores"),
    mem_mb: int = typer.Option(16000, "--mem-mb", "-m", help="Memory in MB"),

    cleaning_fasta_memory: int = typer.Option(4000, "--cleaning-fasta-memory", help="memory (MB) for cleaning_fasta"),
    making_PAF_threads: int = typer.Option(4, "--making-PAF-threads", help="threads for making_PAF"),
    making_PAF_memory: int = typer.Option(16000, "--making-PAF-memory", help="memory (MB) for making_PAF"),
    making_chain_file_memory: int = typer.Option(16000, "--making-chain-file-memory", help="memory (MB) for making_chain_file")
):
    """
    Generate a chain file to lift over genomic coordinates between two reference genomes.
    """
    validate_chain_inputs(
        old_ref=old_ref,
        new_ref=new_ref,
    )
    overrides = {
        "paths": {
            "output_dir": str(output_dir),
            "old_ref": str(old_ref),
            "new_ref": str(new_ref),
        },
        "resources": set_resources({
            "cleaning_fasta": {"memory": cleaning_fasta_memory},
            "making_PAF": {"threads": making_PAF_threads, "memory": making_PAF_memory},
            "making_chain_file": {"memory": making_chain_file_memory},
        }),
    }
    config_path = load_and_write_config("generate_chain_file_config_template.yaml", overrides, output_dir)
    os.environ["TMPDIR"] = str(Path(config_path).parent / "tmp")
    snakefile = Path(__file__).parent / "snakefiles/make_chain_file.smk"
    run_snakemake(snakefile, config_path, cores, mem_mb, dry_run)

