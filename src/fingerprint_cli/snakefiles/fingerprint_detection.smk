from pathlib import Path
import os
from fingerprint_cli.fingerprint_pipeline import detect_unique_snps, bed_maker, chromosome_name_conversion, validation_fingerprints_alt, validation_fingerprints_ref, barcodes_positions_new_ref, chromosome_arms
from fingerprint_cli.io import build_fq_by_sample, generate_report_html

root = Path(config["paths"]["output_dir"])

dir_samples = Path(config["paths"].get("fastq_dir")) if config["paths"].get("fastq_dir") else None

chain_file_old_nip_to_new_nip_ref = config["paths"]["chain_file"]

chromosome_mapping_file = config["paths"]["chromosome_mapping"]

filtered_barcodes = config["paths"]["barcodes"] #The final set of 2-SNPs barcodes

vcf_file = config["paths"]["ref_vcf"] #The strict filtered vcf file that was used for generating the fingerprint

centromere = config["paths"]["centromere"]

query_vcf = config["paths"].get("query_vcf", None)

bam = config["paths"].get("bam", None)

fq_by_sample = build_fq_by_sample(dir_samples)

config["fq_by_sample"] = fq_by_sample

rule all:
    input:
        HTML = expand(str(root / '{sample}' / 'results.html'),sample=config["fq_by_sample"].keys())


rule unique_SNPs:
    """
    Detecting the unique SNPs in the combinations of 2-SNV barcodes
    """
    input:
        Fingerprints = filtered_barcodes,
    output:
        Unique_SNPs = root / 'Fingerprint_Detection' / 'unique_snps.tsv',
    threads: config["resources"]["unique_SNPs"]["threads"]
    resources:
        mem_mb = config["resources"]["unique_SNPs"]["memory"]
    benchmark: root / 'Fingerprint_Detection' / 'benchmarks' / 'unique_SNPs.txt'
    run:
        detect_unique_snps(input.Fingerprints,output.Unique_SNPs)


rule bed_maker:
    """
    Making a bed file out of the unique SNPs 
    """
    input:
        FINGERPRINTS = rules.unique_SNPs.output.Unique_SNPs,
    output:
        FINGER_BED = root / 'Fingerprint_Detection' / 'fingerprints.bed'
    threads: config["resources"]["bed_maker"]["threads"]
    resources:
        mem_mb = config["resources"]["bed_maker"]["memory"]
    benchmark: root / 'Fingerprint_Detection' / 'benchmarks' / 'bed_maker.txt'
    run:
        bed_maker(input.FINGERPRINTS,output.FINGER_BED)

rule liftover:
    """
    Executing liftOver to get the SNP positions in the new reference genome
    """
    input:
        FINGER_BED = rules.bed_maker.output.FINGER_BED,
        CHAINFILE = chain_file_old_nip_to_new_nip_ref
    output:
        FINGER_BED_LIFTOVER = root / 'Fingerprint_Detection' / 'fingerprints_liftover.bed',
        UNMAPPED = root / 'Fingerprint_Detection' / 'unmapped.bed'
    threads: config["resources"]["liftover"]["threads"]
    resources:
        mem_mb = config["resources"]["liftover"]["memory"]
    benchmark: root / 'Fingerprint_Detection' / 'benchmarks' / 'liftover.txt'
    shell:
        """
        liftOver {input.FINGER_BED} {input.CHAINFILE} {output.FINGER_BED_LIFTOVER} {output.UNMAPPED}
        """

rule chromosome_name_conversion:
    """
    Converting the chromosome names to match the new reference genome
    """
    input:
        CHROM_MAP = chromosome_mapping_file,
        FINGER_BED_LIFTOVER = rules.liftover.output.FINGER_BED_LIFTOVER
    output:
        FINGER_BED_CONVERT_CHR = root / 'Fingerprint_Detection' / 'fingerprints_liftover_converted_chr_names.bed'
    threads: config["resources"]["chromosome_name_conversion"]["threads"]
    resources:
        mem_mb = config["resources"]["chromosome_name_conversion"]["memory"]
    benchmark: root / 'Fingerprint_Detection' / 'benchmarks' / 'chromosome_name_conversion.txt'
    run:
        chromosome_name_conversion(input.CHROM_MAP, input.FINGER_BED_LIFTOVER, output.FINGER_BED_CONVERT_CHR)


rule filter_fingerprints:
    """
    Making two separate bed files for reference and alternative alleles
    Note: The decision of a SNP being reference or alternative allele is made based on the new reference genome 
    """
    input:
        FINGER_BED_CONVERT_CHR = rules.chromosome_name_conversion.output.FINGER_BED_CONVERT_CHR,
        VCF_STRICT = vcf_file,
    output:
        REF_BED = root / 'Fingerprint_Detection' / 'reference_allele_fingerprints.bed',
        ALT_BED = root / 'Fingerprint_Detection' / 'alternative_allele_fingerprints.bed',
    threads: config["resources"]["filter_fingerprints"]["threads"]
    resources:
        mem_mb = config["resources"]["filter_fingerprints"]["memory"]
    benchmark: root / 'Fingerprint_Detection' / 'benchmarks' / 'filter_fingerprints.txt'
    shell:
        """
        bedtools subtract -a {input.FINGER_BED_CONVERT_CHR} -b {input.VCF_STRICT} > {output.REF_BED}
        bedtools intersect -a {input.FINGER_BED_CONVERT_CHR} -b {input.VCF_STRICT} > {output.ALT_BED}
        """

rule validation_fingerprints_alt:
    """
    Checking if we can find the alternative SNPs in the samples
    """
    input:
        FINGER_BED_CONVERT_CHR = rules.filter_fingerprints.output.ALT_BED,
        VCF = lambda wildcards: query_vcf or root / wildcards.sample / "deep_variant" / "Filtering" / "filtered_variants.vcf.gz",
        Barcodes = filtered_barcodes,
    output:
        Detected_snps = root / 'Fingerprint_Detection' / 'Alternative_Alleles' / '{sample}' / 'detected_snps.txt',
        Undetected_snps = root / 'Fingerprint_Detection' / 'Alternative_Alleles' / '{sample}' / 'undetected_snps.txt',
        Detected_Barcodes = root / 'Fingerprint_Detection' / 'Alternative_Alleles' / '{sample}' / 'tmp_detected_barcodes.tsv',
    params:
        config = lambda wildcards: root / wildcards.sample / 'config_results.yaml',
        config_dir = lambda wildcards: root / wildcards.sample
    threads: config["resources"]["validation_fingerprints_alt"]["threads"]
    resources:
        mem_mb = config["resources"]["validation_fingerprints_alt"]["memory"]
    benchmark: root / 'Fingerprint_Detection' / 'Alternative_Alleles' / '{sample}' / 'benchmarks' / 'validation_fingerprints_alt.txt'
    run:
        os.makedirs(str(params.config_dir), exist_ok=True)
        validation_fingerprints_alt(input.FINGER_BED_CONVERT_CHR, input.VCF, input.Barcodes, output.Detected_snps, output.Undetected_snps, output.Detected_Barcodes, params.config)


rule validation_fingerprints_ref:
    """
    Checking if we can find the reference SNPs in the samples
    """
    input:
        FINGER_BED_CONVERT_CHR = rules.filter_fingerprints.output.REF_BED,
        VCF = lambda wildcards: query_vcf or root / wildcards.sample / "deep_variant" / "Filtering" / "filtered_variants.vcf.gz",
        BAM = lambda wildcards: bam or root / wildcards.sample / 'mapping' / 'mapped_reads_sorted.bam',
        Barcodes = rules.validation_fingerprints_alt.output.Detected_Barcodes
    output:
        Detected_snps = root / 'Fingerprint_Detection' / 'Reference_Alleles' / '{sample}' / 'detected_snps.txt',
        Undetected_snps = root / 'Fingerprint_Detection' / 'Reference_Alleles' / '{sample}' / 'undetected_snps.txt',
        Detected_Barcodes = root / 'Fingerprint_Detection' / 'Final_Barcodes' / '{sample}' / 'detected_barcodes.tsv'
    params:
        config = lambda wildcards: root / wildcards.sample / 'config_results.yaml'
    threads: config["resources"]["validation_fingerprints_ref"]["threads"]
    resources:
        mem_mb = config["resources"]["validation_fingerprints_ref"]["memory"]
    benchmark: root / 'Fingerprint_Detection' / 'Reference_Alleles' / '{sample}' / 'benchmarks' / 'validation_fingerprints_ref.txt'
    run:
        validation_fingerprints_ref(input.FINGER_BED_CONVERT_CHR, input.VCF, input.BAM, input.Barcodes, output.Detected_snps, output.Undetected_snps, output.Detected_Barcodes, params.config)


rule barcodes_positions_new_ref:
    """
    Generating the final set of fingerprint that could be detected in the samples with their positions in the new reference genome
    """
    input:
        BARCODES = rules.validation_fingerprints_ref.output.Detected_Barcodes,
        SNPs = rules.liftover.output.FINGER_BED_LIFTOVER
    output:
        BARCODES_LIFTOVER = root / 'Fingerprint_Detection' / 'Final_Barcodes' / '{sample}' / 'detected_barcodes.liftover.tsv'
    threads: config["resources"]["barcodes_positions_new_ref"]["threads"]
    resources:
        mem_mb = config["resources"]["barcodes_positions_new_ref"]["memory"]
    benchmark: root / 'Fingerprint_Detection' / 'Final_Barcodes' / '{sample}' / 'benchmarks' / 'barcodes_positions_new_ref.txt'
    run:
        barcodes_positions_new_ref(input.BARCODES, input.SNPs, output.BARCODES_LIFTOVER)


rule chromosome_arms:
    """
    Checking which chromosome arms the SNPs are located in
    """
    input:
        CENTROMERE = centromere,
        BARCODES = rules.barcodes_positions_new_ref.output.BARCODES_LIFTOVER,
    output:
        OUT_TXT = root / 'Fingerprint_Detection' / 'Final_Barcodes' / '{sample}' / 'chr_arms.txt'
    params:
        config = lambda wildcards: root / wildcards.sample / 'config_results.yaml'
    threads: config["resources"]["chromosome_arms"]["threads"]
    resources:
        mem_mb = config["resources"]["chromosome_arms"]["memory"]
    benchmark: root / 'Fingerprint_Detection' / 'Final_Barcodes' / '{sample}' / 'benchmarks' / 'chromosome_arms.txt'
    run:
        chromosome_arms(input.CENTROMERE, input.BARCODES, output.OUT_TXT, params.config)

rule report_html:
    """
    Generating the final HTML report
    """
    input:
        CENTROMERE = rules.chromosome_arms.output.OUT_TXT, # To ensure this rule runs after chromosome_arms
    output:
        HTML = root / '{sample}' / 'results.html'
    params:
        config = lambda wildcards: root / wildcards.sample / 'config_results.yaml'
    threads: config["resources"]["report_html"]["threads"]
    resources:
        mem_mb = config["resources"]["report_html"]["memory"]
    run:
        generate_report_html(output.HTML, params.config)
