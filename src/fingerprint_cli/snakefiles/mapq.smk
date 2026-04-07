from pathlib import Path
from fingerprint_cli.qc import calculate_mean_and_median_depth, calculate_ref_length, calculate_breadth_of_coverage, calculate_fold80
from fingerprint_cli.io import build_fq_by_sample

root = Path(config["paths"]["output_dir"])

dir_samples = Path(config["paths"]["fastq_dir"])

ref = config["paths"]["reference"]

fq_by_sample = build_fq_by_sample(dir_samples)

config["fq_by_sample"] = fq_by_sample

ref_len = calculate_ref_length(ref)

rule all:
    input:
        FOLD = expand(str(root / '{sample}' / 'fold_80' / 'fold_80.txt'),sample=config["fq_by_sample"].keys())

rule qualimap:
    """
    Checking the quality of the alignment 
    """
    input:
        BAM = root / '{sample}' / 'mapping' / 'mapped_reads_sorted.bam'
    output:
        STATS = root / '{sample}' / 'qualimap' / 'qualimapReport.html',
        TSV = root / '{sample}' / 'qualimap' / "raw_data_qualimapReport" / "coverage_histogram.txt"
    params:
        java_mem_size = f"{config['resources']['qualimap']['memory']}M",
        output_dir = lambda wildcards: root / wildcards.sample / 'qualimap'
    threads: config["resources"]["qualimap"]["threads"]
    resources:
        mem_mb = config['resources']['qualimap']['memory']
    benchmark: root / '{sample}' / 'benchmarks' / 'qualimap.txt'
    shell:
        """
        qualimap bamqc -bam {input.BAM} -outdir {params.output_dir} -nt {threads} --java-mem-size={params.java_mem_size}
        """

rule samtools_depth:
    """
    Determining the depth of coverage at each genomic position
    """
    input:
        BAM = root / '{sample}' / 'mapping' / 'mapped_reads_sorted.bam'
    output:
        TSV = root / '{sample}' / 'depth' / 'depth.tsv'
    threads: config["resources"]["samtools_depth"]["threads"]
    resources:
        mem_mb = config['resources']['samtools_depth']['memory']
    benchmark: root / '{sample}' / 'benchmarks' / 'samtools_depth.txt'
    shell:
        """
        samtools depth -@ {threads} -H -a {input.BAM} > {output.TSV}; # -a: output all positions, also those with zero depth
        """

rule calculate_mean_and_median_depth:
    """
    Calculating the mean and median depth
    """
    input:
        TSV = rules.samtools_depth.output.TSV,
    output:
        DEPTH = root / '{sample}' / 'depth' / 'mean_and_median_depth.txt'
    params:
        config = lambda wildcards: root / wildcards.sample / 'config_results.yaml'
    threads: config["resources"]["calculate_mean_and_median_depth"]["threads"]
    resources:
        mem_mb = config['resources']['calculate_mean_and_median_depth']['memory']
    benchmark: root / '{sample}' / 'benchmarks' / 'calculate_mean_and_median_depth.txt'
    run:
        calculate_mean_and_median_depth(input.TSV,output.DEPTH,params.config)

rule breadth_of_coverage:
    """
    Calculating % of reference covered by mapped reads, with the depth cutoff set to 20X
    """
    input:
        TSV = rules.samtools_depth.output.TSV,
        DEPTH = rules.calculate_mean_and_median_depth.output.DEPTH  #To make sure this step is done after calculating mean and median depth
    output:
        BREADTH = root / '{sample}' / 'breadth' / 'breadth.txt',
    params:
        length = ref_len,
        depth_cutoff = config["parameters"]["breadth_of_coverage"]["depth_cutoff"],
        config = lambda wildcards: root / wildcards.sample / 'config_results.yaml'
    threads: config["resources"]["breadth_of_coverage"]["threads"]
    resources:
        mem_mb = config['resources']['breadth_of_coverage']['memory']
    benchmark: root / '{sample}' / 'benchmarks' / 'breadth_of_coverage.txt'
    run:
        calculate_breadth_of_coverage(input.TSV,output.BREADTH,params.depth_cutoff,params.length,params.config)


rule count_bam:
    """
    Counting the number of mapped primary reads
    """
    input:
        BAM = root / '{sample}' / 'mapping' / 'mapped_reads_sorted.bam',
        Breadth = rules.breadth_of_coverage.output.BREADTH #To make sure this step is done after calculating breadth of coverage
    output:
        COUNT = root / '{sample}' / 'count' / 'mapped_count.txt',
    threads: config["resources"]["count_bam"]["threads"]
    resources:
        mem_mb = config['resources']['count_bam']['memory']
    benchmark: root / '{sample}' / 'benchmarks' / 'count_bam.txt'
    shell:
        """
        samtools view -@ {threads} -c -F 260 {input.BAM} > {output.COUNT}
        """

rule fold_80:
    """
    Fold80 is a way to calculate the evenness of coverage. It expresses the amount of additional sequencing needed to have
    80% of all targets covered at the currently observed mean. It is computed as the mean coverage divided by the 20th
    percentile.
    """
    input:
        COV = root / '{sample}' / 'qualimap' / "raw_data_qualimapReport" / "coverage_histogram.txt",
        depth = root / '{sample}' / 'depth' / 'depth.tsv',
        COUNT=  rules.count_bam.output.COUNT
    output:
        FOLD = root / '{sample}' / 'fold_80' / 'fold_80.txt',
    params:
        length = ref_len,
        config = lambda wildcards: root / wildcards.sample / 'config_results.yaml'
    threads: config["resources"]["fold_80"]["threads"]
    resources:
        mem_mb = config['resources']['fold_80']['memory']
    benchmark: root / '{sample}' / 'benchmarks' / 'fold_80.txt'
    run:
        calculate_fold80(input.COV, input.depth, input.COUNT, params.length, output.FOLD, params.config)
