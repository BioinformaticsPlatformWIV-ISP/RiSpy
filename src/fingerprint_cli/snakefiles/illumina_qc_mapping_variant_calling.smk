from pathlib import Path

root = Path(config["paths"]["output_dir"])

dir_samples = Path(config["paths"]["fastq_dir"])

ref = config["paths"]["reference"]

fq_by_sample = {}
for file in dir_samples.iterdir():
    name = file.name.split('_')[0]
    if name not in fq_by_sample.keys():
        fq_by_sample[name] = [next(dir_samples.glob(name + '_R1.fastq.gz')), next(dir_samples.glob(name + '_R2.fastq.gz'))]

rule all:
    input:
        FASTQC_Trimmed_1P = expand(str(root / '{sample}' / 'fastqc_after_trimming' / 'trimmed_reads_1P_fastqc.html'),sample=fq_by_sample.keys()),
        FASTQC_Raw_F = expand(str(root / '{sample}' / 'fastqc_before_trimming' / '{sample}_R1_fastqc.html'),sample=fq_by_sample.keys()),
        Filtered = expand(str(root / '{sample}' / 'deep_variant' / 'Filtering' / 'filtered_variants.vcf.gz'),sample=fq_by_sample.keys()),

rule fastqc_raw:
    """
    Quality control of raw reads
    """
    input:
        FQ_fwd = lambda wildcards: str(fq_by_sample[wildcards.sample][0]),
        FQ_rev = lambda wildcards: str(fq_by_sample[wildcards.sample][1])
    output:
        FASTQC_Raw_F = root / '{sample}' / 'fastqc_before_trimming' / '{sample}_R1_fastqc.html',
        FASTQC_Raw_R = root / '{sample}' / 'fastqc_before_trimming' / '{sample}_R2_fastqc.html'
    params:
        outdir= lambda wildcards: root / wildcards.sample / 'fastqc_before_trimming'
    threads: config["resources"]["fastqc_raw"]["threads"]
    resources:
        mem_mb = config["resources"]["fastqc_raw"]["memory"]
    benchmark: root / '{sample}' / 'benchmarks' / 'fastqc_raw.txt'
    shell:
        """
        fastqc -t {threads} {input.FQ_fwd} {input.FQ_rev} -o {params.outdir}
        """

rule trim_reads:
    """
    Trimming of paired end reads
    """
    input:
        FQ_fwd = lambda wildcards: str(fq_by_sample[wildcards.sample][0]),
        FQ_rev = lambda wildcards: str(fq_by_sample[wildcards.sample][1])
    output:
        FQ_1P = root / '{sample}' / 'trimming' / 'trimmed_reads_1P.fastq.gz',
        FQ_2P = root / '{sample}' / 'trimming' / 'trimmed_reads_2P.fastq.gz',
    params:
        window_size = config["parameters"]["fastp"]["window_size"],
        quality = config["parameters"]["fastp"]["quality"],
        min_len = config["parameters"]["fastp"]["min_len"],
        json = lambda wildcards: root / wildcards.sample / 'trimming' / 'fastp.json',
        html = lambda wildcards: root / wildcards.sample / 'trimming' / 'fastp.html'
    threads: config["resources"]["trim_reads"]["threads"]
    resources:
        mem_mb = config["resources"]["trim_reads"]["memory"]
    benchmark: root / '{sample}' / 'benchmarks' / 'trim_reads.txt'
    shell:
        """
        fastp --cut_right --cut_right_window_size {params.window_size} --cut_right_mean_quality {params.quality} \
        --length_required {params.min_len} --detect_adapter_for_pe --thread {threads} -j {params.json} -h {params.html} \
        -i {input.FQ_fwd} -I {input.FQ_rev} -o {output.FQ_1P} -O {output.FQ_2P}
        """

rule fastqc_trimmed:
    """
    Quality control of trimmed reads
    """
    input:
        FQ_1P = rules.trim_reads.output.FQ_1P,
        FQ_2P = rules.trim_reads.output.FQ_2P,
    output:
        FASTQC_Trimmed_1P = root / '{sample}' / 'fastqc_after_trimming' / 'trimmed_reads_1P_fastqc.html',
        FASTQC_Trimmed_2P = root / '{sample}' / 'fastqc_after_trimming' / 'trimmed_reads_2P_fastqc.html',
    params:
        outdir = lambda wildcards: root / wildcards.sample / 'fastqc_after_trimming',
    threads: config["resources"]["fastqc_trimmed"]["threads"]
    resources:
        mem_mb = config["resources"]["fastqc_trimmed"]["memory"]
    benchmark: root / '{sample}' / 'benchmarks' / 'fastqc_trimmed.txt'
    shell:
        """
        fastqc -t {threads} {input.FQ_1P} {input.FQ_2P} -o {params.outdir}
        """

rule consensus_fasta_index:
    """
    Indexing the reference fasta
    """
    input:
        FASTA = ref,
    output:
        REFERENCE = root / 'reference' / 'reference.fasta'
    threads: config["resources"]["consensus_fasta_index"]["threads"]
    resources:
        mem_mb = config["resources"]["consensus_fasta_index"]["memory"]
    benchmark: root / 'reference' / 'benchmarks' / 'consensus_fasta_index.txt'
    shell:
        """
        cp {input.FASTA} {output.REFERENCE}
        bwa index {output.REFERENCE}
        samtools faidx -@ {threads} {output.REFERENCE}
        """

rule alignment:
    """
    Mapping the short reads to reference and convert to bam, then sorting and indexing the bam file
    """
    input:
        FASTA = rules.consensus_fasta_index.output.REFERENCE,
        FQ_1P = rules.trim_reads.output.FQ_1P,
        FQ_2P = rules.trim_reads.output.FQ_2P,
    output:
        BAM = root / '{sample}' / 'mapping' / 'mapped_reads_sorted.bam'
    threads: config["resources"]["alignment"]["threads"]
    params:
        outdir = root
    resources:
        mem_mb = config["resources"]["alignment"]["memory"]
    benchmark: root / '{sample}' / 'benchmarks' / 'alignment.txt'
    shell:
        """
        bwa mem -K 100000000 -v 3 -t {threads} -k 19 -r 1.5 -M {input.FASTA} {input.FQ_1P} {input.FQ_2P} \
        | samtools view -b - -@ {threads} | samtools sort - -o {output.BAM} -T {params.outdir} -@ {threads}
        samtools index -@ {threads} {output.BAM}
        """

rule deep_variant:
    """
    Performing germline variant calling
    """
    input:
        BAM = rules.alignment.output.BAM
    output:
        VCF = root  / '{sample}' / 'deep_variant' / 'variants.vcf.gz',
        GVCF = root / '{sample}' / 'deep_variant' / 'variants.g.vcf.gz',
    params:
        BIN_VERSION = config["parameters"]["deep_variant"]["BIN_VERSION"],
        input_dir = lambda wildcards: root / wildcards.sample / 'mapping',
        output_dir = lambda wildcards: root / wildcards.sample / 'deep_variant',
        reference = root / 'reference',
        tmpdir = lambda wildcards: root / wildcards.sample / 'tmp'
    threads: config["resources"]["deep_variant"]["threads"]
    resources:
        mem_mb = config["resources"]["deep_variant"]["memory"]
    benchmark: root / '{sample}' / 'deep_variant' / 'benchmarks' / 'deep_variant.txt'
    shell:
        """
        mkdir {params.tmpdir}
        docker run -u $(id -u):$(id -g) --cpus={threads} --rm \
          -v {params.input_dir}:"/input" \
          -v {params.output_dir}:"/output" \
          -v {params.reference}:"/reference" \
          -v {params.tmpdir}:"/tmp_dir" \
          google/deepvariant:{params.BIN_VERSION} \
          /opt/deepvariant/bin/run_deepvariant \
          --model_type=WGS \
          --ref=/reference/reference.fasta \
          --reads=/input/mapped_reads_sorted.bam \
          --output_vcf=/output/variants.vcf.gz \
          --output_gvcf=/output/variants.g.vcf.gz \
          --num_shards={threads} \
          --intermediate_results_dir /tmp_dir \
          --vcf_stats_report=true
        rm -r {params.tmpdir}
        """

rule variant_filtering:
    """
    Filtering the variants with quality cutoff of 20 and depth cutoff of 5, and getting the stats
    """
    input:
        VCF_gz = rules.deep_variant.output.VCF
    output:
        STATS = root / '{sample}' / 'deep_variant' / 'stats' / 'filtered_variants.stats.vchk',
        Filtered = root / '{sample}' / 'deep_variant' / 'Filtering' / 'filtered_variants.vcf.gz',
    params:
        VCF = lambda wildcards: root / wildcards.sample / 'deep_variant' / 'variants.vcf',
        STATS_RAW = lambda wildcards: root / wildcards.sample / 'deep_variant' / 'stats' / 'raw_variants.stats.vchk',
        QUAL = config["parameters"]["variant_filtering"]["min_QUAL"],
        DP = config["parameters"]["variant_filtering"]["min_DP"],
    threads: config["resources"]["variant_filtering"]["threads"]
    resources:
        mem_mb = config["resources"]["variant_filtering"]["memory"]
    benchmark: root / '{sample}' / 'deep_variant' / 'benchmarks' / 'variant_filtering.txt'
    shell:
        """
        gunzip -k {input.VCF_gz}
        bcftools stats {params.VCF} > {params.STATS_RAW}
        bcftools filter -i 'QUAL>={params.QUAL} & FORMAT/DP>={params.DP}' {params.VCF} -o {output.Filtered} -O z --write-index="tbi"
        gunzip -k {output.Filtered}
        bcftools stats {output.Filtered} > {output.STATS}
        """

