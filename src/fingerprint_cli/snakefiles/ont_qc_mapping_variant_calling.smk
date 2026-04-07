from pathlib import Path

root = Path(config["paths"]["output_dir"])

dir_samples = Path(config["paths"]["fastq_dir"])

ref = config["paths"]["reference"]

fq_by_sample = {}
for file in dir_samples.iterdir():
    name = file.name.split('.')[0]
    fq_by_sample[name] = [next(dir_samples.glob(name + '.fastq.gz'))]

rule all:
    input:
        NanoStats_UNIQUE = expand (str(root / '{sample}' / 'nanoplot_after_filtering' / 'NanoStats.txt'), sample = fq_by_sample.keys()),
        NanoStats = expand (str(root / '{sample}' / 'nanoplot_raw' / 'NanoStats.txt'), sample = fq_by_sample.keys()),
        Filtered = root / 'merged_flowcells' / 'deep_variant' / 'Filtering' / 'filtered_variants.vcf.gz'

rule nanoplot_raw:
    """
    Quality control of raw reads 
    """
    input:
        FQ = lambda wildcards: str(fq_by_sample[wildcards.sample][0]),
    output:
        NanoStats = root / '{sample}' / 'nanoplot_raw' / 'NanoStats.txt',
    params:
        NanoStats_dir = lambda wildcards: root / wildcards.sample / 'nanoplot_raw',
    threads: config["resources"]["nanoplot_raw"]["threads"]
    resources:
        mem_mb = config["resources"]["nanoplot_raw"]["memory"]
    benchmark: root / '{sample}' / 'benchmarks' / 'nanoplot_raw.txt'
    shell:
        """
        NanoPlot --fastq {input.FQ} -t {threads} -o {params.NanoStats_dir}
        """

rule extract_headers:
    """
    Extract the read headers of simplex and duplex reads
    """
    input:
        FQ = lambda wildcards: str(fq_by_sample[wildcards.sample][0]),
    output:
        TSV_ALL = root / '{sample}' / 'read-ids' / 'all.tsv', #tsv of headers of all reads
        TSV_ids = root / '{sample}' / 'read-ids' / 'idstofilter.tsv', #tsv of read identifiers of simplex reads used to generate duplex reads
        TSV_filter = root / '{sample}' / 'read-ids' / 'filter.tsv' #tsv of full headers of simplex reads used to generate duplex reads
    threads: config["resources"]["extract_headers"]["threads"]
    resources:
        mem_mb = config["resources"]["extract_headers"]["memory"]
    benchmark: root / '{sample}' / 'benchmarks' / 'extract_headers.txt'
    shell:
        """
        seqkit fx2tab -n {input.FQ} > {output.TSV_ALL}
        seqkit fx2tab -n {input.FQ} | cut -f 1 | grep ';' | tr ';' '\n' > {output.TSV_ids} || true
        cat {output.TSV_ALL} | grep -f {output.TSV_ids} | grep -v ';' > {output.TSV_filter} || true
        """

rule filter_reads:
    """
    Filter simplex reads that were used to generate duplex reads from fastq 
    """
    input:
        FQ = lambda wildcards: str(fq_by_sample[wildcards.sample][0]),
        TSV = root / '{sample}' / 'read-ids' / 'filter.tsv'
    output:
        UNIQUE = root / '{sample}' / 'unique_fastq' / 'unique.fastq.gz',
    threads: config["resources"]["filter_reads"]["threads"]
    resources:
        mem_mb = config["resources"]["filter_reads"]["memory"]
    benchmark: root / '{sample}' / 'benchmarks' / 'filter_reads.txt'
    shell:
        """
        seqkit grep -n -v -f {input.TSV} {input.FQ} -o {output.UNIQUE} || true
        """

rule trim_reads:
    """
    Trimming of ONT reads
    """
    input:
        FQ = root / '{sample}' / 'unique_fastq' / 'unique.fastq.gz',
    output:
        ONT_FILT = root / '{sample}' / 'chopper' / 'unique_filter.fastq.gz',
    params:
        q = config["parameters"]["chopper"]["quality"],
        l = config["parameters"]["chopper"]["min_len"]
    threads: config["resources"]["trim_reads"]["threads"]
    resources:
        mem_mb = config["resources"]["trim_reads"]["memory"]
    benchmark: root / '{sample}' / 'benchmarks' / 'trim_reads.txt'
    shell:
        """
        gunzip -c {input.FQ} | chopper -q {params.q} -l {params.l} --threads {threads} | gzip > {output.ONT_FILT}
        """

rule nanoplot_after_filtering:
    """
    Quality control of trimmed reads
    """
    input:
        ONT_FILT = root / '{sample}' / 'chopper' / 'unique_filter.fastq.gz',
    output:
        NanoStats_UNIQUE = root / '{sample}' / 'nanoplot_after_filtering' / 'NanoStats.txt',
    params:
        NanoStats_dir = lambda wildcards: root / wildcards.sample / 'nanoplot_after_filtering' ,
    threads: config["resources"]["nanoplot_after_filtering"]["threads"]
    resources:
        mem_mb = config["resources"]["nanoplot_after_filtering"]["memory"]
    benchmark: root / '{sample}' / 'benchmarks' / 'nanoplot_after_filtering.txt'
    shell:
        """
        NanoPlot --fastq {input.ONT_FILT} -t {threads} -o {params.NanoStats_dir}
        """

rule index_reference:
    """
    Indexing the reference fasta
    """
    input:
        FASTA = ref
    output:
        MMI = root / 'reference' / 'refgenomes.mmi'
    params:
        out_ref = root / 'reference' / 'reference.fna',
        minimap2_preset = 'map-ont'
    threads: config["resources"]["index_reference"]["threads"]
    resources:
        mem_mb = config["resources"]["index_reference"]["memory"]
    benchmark: root / 'reference' / 'benchmarks' / 'index_reference.txt'
    shell:
        """
        cp {input.FASTA} {params.out_ref}
        minimap2 -I4G -t {threads} -x {params.minimap2_preset} -d {output.MMI} {params.out_ref}
        samtools faidx {params.out_ref}
        """

rule alignment:
    """
    Mapping the long reads to reference and convert to bam, then sorting and indexing the bam file
    """
    input:
        ONT_FQ = root / '{sample}' / 'chopper' / 'unique_filter.fastq.gz',
        MMI = root / 'reference' / 'refgenomes.mmi'
    output:
        BAM = root / '{sample}' / 'mapping' / 'mapped_reads_sorted.bam',
    params:
        minimap2_preset = 'map-ont',
        tmp_dir= root / 'tmp'
    threads: config["resources"]["alignment"]["threads"]
    resources:
        mem_mb = config["resources"]["alignment"]["memory"]
    benchmark: root / '{sample}' / 'benchmarks' / 'alignment.txt'
    shell:
        """
        minimap2 -a -x {params.minimap2_preset} {input.MMI} {input.ONT_FQ} -t {threads} \
        | samtools view -b - -@ {threads} | samtools sort - -o {output.BAM} -T {params.tmp_dir} -@ {threads}
        samtools index {output.BAM}
        """

rule merge_bams:
    """
    Merging the bam files of all flowcells
    """
    input:
        BAM = expand(str(root / "{sample}" / "mapping" / "mapped_reads_sorted.bam"), sample=fq_by_sample.keys())
    output:
        Merged_BAM = root / 'merged_flowcells' / 'mapping' / 'mapped_reads_sorted.bam'
    threads: config["resources"]["merge_bams"]["threads"]
    resources:
        mem_mb = config['resources']['merge_bams']['memory']
    benchmark: root / 'merged_flowcells' / 'benchmarks' / 'merge_bams.txt'
    shell:
        """
        samtools merge -o {output.Merged_BAM} -@ {threads} {input.BAM}
        samtools index {output.Merged_BAM}
        """

rule deep_variant:
    """
    Performing germline variant calling
    """
    input:
        BAM = rules.merge_bams.output.Merged_BAM,
    output:
        VCF = root / 'merged_flowcells' / 'deep_variant' / 'variants.vcf.gz',
        GVCF = root / 'merged_flowcells' / 'deep_variant' / 'variants.g.vcf.gz',
    params:
        BIN_VERSION = config["parameters"]["deep_variant"]["BIN_VERSION"],
        input_dir = root / 'merged_flowcells'/ 'mapping',
        output_dir = root / 'merged_flowcells' / 'deep_variant',
        reference = root / 'reference',
        tmpdir = root / 'merged_flowcells' / 'deep_variant' / 'tmp'
    threads: config["resources"]["deep_variant"]["threads"]
    resources:
        mem_mb = config['resources']['deep_variant']['memory']
    benchmark: root / 'merged_flowcells' / 'deep_variant' / 'benchmarks' / 'deep_variant.txt'
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
          --model_type=ONT_R104 \
          --ref=/reference/reference.fna \
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
    Filtering the variants with quality cutoff of 10 and depth cutoff of 5, and getting the stats
    """
    input:
        VCF_gz = rules.deep_variant.output.VCF
    output:
        STATS = root / 'merged_flowcells' / 'deep_variant' / 'stats' / 'filtered_variants.stats.vchk',
        Filtered = root / 'merged_flowcells' / 'deep_variant' / 'Filtering' / 'filtered_variants.vcf.gz',
    params:
        VCF = root / 'merged_flowcells' / 'deep_variant' / 'variants.vcf',
        STATS_RAW = root / 'merged_flowcells' / 'deep_variant' / 'stats' / 'raw_variants.stats.vchk',
        QUAL = config["parameters"]["variant_filtering"]["min_QUAL"],
        DP = config["parameters"]["variant_filtering"]["min_DP"],
    threads: config["resources"]["variant_filtering"]["threads"]
    resources:
        mem_mb = config['resources']['variant_filtering']['memory']
    benchmark: root / 'merged_flowcells' / 'deep_variant' / 'benchmarks' / 'variant_filtering.txt'
    shell:
        """
        gunzip -k {input.VCF_gz}
        bcftools stats {params.VCF} > {params.STATS_RAW}
        bcftools filter -i 'QUAL>={params.QUAL} & FORMAT/DP>={params.DP}' {params.VCF} -o {output.Filtered} -O z --write-index="tbi"
        gunzip -k {output.Filtered}
        bcftools stats {output.Filtered} > {output.STATS}
        """


