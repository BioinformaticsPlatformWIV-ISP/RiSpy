from .common import check_exists, check_is_dir, check_tool

def validate_ont_inputs(
    reference,
    barcodes,
    ref_vcf,
    fastq_dir,
    chain_file,
    chromosome_mapping,
    centromere,
):
    check_exists(reference, "Reference genome")
    check_exists(barcodes, "Barcode TSV")
    check_exists(ref_vcf, "Reference VCF")
    check_exists(chain_file, "Chain file")
    check_exists(chromosome_mapping, "Chromosome mapping file")
    check_exists(centromere, "Centromere TSV file")

    check_is_dir(fastq_dir, "FASTQ directory")

    fastqs = list(fastq_dir.glob("*.fastq.gz"))
    if not fastqs:
        raise ValueError(f"No FASTQ files found in {fastq_dir}")

    for tool in ["NanoPlot", "seqkit", "chopper", "minimap2", "samtools", "bcftools", "qualimap"]:
        check_tool(tool)
