from .common import check_exists, check_tool

def validate_detection_inputs(
    barcodes,
    ref_vcf,
    bam,
    query_vcf,
    chain_file,
    chromosome_mapping,
    centromere,
):
    check_exists(barcodes, "Barcode TSV")
    check_exists(ref_vcf, "Reference VCF")
    check_exists(bam, "BAM file")
    check_exists(query_vcf, "Query VCF")
    check_exists(chain_file, "Chain file")
    check_exists(chromosome_mapping, "Chromosome mapping file")
    check_exists(centromere, "Centromere TSV file")

    for tool in ["liftOver", "bedtools"]:
        check_tool(tool)
