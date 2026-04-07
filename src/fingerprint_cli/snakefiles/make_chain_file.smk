from pathlib import Path
from Bio import SeqIO
import pysam
from fingerprint_cli.io import rename_and_filter_fasta

root = Path(config["paths"]["output_dir"])

ref_3KRG = config["paths"]["old_ref"]

new_ref = config["paths"]["new_ref"]

rule all:
    input:
        Chain_File = root / 'chain_file' / 'CHAINFILE.chain'

rule cleaning_fasta:
    """
    Converting the chromosome names to [chr1:chr12] and removing organelles contigs
    """
    input:
        old_ref = ref_3KRG,
        new_ref = new_ref
    output:
        old_ref_renamed = root / 'chain_file' / 'old_ref_renamed.fasta',
        new_ref_renamed = root / 'chain_file' / 'new_ref_renamed.fasta',
    threads: config["resources"]["cleaning_fasta"]["threads"]
    resources:
        mem_mb = config["resources"]["cleaning_fasta"]["memory"]
    benchmark: root / 'chain_file' / 'benchmarks' / 'cleaning_fasta.txt'
    run:
        rename_and_filter_fasta(input.old_ref,output.old_ref_renamed)
        rename_and_filter_fasta(input.new_ref,output.new_ref_renamed)


rule making_PAF:
    """
    Generating the PAF file using minimap2
    """
    input:
        old_ref_renamed = rules.cleaning_fasta.output.old_ref_renamed,
        new_ref_renamed = rules.cleaning_fasta.output.new_ref_renamed,
    output:
        PAF = root / 'chain_file' / 'PAF_FILE.paf',
    threads: config["resources"]["making_PAF"]["threads"]
    resources:
        mem_mb = config["resources"]["making_PAF"]["memory"]
    benchmark: root / 'chain_file' / 'benchmarks' / 'making_PAF.txt'
    shell:
        """
        minimap2 -t {threads} -cx asm5 --cs {input.new_ref_renamed} {input.old_ref_renamed} > {output.PAF}
        """

rule making_chain_file:
    """
    Generating the chain file 
    """
    input:
        PAF = rules.making_PAF.output.PAF
    output:
        Chain_File = root / 'chain_file' / 'CHAINFILE.chain',
    threads: config["resources"]["making_chain_file"]["threads"]
    resources:
        mem_mb = config["resources"]["making_chain_file"]["memory"]
    benchmark: root / 'chain_file' / 'benchmarks' / 'making_chain_file.txt'
    shell:
        """
        transanno minimap2chain {input.PAF} --output {output.Chain_File}
        """

