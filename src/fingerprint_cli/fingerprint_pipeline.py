import pandas as pd
import pysam
import yaml
from typing import Dict, List
import os

def _check_file_exists(path: str) -> None:
    """
    Raise FileNotFoundError if file does not exist
    """
    if not os.path.exists(path):
        raise FileNotFoundError(f"Required file not found: {path}")

def detect_unique_snps(input: str, output: str) -> None:
    """
    Detect unique SNPs from a 2‑SNV barcode table.

    Parameters
    ----------
    input : str
        Path to a TSV file containing barcode combinations.
        Expected columns:
            index1, CHROMOSOME_x, POSITION_x, REFCALL_x, GENOTYPE_x, Scores_x
            index2, CHROMOSOME_y, POSITION_y, REFCALL_y, GENOTYPE_y, Scores_y

    output : str
        Path to the output TSV file containing unique SNPs.

    Output
    ------
    Writes a TSV file with columns:
        snp_index, chromosome, position, refcall, genotype, score
    """
    _check_file_exists(input)
    with open(input, 'r') as fing:
        first_line = fing.readline().strip()
    #Check if the file contains the specific message
    if first_line.startswith('#No'):
        #Write the message to the output file
        with open(output, 'w') as out_f:
            out_f.write('#No Fingerprints was detected\n')
    else:
        df = pd.read_csv(input, sep='\t')
        df_index1 = df[['index1', 'CHROMOSOME_x', 'POSITION_x', 'REFCALL_x', 'GENOTYPE_x', 'Scores_x']].copy()
        df_index1.columns = ['snp_index', 'chromosome', 'position', 'refcall', 'genotype', 'score']
        df_index2 = df[['index2', 'CHROMOSOME_y', 'POSITION_y', 'REFCALL_y', 'GENOTYPE_y', 'Scores_y']].copy()
        df_index2.columns = ['snp_index', 'chromosome', 'position', 'refcall', 'genotype', 'score']
        df_combined = pd.concat([df_index1, df_index2], ignore_index=True)
        df_unique = df_combined.drop_duplicates(subset='snp_index')
        df_unique = df_unique.sort_values(by='score', ascending=False).reset_index(drop=True)
        df_unique.to_csv(output, sep='\t', index=False)

def bed_maker(input: str, output: str) -> None:
    """
    Convert a unique SNP table into a BED file.

    Parameters
    ----------
    input : str
        Path to a TSV file produced by detect_unique_snps().
        Required columns: chromosome, position, snp_index, refcall, genotype

    output : str
        Path to the BED file to write.

    Output
    ------
    BED file with columns:
        chromosome, start, end, snp_index, refcall, genotype
    """
    _check_file_exists(input)
    df = pd.read_csv(input, sep='\t')
    positions = df[["chromosome", "position", "snp_index", "refcall", "genotype"]].copy()
    positions['chromosome'] = positions['chromosome'].astype(str)
    #Prepend 'chr' to chromosome identifiers
    positions['chromosome'] = 'chr' + positions['chromosome']
    #Convert positions from 1-based to 0-based for BED format
    positions['start'] = positions['position'] - 1
    positions['end'] = positions['position']
    #Define the desired chromosome order
    chrom_order = [f'chr{i}' for i in range(1, 13)]
    #Convert 'chrom' column to a categorical type with the specified order
    positions['chromosome'] = pd.Categorical(positions['chromosome'], categories=chrom_order, ordered=True)
    #Select and reorder columns for BED format
    bed_df = positions[['chromosome', 'start', 'end', 'snp_index', 'refcall', 'genotype']]
    #Sort by chromosome and start position
    bed_df = bed_df.sort_values(['chromosome', 'start'])
    #Write to BED file without header and index
    bed_df.to_csv(output, sep='\t', header=False, index=False)

def chromosome_name_conversion(mapping_file: str, bed_file: str, output: str) -> None:
    """
    Convert chromosome names in a BED file using a mapping table.

    Parameters
    ----------
    mapping_file : str
        Two-column whitespace-separated file: old_name new_name

    bed_file : str
        BED file with chromosome names in column 1.

    output : str
        Path to write the converted BED file.
    """
    _check_file_exists(mapping_file)
    _check_file_exists(bed_file)
    chromosome_mapping = {}
    with open(mapping_file, 'r') as f:
        for line in f:
            parts = line.strip().split()
            chromosome_mapping[parts[0]] = parts[1]
    with open(bed_file, 'r') as f, open(output, 'w') as out:
        for line in f:
            parts = line.strip().split()
            if parts[0] in chromosome_mapping:
                parts[0] = chromosome_mapping[parts[0]]
            out.write('\t'.join(parts) + '\n')

def prepend_summary(file_path: str, summary_line: str) -> None:
    """
    Prepend a summary line to a file.
    """
    _check_file_exists(file_path)
    with open(file_path, 'r+') as file:
        content = file.read()
        file.seek(0)
        file.write(summary_line + '\n' + '\n' + content)

def validation_fingerprints_alt(bed_file: str, vcf_file: str, barcode_file: str, detected_snps: str, undetected_snps: str, detected_barcodes: str, config: str) -> None:
    """
    Validate alternative SNPs by checking whether the ALT allele appears in the VCF.

    Parameters
    ----------
    bed_file : str
        BED file with SNP positions.
    vcf_file : str
        VCF file to query.
    barcode_file : str
        TSV file with barcode pairs.
    detected_snps : str
        Output file for detected SNPs.
    undetected_snps : str
        Output file for undetected SNPs.
    detected_barcodes : str
        Output file for filtered barcodes.
    config : str
        YAML file to append summary results.
    """
    for f in [bed_file, vcf_file, barcode_file]:
        _check_file_exists(f)
    alt = pd.read_csv(bed_file, sep='\t', names=["chr", "start", "end", "index", "ref", "genotype"])
    barcodes = pd.read_csv(barcode_file, sep='\t')
    vcf = pysam.VariantFile(vcf_file)
    with open(detected_snps, "w") as detected_file, open(undetected_snps, "w") as undetected_file:
        detected_count = 0
        undetected_count = 0
        for index, row in alt.iterrows():
            chrom = str(row["chr"])
            pos = int(row["end"])
            genotype = row["genotype"]
            ind = row["index"]
            found = False
            #fetch is [start, end) and it's 0-based
            for rec in vcf.fetch(chrom, pos - 1, pos):
                if rec.chrom == chrom and rec.pos == pos and genotype in rec.alts:
                    detected_file.write(str(rec) + "\n")
                    detected_count += 1
                    found = True
                    break
            if not found:
                #If no matching record is found, write the SNP information to the undetected file
                undetected_file.write(f"{chrom}\t{pos}\t{genotype}\n")
                undetected_count += 1
                #Filter out barcodes
                barcodes = barcodes[~((barcodes["index1"] == ind) | (barcodes["index2"] == ind))]
    barcodes.to_csv(detected_barcodes, sep="\t", index=False)
    #Prepend summary lines to the respective files
    prepend_summary(detected_snps, f"Number of Detected SNPs: {detected_count}")
    prepend_summary(undetected_snps, f"Number of Undetected SNPs: {undetected_count}")
    results = {"Number of Detected Alternative SNVs": str(detected_count)}
    with open(config, "a") as f:
        yaml.dump(results, f, sort_keys=False)

def validation_fingerprints_ref(bed_file: str, vcf_file: str, bam_file: str, barcodes_file: str, detected_snps: str, undetected_snps: str, detected_barcodes: str, config: str) -> None:
    """
    Validate reference SNPs by checking coverage and REF/ALT consistency.

    Parameters
    ----------
    bed_file : str
        BED file with SNP positions.
    vcf_file : str
        VCF file to query.
    bam_file : str
        BAM file with alignments.
    barcodes_file : str
        TSV file with barcode pairs.
    detected_snps : str
        Output file for detected SNPs.
    undetected_snps : str
        Output file for undetected SNPs.
    detected_barcodes : str
        Output file for filtered barcodes.
    config : str
        YAML file to append summary results.
    """
    for f in [bed_file, vcf_file, bam_file, barcodes_file]:
        _check_file_exists(f)
    ref = pd.read_csv(bed_file, sep='\t', names=["chr", "start", "end", "index", "ref", "genotype"])
    barcodes = pd.read_csv(barcodes_file, sep='\t')
    vcf = pysam.VariantFile(vcf_file)
    bamfile = pysam.AlignmentFile(bam_file, "rb")
    with open(detected_snps, "w") as detected_file, open(undetected_snps, "w") as undetected_file:
        detected_count = 0
        undetected_count = 0
        for index, row in ref.iterrows():
            chrom = str(row["chr"])
            pos = int(row["end"])
            genotype = row["genotype"]
            ind = row["index"]
            found = False
            #Get depth at this position (0-based start)
            a, c, g, t = bamfile.count_coverage(
                chrom,
                start=pos - 1,
                end=pos,
                quality_threshold=0
            )
            depth = a[0] + c[0] + g[0] + t[0]
            if depth < 5:
                undetected_file.write(f"Coverage below 5 at: {chrom}:{pos} \n")
                undetected_count += 1
                barcodes = barcodes[~((barcodes["index1"] == ind) | (barcodes["index2"] == ind))]
                found = True
            else:
                #fetch is [start, end) and it's 0-based
                for rec in vcf.fetch(chrom, pos - 1, pos):
                    if rec.chrom == chrom and rec.pos == pos and genotype not in rec.alts:
                        undetected_file.write(str(rec) + "\n")
                        undetected_count += 1
                        """
                        Filter out barcodes
                        """
                        barcodes = barcodes[~((barcodes["index1"] == ind) | (barcodes["index2"] == ind))]
                        found = True
                        break
            if not found:
                #If no matching record is found in the vcf file and the depth is at least 5X, write the SNP information to the detected file
                detected_file.write(f"{chrom}\t{pos}\t{genotype}\n")
                detected_count += 1
    num_barcodes = barcodes.shape[0]
    barcodes.to_csv(detected_barcodes, sep="\t", index=False)
    #Prepend summary lines to the respective files
    prepend_summary(detected_snps, f"Number of Detected SNPs: {detected_count}")
    prepend_summary(undetected_snps, f"Number of Undetected SNPs: {undetected_count}")
    results = {"Number of Detected Reference SNVs": str(detected_count), "Number of Detected Barcodes": str(num_barcodes)}
    with open(config, "a") as f:
        yaml.dump(results, f, sort_keys=False)


def barcodes_positions_new_ref(barcodes_file: str, bed_file: str, output: str) -> None:
    """
    Update barcode SNP positions using a BED file mapped to a new reference genome.

    Parameters
    ----------
    barcodes_file : str
        TSV file with barcode pairs and SNP positions.

    bed_file : str
        BED file with updated SNP positions.

    output : str
        Path to write the updated barcode table.
    """
    for f in [barcodes_file, bed_file]:
        _check_file_exists(f)
    barcodes = pd.read_csv(barcodes_file, sep="\t")
    snps = pd.read_csv(bed_file, sep="\t", header=None, names=['chromosome', 'start', 'end', 'snp_index', "refcall", "genotype"])
    snps_start = snps[['snp_index', 'end']].rename(columns={'snp_index': 'INDEX', 'end': 'NEW_POS'})
    barcodes = barcodes.merge(
        snps_start,
        how='left',
        left_on='index1',
        right_on='INDEX',
        suffixes=('', '_1_new')
    )
    barcodes.loc[~barcodes['NEW_POS'].isnull(), 'POSITION_x'] = barcodes.loc[~barcodes['NEW_POS'].isnull(), 'NEW_POS']
    barcodes.drop(['INDEX', 'NEW_POS'], axis=1, inplace=True)
    snps_start2 = snps_start.rename(columns={'INDEX': 'INDEX2', 'NEW_POS': 'NEW_POS2'})
    barcodes = barcodes.merge(
        snps_start2,
        how='left',
        left_on='index2',
        right_on='INDEX2'
    )
    barcodes.loc[~barcodes['NEW_POS2'].isnull(), 'POSITION_y'] = barcodes.loc[~barcodes['NEW_POS2'].isnull(), 'NEW_POS2']
    barcodes.drop(['INDEX2', 'NEW_POS2', 'Unnamed: 0', 'REFCALL_x', 'REFCALL_y'], axis=1, inplace=True)
    barcodes.to_csv(output, sep='\t', header=True, index=False)

def classify(chrom: str, pos: int, bounds: dict[str, dict[str, int]]) -> str:
    """
    Classify a SNP position relative to the centromere.

    Returns
    -------
    str
        One of:
            '<chrom>_1'  (short arm)
            'centromere'
            '<chrom>_2'  (long arm)
            None         (if chromosome not in bounds)
    """
    if chrom not in bounds:
        return None
    s, e = bounds[chrom]['Start'], bounds[chrom]['End']
    if pos < s:
        return f"{chrom}_1"
    elif pos > e:
        return f"{chrom}_2"
    else:
        return "centromere"

def chromosome_arms(centromere_file: str, barcodes_file: str, output: str, config: str) -> None:
    """
    Determine which chromosome arms contain SNPs from barcode pairs.

    Parameters
    ----------
    centromere_file : str
        TSV file with columns: Chromosome, Start, End

    barcodes_file : str
        TSV file with barcode pairs and SNP positions.

    output : str
        Path to write the chromosome arm classification.

    config : str
        YAML file to append summary results.
    """
    for f in [centromere_file, barcodes_file]:
        _check_file_exists(f)
    df = pd.read_csv(centromere_file, sep="\t")
    df = df[["Chromosome", "Start", "End"]]
    barcodes = pd.read_csv(barcodes_file, sep="\t")
    #Prepare chromosome boundaries as a dict
    bounds = df.set_index('Chromosome')[['Start', 'End']].to_dict('index')
    #Apply classifier to both SNPs in barcodes
    barcodes['snp1'] = barcodes.apply(lambda r: classify(f"Chr{r['CHROMOSOME_x']}", r['POSITION_x'], bounds), axis=1)
    barcodes['snp2'] = barcodes.apply(lambda r: classify(f"Chr{r['CHROMOSOME_y']}", r['POSITION_y'], bounds), axis=1)
    #Select only the desired columns for output
    result = barcodes[['index1', 'index2', 'snp1', 'snp2']]
    #Count centromere SNPs and collect their row indices
    stacked = result[['snp1', 'snp2']].stack().reset_index()
    stacked.columns = ['barcode_row', 'snp_col', 'snp_val']
    centromeres = stacked[stacked['snp_val'] == 'centromere']
    #Build detailed info
    centromere_info = []
    for _, row in centromeres.iterrows():
        br = int(row['barcode_row'])
        row_vals = result.loc[br, ['index1', 'index2', 'snp1', 'snp2']]
        info = (
            f"index1: {row_vals['index1']}, "
            f"index2: {row_vals['index2']}, "
            f"snp1: {row_vals['snp1']}, "
            f"snp2: {row_vals['snp2']}, "
            f"col: {row['snp_col']}"
        )
        centromere_info.append(info)
    n_cent = len(centromere_info)
    #Determine covered chromosome arms
    arms = set(result['snp1'].dropna()) | set(result['snp2'].dropna())
    arms.discard('centromere')
    full_arms = {f"Chr{i}_{j}" for i in range(1, 13) for j in (1, 2)}
    seen = sorted(arms)
    missing = sorted(full_arms - arms)
    #Summary lines
    lines = [
        f"# Covered {len(seen)}/24 chromosome arms: {', '.join(seen) or 'none'}",
        f"# Missing arms: {', '.join(missing) or 'none'}",
        f"# SNPs within centromere: {n_cent}",
    ]
    if n_cent:
        lines.append("# Centromere details:")
        lines.extend(f"#  • {info}" for info in centromere_info)
    #Write output
    with open(output, 'w') as f:
        f.write("\n".join(lines) + "\n")
        result.to_csv(f, sep="\t", index=False)
    results = {"Number of Chromosome Arms Covered with Barcodes": f"{len(seen)} / 24"}
    with open(config, "a") as f:
        yaml.dump(results, f, sort_keys=False)
