from Bio import SeqIO
import yaml
from datetime import datetime
from pathlib import Path

def build_fq_by_sample(dir_samples):
    fq_by_sample = {}

    # Case 1: fastq_dir is None → return {"sample": None}
    if dir_samples is None:
        fq_by_sample["sample"] = None
        return fq_by_sample

    dir_samples = Path(dir_samples)

    # Detect ONT
    ont_mode = any(
        f.name.endswith(".fastq.gz")
        and "_R1" not in f.name
        and "_R2" not in f.name
        for f in dir_samples.iterdir()
    )

    if ont_mode:
        fq_by_sample["merged_flowcells"] = []
        return fq_by_sample

    # Illumina logic
    for file in dir_samples.iterdir():
        if not file.name.endswith(".fastq.gz"):
            continue
        if "_R1.fastq.gz" in file.name:
            name = file.name.split("_R1.fastq.gz")[0]
            if name not in fq_by_sample:
                r1 = list(dir_samples.glob(f"{name}_R1.fastq.gz"))
                r2 = list(dir_samples.glob(f"{name}_R2.fastq.gz"))
                if r1 and r2:
                    fq_by_sample[name] = [str(r1[0]), str(r2[0])]
                else:
                    raise FileNotFoundError(
                        f"Missing R1 or R2 fastq.gz for sample {name}"
                    )

    return fq_by_sample


def rename_and_filter_fasta(input_fasta: str, output_fasta: str) -> None:
    """
    Define new chromosome names for the first 12 sequences
    """
    chr_names = [f"chr{i}" for i in range(1, 13)]
    with open(input_fasta, "r") as in_handle, open(output_fasta, "w") as out_handle:
        records = list(SeqIO.parse(in_handle, "fasta"))
        #Take only first 12 records (chromosomes)
        selected_records = records[:12]
        for rec, new_name in zip(selected_records, chr_names):
            rec.id = new_name
            #remove description to keep header clean
            rec.description = ""
        SeqIO.write(selected_records, out_handle, "fasta")

def write_table(f: str, title: str, metrics: str) -> None:
    """
    Writes a single table to an HTML file.
    f: file object (already open)
    title: string, table title
    metrics: list of (key, value) tuples
    """
    f.write(f"<h2>{title}</h2>\n")
    f.write("<table>\n")
    f.write("<tr><th>Metric</th><th>Value</th></tr>\n")
    for key, value in metrics:
        f.write(f"<tr><td class='metric'>{key}</td><td>{value}</td></tr>\n")
    f.write("</table>\n\n")

def generate_report_html(html_file: str, config: str) -> None:
    """
    Generate an HTML report with optional QC and Fingerprint metrics.
    """
    qc_metrics = []
    fingerprint_metrics = []
    with open(config, "r") as f:
        result = yaml.safe_load(f)
    if list(result.keys())[0] == "Mean Depth of Coverage":
        qc_metrics.append(("Mean Depth of Coverage", result["Mean Depth of Coverage"]))
        qc_metrics.append(("Median Depth of Coverage", result["Median Depth of Coverage"]))
        qc_metrics.append(("Breadth of Coverage", result["Breadth of Coverage"]))
        qc_metrics.append(("Number of Primary Mapped Reads", result["Number of Primary Mapped Reads"]))
        qc_metrics.append(("Fold80", result["Fold80"]))
        fingerprint_metrics.append(("Number of Detected Alternative SNVs", result["Number of Detected Alternative SNVs"]))
        fingerprint_metrics.append(("Number of Detected Reference SNVs", result["Number of Detected Reference SNVs"]))
        fingerprint_metrics.append(("Number of Detected Barcodes", result["Number of Detected Barcodes"]))
        fingerprint_metrics.append(("Number of Chromosome Arms Covered with Barcodes", result["Number of Chromosome Arms Covered with Barcodes"]))
    elif list(result.keys())[0] == "Number of Detected Alternative SNVs":
        fingerprint_metrics.append(("Number of Detected Alternative SNVs", result["Number of Detected Alternative SNVs"]))
        fingerprint_metrics.append(("Number of Detected Reference SNVs", result["Number of Detected Reference SNVs"]))
        fingerprint_metrics.append(("Number of Detected Barcodes", result["Number of Detected Barcodes"]))
        fingerprint_metrics.append(("Number of Chromosome Arms Covered with Barcodes",result["Number of Chromosome Arms Covered with Barcodes"]))

    # Get current datetime
    report_datetime = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    qc_metrics = qc_metrics or []
    fingerprint_metrics = fingerprint_metrics or []
    with open(html_file, "w") as f:
        f.write("""<!DOCTYPE html>
        <html lang="en">
        <head>
        <meta charset="UTF-8">
        <title>Final Report</title>
        <style>
        body { font-family: Arial, sans-serif; margin: 30px; background-color: #f9f9f9; }
        h1 { text-align: center; color: #2c3e50; }
        h2 { color: #34495e; margin-top: 30px; text-align: center; }
        table { width: 60%; margin: 20px auto; border-collapse: collapse; box-shadow: 0 2px 8px rgba(0,0,0,0.1); }
        th, td { padding: 12px 15px; text-align: left; }
        th { background-color: #3498db; color: white; }
        tr:nth-child(even) { background-color: #ecf0f1; }
        tr:hover { background-color: #d1ecf1; }
        .metric { font-weight: bold; color: #2c3e50; }
        </style>
        </head>
        <body>
        <h1>Final Report</h1>
        """)
        f.write(f"""<p style="text-align:center; font-style: italic; color: #7f8c8d;">
        Generated on: {report_datetime}</p>
        """)
        if qc_metrics:
            write_table(f, "QC Report", qc_metrics)
        if fingerprint_metrics:
            write_table(f, "Fingerprint Detection Report", fingerprint_metrics)
        f.write("</body></html>")