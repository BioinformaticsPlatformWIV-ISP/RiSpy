from pathlib import Path
import importlib.resources as resources
import yaml
import typer
import subprocess
import fingerprint_cli


def get_default_centromere_path() -> Path:
    return resources.files("fingerprint_cli.data") / "Centromere.tsv"

def get_default_chain_file_path() -> Path:
    return resources.files("fingerprint_cli.data") / "CHAINFILE.chain"

def get_default_chromosome_mapping_path() -> Path:
    return resources.files("fingerprint_cli.data") / "chromosome_mapping.txt"

def get_default_Nip_barcodes_path() -> Path:
    return resources.files("fingerprint_cli.data") / "Nipponbare_Barcodes.tsv"

def get_default_Nip_ref_vcf_path() -> Path:
    return resources.files("fingerprint_cli.data") / "Nipponbare.vcf"

def merge_config(template: dict, overrides: dict) -> dict:
    """
    Recursively merge overrides into template.
    Only assign override if the template key is missing or None.
    """
    for k, v in overrides.items():
        if k in template:
            if isinstance(v, dict) and isinstance(template[k], dict):
                merge_config(template[k], v)
            elif template[k] is None:
                template[k] = v
        else:
            template[k] = v
    return template

def load_and_write_config(template_name: str, overrides: dict, output_dir = str) -> Path:
    """
    Load template config, apply overrides only if the key is missing or has no value, and write to config.yaml.
    """
    with resources.files(fingerprint_cli).joinpath(f"templates/{template_name}").open("r") as f:
        config = yaml.safe_load(f)

    config = merge_config(config, overrides)
    output_dir.mkdir(exist_ok=True)
    config_path = output_dir / "config.yaml"
    with open(config_path, "w") as f:
        yaml.safe_dump(config, f, sort_keys=False)
    typer.echo(f"Config written to {config_path}")
    return config_path


def run_snakemake(snakefile: Path, config_path: Path, cores: int, mem_mb: int, dry_run: bool) -> None:
    """
    Run Snakemake with given parameters.
    """
    cmd = [
        "snakemake",
        "--snakefile", str(snakefile),
        "--configfile", str(config_path),
        "--cores", str(cores),
        "--resources", f"mem_mb={mem_mb}",
        "--rerun-triggers", "mtime",
    ]
    if dry_run:
        cmd.extend(["--dry-run", "--printshellcmds"])
    subprocess.run(cmd, check=True)

def set_paths(reference: str, barcodes: str, ref_vcf: str, output_dir: str, fastq_dir: str, chain_file: str, chromosome_mapping: str, centromere: str) -> dict:
    """
    Return dictionary of paths for config overrides.
    """
    return {
        "reference": str(reference),
        "barcodes": str(barcodes),
        "ref_vcf": str(ref_vcf),
        "output_dir": str(output_dir),
        "fastq_dir": str(fastq_dir),
        "chain_file": str(chain_file),
        "chromosome_mapping": str(chromosome_mapping),
        "centromere": str(centromere),
    }

def set_resources(resource_args: dict) -> dict:
    """
    Return dictionary of resources for config overrides.
    """
    resources_dict = {}
    for key, values in resource_args.items():
        resources_dict[key] = {k: int(v) for k, v in values.items()}
    return resources_dict

