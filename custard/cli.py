"""CUSTARD command-line interface."""

import os
import sys
import shutil
import subprocess

import click

from custard import __version__

# Import from scripts/ directory (alongside Snakefile)
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "scripts"))
from config_builder import (
    build_config, write_config, find_pipeline_dir,
    detect_read_suffixes, count_samples, find_metadata,
    generate_metadata_template, validate_metadata
)


def _pipeline_dir():
    try:
        return find_pipeline_dir()
    except FileNotFoundError:
        click.echo("Error: Cannot find CUSTARD Snakefile.", err=True)
        click.echo("Run this command from the CUSTARD directory, or install CUSTARD properly.", err=True)
        sys.exit(1)


@click.group()
@click.version_option(__version__, prog_name="custard")
def cli():
    """CUSTARD - rpoC amplicon sequencing pipeline.

    Bacterial community profiling using the rpoC gene as a phylogenetic
    marker. Takes paired-end Illumina reads and produces taxonomy tables,
    species-level classifications, KO functional predictions, and
    interactive visualizations.
    """
    pass


@cli.command()
def setup():
    """Download reference databases (one-time, ~500 MB)."""
    pdir = _pipeline_dir()
    click.echo("Downloading CUSTARD reference databases...")
    result = subprocess.run(
        ["snakemake", "--cores", "1", "setup"],
        cwd=pdir
    )
    if result.returncode == 0:
        click.echo("Databases ready.")
    else:
        click.echo("Database setup failed. Check your internet connection and try again.", err=True)
        sys.exit(1)


@cli.command()
@click.option("--input", "input_dir", required=True, type=click.Path(exists=True),
              help="Directory containing paired-end FASTQ files.")
@click.option("--output", "output_dir", default="results", type=click.Path(),
              help="Output directory (default: results/).")
@click.option("--metadata", default=None, type=click.Path(),
              help="Sample metadata file (TSV with SampleID column). Auto-generated if not found.")
@click.option("--threads", default=4, type=int, show_default=True,
              help="Number of CPU threads.")
@click.option("--method", default="dada2", type=click.Choice(["dada2", "kraken2"]),
              show_default=True, help="Taxonomy classification method.")
@click.option("--dry-run", is_flag=True, help="Validate inputs without running the pipeline.")
@click.option("--verbose", is_flag=True, help="Show full Snakemake output.")
def run(input_dir, output_dir, metadata, threads, method, dry_run, verbose):
    """Run the CUSTARD pipeline on your FASTQ files.

    \b
    Minimal usage:
        custard run --input my_fastqs/

    \b
    With metadata and more threads:
        custard run --input my_fastqs/ --metadata samples.tsv --threads 16
    """
    pdir = _pipeline_dir()

    db_sentinel = os.path.join(pdir, "data", "db", ".setup_complete")
    taxa_db = os.path.join(pdir, "data", "db", "rpocdb_v3.1_dada2_toGenus.fasta")
    if not os.path.exists(db_sentinel) and not os.path.exists(taxa_db):
        click.echo("Databases not found. Running 'custard setup' first...")
        ctx = click.get_current_context()
        ctx.invoke(setup)

    try:
        r1, r2 = detect_read_suffixes(input_dir)
    except (FileNotFoundError, ValueError) as e:
        click.echo(f"Error: {e}", err=True)
        sys.exit(1)

    n_samples = count_samples(input_dir, r1)

    if metadata and metadata != "none" and not os.path.exists(metadata):
        metadata = None

    if metadata is None:
        metadata = find_metadata(input_dir)

    if metadata is None:
        template_path = os.path.join(os.getcwd(), "metadata.tsv")
        generated = generate_metadata_template(input_dir, r1, template_path)
        if generated:
            click.echo()
            click.echo("No metadata file found. A template has been created:")
            click.echo(f"  {generated}")
            click.echo()
            click.echo("Please fill in the 'Group' column for each sample, then re-run:")
            click.echo(f"  custard run --input {input_dir} --metadata {generated}")
            click.echo()
            click.echo("Or run without metadata (skips diversity analysis):")
            click.echo(f"  custard run --input {input_dir} --metadata none")
            sys.exit(0)

    if metadata and metadata != "none" and os.path.exists(metadata):
        ok, msg = validate_metadata(metadata)
        if not ok:
            click.echo()
            click.echo(f"Metadata issue: {msg}")
            click.echo()
            click.echo("Fill in the Group column and re-run, or skip metadata:")
            click.echo(f"  custard run --input {input_dir} --metadata none")
            sys.exit(1)

    if metadata == "none":
        metadata = None

    click.echo(f"CUSTARD v{__version__}")
    click.echo(f"  Input:    {os.path.abspath(input_dir)}")
    click.echo(f"  Samples:  {n_samples} (detected suffixes: {r1})")
    click.echo(f"  Output:   {os.path.abspath(output_dir)}")
    click.echo(f"  Metadata: {metadata if metadata else 'none (diversity plots will be skipped)'}")
    click.echo(f"  Method:   {method}")
    click.echo(f"  Threads:  {threads}")
    click.echo()

    config, _ = build_config(
        input_dir=input_dir,
        output_dir=output_dir,
        metadata=metadata,
        threads=threads,
        method=method,
        pipeline_dir=pdir,
    )

    config_path = os.path.join(pdir, "config.yaml")
    write_config(config, config_path)

    cmd = ["snakemake", "--cores", str(threads)]

    if dry_run:
        cmd.append("--dry-run")
        click.echo("Dry run — validating inputs...")
    else:
        click.echo("Starting pipeline...")

    if not verbose:
        cmd.append("--quiet")

    result = subprocess.run(cmd, cwd=pdir)

    if result.returncode == 0:
        if dry_run:
            click.echo("Validation passed. Remove --dry-run to execute.")
        else:
            click.echo()
            click.echo("Pipeline complete!")
            click.echo(f"  Results: {os.path.abspath(output_dir)}/")
            click.echo(f"  Report:  {os.path.abspath(output_dir)}/report.html")
            if os.path.exists(os.path.join(output_dir, "taxonomy", "krona.html")):
                click.echo(f"  Krona:   {os.path.abspath(output_dir)}/taxonomy/krona.html")
    else:
        click.echo("Pipeline failed. Run with --verbose to see details.", err=True)
        sys.exit(1)


@cli.command()
@click.option("--output", "output_dir", default="results", type=click.Path(),
              help="Output directory to clean.")
@click.confirmation_option(prompt="This will delete all results. Continue?")
def clean(output_dir):
    """Remove all pipeline outputs."""
    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)
        click.echo(f"Removed {output_dir}/")
    else:
        click.echo("Nothing to clean.")


def main():
    cli()


if __name__ == "__main__":
    main()
