"""Snakemake adapter for generation scope discovery."""

import os
from pathlib import Path
import sys
from types import SimpleNamespace

if "snakemake" in globals():
    os.environ.setdefault("POLARS_MAX_THREADS", str(snakemake.threads))

SCRIPTS_ROOT = Path(snakemake.scriptdir).parent
sys.path.insert(0, str(SCRIPTS_ROOT))
from _library import add_library_path

add_library_path(snakemake.scriptdir)

from generation_scope import prepare_generation_scope


def main() -> None:
    """Decode Snakemake bindings and delegate domain work to the library."""
    prepare_generation_scope(
        SimpleNamespace(
            population_vcf=str(snakemake.input.population_vcf),
            score_snv=str(snakemake.input.score_snv),
            coverage=str(snakemake.input.coverage),
            testing_units=str(snakemake.input.testing_units),
            par_regions=str(snakemake.input.par_regions),
            exclusions=[str(path) for path in snakemake.input.exclusions],
            fasta=str(snakemake.input.fasta),
            fasta_fai=str(snakemake.input.fasta_fai),
            mutation_rate=[str(path) for path in snakemake.input.mutation_rate],
        ),
        SimpleNamespace(**{name: str(snakemake.output[name]) for name in snakemake.output.keys()}),
        SimpleNamespace(
            coverage_minimum=float(snakemake.params.coverage_minimum),
            partition_size=int(snakemake.params.partition_size),
            unit_bucket_size=int(snakemake.params.unit_bucket_size),
            compression=str(snakemake.params.compression),
        ),
    )


if __name__ == "__main__":
    main()
