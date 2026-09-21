"""Snakemake adapter for one bounded generation partition."""

import os
from pathlib import Path
import sys
from types import SimpleNamespace

if "snakemake" in globals():
    os.environ.setdefault("POLARS_MAX_THREADS", str(snakemake.threads))

SHARED_SCRIPT_DIR = Path(snakemake.scriptdir).parents[1] / "shared"
sys.path.insert(0, str(SHARED_SCRIPT_DIR))
from _library import add_library_path

add_library_path(snakemake.scriptdir)

from generation_partition import prepare_generation_partition


def main() -> None:
    """Decode Snakemake bindings and delegate partition joins to the library."""
    prepare_generation_partition(
        SimpleNamespace(
            partitions=str(snakemake.input.partitions),
            units=str(snakemake.input.units),
            unit_intervals=str(snakemake.input.unit_intervals),
            mask=str(snakemake.input.mask),
            population_vcf=str(snakemake.input.population_vcf),
            population_vcf_index=(str(snakemake.input.population_vcf_index[0]) if snakemake.input.population_vcf_index else None),
            score_snv=str(snakemake.input.score_snv),
            score_snv_index=(str(snakemake.input.score_snv_index[0]) if snakemake.input.score_snv_index else None),
            fasta=str(snakemake.input.fasta),
            fasta_fai=str(snakemake.input.fasta_fai),
        ),
        SimpleNamespace(**{name: str(snakemake.output[name]) for name in snakemake.output.keys()}),
        SimpleNamespace(
            af_tags=dict(snakemake.params.af_tags),
            x_af_policy=str(snakemake.params.x_af_policy),
            hemizygote=bool(snakemake.params.hemizygote),
            compression=str(snakemake.params.compression),
        ),
        str(snakemake.wildcards.partition),
    )


if __name__ == "__main__":
    main()
