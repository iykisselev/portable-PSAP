"""Snakemake adapter for generation scope discovery."""

import os
from pathlib import Path
import sys

if "snakemake" in globals():
    os.environ.setdefault("POLARS_MAX_THREADS", str(snakemake.threads))

SHARED_SCRIPT_DIR = Path(snakemake.scriptdir).parents[1] / "shared"
sys.path.insert(0, str(SHARED_SCRIPT_DIR))
from _library import add_library_path

add_library_path(snakemake.scriptdir)

from generation_scope import prepare_generation_scope
from generation_bindings import (
    GenerationScopeInputs,
    GenerationScopeOutputs,
    GenerationScopeParameters,
)


def main() -> None:
    """Decode Snakemake bindings and delegate domain work to the library."""
    prepare_generation_scope(
        GenerationScopeInputs(
            population_vcf=str(snakemake.input.population_vcf),
            population_vcf_index=(
                str(snakemake.input.population_vcf_index[0])
                if snakemake.input.population_vcf_index
                else None
            ),
            score_snv=str(snakemake.input.score_snv),
            score_snv_index=(
                str(snakemake.input.score_snv_index[0])
                if snakemake.input.score_snv_index
                else None
            ),
            exclusion=str(snakemake.input.exclusion[0]) if snakemake.input.exclusion else None,
            testing_units=str(snakemake.input.testing_units),
            par_regions=str(snakemake.input.par_regions),
            exclusions={
                str(identifier): str(path)
                for identifier, path in snakemake.params.exclusions.items()
            },
            fasta=str(snakemake.input.fasta),
            fasta_fai=str(snakemake.input.fasta_fai),
            mutation_rate=(
                str(snakemake.input.mutation_rate[0])
                if snakemake.input.mutation_rate
                else None
            ),
            mutation_rate_index=(
                str(snakemake.input.mutation_rate_index[0])
                if snakemake.input.mutation_rate_index
                else None
            ),
            score_indel=(str(snakemake.input.score_indel[0]) if snakemake.input.score_indel else None),
            score_indel_index=(
                str(snakemake.input.score_indel_index[0])
                if snakemake.input.score_indel_index else None
            ),
        ),
        GenerationScopeOutputs(
            partitions=str(snakemake.output.partitions),
            buckets=str(snakemake.output.buckets),
            units=str(snakemake.output.units),
            unit_intervals=str(snakemake.output.unit_intervals),
            masked_unit_intervals=str(snakemake.output.masked_unit_intervals),
            mask=str(snakemake.output.mask),
            resource_identities=str(snakemake.output.resource_identities),
        ),
        GenerationScopeParameters(
            minimum_retained_fraction=float(snakemake.params.minimum_retained_fraction),
            partition_size=int(snakemake.params.partition_size),
            unit_bucket_size=int(snakemake.params.unit_bucket_size),
            compression=str(snakemake.params.compression),
            af_tags=dict(snakemake.params.af_tags),
            include_indels=bool(snakemake.params.include_indels),
            indel_max_reference_length=int(snakemake.params.indel_max_reference_length),
        ),
    )


if __name__ == "__main__":
    main()
