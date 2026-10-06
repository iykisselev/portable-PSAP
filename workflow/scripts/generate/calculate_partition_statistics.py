"""Snakemake adapter for bounded probability aggregation."""

import os
from pathlib import Path
import sys

if "snakemake" in globals():
    os.environ.setdefault("POLARS_MAX_THREADS", str(snakemake.threads))

SHARED_SCRIPT_DIR = Path(snakemake.scriptdir).parents[1] / "shared"
sys.path.insert(0, str(SHARED_SCRIPT_DIR))
from _library import add_library_path

add_library_path(snakemake.scriptdir)

from generation_statistics import calculate_partition_statistics
from generation_bindings import (
    PartitionStatisticsInputs,
    PartitionStatisticsOutputs,
    PartitionStatisticsParameters,
)
from generation_probability import statistics_for_role


def main() -> None:
    """Decode Snakemake bindings and delegate probability aggregation."""
    calculate_partition_statistics(
        PartitionStatisticsInputs(
            partitions=str(snakemake.input.partitions),
            buckets=str(snakemake.input.buckets),
            allele_states=str(snakemake.input.allele_states),
            allele_membership=str(snakemake.input.allele_membership),
            allele_scores=str(snakemake.input.allele_scores),
            mutation_rate=(str(snakemake.input.mutation_rate[0]) if snakemake.input.mutation_rate else None),
            mutation_rate_index=(str(snakemake.input.mutation_rate_index[0]) if snakemake.input.mutation_rate_index else None),
            mutation_rate_indel=(str(snakemake.input.mutation_rate_indel[0]) if snakemake.input.mutation_rate_indel else None),
            mutation_rate_indel_index=(str(snakemake.input.mutation_rate_indel_index[0]) if snakemake.input.mutation_rate_indel_index else None),
        ),
        PartitionStatisticsOutputs(statistics=str(snakemake.output.statistics)),
        PartitionStatisticsParameters(
            x_af_policy=str(snakemake.params.x_af_policy),
            unobserved_policy=str(snakemake.params.unobserved_policy),
            x_mutation_rate_policy=str(snakemake.params.x_mutation_rate_policy),
            score_max=float(snakemake.params.score_max),
            bin_width=float(snakemake.params.bin_width),
            roles=tuple(snakemake.params.roles),
            compression=str(snakemake.params.compression),
        ),
        str(snakemake.wildcards.partition),
        str(snakemake.wildcards.bucket),
    )


if __name__ == "__main__":
    main()
