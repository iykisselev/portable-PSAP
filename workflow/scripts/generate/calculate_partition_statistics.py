"""Snakemake adapter for bounded probability aggregation."""

import sys

sys.path.insert(0, snakemake.scriptdir)
from _bootstrap import bootstrap

bootstrap(snakemake, configure_polars=True)

from _library import bind

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
        bind(PartitionStatisticsInputs, snakemake.input),
        bind(PartitionStatisticsOutputs, snakemake.output),
        bind(PartitionStatisticsParameters, snakemake.params),
        str(snakemake.wildcards.partition),
        str(snakemake.wildcards.bucket),
    )


if __name__ == "__main__":
    main()
