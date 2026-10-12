"""Snakemake adapter for bounded probability aggregation."""

from pathlib import Path
import sys

sys.path.insert(0, str(Path(snakemake.scriptdir).parents[1] / "shared"))
from _bootstrap import bootstrap

bootstrap(snakemake, configure_polars=True)

from _library import bind

from generation_statistics import calculate_partition_statistics
from generation_bindings import (
    PartitionStatisticsInputs,
    PartitionStatisticsOutputs,
    PartitionStatisticsParameters,
)


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
