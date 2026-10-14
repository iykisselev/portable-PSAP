"""Snakemake adapter for bounded probability aggregation."""

from pathlib import Path
import sys

sys.path.insert(0, str(Path(snakemake.scriptdir).parents[1] / "shared"))
from _bootstrap import bootstrap

bootstrap(snakemake, configure_polars=True)

from _library import bind

from generation_statistics import calculate_window_statistics
from generation_bindings import (
    WindowStatisticsInputs,
    WindowStatisticsOutputs,
    WindowStatisticsParameters,
)


def main() -> None:
    """Decode Snakemake bindings and delegate probability aggregation."""
    calculate_window_statistics(
        bind(WindowStatisticsInputs, snakemake.input),
        bind(WindowStatisticsOutputs, snakemake.output),
        bind(WindowStatisticsParameters, snakemake.params),
        str(snakemake.wildcards.window),
        str(snakemake.wildcards.bucket),
    )


if __name__ == "__main__":
    main()
