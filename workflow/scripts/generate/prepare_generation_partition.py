"""Snakemake adapter for one bounded generation partition."""

from pathlib import Path
import sys

sys.path.insert(0, str(Path(snakemake.scriptdir).parents[1] / "shared"))
from _bootstrap import bootstrap

bootstrap(snakemake, configure_polars=True)

from _library import bind

from generation_partition import prepare_generation_partition
from generation_bindings import (
    GenerationPartitionInputs,
    GenerationPartitionOutputs,
    GenerationPartitionParameters,
)


def main() -> None:
    """Decode Snakemake bindings and delegate partition joins to the library."""
    prepare_generation_partition(
        bind(GenerationPartitionInputs, snakemake.input),
        bind(GenerationPartitionOutputs, snakemake.output),
        bind(GenerationPartitionParameters, snakemake.params),
        str(snakemake.wildcards.partition),
    )


if __name__ == "__main__":
    main()
