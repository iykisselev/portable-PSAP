"""Snakemake adapter for one bounded generation window."""

from pathlib import Path
import sys

sys.path.insert(0, str(Path(snakemake.scriptdir).parents[1] / "shared"))
from _bootstrap import bootstrap

bootstrap(snakemake, configure_polars=True)

from _library import bind

from generation_window import prepare_generation_window
from generation_bindings import (
    GenerationWindowInputs,
    GenerationWindowOutputs,
    GenerationWindowParameters,
)


def main() -> None:
    """Decode Snakemake bindings and delegate window joins to the library."""
    prepare_generation_window(
        bind(GenerationWindowInputs, snakemake.input),
        bind(GenerationWindowOutputs, snakemake.output),
        bind(GenerationWindowParameters, snakemake.params),
        str(snakemake.wildcards.window),
    )


if __name__ == "__main__":
    main()
