"""Snakemake adapter for generation scope discovery."""

import sys

sys.path.insert(0, snakemake.scriptdir)
from _bootstrap import bootstrap

bootstrap(snakemake, configure_polars=True)

from _library import bind

from generation_scope import prepare_generation_scope
from generation_bindings import (
    GenerationScopeInputs,
    GenerationScopeOutputs,
    GenerationScopeParameters,
)


def main() -> None:
    """Decode Snakemake bindings and delegate domain work to the library."""
    prepare_generation_scope(
        bind(GenerationScopeInputs, snakemake.input),
        bind(GenerationScopeOutputs, snakemake.output),
        bind(GenerationScopeParameters, snakemake.params),
    )


if __name__ == "__main__":
    main()
