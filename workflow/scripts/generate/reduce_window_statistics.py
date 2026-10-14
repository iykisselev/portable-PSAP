"""Snakemake adapter for lookup-bucket rendering."""

from pathlib import Path
import sys

sys.path.insert(0, str(Path(snakemake.scriptdir).parents[1] / "shared"))
from _bootstrap import bootstrap

bootstrap(snakemake, configure_polars=True)

from _library import bind

from generation_lookup import reduce_window_statistics
from generation_bindings import (
    LookupReductionInputs,
    LookupReductionOutputs,
    LookupReductionParameters,
)


def main() -> None:
    """Decode Snakemake bindings and delegate lookup rendering."""
    reduce_window_statistics(
        bind(LookupReductionInputs, snakemake.input),
        bind(LookupReductionOutputs, snakemake.output),
        bind(LookupReductionParameters, snakemake.params),
        str(snakemake.wildcards.bucket),
    )


if __name__ == "__main__":
    main()
