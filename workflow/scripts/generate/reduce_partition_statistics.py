"""Snakemake adapter for lookup-bucket rendering."""

import sys

sys.path.insert(0, snakemake.scriptdir)
from _bootstrap import bootstrap

bootstrap(snakemake, configure_polars=True)

from _library import bind

from generation_lookup import reduce_partition_statistics
from generation_bindings import (
    LookupReductionInputs,
    LookupReductionOutputs,
    LookupReductionParameters,
)


def main() -> None:
    """Decode Snakemake bindings and delegate lookup rendering."""
    reduce_partition_statistics(
        bind(LookupReductionInputs, snakemake.input),
        bind(LookupReductionOutputs, snakemake.output),
        bind(LookupReductionParameters, snakemake.params),
        str(snakemake.wildcards.bucket),
    )


if __name__ == "__main__":
    main()
