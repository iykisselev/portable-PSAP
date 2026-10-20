"""Snakemake adapter for generated-bundle finalization."""

from pathlib import Path
import sys

sys.path.insert(0, str(Path(snakemake.scriptdir).parents[1] / "shared"))
from _bootstrap import bootstrap

bootstrap(snakemake, configure_polars=True)

from _library import bind

from bundle_finalization import finalize_partitioned_bundle
from generation_bindings import (
    BundleFinalizationInputs,
    BundleFinalizationOutputs,
    BundleFinalizationParameters,
)


def main() -> None:
    """Decode Snakemake bindings and publish the bundle through the library."""
    finalize_partitioned_bundle(
        bind(BundleFinalizationInputs, snakemake.input),
        bind(BundleFinalizationOutputs, snakemake.output),
        bind(BundleFinalizationParameters, snakemake.params),
    )


if __name__ == "__main__":
    main()
