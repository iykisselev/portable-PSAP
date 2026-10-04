"""Snakemake adapter for lookup-bucket rendering."""

import os
from pathlib import Path
import sys

if "snakemake" in globals():
    os.environ.setdefault("POLARS_MAX_THREADS", str(snakemake.threads))

SHARED_SCRIPT_DIR = Path(snakemake.scriptdir).parents[1] / "shared"
sys.path.insert(0, str(SHARED_SCRIPT_DIR))
from _library import add_library_path

add_library_path(snakemake.scriptdir)

from generation_lookup import reduce_partition_statistics
from generation_bindings import (
    LookupReductionInputs,
    LookupReductionOutputs,
    LookupReductionParameters,
)


def main() -> None:
    """Decode Snakemake bindings and delegate lookup rendering."""
    hem_output = getattr(snakemake.output, "lookup_hem", [])
    reduce_partition_statistics(
        LookupReductionInputs(
            buckets=str(snakemake.input.buckets),
            partition_statistics=tuple(
                str(path) for path in snakemake.input.partition_statistics
            ),
        ),
        LookupReductionOutputs(
            lookup_het=str(snakemake.output.lookup_het),
            lookup_hom=str(snakemake.output.lookup_hom),
            lookup_hem=(str(hem_output) if isinstance(hem_output, (str, Path)) else str(hem_output[0])) if hem_output else None,
            diagnostics=str(snakemake.output.diagnostics),
        ),
        LookupReductionParameters(
            bin_count=int(snakemake.params.bin_count),
            roles=tuple(snakemake.params.roles),
        ),
        str(snakemake.wildcards.bucket),
    )


if __name__ == "__main__":
    main()
