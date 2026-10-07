"""Shared runtime setup for standalone generation Snakemake adapters."""

import os
from pathlib import Path
import sys


def bootstrap(snakemake: object, *, configure_polars: bool = False) -> None:
    """Expose shared and workflow-library modules before adapter imports."""
    scriptdir = Path(snakemake.scriptdir)  # type: ignore[attr-defined]
    shared_dir = scriptdir.parents[1] / "shared"
    shared = str(shared_dir)
    if shared not in sys.path:
        sys.path.insert(0, shared)

    if configure_polars:
        os.environ.setdefault(
            "POLARS_MAX_THREADS", str(snakemake.threads)  # type: ignore[attr-defined]
        )

    from _library import add_library_path

    add_library_path(scriptdir)
