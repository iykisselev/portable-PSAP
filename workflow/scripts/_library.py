"""Make the workflow's reusable Python library importable to script adapters."""

from pathlib import Path
import sys


def add_library_path(scriptdir: str | Path) -> None:
    """Prepend the workflow's shared Python library to the import path.

    Args:
        scriptdir: Directory containing the executing Snakemake script.

    Side Effects:
        Inserts ``workflow/lib/python`` at the front of ``sys.path`` when it is
        not already present, so local workflow modules take precedence.
    """
    library_dir = Path(scriptdir).parent / "lib" / "python"
    library = str(library_dir)
    if library not in sys.path:
        sys.path.insert(0, library)
