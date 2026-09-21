"""Make the workflow's Python library available to script adapters."""

from pathlib import Path
import sys


def workflow_root(scriptdir: str | Path) -> Path:
    """Find the workflow root from either a flat or nested script directory."""
    start = Path(scriptdir).resolve()
    for candidate in (start, *start.parents):
        if (candidate / "lib" / "python").is_dir():
            return candidate
    raise ValueError(f"Could not find workflow/lib/python above script directory {start}")


def add_library_path(scriptdir: str | Path) -> None:
    """Prepend the workflow's shared Python library to the import path.

    Args:
        scriptdir: Directory containing the executing Snakemake script.

    Side Effects:
        Inserts ``workflow/lib/python`` at the front of ``sys.path`` when it is
        not already present, so local workflow modules take precedence.
    """
    library_dir = workflow_root(scriptdir) / "lib" / "python"
    library = str(library_dir)
    if library not in sys.path:
        sys.path.insert(0, library)
