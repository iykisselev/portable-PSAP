"""Make the workflow's reusable Python library importable to script adapters."""

from pathlib import Path
import sys


def add_library_path(scriptdir: str | Path) -> None:
    library_dir = Path(scriptdir).parent / "lib" / "python"
    library = str(library_dir)
    if library not in sys.path:
        sys.path.insert(0, library)
