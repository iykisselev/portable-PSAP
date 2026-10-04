"""Small, shell-safe argument builders used by workflow rules."""

from __future__ import annotations

from pathlib import Path
from shlex import quote


def optional_regions_argument(path: str | Path | None) -> str:
    """Return a quoted bcftools region-file argument when a path is configured."""
    if path is None:
        return ""
    return f"--regions-file {quote(str(path))}"
