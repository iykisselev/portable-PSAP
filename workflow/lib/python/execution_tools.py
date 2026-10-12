"""Executable discovery shared by workflow libraries and script adapters."""

import shutil


def require_tool(name: str) -> str:
    executable = shutil.which(name)
    if executable is None:
        raise RuntimeError(f"Required executable is unavailable: {name}")
    return executable
