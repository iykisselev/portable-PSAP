"""Dependency-free helpers for declaring optional workflow file inputs."""

from pathlib import Path


def as_inputs(value: str | Path | None) -> list[str]:
    """Represent an optional scalar path as zero or one declared inputs."""
    return [str(value)] if value else []
