"""Prevent repeated whole-file scans of production generation resources."""

from pathlib import Path


MAX_UNINDEXED_FIXTURE_BYTES = 1024 * 1024
COMPRESSED_SUFFIXES = {".gz", ".bgz", ".bgzf", ".bz2", ".xz", ".zip", ".bcf"}


def check_generation_index_binding(
    path: str | Path,
    index: str | None,
    *,
    key: str,
    allow_unindexed_fixtures: bool = False,
) -> None:
    """Validate configuration without reading resources during DAG construction."""
    if index:
        return
    if not allow_unindexed_fixtures or Path(path).suffix.lower() in COMPRESSED_SUFFIXES:
        raise ValueError(
            f"{key}.index is required for generation resource {path!s}. "
            "Declare the matching .tbi or .csi index explicitly; adjacent indexes "
            "are not auto-discovered. Without it, each partition scans the entire file. "
            "Only uncompressed fixtures up to 1 MiB may use "
            "generation.execution.allow_unindexed_fixture_inputs: true."
        )


def check_unindexed_fixture(path: str | Path) -> None:
    """Bound the fixture fallback before opening a sequential resource stream."""
    resource = Path(path)
    if (
        resource.suffix.lower() in COMPRESSED_SUFFIXES
        or resource.stat().st_size > MAX_UNINDEXED_FIXTURE_BYTES
    ):
        raise ValueError(
            f"Generation resource {resource!s} requires an explicit .tbi or .csi index. "
            "The unindexed fallback is limited to uncompressed fixtures up to 1 MiB."
        )
    with resource.open("rb") as handle:
        if handle.read(2) == b"\x1f\x8b":
            raise ValueError(
                f"Compressed generation resource {resource!s} requires an explicit .tbi or .csi index."
            )
