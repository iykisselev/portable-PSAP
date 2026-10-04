"""Atomic artifact publication and content identity helpers.

Writers in this module publish a complete sibling temporary file before
replacing the destination.  Generation jobs can therefore use the same
failure-safe semantics for manifests, Parquet intermediates, and compressed
bundle components.
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import io
import json
import os
from pathlib import Path
from typing import Iterable, Mapping, Sequence

import polars as pl


def sha256_file(path: str | Path) -> str:
    """Hash one completed artifact without loading it into memory."""
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _temporary_path(destination: Path, *, create_parents: bool = True) -> Path:
    if create_parents:
        destination.parent.mkdir(parents=True, exist_ok=True)
    return destination.with_suffix(destination.suffix + ".tmp")


def write_json_atomic(path: str | Path, value: object) -> None:
    """Write a stable JSON document only after its complete content is ready."""
    destination = Path(path)
    temporary = _temporary_path(destination)
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    os.replace(temporary, destination)


def write_text_atomic(
    path: str | Path, value: str, *, create_parents: bool = True
) -> Path:
    """Atomically publish text and return its destination path.

    Set ``create_parents`` to false when a missing output directory should be
    treated as a caller error rather than created implicitly.
    """
    destination = Path(path)
    temporary = _temporary_path(destination, create_parents=create_parents)
    temporary.write_text(value, encoding="utf-8")
    os.replace(temporary, destination)
    return destination


def write_parquet_atomic(frame: pl.DataFrame, path: str | Path, compression: str) -> None:
    """Write one Parquet artifact atomically with a declared compression codec."""
    destination = Path(path)
    temporary = _temporary_path(destination)
    frame.write_parquet(temporary, compression=compression)
    os.replace(temporary, destination)


def write_tsv_atomic(
    path: str | Path,
    header: Sequence[object],
    rows: Iterable[Sequence[object]],
) -> None:
    """Atomically write a canonical newline-delimited TSV intermediate."""
    destination = Path(path)
    temporary = _temporary_path(destination)
    with temporary.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)
    os.replace(temporary, destination)


class _DigestWriter:
    """Binary-file adapter that records bytes emitted by ``gzip``."""

    def __init__(self, handle, digest: "hashlib._Hash") -> None:
        self.handle = handle
        self.digest = digest

    def write(self, data: bytes) -> int:
        self.digest.update(data)
        return self.handle.write(data)

    def flush(self) -> None:
        self.handle.flush()

    def tell(self) -> int:
        return self.handle.tell()


def write_gzip_tsv_atomic(
    path: str | Path,
    header: Sequence[object],
    rows: Iterable[Sequence[object]],
) -> tuple[int, dict[str, object]]:
    """Atomically write a deterministic gzip TSV and its component inventory.

    Gzip's fixed ``mtime`` is required: bundle component checksums are part of
    the public certificate and must not vary between equivalent reruns.
    """
    destination = Path(path)
    temporary = _temporary_path(destination)
    digest = hashlib.sha256()
    count = 0
    with temporary.open("wb") as raw:
        digest_writer = _DigestWriter(raw, digest)
        with gzip.GzipFile(fileobj=digest_writer, mode="wb", mtime=0) as compressed:
            with io.TextIOWrapper(compressed, encoding="utf-8", newline="") as text:
                writer = csv.writer(text, delimiter="\t", lineterminator="\n")
                writer.writerow(header)
                for row in rows:
                    writer.writerow(row)
                    count += 1
    os.replace(temporary, destination)
    return count, component_inventory(destination, digest.hexdigest())


def component_inventory(path: str | Path, sha256: str | None = None) -> dict[str, object]:
    """Return the filename, compressed size, and SHA-256 used in bundle metadata."""
    component = Path(path)
    return {
        "path": component.name,
        "sha256": sha256 or sha256_file(component),
        "bytes": component.stat().st_size,
    }


def resource_inventory(paths: Mapping[str, str | Path]) -> dict[str, dict[str, str]]:
    """Build a provenance record for raw resources used by a generation run."""
    return {
        name: {"path": str(path), "sha256": sha256_file(path)}
        for name, path in paths.items()
    }
