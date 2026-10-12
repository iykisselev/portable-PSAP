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
from dataclasses import dataclass
from typing import Callable, Iterable, Iterator, Sequence, TextIO

import polars as pl


MISSING_TOKENS = frozenset({"", ".", "-", "na", "nan", "none", "null"})
"""Case-insensitive tokens denoting missing identifiers or table cells."""

PED_MISSING_SEX_TOKENS = MISSING_TOKENS | frozenset({"-9", "0"})
"""Additional conventional PED sex sentinels; not generic missing IDs."""


def is_missing_token(value: object) -> bool:
    """Return whether a table-cell or identifier value is contractually absent."""
    return value is None or str(value).strip().lower() in MISSING_TOKENS


@dataclass(frozen=True)
class BedRow:
    """One structurally valid zero-based, half-open BED record."""

    line_number: int
    fields: tuple[str, ...]
    chromosome: str
    start: int
    end: int


class BedFormatError(ValueError):
    """A malformed BED row, retaining its source line for diagnostics."""


def open_text(
    path: str | Path,
    *,
    errors: str = "strict",
    newline: str | None = None,
    require_gzip: bool = False,
) -> TextIO:
    """Open plain or gzip text using content detection rather than a suffix.

    ``require_gzip`` is for public bundle components whose compression is part
    of their file-format contract, rather than merely a transport detail.
    """
    resource = Path(path)
    with resource.open("rb") as handle:
        compressed = handle.read(2) == b"\x1f\x8b"
    if require_gzip and not compressed:
        raise ValueError(f"Expected a gzip-compressed text file: {resource}")
    if compressed:
        return gzip.open(
            resource, "rt", encoding="utf-8", errors=errors, newline=newline
        )
    return resource.open("rt", encoding="utf-8", errors=errors, newline=newline)


def iter_bed(
    path: str | Path,
    min_columns: int = 3,
    *,
    errors: str = "strict",
    on_error: Callable[[BedFormatError], None] | None = None,
) -> Iterator[BedRow]:
    """Yield valid BED rows after standard comment and coordinate checks.

    By default a malformed row stops the caller.  Validators may pass
    ``on_error`` to collect diagnostics and continue with later records.
    """
    if min_columns < 3:
        raise ValueError("BED parsing requires at least three columns")
    with open_text(path, errors=errors) as handle:
        for line_number, raw in enumerate(handle, start=1):
            stripped = raw.strip()
            if not stripped or stripped.startswith(("#", "track", "browser")):
                continue
            fields = tuple(stripped.split("\t"))
            try:
                if len(fields) < min_columns:
                    raise BedFormatError(
                        f"BED line {line_number} needs at least {min_columns} columns"
                    )
                chromosome = fields[0].strip()
                if not chromosome:
                    raise BedFormatError(
                        f"BED line {line_number} has an empty chromosome"
                    )
                try:
                    start, end = int(fields[1]), int(fields[2])
                except ValueError as exc:
                    raise BedFormatError(
                        f"BED line {line_number} has non-integer start or end coordinates"
                    ) from exc
                if start < 0 or end <= start:
                    raise BedFormatError(
                        f"BED line {line_number} must satisfy 0 <= start < end; "
                        f"received {start}, {end}"
                    )
            except BedFormatError as exc:
                if on_error is None:
                    raise
                on_error(exc)
                continue
            yield BedRow(line_number, fields, chromosome, start, end)


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
    temporary.write_text(
        json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
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


def write_parquet_atomic(
    frame: pl.DataFrame,
    path: str | Path,
    compression: str,
    *,
    row_group_size: int | None = None,
) -> None:
    """Write one Parquet artifact atomically with a declared compression codec."""
    destination = Path(path)
    temporary = _temporary_path(destination)
    frame.write_parquet(
        temporary, compression=compression, row_group_size=row_group_size
    )
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


class GzipTsvAtomicWriter:
    """Incrementally write one deterministic gzip TSV, then publish it atomically."""

    def __init__(self, path: str | Path, header: Sequence[object]):
        self.destination = Path(path)
        self.temporary = _temporary_path(self.destination)
        self.header = header
        self.count = 0
        self.inventory: dict[str, object] | None = None
        self._raw = None
        self._text = None
        self._writer = None
        self._digest = hashlib.sha256()

    def __enter__(self) -> "GzipTsvAtomicWriter":
        self._raw = self.temporary.open("wb")
        digest_writer = _DigestWriter(self._raw, self._digest)
        compressed = gzip.GzipFile(fileobj=digest_writer, mode="wb", mtime=0)
        self._text = io.TextIOWrapper(compressed, encoding="utf-8", newline="")
        self._writer = csv.writer(self._text, delimiter="\t", lineterminator="\n")
        self._writer.writerow(self.header)
        return self

    def write_row(self, row: Sequence[object]) -> None:
        """Append one row while keeping only the active row in memory."""
        assert self._writer is not None
        self._writer.writerow(row)
        self.count += 1

    def __exit__(self, exc_type, exc_value, traceback) -> bool:
        assert self._text is not None and self._raw is not None
        self._text.close()
        self._raw.close()
        if exc_type is not None:
            self.temporary.unlink(missing_ok=True)
            return False
        os.replace(self.temporary, self.destination)
        self.inventory = component_inventory(self.destination, self._digest.hexdigest())
        return False


def write_gzip_tsv_atomic(
    path: str | Path,
    header: Sequence[object],
    rows: Iterable[Sequence[object]],
) -> tuple[int, dict[str, object]]:
    """Atomically write a deterministic gzip TSV and its component inventory.

    Gzip's fixed ``mtime`` is required: bundle component checksums are part of
    the public certificate and must not vary between equivalent reruns.
    """
    with GzipTsvAtomicWriter(path, header) as writer:
        for row in rows:
            writer.write_row(row)
    assert writer.inventory is not None
    return writer.count, writer.inventory


def component_inventory(
    path: str | Path, sha256: str | None = None
) -> dict[str, object]:
    """Return the filename, compressed size, and SHA-256 used in bundle metadata."""
    component = Path(path)
    return {
        "path": component.name,
        "sha256": sha256 or sha256_file(component),
        "bytes": component.stat().st_size,
    }
