"""Versioned, validated manifests for data-dependent generation topology."""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
import json
from pathlib import Path


PARTITION_MANIFEST_SCHEMA_VERSION = 1
BUCKET_MANIFEST_SCHEMA_VERSION = 2


class GenerationManifestError(ValueError):
    """Raised when a generated topology manifest is malformed or stale."""


@dataclass(frozen=True)
class PartitionManifestEntry:
    """One canonical half-open genomic window published by scope preparation."""

    identifier: str
    chromosome: str
    start: int
    end: int
    ordinal: int


@dataclass(frozen=True)
class BucketManifestEntry:
    """One deterministic, sorted range of eligible testing-unit identifiers."""

    identifier: str
    ordinal: int
    unit_ids: tuple[str, ...]
    partition_ids: tuple[str, ...]


def _load_mapping(path: str | Path, label: str) -> dict:
    try:
        value = json.loads(Path(path).read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise GenerationManifestError(
            f"Unable to read {label} manifest {path!s}"
        ) from exc
    if not isinstance(value, dict):
        raise GenerationManifestError(
            f"{label.capitalize()} manifest must be a JSON object"
        )
    return value


def _partition_identifier_sort_key(identifier: str) -> tuple[str, int, int]:
    """Return the ordering used by the partition producer for one identifier."""
    fields = identifier.rsplit("_", 2)
    if len(fields) != 3:
        raise ValueError("Partition identifier has the wrong shape")
    chromosome, raw_start, raw_end = fields
    start, end = int(raw_start), int(raw_end)
    if (
        not chromosome
        or start < 0
        or end <= start
        or identifier != f"{chromosome}_{start:012d}_{end:012d}"
    ):
        raise ValueError("Partition identifier is noncanonical")
    return chromosome, start, end


def _read_partition_manifest_uncached(
    path: str | Path,
) -> tuple[PartitionManifestEntry, ...]:
    """Read and validate the compact partition-discovery manifest."""
    manifest = _load_mapping(path, "partition")
    if manifest.get("schema_version") != PARTITION_MANIFEST_SCHEMA_VERSION:
        raise GenerationManifestError(
            "Unsupported generation partition manifest schema"
        )
    rows = manifest.get("partitions")
    if not isinstance(rows, list):
        raise GenerationManifestError("Partition manifest needs a partitions array")

    entries: list[PartitionManifestEntry] = []
    identifiers: set[str] = set()
    previous: PartitionManifestEntry | None = None
    for ordinal, row in enumerate(rows):
        if not isinstance(row, dict):
            raise GenerationManifestError("Partition manifest entries must be objects")
        try:
            entry = PartitionManifestEntry(
                identifier=str(row["partition_id"]),
                chromosome=str(row["chromosome"]),
                start=int(row["start"]),
                end=int(row["end"]),
                ordinal=int(row["ordinal"]),
            )
        except (KeyError, TypeError, ValueError) as exc:
            raise GenerationManifestError(
                "Partition manifest has an invalid entry"
            ) from exc
        expected_identifier = f"{entry.chromosome}_{entry.start:012d}_{entry.end:012d}"
        if (
            entry.ordinal != ordinal
            or entry.identifier != expected_identifier
            or entry.start < 0
            or entry.end <= entry.start
            or entry.identifier in identifiers
        ):
            raise GenerationManifestError("Partition manifest has noncanonical entries")
        if (
            previous
            and entry.chromosome == previous.chromosome
            and entry.start < previous.end
        ):
            raise GenerationManifestError(
                "Partition manifest contains overlapping windows"
            )
        identifiers.add(entry.identifier)
        entries.append(entry)
        previous = entry
    return tuple(entries)


def _read_bucket_manifest_uncached(path: str | Path) -> tuple[BucketManifestEntry, ...]:
    """Read and validate the deterministic eligible-unit bucket manifest."""
    manifest = _load_mapping(path, "bucket")
    if manifest.get("schema_version") != BUCKET_MANIFEST_SCHEMA_VERSION:
        raise GenerationManifestError("Unsupported generation bucket manifest schema")
    rows = manifest.get("buckets")
    if not isinstance(rows, list):
        raise GenerationManifestError("Bucket manifest needs a buckets array")

    entries: list[BucketManifestEntry] = []
    seen_units: set[str] = set()
    for ordinal, row in enumerate(rows):
        if not isinstance(row, dict):
            raise GenerationManifestError("Bucket manifest entries must be objects")
        try:
            identifier = str(row["bucket_id"])
            entry_ordinal = int(row["ordinal"])
            unit_ids = tuple(str(value) for value in row["unit_ids"])
            partition_ids = tuple(str(value) for value in row["partition_ids"])
            partition_keys = tuple(
                _partition_identifier_sort_key(value) for value in partition_ids
            )
        except (KeyError, TypeError, ValueError) as exc:
            raise GenerationManifestError(
                "Bucket manifest has an invalid entry"
            ) from exc
        if (
            entry_ordinal != ordinal
            or identifier != f"{ordinal:05d}"
            or not unit_ids
            or tuple(sorted(unit_ids)) != unit_ids
            or len(set(unit_ids)) != len(unit_ids)
            or seen_units.intersection(unit_ids)
            or not partition_ids
            or tuple(sorted(partition_keys)) != partition_keys
            or len(set(partition_ids)) != len(partition_ids)
        ):
            raise GenerationManifestError("Bucket manifest has noncanonical entries")
        seen_units.update(unit_ids)
        entries.append(
            BucketManifestEntry(
                identifier=identifier,
                ordinal=entry_ordinal,
                unit_ids=unit_ids,
                partition_ids=partition_ids,
            )
        )
    return tuple(entries)


def _manifest_cache_key(path: str | Path) -> tuple[str, int | None, int | None]:
    """Identify a published manifest while invalidating after atomic replacement."""
    manifest_path = Path(path).resolve()
    try:
        status = manifest_path.stat()
    except OSError:
        # Preserve the public reader's GenerationManifestError for unreadable
        # paths by letting the cached parser delegate to ``_load_mapping``.
        return str(manifest_path), None, None
    return str(manifest_path), status.st_mtime_ns, status.st_size


@lru_cache(maxsize=8)
def _cached_partition_manifest(
    path: str,
    mtime_ns: int | None,
    size: int | None,
) -> tuple[PartitionManifestEntry, ...]:
    return _read_partition_manifest_uncached(path)


@lru_cache(maxsize=8)
def _cached_bucket_manifest(
    path: str,
    mtime_ns: int | None,
    size: int | None,
) -> tuple[BucketManifestEntry, ...]:
    return _read_bucket_manifest_uncached(path)


@lru_cache(maxsize=8)
def _cached_bucket_index(
    path: str,
    mtime_ns: int | None,
    size: int | None,
) -> dict[str, BucketManifestEntry]:
    return {
        entry.identifier: entry
        for entry in _cached_bucket_manifest(path, mtime_ns, size)
    }


def read_partition_manifest(path: str | Path) -> tuple[PartitionManifestEntry, ...]:
    """Read and validate the compact partition-discovery manifest once per version."""
    return _cached_partition_manifest(*_manifest_cache_key(path))


def read_bucket_manifest(path: str | Path) -> tuple[BucketManifestEntry, ...]:
    """Read and validate the eligible-unit bucket manifest once per version."""
    return _cached_bucket_manifest(*_manifest_cache_key(path))


def partition_ids(path: str | Path) -> tuple[str, ...]:
    """Return the stable wildcard values from a validated partition manifest."""
    return tuple(entry.identifier for entry in read_partition_manifest(path))


def bucket_ids(path: str | Path) -> tuple[str, ...]:
    """Return the stable wildcard values from a validated bucket manifest."""
    return tuple(entry.identifier for entry in read_bucket_manifest(path))


def bucket_from_manifest(path: str | Path, identifier: str) -> BucketManifestEntry:
    """Resolve one bucket wildcard only after validating the full manifest."""
    try:
        return _cached_bucket_index(*_manifest_cache_key(path))[str(identifier)]
    except KeyError as exc:
        raise GenerationManifestError(
            f"Unknown generation bucket {identifier!r}"
        ) from exc
