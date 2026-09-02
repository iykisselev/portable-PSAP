"""Canonical, immutable sample indexes for apply-rule expansion and validation."""

from __future__ import annotations

import csv
from dataclasses import dataclass
import io
from pathlib import Path
from types import MappingProxyType
from typing import Iterable, Mapping

from validation_reporting import write_atomic_text


SAMPLE_MANIFEST_COLUMNS = ("family_id", "pid", "sex", "phenotype")


@dataclass(frozen=True)
class SampleRecord:
    """One normalized, processable sample from the canonical manifest."""

    family_id: str
    pid: str
    sex: str
    phenotype: str


@dataclass(frozen=True)
class SampleIndex:
    """One immutable, order-preserving view of canonical sample metadata."""

    records: tuple[SampleRecord, ...]
    by_pid: Mapping[str, SampleRecord]

    @classmethod
    def from_rows(cls, rows: Iterable[Mapping[str, str]]) -> "SampleIndex":
        records = tuple(
            SampleRecord(
                family_id=str(row["family_id"]),
                pid=str(row["pid"]),
                sex=str(row["sex"]),
                phenotype=str(row["phenotype"]),
            )
            for row in rows
        )
        pids = [record.pid for record in records]
        if any(not pid for pid in pids):
            raise ValueError("Sample manifest contains an empty pid value")
        if len(set(pids)) != len(pids):
            raise ValueError("Sample manifest contains duplicate pid values")
        return cls(
            records=records,
            by_pid=MappingProxyType({record.pid: record for record in records}),
        )

    def sample_ids(self, phenotype: str) -> list[str]:
        """Return manifest-order PIDs belonging to one normalized phenotype."""
        return [record.pid for record in self.records if record.phenotype == phenotype]

    def sample_sexes(self, phenotype: str) -> list[str]:
        """Return manifest-order sex values belonging to one normalized phenotype."""
        return [record.sex for record in self.records if record.phenotype == phenotype]

    def sample(self, pid: str, expected_phenotype: str | None = None) -> SampleRecord:
        """Return one sample record, optionally enforcing its phenotype role."""
        try:
            record = self.by_pid[str(pid)]
        except KeyError as exc:
            raise ValueError(f"Sample {pid!r} is not uniquely present in the manifest") from exc
        if expected_phenotype is not None and record.phenotype != expected_phenotype:
            raise ValueError(
                f"Sample {pid!r} has phenotype {record.phenotype!r}; "
                f"expected {expected_phenotype!r}"
            )
        return record


def read_sample_manifest(path: str | Path) -> "pd.DataFrame":
    """Load the canonical manifest and reject malformed or duplicate samples."""
    import pandas as pd

    manifest = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    missing = set(SAMPLE_MANIFEST_COLUMNS) - set(manifest.columns)
    if missing:
        raise ValueError(
            f"Sample manifest {path!s} has columns {manifest.columns.tolist()!r}; "
            "missing: " + ", ".join(sorted(missing))
        )
    manifest = manifest.loc[:, SAMPLE_MANIFEST_COLUMNS]
    SampleIndex.from_rows(manifest.to_dict("records"))
    return manifest


def read_sample_manifest_index(path: str | Path) -> SampleIndex:
    """Read and index the physical canonical manifest exactly once."""
    with Path(path).open("rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        headers = tuple(reader.fieldnames or ())
        missing = set(SAMPLE_MANIFEST_COLUMNS) - set(headers)
        if missing:
            raise ValueError(
                f"Sample manifest {path!s} has columns {list(headers)!r}; missing: "
                + ", ".join(sorted(missing))
            )
        return SampleIndex.from_rows(
            {
                column: (row.get(column) or "")
                for column in SAMPLE_MANIFEST_COLUMNS
            }
            for row in reader
        )


def write_sample_manifest(
    rows: Iterable[Mapping[str, str]],
    path: str | Path,
    *,
    create_parents: bool = True,
) -> Path:
    """Atomically publish the canonical manifest consumed by cohort validation."""
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(
        stream,
        fieldnames=SAMPLE_MANIFEST_COLUMNS,
        delimiter="\t",
        lineterminator="\n",
    )
    writer.writeheader()
    writer.writerows(rows)
    return write_atomic_text(path, stream.getvalue(), create_parents=create_parents)


def configured_sample_index(path: str | Path, metadata_kind: str) -> SampleIndex:
    """Resolve safe, processable samples from a configured PED or sample TSV.

    This intentionally tolerates invalid source rows: the normal validation
    marker remains responsible for reporting malformed metadata and prevents
    every expensive descendant from starting.  The planning index only needs a
    deterministic, duplicate-free set of safe wildcard values.
    """
    records: list[dict[str, str]] = []
    seen: set[str] = set()
    source = Path(path)

    try:
        if metadata_kind == "ped":
            with source.open("rt", encoding="utf-8", errors="replace", newline="") as handle:
                for row in csv.reader(handle, delimiter="\t"):
                    if len(row) < 6:
                        continue
                    _append_configured_record(records, seen, row[0], row[1], row[4], row[5])
        elif metadata_kind == "samples":
            with source.open("rt", encoding="utf-8", errors="replace", newline="") as handle:
                for row in csv.DictReader(handle, delimiter="\t"):
                    _append_configured_record(
                        records,
                        seen,
                        row.get("family_id", ""),
                        row.get("pid", ""),
                        row.get("sex", ""),
                        row.get("phenotype", ""),
                    )
        else:
            raise ValueError(f"Unsupported sample metadata kind: {metadata_kind!r}")
    except OSError:
        # The validation rule emits the actionable error. An empty planning
        # index keeps Snakemake from expanding an unsafe partial workflow.
        return SampleIndex.from_rows(())

    return SampleIndex.from_rows(records)


def _append_configured_record(
    records: list[dict[str, str]],
    seen: set[str],
    family_id: str | None,
    pid: str | None,
    sex: str | None,
    phenotype: str | None,
) -> None:
    """Append one safe processing record, retaining first source occurrence."""
    sample_id = (pid or "").strip()
    label = {"1": "unaffected", "2": "affected"}.get((phenotype or "").strip())
    if (
        not sample_id
        or sample_id in seen
        or "/" in sample_id
        or "\\" in sample_id
        or label is None
    ):
        return
    seen.add(sample_id)
    records.append(
        {
            "family_id": (family_id or "").strip() or "NA",
            "pid": sample_id,
            "sex": (sex or "").strip() if (sex or "").strip() in {"1", "2"} else "NA",
            "phenotype": label,
        }
    )


def _as_index(manifest: "pd.DataFrame" | SampleIndex) -> SampleIndex:
    if isinstance(manifest, SampleIndex):
        return manifest
    return SampleIndex.from_rows(manifest.to_dict("records"))


def sample_metadata(
    manifest: "pd.DataFrame" | SampleIndex,
    pid: str,
    expected_phenotype: str | None = None,
) -> "pd.Series" | SampleRecord:
    """Return one sample record, optionally checking its phenotype role."""
    if isinstance(manifest, SampleIndex):
        return manifest.sample(pid, expected_phenotype)
    _as_index(manifest).sample(pid, expected_phenotype)
    return manifest.loc[manifest["pid"] == str(pid)].iloc[0]


def sample_ids_by_phenotype(
    manifest: "pd.DataFrame" | SampleIndex, phenotype: str
) -> list[str]:
    """Return stable sample IDs for one workflow phenotype role."""
    return _as_index(manifest).sample_ids(phenotype)


def sample_sexes_by_phenotype(
    manifest: "pd.DataFrame" | SampleIndex, phenotype: str
) -> list[str]:
    """Return stable sex values for one workflow phenotype role."""
    return _as_index(manifest).sample_sexes(phenotype)
