"""Shared parsing and normalization for PED and sample-sheet metadata."""

from __future__ import annotations

import csv
from collections.abc import Callable, Iterator
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class SampleRow:
    """One decoded source row using the workflow's four metadata fields."""

    line_number: int
    family_id: str | None
    pid: str | None
    sex: str | None
    phenotype: str | None


@dataclass(frozen=True)
class NormalizedSample:
    """Trimmed metadata and its safe, processable manifest representation."""

    family_id: str
    pid: str
    sex: str
    phenotype: str
    phenotype_label: str | None

    @property
    def is_safe(self) -> bool:
        return (
            bool(self.pid)
            and "/" not in self.pid
            and "\\" not in self.pid
            and self.phenotype_label is not None
        )

    def manifest_row(self) -> dict[str, str]:
        if not self.is_safe or self.phenotype_label is None:
            raise ValueError("Only safe normalized samples can form manifest rows")
        return {
            "family_id": self.family_id or "NA",
            "pid": self.pid,
            "sex": self.sex if self.sex in {"1", "2"} else "NA",
            "phenotype": self.phenotype_label,
        }


class SampleMetadataFormatError(ValueError):
    """A format-level metadata issue with an optional physical line number."""

    def __init__(self, message: str, line_number: int | None = None) -> None:
        super().__init__(message)
        self.line_number = line_number


def normalize_sample(row: SampleRow) -> NormalizedSample:
    """Trim source cells and classify the workflow phenotype coding."""
    family_id = (row.family_id or "").strip()
    pid = (row.pid or "").strip()
    sex = (row.sex or "").strip()
    phenotype = (row.phenotype or "").strip()
    return NormalizedSample(
        family_id=family_id,
        pid=pid,
        sex=sex,
        phenotype=phenotype,
        phenotype_label={"1": "unaffected", "2": "affected"}.get(phenotype),
    )


def iter_sample_rows(
    path: str | Path,
    metadata_kind: str,
    *,
    errors: str = "strict",
    on_error: Callable[[SampleMetadataFormatError], None] | None = None,
) -> Iterator[SampleRow]:
    """Yield decoded PED or sample-sheet rows, optionally reporting bad rows.

    An error callback preserves validation's collect-and-continue behavior;
    callers without one receive the first malformed-record exception.
    """
    source = Path(path)

    def report(message: str, line_number: int | None = None) -> bool:
        issue = SampleMetadataFormatError(message, line_number)
        if on_error is None:
            raise issue
        on_error(issue)
        return False

    if metadata_kind == "ped":
        with source.open("rt", encoding="utf-8", errors=errors, newline="") as handle:
            for line_number, fields in enumerate(csv.reader(handle, delimiter="\t"), start=1):
                if not fields or all(not value.strip() for value in fields):
                    continue
                if len(fields) < 6:
                    report(f"has {len(fields)} columns; at least 6 are required", line_number)
                    continue
                yield SampleRow(line_number, fields[0], fields[1], fields[4], fields[5])
        return

    if metadata_kind == "samples":
        with source.open("rt", encoding="utf-8", errors=errors, newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            required = ("family_id", "pid", "sex", "phenotype")
            missing = [name for name in required if name not in (reader.fieldnames or [])]
            if missing:
                report("is missing required column(s): " + ", ".join(missing))
                return
            for line_number, fields in enumerate(reader, start=2):
                if not fields or not any((value or "").strip() for value in fields.values()):
                    continue
                yield SampleRow(
                    line_number,
                    fields.get("family_id"),
                    fields.get("pid"),
                    fields.get("sex"),
                    fields.get("phenotype"),
                )
        return

    raise ValueError(f"Unsupported sample metadata kind: {metadata_kind!r}")
