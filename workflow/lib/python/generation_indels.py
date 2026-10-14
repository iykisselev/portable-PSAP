"""Bounded, normalized indel catalog readers for generation windows."""

from __future__ import annotations

import math
from pathlib import Path
from typing import Protocol

import polars as pl

from generation_population import window_records
from genomic_contracts import (
    AlleleContractError,
    allele_key_text,
    is_left_aligned_allele,
    is_minimal_allele,
    reference_span,
    validated_allele_key,
    variant_class,
)


class ReferenceSequence(Protocol):
    """The bounded FASTA access required to validate one candidate row."""

    def fetch(self, chromosome: str, start: int, end: int) -> str:
        """Return uppercase reference bases for one zero-based half-open span."""


INDEL_SCORE_SCHEMA = {
    "CHROM": pl.String,
    "POS": pl.Int64,
    "REF": pl.String,
    "ALT": pl.String,
    "SPAN_START": pl.Int64,
    "SPAN_END": pl.Int64,
    "score": pl.Float64,
}

INDEL_MEMBERSHIP_SCHEMA = {
    "testing_unit_id": pl.String,
    "CHROM": pl.String,
    "POS": pl.Int64,
    "REF": pl.String,
    "ALT": pl.String,
    "SPAN_START": pl.Int64,
    "SPAN_END": pl.Int64,
    "variant_class": pl.String,
}

INDEL_REJECTION_SCHEMA = {
    "CHROM": pl.String,
    "POS": pl.Int64,
    "REF": pl.String,
    "ALT": pl.String,
    "SPAN_START": pl.Int64,
    "SPAN_END": pl.Int64,
    "rejection_reason": pl.String,
}


class IndelCandidateError(ValueError):
    """Raised when an owner-window indel catalog row is not usable."""


def _candidate_error(
    message: str, *, line_number: int | None = None
) -> IndelCandidateError:
    prefix = (
        f"Indel score resource line {line_number}: "
        if line_number
        else "Indel score resource: "
    )
    return IndelCandidateError(prefix + message)


def indel_score_rows(
    path: str | Path,
    index: str | Path | None,
    chromosome: str,
    start: int,
    end: int,
    reference: ReferenceSequence,
    *,
    max_reference_length: int,
    max_alternate_length: int,
) -> pl.DataFrame:
    """Read and validate catalog indels owned by one anchor-position window.

    The score catalog uses the same six-column CADD layout as the existing
    SNV reader.  Tabix is queried only for the owner window ``[start, end)``
    and each accepted row satisfies ``start < POS <= end``.  Genuine duplicate
    catalog rows are intentionally retained so the later window join can
    report them as a structural source error.
    """
    if max_reference_length < 1 or max_alternate_length < 1:
        raise ValueError("Indel length limits must be positive")
    records: list[dict[str, object]] = []
    for line_number, fields, key in window_records(
        path, index, chromosome, start, end, min_columns=6
    ):
        try:
            key = validated_allele_key(*key)
        except AlleleContractError as exc:
            raise _candidate_error(str(exc), line_number=line_number) from exc
        key_text = allele_key_text(key)
        if variant_class(key[2], key[3]) != "indel":
            raise _candidate_error(
                f"{key_text} is not a length-changing A/C/G/T indel",
                line_number=line_number,
            )
        if len(key[2]) > max_reference_length or len(key[3]) > max_alternate_length:
            raise _candidate_error(
                f"{key_text} exceeds configured indel length limits "
                f"REF<={max_reference_length}, ALT<={max_alternate_length}",
                line_number=line_number,
            )
        if not is_minimal_allele(key[2], key[3]):
            raise _candidate_error(
                f"{key_text} is not minimal", line_number=line_number
            )
        if not is_left_aligned_allele(key[1], key[2], key[3]):
            raise _candidate_error(
                f"{key_text} is not left-aligned", line_number=line_number
            )
        span_start, span_end = reference_span(key[1], key[2])
        try:
            observed_reference = reference.fetch(key[0], span_start, span_end)
        except (KeyError, ValueError) as exc:
            raise _candidate_error(
                f"{key_text} has a reference span outside the configured FASTA",
                line_number=line_number,
            ) from exc
        if observed_reference != key[2]:
            raise _candidate_error(
                f"{key_text} REF does not match the configured FASTA",
                line_number=line_number,
            )
        try:
            score = float(fields[5])
        except ValueError as exc:
            raise _candidate_error(
                f"{key_text} has a nonnumeric score", line_number=line_number
            ) from exc
        if not math.isfinite(score):
            raise _candidate_error(
                f"{key_text} has a non-finite score", line_number=line_number
            )
        records.append(
            {
                "CHROM": key[0],
                "POS": key[1],
                "REF": key[2],
                "ALT": key[3],
                "SPAN_START": span_start,
                "SPAN_END": span_end,
                "score": score,
            }
        )
    return pl.DataFrame(records, schema=INDEL_SCORE_SCHEMA, strict=False)


def indel_candidate_memberships(
    candidates: pl.DataFrame,
    mask: list[tuple[int, int]],
    unit_intervals: list[tuple[str, int, int]],
) -> tuple[pl.DataFrame, pl.DataFrame]:
    """Assign accepted catalog candidates to mask-contained overlapping units.

    ``mask`` and ``unit_intervals`` must already be restricted to the
    candidate chromosome.  A candidate's full REF span must fit in one mask
    interval; membership uses any reference-span overlap, allowing a deletion
    anchored just before a testing unit to belong to that unit.  Candidates
    not retained by either rule are returned once in the audit frame.
    """
    required = set(INDEL_SCORE_SCHEMA)
    missing = sorted(required - set(candidates.columns))
    if missing:
        raise ValueError("Indel candidates are missing columns: " + ", ".join(missing))
    memberships: list[dict[str, object]] = []
    rejected: list[dict[str, object]] = []
    for candidate in candidates.select(*INDEL_SCORE_SCHEMA).iter_rows(named=True):
        try:
            key = validated_allele_key(
                candidate["CHROM"], candidate["POS"], candidate["REF"], candidate["ALT"]
            )
        except AlleleContractError as exc:
            raise ValueError("Indel candidate has an invalid allele key") from exc
        if variant_class(key[2], key[3]) != "indel":
            raise ValueError("Indel candidate has a non-indel allele class")
        span_start = int(candidate["SPAN_START"])
        span_end = int(candidate["SPAN_END"])
        if (span_start, span_end) != reference_span(key[1], key[2]):
            raise ValueError(
                "Indel candidate has a reference span inconsistent with POS and REF"
            )
        audit_row = {
            "CHROM": key[0],
            "POS": key[1],
            "REF": key[2],
            "ALT": key[3],
            "SPAN_START": span_start,
            "SPAN_END": span_end,
        }
        if not any(
            mask_start <= span_start and span_end <= mask_end
            for mask_start, mask_end in mask
        ):
            rejected.append(
                {
                    **audit_row,
                    "rejection_reason": "span_not_contained_in_calibration_mask",
                }
            )
            continue
        matched_units = [
            unit_id
            for unit_id, unit_start, unit_end in unit_intervals
            if span_start < unit_end and span_end > unit_start
        ]
        if not matched_units:
            rejected.append(
                {**audit_row, "rejection_reason": "no_eligible_testing_unit_overlap"}
            )
            continue
        for unit_id in matched_units:
            memberships.append(
                {
                    "testing_unit_id": unit_id,
                    **audit_row,
                    "variant_class": "indel",
                }
            )
    return (
        pl.DataFrame(memberships, schema=INDEL_MEMBERSHIP_SCHEMA, strict=False).unique(
            maintain_order=True
        ),
        pl.DataFrame(rejected, schema=INDEL_REJECTION_SCHEMA, strict=False).unique(
            maintain_order=True
        ),
    )
