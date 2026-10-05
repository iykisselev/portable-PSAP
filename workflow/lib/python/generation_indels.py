"""Bounded, normalized indel catalog readers for future null generation.

This module is intentionally not wired into the active generator yet.  It
provides the candidate contract that later preparation stages will consume
once the generated-bundle and apply contracts support indels end-to-end.
"""

from __future__ import annotations

import math
from pathlib import Path
from typing import Protocol

import polars as pl

from generation_population import POPULATION_SCHEMA, indexed_lines
from genomic_contracts import (
    AlleleContractError,
    allele_key_text,
    canonical_chromosome,
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

INDEL_MUTATION_RATE_SCHEMA = {
    "CHROM": pl.String,
    "POS": pl.Int64,
    "REF": pl.String,
    "ALT": pl.String,
    "mu": pl.Float64,
}

INDEL_MUTATION_ASSIGNMENT_SCHEMA = {
    "CHROM": pl.String,
    "POS": pl.Int64,
    "REF": pl.String,
    "ALT": pl.String,
    "mu": pl.Float64,
    "assignment_source": pl.String,
}


class IndelCandidateError(ValueError):
    """Raised when an owner-partition indel catalog row is not usable."""


def _candidate_error(message: str, *, line_number: int | None = None) -> IndelCandidateError:
    prefix = f"Indel score resource line {line_number}: " if line_number else "Indel score resource: "
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
    SNV reader.  Tabix is queried only for the owner partition ``[start, end)``
    and each accepted row satisfies ``start < POS <= end``.  Genuine duplicate
    catalog rows are intentionally retained so the later partition join can
    report them as a structural source error.
    """
    if max_reference_length < 1 or max_alternate_length < 1:
        raise ValueError("Indel length limits must be positive")
    chromosome = canonical_chromosome(chromosome)
    records: list[dict[str, object]] = []
    for line_number, raw in enumerate(indexed_lines(path, index, chromosome, start, end), start=1):
        if not raw.strip() or raw.startswith("#"):
            continue
        fields = raw.rstrip("\n").split("\t")
        if len(fields) < 6:
            raise _candidate_error("needs at least six CADD columns", line_number=line_number)
        record_chromosome = canonical_chromosome(fields[0])
        try:
            position = int(fields[1])
        except ValueError as exc:
            if record_chromosome == chromosome:
                raise _candidate_error("has a non-integer POS", line_number=line_number) from exc
            continue
        # Plain fixture reads scan the small whole resource, whereas indexed
        # production reads already return this owner window.
        if record_chromosome != chromosome or not start < position <= end:
            continue
        try:
            key = validated_allele_key(record_chromosome, position, fields[2], fields[3])
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
            raise _candidate_error(f"{key_text} is not minimal", line_number=line_number)
        if not is_left_aligned_allele(key[1], key[2], key[3]):
            raise _candidate_error(f"{key_text} is not left-aligned", line_number=line_number)
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
            raise _candidate_error(f"{key_text} has a nonnumeric score", line_number=line_number) from exc
        if not math.isfinite(score):
            raise _candidate_error(f"{key_text} has a non-finite score", line_number=line_number)
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
            raise ValueError("Indel candidate has a reference span inconsistent with POS and REF")
        audit_row = {
            "CHROM": key[0],
            "POS": key[1],
            "REF": key[2],
            "ALT": key[3],
            "SPAN_START": span_start,
            "SPAN_END": span_end,
        }
        if not any(mask_start <= span_start and span_end <= mask_end for mask_start, mask_end in mask):
            rejected.append({**audit_row, "rejection_reason": "span_not_contained_in_calibration_mask"})
            continue
        matched_units = [
            unit_id
            for unit_id, unit_start, unit_end in unit_intervals
            if span_start < unit_end and span_end > unit_start
        ]
        if not matched_units:
            rejected.append({**audit_row, "rejection_reason": "no_eligible_testing_unit_overlap"})
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


def indel_population_states(candidates: pl.DataFrame, panel: pl.DataFrame) -> pl.DataFrame:
    """Join validated catalog indels to exact population allele states.

    This preserves the population reader's AF-zero-imputation provenance.  An
    unmatched catalog allele is explicitly ``absent``; a matched invalid AF
    remains a hard error, while non-PASS records remain ``rejected`` for the
    later aggregation stage to exclude.
    """
    allele_columns = ("CHROM", "POS", "REF", "ALT")
    missing_candidates = sorted(set(allele_columns) - set(candidates.columns))
    missing_panel = sorted(set(POPULATION_SCHEMA) - set(panel.columns))
    if missing_candidates:
        raise ValueError("Indel candidates are missing columns: " + ", ".join(missing_candidates))
    if missing_panel:
        raise ValueError("Population states are missing columns: " + ", ".join(missing_panel))
    keys = candidates.select(*allele_columns)
    if keys.is_duplicated().any():
        raise IndelCandidateError("Indel catalog has duplicate candidate allele keys")
    matching_panel = panel.join(keys, on=list(allele_columns), how="inner")
    if matching_panel.select(*allele_columns).is_duplicated().any():
        raise IndelCandidateError("Population VCF has duplicate states for an indel candidate")
    states = (
        keys.join(panel, on=list(allele_columns), how="left")
        .with_columns(
            pl.when(pl.col("allele_state").is_null())
            .then(pl.lit("absent"))
            .otherwise(pl.col("allele_state"))
            .alias("allele_state"),
            pl.when(pl.col("state_reason").is_null())
            .then(pl.lit("not_present_in_population_panel"))
            .otherwise(pl.col("state_reason"))
            .alias("state_reason"),
        )
        .select(*POPULATION_SCHEMA)
    )
    invalid = states.filter(pl.col("allele_state") == "invalid")
    if invalid.height:
        examples = "; ".join(
            f"{chrom}:{position}:{reference}:{alternate}: {reason}"
            for chrom, position, reference, alternate, reason in invalid.head(10).select(
                *allele_columns, "state_reason"
            ).iter_rows()
        )
        raise IndelCandidateError(
            f"Population VCF has invalid selected AF for {invalid.height} indel candidate(s): {examples}"
        )
    return states


def indel_mutation_rate_rows(
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
    """Read normalized, allele-specific indel mutation rates by owner window."""
    if max_reference_length < 1 or max_alternate_length < 1:
        raise ValueError("Indel length limits must be positive")
    chromosome = canonical_chromosome(chromosome)
    records: list[dict[str, object]] = []
    for line_number, raw in enumerate(indexed_lines(path, index, chromosome, start, end), start=1):
        if not raw.strip() or raw.startswith("#"):
            continue
        fields = raw.rstrip("\n").split("\t")
        if len(fields) < 5:
            raise IndelCandidateError(
                f"Indel mutation-rate resource row {line_number} needs CHROM, POS, REF, ALT, and mu"
            )
        record_chromosome = canonical_chromosome(fields[0])
        try:
            position = int(fields[1])
        except ValueError as exc:
            if record_chromosome == chromosome:
                raise IndelCandidateError(
                    f"Indel mutation-rate resource row {line_number} has a non-integer POS"
                ) from exc
            continue
        if record_chromosome != chromosome or not start < position <= end:
            continue
        try:
            key = validated_allele_key(record_chromosome, position, fields[2], fields[3])
        except AlleleContractError as exc:
            raise IndelCandidateError(
                f"Indel mutation-rate resource row {line_number}: {exc}"
            ) from exc
        key_text = allele_key_text(key)
        if variant_class(key[2], key[3]) != "indel":
            raise IndelCandidateError(
                f"Indel mutation-rate resource row {line_number}: {key_text} is not an indel"
            )
        if len(key[2]) > max_reference_length or len(key[3]) > max_alternate_length:
            raise IndelCandidateError(
                f"Indel mutation-rate resource row {line_number}: {key_text} exceeds configured indel length limits"
            )
        if not is_minimal_allele(key[2], key[3]) or not is_left_aligned_allele(
            key[1], key[2], key[3]
        ):
            raise IndelCandidateError(
                f"Indel mutation-rate resource row {line_number}: {key_text} is not normalized"
            )
        span_start, span_end = reference_span(key[1], key[2])
        if reference.fetch(key[0], span_start, span_end) != key[2]:
            raise IndelCandidateError(
                f"Indel mutation-rate resource row {line_number}: {key_text} REF does not match the configured FASTA"
            )
        try:
            rate = float(fields[4])
        except ValueError as exc:
            raise IndelCandidateError(
                f"Indel mutation-rate resource row {line_number}: {key_text} has a nonnumeric mu"
            ) from exc
        if not math.isfinite(rate) or not 0 <= rate <= 1:
            raise IndelCandidateError(
                f"Indel mutation-rate resource row {line_number}: {key_text} has invalid mu"
            )
        records.append({"CHROM": key[0], "POS": key[1], "REF": key[2], "ALT": key[3], "mu": rate})
    return pl.DataFrame(records, schema=INDEL_MUTATION_RATE_SCHEMA, strict=False)


def indel_mutation_assignments(
    states: pl.DataFrame,
    rates: pl.DataFrame | None,
    policy: str,
) -> pl.DataFrame:
    """Assign zero or a class-specific mutation rate to absent indels only."""
    allele_columns = ("CHROM", "POS", "REF", "ALT")
    required_states = set(allele_columns) | {"allele_state"}
    missing_states = sorted(required_states - set(states.columns))
    if missing_states:
        raise ValueError("Indel states are missing columns: " + ", ".join(missing_states))
    if policy not in {"zero", "mutation_rate"}:
        raise ValueError(f"Unsupported indel unobserved policy {policy!r}")
    if states.select(*allele_columns).is_duplicated().any():
        raise ValueError("Indel states contain duplicate allele keys")
    if policy == "mutation_rate":
        if rates is None:
            raise IndelCandidateError("Indel mutation-rate fallback requires a separate mutation-rate resource")
        missing_rates = sorted(set(INDEL_MUTATION_RATE_SCHEMA) - set(rates.columns))
        if missing_rates:
            raise ValueError("Indel mutation rates are missing columns: " + ", ".join(missing_rates))
        if rates.select(*allele_columns).is_duplicated().any():
            raise IndelCandidateError("Indel mutation-rate resource has duplicate allele keys")
        joined = states.join(rates.select(*allele_columns, "mu"), on=list(allele_columns), how="left")
        absent_missing = joined.filter((pl.col("allele_state") == "absent") & pl.col("mu").is_null())
        if absent_missing.height:
            examples = ", ".join(
                f"{chrom}:{position}:{reference}:{alternate}"
                for chrom, position, reference, alternate in absent_missing.head(10).select(*allele_columns).iter_rows()
            )
            raise IndelCandidateError(
                f"Indel mutation-rate resource is missing {absent_missing.height} absent candidate(s): {examples}"
            )
        return joined.select(
            *allele_columns,
            pl.when(pl.col("allele_state") == "absent").then(pl.col("mu")).otherwise(pl.lit(None)).alias("mu"),
            pl.when(pl.col("allele_state") == "absent")
            .then(pl.lit("absent_mutation_rate"))
            .otherwise(pl.lit("observed_or_rejected"))
            .alias("assignment_source"),
        ).cast(INDEL_MUTATION_ASSIGNMENT_SCHEMA)
    return states.select(
        *allele_columns,
        pl.when(pl.col("allele_state") == "absent").then(pl.lit(0.0)).otherwise(pl.lit(None)).alias("mu"),
        pl.when(pl.col("allele_state") == "absent")
        .then(pl.lit("absent_zero"))
        .otherwise(pl.lit("observed_or_rejected"))
        .alias("assignment_source"),
    ).cast(INDEL_MUTATION_ASSIGNMENT_SCHEMA)
