"""Class-aware event assembly for future combined SNV and indel nulls."""

from __future__ import annotations

import polars as pl

from generation_probability import role_probability_expression
from generation_scores import ALLELE_COLUMNS


EVENT_SCHEMA = {
    "testing_unit_id": pl.String,
    "CHROM": pl.String,
    "POS": pl.Int64,
    "REF": pl.String,
    "ALT": pl.String,
    "variant_class": pl.String,
    "score": pl.Float64,
    "allele_state": pl.String,
    "assignment_source": pl.String,
    "effective_af_standard": pl.Float64,
    "effective_af_hem": pl.Float64,
}

CLASS_DIAGNOSTIC_SCHEMA = {
    "testing_unit_id": pl.String,
    "role": pl.String,
    "variant_class": pl.String,
    "event_count": pl.Int64,
}


class EventAssemblyError(ValueError):
    """Raised when class-specific event artifacts cannot be joined safely."""


def _require_columns(frame: pl.DataFrame, columns: set[str], label: str) -> None:
    missing = sorted(columns - set(frame.columns))
    if missing:
        raise EventAssemblyError(f"{label} is missing columns: " + ", ".join(missing))


def _require_unique(frame: pl.DataFrame, columns: list[str], label: str) -> None:
    if frame.select(*columns).is_duplicated().any():
        raise EventAssemblyError(f"{label} contains duplicate {'/'.join(columns)} values")


def _join_required(
    events: pl.DataFrame,
    source: pl.DataFrame,
    source_columns: list[str],
    label: str,
    *,
    presence_column: str | None = None,
) -> pl.DataFrame:
    """Attach a unique allele-keyed source and fail for unmatched events."""
    joined = events.join(source.select(*ALLELE_COLUMNS, *source_columns), on=ALLELE_COLUMNS, how="left")
    missing = joined.filter(pl.col(presence_column or source_columns[0]).is_null())
    if missing.height:
        examples = ", ".join(
            f"{chrom}:{position}:{reference}:{alternate}"
            for chrom, position, reference, alternate in missing.head(10).select(*ALLELE_COLUMNS).iter_rows()
        )
        raise EventAssemblyError(
            f"{label} is missing {missing.height} indel event allele(s): {examples}"
        )
    return joined


def assemble_indel_events(
    memberships: pl.DataFrame,
    scores: pl.DataFrame,
    states: pl.DataFrame,
    assignments: pl.DataFrame,
    *,
    x_af_policy: str,
    x_mutation_rate_policy: str,
) -> pl.DataFrame:
    """Build one event per retained indel/testing-unit membership.

    Population states and indel mutation assignments are joined only by the
    canonical biallelic key.  Rejected population records are excluded before
    probability construction.  The returned table retains ``variant_class``
    so a future combined reducer can publish per-class diagnostics while using
    the same independent-event probability model for all rows.
    """
    if x_af_policy not in {"sex_specific", "overall_proxy"}:
        raise EventAssemblyError(f"Unsupported X chromosome AF policy {x_af_policy!r}")
    if x_mutation_rate_policy not in {"reject", "sex_averaged_proxy"}:
        raise EventAssemblyError(
            f"Unsupported X chromosome mutation-rate policy {x_mutation_rate_policy!r}"
        )
    _require_columns(memberships, {"testing_unit_id", "variant_class", *ALLELE_COLUMNS}, "Indel memberships")
    _require_columns(scores, {*ALLELE_COLUMNS, "score"}, "Indel scores")
    _require_columns(
        states,
        {*ALLELE_COLUMNS, "allele_state", "af_overall", "af_female", "af_male"},
        "Indel states",
    )
    _require_columns(assignments, {*ALLELE_COLUMNS, "mu", "assignment_source"}, "Indel mutation assignments")
    _require_unique(memberships, ["testing_unit_id", *ALLELE_COLUMNS], "Indel memberships")
    for frame, label in ((scores, "Indel scores"), (states, "Indel states"), (assignments, "Indel mutation assignments")):
        _require_unique(frame, ALLELE_COLUMNS, label)
    if not memberships.filter(pl.col("variant_class") != "indel").is_empty():
        raise EventAssemblyError("Indel memberships contain a non-indel variant class")

    events = memberships.select("testing_unit_id", "variant_class", *ALLELE_COLUMNS)
    events = _join_required(events, scores, ["score"], "Indel scores")
    events = _join_required(
        events,
        states,
        ["allele_state", "af_overall", "af_female", "af_male"],
        "Indel states",
    )
    events = _join_required(
        events,
        assignments,
        ["mu", "assignment_source"],
        "Indel mutation assignments",
        presence_column="assignment_source",
    )
    events = events.filter(pl.col("allele_state") != "rejected")
    if (
        x_mutation_rate_policy == "reject"
        and events.filter(
            (pl.col("CHROM") == "X")
            & (pl.col("allele_state") == "absent")
            & (pl.col("assignment_source") == "absent_mutation_rate")
        ).height
    ):
        raise EventAssemblyError(
            "Indel mutation-rate fallback encountered absent chrX alleles, but the configured policy is reject"
        )

    standard_observed = pl.when((pl.col("CHROM") == "X") & (x_af_policy == "sex_specific")).then(
        pl.col("af_female")
    ).otherwise(pl.col("af_overall"))
    hem_observed = pl.when((pl.col("CHROM") == "X") & (x_af_policy == "sex_specific")).then(
        pl.col("af_male")
    ).otherwise(pl.col("af_overall"))
    events = events.with_columns(
        pl.when(pl.col("allele_state") == "absent")
        .then(pl.col("mu"))
        .otherwise(standard_observed)
        .alias("effective_af_standard"),
        pl.when(pl.col("allele_state") == "absent")
        .then(pl.col("mu"))
        .otherwise(hem_observed)
        .alias("effective_af_hem"),
    )
    invalid_frequency = (
        pl.col("effective_af_standard").is_null()
        | ~pl.col("effective_af_standard").is_finite()
        | (pl.col("effective_af_standard") < 0)
        | (pl.col("effective_af_standard") > 1)
        | pl.col("effective_af_hem").is_null()
        | ~pl.col("effective_af_hem").is_finite()
        | (pl.col("effective_af_hem") < 0)
        | (pl.col("effective_af_hem") > 1)
    )
    if events.filter(invalid_frequency).height:
        raise EventAssemblyError("Indel events have missing or invalid effective allele frequencies")
    return events.select(*EVENT_SCHEMA).cast(EVENT_SCHEMA)


def role_event_frames(events: pl.DataFrame, roles: tuple[str, ...]) -> dict[str, pl.DataFrame]:
    """Adapt assembled events to the explicit HET/HOM/HEM probability roles."""
    _require_columns(events, set(EVENT_SCHEMA), "Assembled events")
    frames: dict[str, pl.DataFrame] = {}
    for role in roles:
        if role in {"het", "hom"}:
            frames[role] = events.with_columns(pl.col("effective_af_standard").alias("effective_af"))
        elif role == "hem_chrX":
            frames[role] = events.filter(pl.col("CHROM") == "X").with_columns(
                pl.col("effective_af_hem").alias("effective_af")
            )
        else:
            raise EventAssemblyError(f"Unsupported generated null-distribution role {role!r}")
    return frames


def event_class_diagnostics(events: pl.DataFrame, roles: tuple[str, ...]) -> pl.DataFrame:
    """Count non-zero independent events by unit, role, and variant class."""
    records: list[pl.DataFrame] = []
    for role, frame in role_event_frames(events, roles).items():
        records.append(
            frame.with_columns(role_probability_expression(role).alias("q"))
            .filter(pl.col("q") > 0)
            .group_by("testing_unit_id", "variant_class")
            .agg(pl.len().alias("event_count"))
            .with_columns(pl.lit(role).alias("role"))
            .select(*CLASS_DIAGNOSTIC_SCHEMA)
        )
    if not records:
        return pl.DataFrame(schema=CLASS_DIAGNOSTIC_SCHEMA)
    return pl.concat(records, how="vertical").cast(CLASS_DIAGNOSTIC_SCHEMA)
