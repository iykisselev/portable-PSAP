"""Canonical effective-frequency and role assembly for generation events."""

from __future__ import annotations

import polars as pl
from bundle_contract import GENERATED_ROLES


class EventAssemblyError(ValueError):
    """Raised when generation event probabilities cannot be assembled safely."""


def _require_columns(frame: pl.DataFrame, columns: set[str], label: str) -> None:
    missing = sorted(columns - set(frame.columns))
    if missing:
        raise EventAssemblyError(f"{label} is missing columns: " + ", ".join(missing))


def assign_effective_frequencies(
    events: pl.DataFrame,
    *,
    unobserved_policy: str,
    x_af_policy: str,
    x_mutation_rate_policy: str,
    roles: tuple[str, ...],
) -> pl.DataFrame:
    """Assign allele frequencies for retained SNV and indel membership events.

    ``mu`` holds SNV fallback rates and ``mu_indel`` holds optional indel
    fallback rates.  Missing indel rates deliberately resolve to zero, which
    preserves the finite score-catalog universe when no indel rate resource is
    configured.  Missing SNV rates remain a hard error.
    """
    if not roles or set(roles) - set(GENERATED_ROLES):
        raise EventAssemblyError(f"Unsupported generated roles: {roles!r}")
    if x_af_policy not in {"sex_specific", "overall_proxy"}:
        raise EventAssemblyError(f"Unsupported X chromosome AF policy {x_af_policy!r}")
    if x_mutation_rate_policy not in {"reject", "sex_averaged_proxy"}:
        raise EventAssemblyError(
            f"Unsupported X chromosome mutation-rate policy {x_mutation_rate_policy!r}"
        )
    _require_columns(
        events,
        {
            "CHROM",
            "variant_class",
            "allele_state",
            "af_overall",
            "af_female",
            "af_male",
        },
        "Generation events",
    )
    events = events.filter(pl.col("allele_state") != "rejected")
    if unobserved_policy not in {"zero", "mutation_rate"}:
        raise EventAssemblyError(f"Unsupported unobserved policy {unobserved_policy!r}")
    standard_observed = (
        pl.when((pl.col("CHROM") == "X") & (x_af_policy == "sex_specific"))
        .then(pl.col("af_female"))
        .otherwise(pl.col("af_overall"))
    )
    hem_observed = (
        pl.when((pl.col("CHROM") == "X") & (x_af_policy == "sex_specific"))
        .then(pl.col("af_male"))
        .otherwise(pl.col("af_overall"))
    )
    if unobserved_policy == "zero":
        absent_frequency = pl.lit(0.0)
        assignment_source = pl.lit("absent_zero")
    else:
        _require_columns(events, {"mu", "mu_indel"}, "Mutation-rate events")
        needs_snv_rate = (pl.col("allele_state") == "absent") & (
            pl.col("variant_class") != "indel"
        )
        if events.filter(needs_snv_rate & pl.col("mu").is_null()).height:
            raise EventAssemblyError(
                "Mutation-rate resource is missing absent SNV alleles"
            )
        if (
            x_mutation_rate_policy == "reject"
            and events.filter(
                (pl.col("CHROM") == "X") & (pl.col("allele_state") == "absent")
            ).height
        ):
            raise EventAssemblyError(
                "mutation-rate fallback encountered absent chrX alleles, but the configured policy is reject"
            )
        absent_frequency = (
            pl.when(pl.col("variant_class") == "indel")
            .then(pl.col("mu_indel").fill_null(0.0))
            .otherwise(pl.col("mu"))
        )
        assignment_source = (
            pl.when(pl.col("variant_class") == "indel")
            .then(pl.lit("absent_indel_mutation_rate"))
            .otherwise(pl.lit("absent_mutation_rate"))
        )
    events = events.with_columns(
        pl.when(pl.col("allele_state") == "absent")
        .then(absent_frequency)
        .otherwise(standard_observed)
        .alias("effective_af_standard"),
        pl.when(pl.col("allele_state") == "absent")
        .then(absent_frequency)
        .otherwise(hem_observed)
        .alias("effective_af_hem"),
        pl.when(pl.col("allele_state") == "absent")
        .then(assignment_source)
        .otherwise(pl.lit("observed_af"))
        .alias("assignment_source"),
    )
    invalid_frequency = pl.lit(False)
    for column, applicable in (
        ("effective_af_standard", pl.lit(bool(set(roles) & {"het", "hom"}))),
        ("effective_af_hem", (pl.col("CHROM") == "X") & ("hem_chrX" in roles)),
    ):
        frequency = pl.col(column)
        invalid_frequency |= applicable & (
            frequency.is_null()
            | ~frequency.is_finite()
            | (frequency < 0)
            | (frequency > 1)
        )
    if events.filter(invalid_frequency).height:
        raise EventAssemblyError(
            "Generation events have missing or invalid effective allele frequencies"
        )
    return events


def role_event_frames(
    events: pl.DataFrame, roles: tuple[str, ...]
) -> dict[str, pl.DataFrame]:
    """Adapt assembled events to the explicit HET/HOM/HEM probability roles."""
    _require_columns(
        events,
        {"CHROM", "effective_af_standard", "effective_af_hem"},
        "Assembled events",
    )
    frames: dict[str, pl.DataFrame] = {}
    for role in roles:
        if role in {"het", "hom"}:
            frames[role] = events.with_columns(
                pl.col("effective_af_standard").alias("effective_af")
            )
        elif role == "hem_chrX":
            frames[role] = events.filter(pl.col("CHROM") == "X").with_columns(
                pl.col("effective_af_hem").alias("effective_af")
            )
        else:
            raise EventAssemblyError(
                f"Unsupported generated null-distribution role {role!r}"
            )
    return frames
