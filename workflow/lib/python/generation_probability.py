"""Probability assignment, binning, and additive null-distribution reduction."""

from __future__ import annotations

from collections import defaultdict
from decimal import Decimal, ROUND_CEILING
import math
import sys
from pathlib import Path

import polars as pl


def fixed_width_bin_expression(score_max: float, bin_width: float) -> pl.Expr:
    """Return the producer-side form of the canonical decimal boundary rule."""
    if score_max <= 0 or bin_width <= 0:
        raise ValueError("score_max and bin_width must be positive")
    bin_count = int((Decimal(str(score_max)) / Decimal(str(bin_width))).to_integral_value(rounding=ROUND_CEILING))
    score = pl.col("score").clip(0, score_max)
    scaled = score / bin_width
    nearest = scaled.round(0)
    tolerance = sys.float_info.epsilon * 64 * pl.max_horizontal(pl.lit(1.0), scaled.abs())
    corrected = pl.when((scaled - nearest).abs() <= tolerance).then(nearest).otherwise(scaled)
    return corrected.floor().cast(pl.Int64).add(1).clip(1, bin_count).alias("bin")


def role_probability_expression(role: str, af_column: str = "effective_af") -> pl.Expr:
    """Return the explicit event-probability adapter for one generated role."""
    frequency = pl.col(af_column)
    if role == "het":
        return 2 * frequency * (1 - frequency)
    if role == "hom":
        return frequency**2
    if role == "hem_chrX":
        return frequency
    raise ValueError(f"Unsupported generated null-distribution role {role!r}")


def unit_bucket_frame(
    units: pl.DataFrame,
    bucket_size: int,
    unit_intervals: pl.DataFrame,
) -> pl.DataFrame:
    """Assign eligible units to deterministic, genomically compact buckets."""
    if bucket_size <= 0:
        raise ValueError("generation.execution.unit_bucket_size must be positive")
    eligible = units.filter(pl.col("status") == "eligible").select("testing_unit_id")
    first_locations = unit_intervals.group_by("testing_unit_id").agg(
        pl.col("CHROM").first().alias("_chromosome"),
        pl.col("START").min().alias("_first_start"),
    )
    ordered = eligible.join(first_locations, on="testing_unit_id", how="inner")
    if ordered.height != eligible.height:
        raise ValueError("Eligible testing units require genomic intervals")
    # Scope preparation accepts only canonical autosomes and chrX, so numeric
    # contig order followed by X is stable across platforms and input BED order.
    return ordered.with_columns(
        pl.col("_chromosome").cast(pl.Int64, strict=False).fill_null(23).alias("_chromosome_order")
    ).sort("_chromosome_order", "_first_start", "testing_unit_id").select("testing_unit_id").with_row_index("unit_ordinal").with_columns(
        (pl.col("unit_ordinal") // bucket_size).cast(pl.Int64).alias("unit_bucket")
    ).drop("unit_ordinal")


def statistics_for_role(events: pl.DataFrame, role: str, score_max: float, bin_width: float) -> pl.DataFrame:
    """Emit additive no-event statistics once per unit, role, and score bin."""
    probability = role_probability_expression(role)
    return (
        events.with_columns(probability.alias("q"), fixed_width_bin_expression(score_max, bin_width))
        .filter(pl.col("q") > 0)
        .group_by("testing_unit_id", "unit_bucket", "bin")
        .agg(
            (pl.col("q") == 1).cast(pl.Int64).sum().alias("certain_count"),
            pl.when(pl.col("q") < 1).then((-pl.col("q")).log1p()).otherwise(pl.lit(0.0)).sum().alias("log_none_sum"),
            pl.len().alias("event_count"),
            (pl.col("variant_class") == "snv").cast(pl.Int64).sum().alias("snv_event_count"),
            (pl.col("variant_class") == "indel").cast(pl.Int64).sum().alias("indel_event_count"),
        )
        .with_columns(pl.lit(role).alias("role"))
        .select(
            "testing_unit_id", "unit_bucket", "role", "bin", "certain_count",
            "log_none_sum", "event_count", "snv_event_count", "indel_event_count",
        )
    )


def distribution(statistics: dict[int, tuple[int, float, int]], bin_count: int) -> tuple[str, float, int, list[float] | None]:
    """Calculate a conditional survival row from additive per-bin statistics."""
    certain = [0] * bin_count
    log_none = [0.0] * bin_count
    event_count = 0
    for bin_number, (certain_count, log_none_sum, count) in statistics.items():
        certain[bin_number - 1] = certain_count
        log_none[bin_number - 1] = log_none_sum
        event_count += count
    if event_count == 0:
        return "zero_event_mass", 0.0, 0, None
    total_log_none = math.fsum(log_none)
    event_mass = 1.0 if sum(certain) else -math.expm1(total_log_none)
    if not math.isfinite(event_mass) or event_mass <= 0 or event_mass > 1:
        return "numerical_failure", event_mass, event_count, None
    survival = [0.0] * bin_count
    cumulative_certain = 0
    cumulative_log_none = 0.0
    for index in range(bin_count - 1, -1, -1):
        cumulative_certain += certain[index]
        cumulative_log_none += log_none[index]
        tail_mass = 1.0 if cumulative_certain else -math.expm1(cumulative_log_none)
        survival[index] = min(1.0, max(0.0, tail_mass / event_mass))
    return "ok", event_mass, event_count, survival


def merge_statistics(
    files: list[Path],
) -> dict[tuple[str, str], dict[int, tuple[int, float, int, int, int]]]:
    """Reduce partition files in manifest order with compensated log summation."""
    accumulated: dict[tuple[str, str, int], list[float | int]] = {}
    for path in files:
        for (
            unit, role, bin_number, certain_count, log_none_sum, event_count,
            snv_event_count, indel_event_count,
        ) in pl.read_parquet(path).iter_rows():
            key = (unit, role, bin_number)
            if key not in accumulated:
                accumulated[key] = [
                    certain_count, log_none_sum, 0.0, event_count,
                    snv_event_count, indel_event_count,
                ]
                continue
            current = accumulated[key]
            adjusted = log_none_sum - current[2]
            total = current[1] + adjusted
            current[2] = (total - current[1]) - adjusted
            current[1] = total
            current[0] += certain_count
            current[3] += event_count
            current[4] += snv_event_count
            current[5] += indel_event_count
    grouped: dict[
        tuple[str, str], dict[int, tuple[int, float, int, int, int]]
    ] = defaultdict(dict)
    for (
        unit, role, bin_number,
    ), (
        certain_count, log_none_sum, _, event_count, snv_event_count,
        indel_event_count,
    ) in accumulated.items():
        grouped[(unit, role)][bin_number] = (
            int(certain_count), float(log_none_sum), int(event_count),
            int(snv_event_count), int(indel_event_count),
        )
    return grouped
