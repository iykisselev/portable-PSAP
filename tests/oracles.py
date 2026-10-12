"""Independent scalar reference implementations for workflow property tests."""

from __future__ import annotations

from collections import defaultdict
import math
from pathlib import Path

import polars as pl


def distribution(
    statistics: dict[int, tuple[int, float, int]], bin_count: int
) -> tuple[str, float, int, list[float] | None]:
    """Reference scalar implementation of conditional survival calculation."""
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
    """Reference scalar reduction using manifest-ordered Kahan summation."""
    accumulated: dict[tuple[str, str, int], list[float | int]] = {}
    for path in files:
        for (
            unit,
            role,
            bin_number,
            certain_count,
            log_none_sum,
            event_count,
            snv_event_count,
            indel_event_count,
        ) in pl.read_parquet(path).iter_rows():
            key = (unit, role, bin_number)
            if key not in accumulated:
                accumulated[key] = [
                    certain_count,
                    log_none_sum,
                    0.0,
                    event_count,
                    snv_event_count,
                    indel_event_count,
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
    grouped: dict[tuple[str, str], dict[int, tuple[int, float, int, int, int]]] = (
        defaultdict(dict)
    )
    for (
        unit,
        role,
        bin_number,
    ), (
        certain_count,
        log_none_sum,
        _,
        event_count,
        snv_event_count,
        indel_event_count,
    ) in accumulated.items():
        grouped[(unit, role)][bin_number] = (
            int(certain_count),
            float(log_none_sum),
            int(event_count),
            int(snv_event_count),
            int(indel_event_count),
        )
    return grouped
