"""Pair-level biological partition ownership and eligibility helpers."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping

from generation_prep import merge_intervals, overlap_sorted, resolve_partition_ownership


@dataclass(frozen=True)
class PartitionScope:
    """Resolved ownership intervals and one auditable row per unit/partition."""

    resolved: Mapping[str, Mapping[str, list[tuple[str, int, int]]]]
    rows: tuple[dict[str, object], ...]


def _bases(intervals: list[tuple[int, int]]) -> int:
    return sum(end - start for start, end in intervals)


def build_partition_scope(
    partitions: list[tuple[str, dict[str, list[tuple[str, int, int]]]]],
    supported_scope: Mapping[str, list[tuple[int, int]]],
    calibration_mask: Mapping[str, list[tuple[int, int]]],
    minimum_retained_fraction: float,
) -> PartitionScope:
    """Resolve precedence and calculate scope eligibility for every declared pair."""
    if not 0 <= minimum_retained_fraction <= 1:
        raise ValueError("minimum_retained_fraction must be in [0, 1]")
    resolved = resolve_partition_ownership(partitions)
    rows: list[dict[str, object]] = []
    for partition_id, units in partitions:
        for unit_id, raw_intervals in sorted(units.items()):
            chromosome = raw_intervals[0][0] if raw_intervals else ""
            raw = merge_intervals([(start, end) for _, start, end in raw_intervals])
            owned = resolved[partition_id][unit_id]
            owned_simple = [(start, end) for _, start, end in owned]
            supported = overlap_sorted(
                owned_simple, supported_scope.get(chromosome, [])
            )
            retained = overlap_sorted(
                owned_simple, calibration_mask.get(chromosome, [])
            )
            raw_bases, resolved_bases = _bases(raw), _bases(owned_simple)
            supported_bases, retained_bases = _bases(supported), _bases(retained)
            fraction = retained_bases / supported_bases if supported_bases else 0.0
            if resolved_bases == 0:
                status = "empty_after_priority"
            elif supported_bases == 0:
                status = "no_calibration_scope"
            elif retained_bases == 0:
                status = "fully_excluded"
            elif fraction < minimum_retained_fraction:
                status = "retained_fraction_below_minimum"
            else:
                status = "eligible"
            rows.append(
                {
                    "testing_unit_id": unit_id,
                    "partition_id": partition_id,
                    "CHROM": chromosome,
                    "raw_bases": raw_bases,
                    "priority_removed_bases": raw_bases - resolved_bases,
                    "resolved_bases": resolved_bases,
                    "supported_bases": supported_bases,
                    "retained_bases": retained_bases,
                    "retained_fraction": fraction,
                    "status": status,
                }
            )
    return PartitionScope(resolved=resolved, rows=tuple(rows))
