"""Pure preparation helpers for the null-distribution generator.

These functions deliberately do not use cohort resources.  They retain enough
state to distinguish an absent panel allele from a malformed, filtered, or
otherwise rejected observed allele before probability assignment.
"""

from __future__ import annotations

from bisect import bisect_right
from collections import defaultdict
from pathlib import Path
import gzip

from genomic_contracts import canonical_chromosome


def _open(path: str | Path):
    file_path = Path(path)
    return gzip.open(file_path, "rt", encoding="utf-8") if file_path.suffix == ".gz" else file_path.open(encoding="utf-8")


def read_bed_units(path: str | Path):
    units: dict[str, list[tuple[str, int, int]]] = defaultdict(list)
    unit_chromosomes: dict[str, str] = {}
    with _open(path) as handle:
        for number, raw in enumerate(handle, start=1):
            if not raw.strip() or raw.startswith(("#", "track", "browser")):
                continue
            fields = raw.rstrip().split("\t")
            if len(fields) < 4:
                raise ValueError(f"Testing-unit BED line {number} needs four columns")
            chrom, start_text, end_text, unit = fields[:4]
            start, end = int(start_text), int(end_text)
            if start < 0 or end <= start or not unit:
                raise ValueError(f"Invalid testing-unit BED line {number}")
            canonical = canonical_chromosome(chrom)
            previous = unit_chromosomes.setdefault(unit, canonical)
            if previous != canonical:
                raise ValueError(f"Testing unit {unit!r} spans chromosomes")
            units[unit].append((canonical, start, end))
    return units


def read_intervals(path: str | Path):
    intervals: dict[str, list[tuple[int, int]]] = defaultdict(list)
    with _open(path) as handle:
        for number, raw in enumerate(handle, start=1):
            if not raw.strip() or raw.startswith(("#", "track", "browser")):
                continue
            fields = raw.rstrip().split("\t")
            if len(fields) < 3:
                raise ValueError(f"Interval BED line {number} needs three columns")
            start, end = int(fields[1]), int(fields[2])
            if start < 0 or end <= start:
                raise ValueError(f"Invalid interval BED line {number}")
            intervals[canonical_chromosome(fields[0])].append((start, end))
    return intervals


def merge_intervals(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
    """Return sorted, non-overlapping half-open intervals."""
    merged: list[tuple[int, int]] = []
    for start, end in sorted(intervals):
        if not merged or start > merged[-1][1]:
            merged.append((start, end))
        else:
            merged[-1] = (merged[-1][0], max(merged[-1][1], end))
    return merged


def intersect_intervals(
    left: list[tuple[int, int]], right: list[tuple[int, int]]
) -> list[tuple[int, int]]:
    """Intersect two sets of half-open intervals."""
    left_merged, right_merged = merge_intervals(left), merge_intervals(right)
    result: list[tuple[int, int]] = []
    left_index = right_index = 0
    while left_index < len(left_merged) and right_index < len(right_merged):
        left_start, left_end = left_merged[left_index]
        right_start, right_end = right_merged[right_index]
        start, end = max(left_start, right_start), min(left_end, right_end)
        if start < end:
            result.append((start, end))
        if left_end <= right_end:
            left_index += 1
        else:
            right_index += 1
    return result


def subtract_intervals(
    intervals: list[tuple[int, int]], excluded: list[tuple[int, int]]
) -> list[tuple[int, int]]:
    """Subtract exclusions in a linear pass after sorting and merging both sets."""
    exclusions = merge_intervals(excluded)
    remaining: list[tuple[int, int]] = []
    exclusion_index = 0
    for start, end in merge_intervals(intervals):
        while exclusion_index < len(exclusions) and exclusions[exclusion_index][1] <= start:
            exclusion_index += 1
        position = start
        while exclusion_index < len(exclusions) and exclusions[exclusion_index][0] < end:
            excluded_start, excluded_end = exclusions[exclusion_index]
            if excluded_start > position:
                remaining.append((position, excluded_start))
            position = max(position, excluded_end)
            if position >= end:
                # This exclusion may also overlap the next input interval.
                break
            exclusion_index += 1
        if position < end:
            remaining.append((position, end))
    return remaining


def supported_scope_mask(
    contig_lengths: dict[str, int],
    excluded: dict[str, list[tuple[int, int]]],
) -> dict[str, list[tuple[int, int]]]:
    """Build supported contig intervals after all static interval exclusions."""
    return {
        chromosome: subtract_intervals([(0, length)], excluded.get(chromosome, []))
        for chromosome, length in sorted(contig_lengths.items())
        if length > 0
    }


def effective_calibration_mask(
    exclusions: dict[str, list[tuple[int, int]]],
    supported_scope: dict[str, list[tuple[int, int]]],
) -> dict[str, list[tuple[int, int]]]:
    """Subtract interval exclusions from supported, non-PAR sequence."""
    return {
        chromosome: subtract_intervals(
            supported_intervals, exclusions.get(chromosome, [])
        )
        for chromosome, supported_intervals in sorted(supported_scope.items())
    }


def merged_length(intervals: list[tuple[int, int]]) -> int:
    return sum(end - start for start, end in merge_intervals(intervals))


def overlap_sorted(
    intervals: list[tuple[int, int]],
    merged_regions: list[tuple[int, int]],
    region_starts: list[int] | None = None,
) -> list[tuple[int, int]]:
    """Return portions of ``intervals`` covered by sorted, merged regions.

    A unit has only a few intervals while a chromosome-scale mask can have many
    thousands.  Binary-searching the mask prevents repeatedly scanning its
    earlier regions for every unit.
    """
    if not merged_regions:
        return []
    if region_starts is None:
        region_starts = [start for start, _ in merged_regions]
    pieces: list[tuple[int, int]] = []
    for start, end in merge_intervals(intervals):
        # The preceding region is the only one that can start before ``start``
        # and still overlap it because regions are non-overlapping.
        index = max(0, bisect_right(region_starts, start) - 1)
        while index < len(merged_regions) and merged_regions[index][0] < end:
            region_start, region_end = merged_regions[index]
            left, right = max(start, region_start), min(end, region_end)
            if left < right:
                pieces.append((left, right))
            index += 1
    return pieces


def overlap_length(intervals: list[tuple[int, int]], regions: list[tuple[int, int]]) -> int:
    """Return covered bases for arbitrary interval collections."""
    merged_regions = merge_intervals(regions)
    return sum(
        end - start
        for start, end in overlap_sorted(
            intervals,
            merged_regions,
            [start for start, _ in merged_regions],
        )
    )


def retention_rows(
    units,
    retained_regions,
    minimum: float,
    denominator_scope=None,
    *,
    retained_region_starts=None,
    denominator_region_starts=None,
    callable_intervals_by_unit=None,
):
    """Calculate per-unit retention and optionally retain callable pieces.

    Callers that already hold chromosome masks can supply their start offsets
    once and receive the retained pieces used downstream, avoiding a second
    mask intersection for every unit.
    """
    if retained_region_starts is None:
        retained_region_starts = {
            chrom: [start for start, _ in regions]
            for chrom, regions in retained_regions.items()
        }
    if denominator_scope is not None and denominator_region_starts is None:
        denominator_region_starts = {
            chrom: [start for start, _ in regions]
            for chrom, regions in denominator_scope.items()
        }
    rows = []
    for unit, intervals in sorted(units.items()):
        chrom = intervals[0][0]
        simple_intervals = [(start, end) for _, start, end in intervals]
        denominator_pieces = (
            overlap_sorted(
                simple_intervals,
                denominator_scope.get(chrom, []),
                denominator_region_starts.get(chrom, []),
            )
            if denominator_scope is not None
            else merge_intervals(simple_intervals)
        )
        retained_pieces = overlap_sorted(
            simple_intervals,
            retained_regions.get(chrom, []),
            retained_region_starts.get(chrom, []),
        )
        if callable_intervals_by_unit is not None:
            callable_intervals_by_unit[unit] = retained_pieces
        total = sum(end - start for start, end in denominator_pieces)
        retained = sum(end - start for start, end in retained_pieces)
        fraction = retained / total if total else 0.0
        if total == 0:
            status = "no_calibration_scope"
        elif retained == 0:
            status = "fully_excluded"
        else:
            status = "eligible" if fraction >= minimum else "retained_fraction_below_minimum"
        rows.append((unit, chrom, total, retained, fraction, status))
    return rows
