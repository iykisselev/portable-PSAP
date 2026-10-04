"""Pure preparation helpers for the null-distribution generator.

These functions deliberately do not use cohort resources.  They retain enough
state to distinguish an absent panel allele from a malformed, filtered, or
otherwise rejected observed allele before probability assignment.
"""

from __future__ import annotations

from collections import defaultdict
from pathlib import Path
import gzip
import math

try:
    from genomic_contracts import allele_key, allele_key_text as key_text, canonical_chromosome
except ModuleNotFoundError:
    from .genomic_contracts import allele_key, allele_key_text as key_text, canonical_chromosome


def _open(path: str | Path):
    file_path = Path(path)
    return gzip.open(file_path, "rt", encoding="utf-8") if file_path.suffix == ".gz" else file_path.open(encoding="utf-8")


def read_population_af(path: str | Path, info_tag: str):
    """Read valid observed AF and exact rejected allele keys from a VCF.

    ``Number=1`` is accepted only for records originally biallelic; splitting
    a multiallelic record must never manufacture per-ALT scalar frequencies.
    """
    number: str | None = None
    observed: dict[tuple[str, int, str, str], float] = {}
    rejected: dict[tuple[str, int, str, str], str] = {}
    with _open(path) as handle:
        for raw in handle:
            if raw.startswith(f"##INFO=<ID={info_tag},"):
                fields = raw.rstrip().split(",")
                for field in fields:
                    if field.startswith("Number="):
                        number = field.split("=", 1)[1].rstrip(">")
            if raw.startswith("#"):
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) < 8:
                raise ValueError("Population VCF record has fewer than eight columns")
            chrom, pos_text, _, ref, alt_text, _, filter_value, info = fields[:8]
            alts = alt_text.split(",")
            keys = [allele_key(chrom, int(pos_text), ref, alt) for alt in alts]
            if filter_value not in {"PASS", "."}:
                for key in keys:
                    rejected[key] = "non_pass_filter"
                continue
            info_values = dict(
                item.split("=", 1) for item in info.split(";") if "=" in item
            )
            raw_value = info_values.get(info_tag)
            if number not in {"A", "R", "1"}:
                raise ValueError(f"INFO tag {info_tag!r} must declare Number=A, R, or 1")
            if raw_value in {None, "", "."}:
                reason = "missing_configured_af"
                for key in keys:
                    rejected[key] = reason
                continue
            values = raw_value.split(",")
            if number == "A":
                selected = values if len(values) == len(alts) else None
            elif number == "R":
                selected = values[1:] if len(values) == len(alts) + 1 else None
            else:
                selected = values if len(alts) == 1 and len(values) == 1 else None
            if selected is None:
                reason = "scalar_af_on_multiallelic_record" if number == "1" and len(alts) > 1 else "invalid_af_cardinality"
                for key in keys:
                    rejected[key] = reason
                continue
            for key, value in zip(keys, selected):
                try:
                    af = float(value)
                except ValueError:
                    af = math.nan
                if not math.isfinite(af) or not 0 <= af <= 1:
                    rejected[key] = "invalid_af"
                elif key in observed or key in rejected:
                    raise ValueError(f"Duplicate/conflicting population allele: {key_text(key)}")
                else:
                    observed[key] = af
    if number is None:
        raise ValueError(f"Population VCF has no INFO header for configured tag {info_tag!r}")
    return observed, rejected


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


def overlap_length(intervals: list[tuple[int, int]], regions: list[tuple[int, int]]) -> int:
    pieces: list[tuple[int, int]] = []
    for start, end in intervals:
        for region_start, region_end in regions:
            left, right = max(start, region_start), min(end, region_end)
            if left < right:
                pieces.append((left, right))
    return merged_length(pieces)


def retention_rows(units, retained_regions, minimum: float, denominator_scope=None):
    rows = []
    for unit, intervals in sorted(units.items()):
        chrom = intervals[0][0]
        simple_intervals = [(start, end) for _, start, end in intervals]
        total = (
            overlap_length(simple_intervals, denominator_scope.get(chrom, []))
            if denominator_scope is not None
            else merged_length(simple_intervals)
        )
        retained = overlap_length(simple_intervals, retained_regions.get(chrom, []))
        fraction = retained / total if total else 0.0
        if total == 0:
            status = "no_calibration_scope"
        elif retained == 0:
            status = "fully_excluded"
        else:
            status = "eligible" if fraction >= minimum else "retained_fraction_below_minimum"
        rows.append((unit, chrom, total, retained, fraction, status))
    return rows
