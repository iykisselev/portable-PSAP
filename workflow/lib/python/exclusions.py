#!/usr/bin/env python3
"""Project mixed exclusions onto full cohort-VCF record spans once."""

from __future__ import annotations

import bisect
import gzip
from collections import defaultdict
from collections.abc import Iterable
from pathlib import Path
from typing import TextIO

from genomic_contracts import canonical_chromosome


def open_text(path: Path):
    with path.open("rb") as handle:
        compressed = handle.read(2) == b"\x1f\x8b"
    if compressed:
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open("rt", encoding="utf-8", errors="replace")


def data_rows(path: Path):
    with open_text(path) as handle:
        for line_number, line in enumerate(handle, start=1):
            stripped = line.strip()
            if (
                not stripped
                or stripped.startswith("#")
                or stripped.startswith("track ")
                or stripped.startswith("browser ")
            ):
                continue
            yield line_number, stripped.split("\t")


def add_bed_interval(
    path: Path,
    line_number: int,
    columns: list[str],
    intervals: dict[str, list[tuple[int, int]]],
) -> None:
    if len(columns) < 3:
        raise ValueError(f"{path}:{line_number}: BED row has fewer than 3 columns")
    start = int(columns[1])
    end = int(columns[2])
    if start < 0 or end <= start:
        raise ValueError(f"{path}:{line_number}: invalid BED interval {start}, {end}")
    intervals[canonical_chromosome(columns[0])].append((start + 1, end))


def add_tsv_target(
    path: Path,
    line_number: int,
    columns: list[str],
    intervals: dict[str, list[tuple[int, int]]],
    variants: set[tuple[str, int, str, str]],
) -> None:
    if len(columns) < 3:
        raise ValueError(f"{path}:{line_number}: TSV row has fewer than 3 columns")
    try:
        start = int(columns[1])
    except ValueError as exc:
        raise ValueError(f"{path}:{line_number}: TSV position is not an integer") from exc

    try:
        end = int(columns[2])
    except ValueError:
        if len(columns) < 4 or start < 1 or not columns[2] or not columns[3]:
            raise ValueError(
                f"{path}:{line_number}: expected CHROM/POS/REF/ALT variant columns"
            )
        for alternate in columns[3].split(","):
            variants.add((canonical_chromosome(columns[0]), start, columns[2], alternate))
        return

    if start < 1 or end < start:
        raise ValueError(
            f"{path}:{line_number}: TSV intervals are 1-based inclusive and must "
            f"satisfy 1 <= start <= end; received {start}, {end}"
        )
    intervals[canonical_chromosome(columns[0])].append((start, end))


def add_vcf_variants(
    path: Path,
    line_number: int,
    columns: list[str],
    variants: set[tuple[str, int, str, str]],
) -> None:
    if len(columns) < 5:
        raise ValueError(f"{path}:{line_number}: VCF row has fewer than 5 columns")
    position = int(columns[1])
    if position < 1:
        raise ValueError(f"{path}:{line_number}: VCF position must be positive")
    for alternate in columns[4].split(","):
        variants.add(
            (canonical_chromosome(columns[0]), position, columns[3], alternate)
        )


def merged_interval_index(
    intervals: dict[str, list[tuple[int, int]]],
) -> dict[str, tuple[list[int], list[int]]]:
    result = {}
    for chromosome, chromosome_intervals in intervals.items():
        merged: list[list[int]] = []
        for start, end in sorted(chromosome_intervals):
            if merged and start <= merged[-1][1] + 1:
                merged[-1][1] = max(merged[-1][1], end)
            else:
                merged.append([start, end])
        result[chromosome] = (
            [interval[0] for interval in merged],
            [interval[1] for interval in merged],
        )
    return result


def interval_overlaps_variant(
    index: dict[str, tuple[list[int], list[int]]],
    chromosome: str,
    variant_start: int,
    variant_end: int,
) -> bool:
    if variant_start < 1 or variant_end < variant_start:
        raise ValueError(
            f"Invalid variant span: {chromosome}:{variant_start}-{variant_end}"
        )
    chromosome_index = index.get(canonical_chromosome(chromosome))
    if chromosome_index is None:
        return False
    starts, ends = chromosome_index
    candidate = bisect.bisect_right(starts, variant_end) - 1
    return candidate >= 0 and ends[candidate] >= variant_start


def vcf_record_span(columns: list[str]) -> tuple[int, int]:
    """Return the full 1-based inclusive span represented by a VCF record."""
    position = int(columns[1])
    if position < 1:
        raise ValueError("VCF position must be positive")

    reference = columns[3]
    if not reference or reference == ".":
        raise ValueError("VCF REF allele must be non-empty")
    end = position + len(reference) - 1

    if len(columns) >= 8 and columns[7] not in {"", "."}:
        for field in columns[7].split(";"):
            if not field.startswith("END="):
                continue
            try:
                info_end = int(field.removeprefix("END="))
            except ValueError as exc:
                raise ValueError(f"VCF INFO/END is not an integer: {field!r}") from exc
            if info_end < position:
                raise ValueError(
                    f"VCF INFO/END precedes POS: {position}, {info_end}"
                )
            end = max(end, info_end)
            break

    return position, end


def collect_exclusions(specifications: list[str]):
    intervals: dict[str, list[tuple[int, int]]] = defaultdict(list)
    variants: set[tuple[str, int, str, str]] = set()
    for specification in specifications:
        kind, separator, path_text = specification.partition(":")
        if not separator or kind not in {"bed", "tsv", "vcf"}:
            raise ValueError(f"Invalid --source value: {specification!r}")
        path = Path(path_text)
        for line_number, columns in data_rows(path):
            if kind == "bed":
                add_bed_interval(path, line_number, columns, intervals)
            elif kind == "tsv":
                add_tsv_target(path, line_number, columns, intervals, variants)
            else:
                add_vcf_variants(path, line_number, columns, variants)
    return merged_interval_index(intervals), variants


def project_vcf_records(
    cohort: Iterable[str],
    output: TextIO,
    intervals: dict[str, tuple[list[int], list[int]]],
    variants: set[tuple[str, int, str, str]],
    source_name: str,
) -> None:
    """Write matching records from an already-filtered, sites-only VCF stream."""
    for line_number, line in enumerate(cohort, start=1):
        if line.startswith("#"):
            output.write(line)
            continue
        columns = line.rstrip("\n").split("\t")
        if len(columns) < 5:
            raise ValueError(
                f"{source_name}:{line_number}: cohort VCF row has fewer than 5 columns"
            )
        try:
            position, record_end = vcf_record_span(columns)
        except ValueError as exc:
            raise ValueError(f"{source_name}:{line_number}: {exc}") from exc
        interval_match = interval_overlaps_variant(
            intervals, columns[0], position, record_end
        )
        chromosome = canonical_chromosome(columns[0])
        variant_match = any(
            (chromosome, position, columns[3], alternate) in variants
            for alternate in columns[4].split(",")
        )
        if interval_match or variant_match:
            output.write(line)


def write_projected_vcf(
    cohort_path: Path,
    output_path: Path,
    intervals: dict[str, tuple[list[int], list[int]]],
    variants: set[tuple[str, int, str, str]],
) -> None:
    """Project a VCF path onto exclusion targets for direct library callers."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open_text(cohort_path) as cohort, output_path.open("w", encoding="utf-8") as output:
        project_vcf_records(cohort, output, intervals, variants, str(cohort_path))
