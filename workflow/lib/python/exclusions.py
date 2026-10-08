#!/usr/bin/env python3
"""Project mixed exclusions onto full cohort-VCF record spans once."""

from __future__ import annotations

import bisect
from collections import defaultdict
from collections.abc import Iterable
from contextlib import ExitStack
from dataclasses import dataclass
from pathlib import Path
from typing import TextIO

from genomic_contracts import canonical_chromosome
from io_contracts import iter_bed, open_text

try:
    import pysam
except ImportError:  # pragma: no cover - production validation environment supplies pysam
    pysam = None


@dataclass(frozen=True)
class ExclusionInterval:
    """One 1-based inclusive exclusion interval."""

    chromosome: str
    start: int
    end: int


@dataclass(frozen=True)
class ExclusionAlleles:
    """One 1-based exclusion anchor and its alternate alleles."""

    chromosome: str
    position: int
    reference: str
    alternates: tuple[str, ...]


def data_rows(path: Path):
    with open_text(path, errors="replace") as handle:
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


def parse_exclusion_tsv_row(
    columns: list[str],
) -> ExclusionInterval | ExclusionAlleles:
    """Parse one interval or allele exclusion TSV row into canonical targets."""
    if len(columns) < 3:
        raise ValueError("must have at least 3 columns")
    chromosome = canonical_chromosome(columns[0])
    if not chromosome:
        raise ValueError("has an empty chromosome")
    try:
        start = int(columns[1])
    except ValueError as exc:
        raise ValueError("has a non-integer position") from exc

    try:
        end = int(columns[2])
    except ValueError:
        if len(columns) < 4 or start < 1 or not columns[2] or not columns[3]:
            raise ValueError("must contain either CHROM,START,END or CHROM,POS,REF,ALT")
        alternates = tuple(columns[3].split(","))
        if any(not alternate for alternate in alternates):
            raise ValueError("has an empty ALT allele")
        return ExclusionAlleles(chromosome, start, columns[2], alternates)

    if start < 1 or end < start:
        raise ValueError(
            "interval must satisfy 1 <= start <= end; "
            f"received {start}, {end}"
        )
    return ExclusionInterval(chromosome, start, end)


def add_tsv_target(
    path: Path,
    line_number: int,
    columns: list[str],
    intervals: dict[str, list[tuple[int, int]]],
    variants: set[tuple[str, int, str, str]],
) -> None:
    """Add one shared-parsed TSV target to staging collections."""
    try:
        target = parse_exclusion_tsv_row(columns)
    except ValueError as exc:
        raise ValueError(f"{path}:{line_number}: TSV row {exc}") from exc
    if isinstance(target, ExclusionInterval):
        intervals[target.chromosome].append((target.start, target.end))
    else:
        variants.update(
            (target.chromosome, target.position, target.reference, alternate)
            for alternate in target.alternates
        )


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


def collect_exclusions(
    specifications: list[str],
    *,
    indexed_tsv_paths: Iterable[str] = (),
):
    """Collect ordinary targets and retain selected large TSVs for tabix lookup."""
    intervals: dict[str, list[tuple[int, int]]] = defaultdict(list)
    variants: set[tuple[str, int, str, str]] = set()
    indexed_paths = {str(Path(path)) for path in indexed_tsv_paths}
    indexed_tsv_sources: list[Path] = []
    for specification in specifications:
        kind, separator, path_text = specification.partition(":")
        if not separator or kind not in {"bed", "tsv", "vcf"}:
            raise ValueError(f"Invalid --source value: {specification!r}")
        path = Path(path_text)
        if kind == "bed":
            for row in iter_bed(path, errors="replace"):
                intervals[canonical_chromosome(row.chromosome)].append(
                    (row.start + 1, row.end)
                )
            continue
        if kind == "tsv" and str(path) in indexed_paths:
            index = Path(f"{path}.tbi")
            if not index.is_file():
                raise ValueError(f"Indexed TSV exclusion has no tabix index: {index}")
            indexed_tsv_sources.append(path)
            continue
        for line_number, columns in data_rows(path):
            if kind == "tsv":
                add_tsv_target(path, line_number, columns, intervals, variants)
            else:
                add_vcf_variants(path, line_number, columns, variants)
    return merged_interval_index(intervals), variants, tuple(indexed_tsv_sources)


def indexed_tsv_matches_variant(
    readers,
    chromosome: str,
    position: int,
    record_end: int,
    reference: str,
    alternates: Iterable[str],
) -> bool:
    """Check indexed TSV exclusions at one VCF record without building a set."""
    if pysam is None:
        raise RuntimeError("Indexed TSV exclusions require pysam")
    canonical = canonical_chromosome(chromosome)
    alternate_set = set(alternates)
    for path, reader in readers:
        try:
            lines = reader.fetch(chromosome, position - 1, record_end)
        except ValueError:
            # This source has no records on the requested contig.
            continue
        for line in lines:
            columns = line.rstrip("\n").split("\t")
            try:
                target = parse_exclusion_tsv_row(columns)
            except ValueError as exc:
                raise ValueError(f"{path}: indexed TSV row {exc}") from exc
            if isinstance(target, ExclusionInterval):
                if (
                    target.chromosome == canonical
                    and target.start <= record_end
                    and target.end >= position
                ):
                    return True
            elif (
                target.chromosome == canonical
                and target.position == position
                and target.reference == reference
                and alternate_set.intersection(target.alternates)
            ):
                return True
    return False


def project_vcf_records(
    cohort: Iterable[str],
    output: TextIO,
    intervals: dict[str, tuple[list[int], list[int]]],
    variants: set[tuple[str, int, str, str]],
    source_name: str,
    indexed_tsv_paths: Iterable[Path] = (),
) -> None:
    """Write matching records from an already-filtered, sites-only VCF stream."""
    if indexed_tsv_paths and pysam is None:
        raise RuntimeError("Indexed TSV exclusions require pysam")
    with ExitStack() as stack:
        readers = [
            (path, stack.enter_context(pysam.TabixFile(str(path))))
            for path in indexed_tsv_paths
        ]
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
            indexed_match = indexed_tsv_matches_variant(
                readers, columns[0], position, record_end, columns[3], columns[4].split(",")
            ) if readers else False
            if interval_match or variant_match or indexed_match:
                output.write(line)


def write_projected_vcf(
    cohort_path: Path,
    output_path: Path,
    intervals: dict[str, tuple[list[int], list[int]]],
    variants: set[tuple[str, int, str, str]],
) -> None:
    """Project a VCF path onto exclusion targets for direct library callers."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open_text(cohort_path, errors="replace") as cohort, output_path.open("w", encoding="utf-8") as output:
        project_vcf_records(cohort, output, intervals, variants, str(cohort_path))
