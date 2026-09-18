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


def _open(path: str | Path):
    file_path = Path(path)
    return gzip.open(file_path, "rt", encoding="utf-8") if file_path.suffix == ".gz" else file_path.open(encoding="utf-8")


def canonical_chromosome(chromosome: str) -> str:
    value = chromosome.strip()
    if value.lower().startswith("chr"):
        value = value[3:]
    value = value.upper()
    return "MT" if value == "M" else value


def allele_key(chromosome: str, position: int, ref: str, alt: str) -> tuple[str, int, str, str]:
    return (canonical_chromosome(chromosome), int(position), ref.upper(), alt.upper())


def key_text(key: tuple[str, int, str, str]) -> str:
    return f"{key[0]}:{key[1]}:{key[2]}:{key[3]}"


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
                raise ValueError(f"Coverage BED line {number} needs three columns")
            start, end = int(fields[1]), int(fields[2])
            if start < 0 or end <= start:
                raise ValueError(f"Invalid coverage BED line {number}")
            intervals[canonical_chromosome(fields[0])].append((start, end))
    return intervals


def merged_length(intervals: list[tuple[int, int]]) -> int:
    total = 0
    end = -1
    for start, stop in sorted(intervals):
        if stop <= end:
            continue
        total += stop - max(start, end)
        end = stop
    return total


def overlap_length(intervals: list[tuple[int, int]], coverage: list[tuple[int, int]]) -> int:
    pieces: list[tuple[int, int]] = []
    for start, end in intervals:
        for covered_start, covered_end in coverage:
            left, right = max(start, covered_start), min(end, covered_end)
            if left < right:
                pieces.append((left, right))
    return merged_length(pieces)


def coverage_rows(units, coverage, minimum: float):
    rows = []
    for unit, intervals in sorted(units.items()):
        chrom = intervals[0][0]
        simple_intervals = [(start, end) for _, start, end in intervals]
        total = merged_length(simple_intervals)
        covered = overlap_length(simple_intervals, coverage.get(chrom, []))
        fraction = covered / total
        rows.append((unit, chrom, total, covered, fraction, "eligible" if fraction >= minimum else "coverage_below_minimum"))
    return rows


def read_fasta(path: str | Path):
    sequences: dict[str, str] = {}
    current: str | None = None
    parts: list[str] = []
    with _open(path) as handle:
        for raw in handle:
            if raw.startswith(">"):
                if current is not None:
                    sequences[current] = "".join(parts).upper()
                current, parts = canonical_chromosome(raw[1:].split()[0]), []
            else:
                parts.append(raw.strip())
    if current is not None:
        sequences[current] = "".join(parts).upper()
    return sequences


def read_scores(path: str | Path):
    scores = {}
    with _open(path) as handle:
        for number, raw in enumerate(handle, start=1):
            if not raw.strip() or raw.startswith("#"):
                continue
            fields = raw.rstrip().split("\t")
            if len(fields) < 6:
                raise ValueError(f"Score line {number} needs six CADD columns")
            key = allele_key(fields[0], int(fields[1]), fields[2], fields[3])
            if len(key[2]) != 1 or len(key[3]) != 1:
                continue
            score = float(fields[5])
            if not math.isfinite(score):
                raise ValueError(f"Score line {number} has non-finite PHRED")
            if key in scores and scores[key] != score:
                raise ValueError(f"Conflicting score for {key_text(key)}")
            scores[key] = score
    return scores
