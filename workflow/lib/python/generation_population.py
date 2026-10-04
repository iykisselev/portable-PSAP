"""Indexed population-panel parsing for partition-local generation work."""

from __future__ import annotations

import gzip
import math
from pathlib import Path
from typing import Iterator

import polars as pl

from genomic_contracts import allele_key, canonical_chromosome
from generation_resource_policy import check_unindexed_fixture

try:
    import pysam
except ImportError:  # The fixture reader supports plain, unindexed text files.
    pysam = None


POPULATION_SCHEMA = {
    "CHROM": pl.String, "POS": pl.Int64, "REF": pl.String, "ALT": pl.String,
    "allele_state": pl.String, "af": pl.Float64, "af_overall": pl.Float64,
    "af_female": pl.Float64, "af_male": pl.Float64, "state_reason": pl.String,
}


def open_text(path: str | Path):
    """Open plain or gzip-compressed text resources with the same iterator API."""
    resource = Path(path)
    if resource.suffix == ".gz":
        return gzip.open(resource, "rt", encoding="utf-8")
    return resource.open(encoding="utf-8")


def indexed_lines(
    path: str | Path,
    index: str | Path | None,
    chromosome: str,
    start: int,
    end: int,
    *,
    intervals: list[tuple[int, int]] | None = None,
) -> Iterator[str]:
    """Read a window or disjoint subintervals using one open indexed reader.

    Callers supplying intervals must merge overlaps first. Unindexed fixtures
    are scanned once; callers remain responsible for filtering their rows.
    """
    regions = [(start, end)] if intervals is None else intervals
    if not regions:
        return
    if index:
        if pysam is None:
            raise RuntimeError("Indexed generation resources require pysam in the generation environment")
        with pysam.TabixFile(str(path), index=str(index)) as reader:
            names = {canonical_chromosome(name): name for name in reader.contigs}
            source_chromosome = names.get(chromosome)
            if source_chromosome is not None:
                for region_start, region_end in regions:
                    yield from reader.fetch(source_chromosome, region_start, region_end)
        return

    check_unindexed_fixture(path)
    with open_text(path) as handle:
        yield from handle


def population_info_definition(path: str | Path, info_tag: str) -> tuple[str, str]:
    """Return the declared cardinality and numeric type for a selected INFO tag."""
    with open_text(path) as handle:
        for raw in handle:
            if raw.startswith(f"##INFO=<ID={info_tag},"):
                fields = dict(
                    field.split("=", 1)
                    for field in raw.rstrip("\n").removeprefix("##INFO=<").removesuffix(">").split(",")
                    if "=" in field
                )
                number, value_type = fields.get("Number"), fields.get("Type")
                if number not in {"A", "R", "1"}:
                    raise ValueError(
                        f"Population VCF INFO {info_tag!r} must declare Number=A, R, or 1; got {number!r}"
                    )
                if value_type not in {"Float", "Integer"}:
                    raise ValueError(
                        f"Population VCF INFO {info_tag!r} must have numeric Type=Float or Integer; got {value_type!r}"
                    )
                return number, value_type
            if raw.startswith("#CHROM"):
                break
    raise ValueError(f"Population VCF has no INFO header for configured tag {info_tag!r}")


def population_info_number(path: str | Path, info_tag: str) -> str:
    """Backward-compatible cardinality accessor for callers outside the workflow."""
    return population_info_definition(path, info_tag)[0]


def selected_af_values(
    values_by_tag: dict[str, str],
    tag: str,
    number: str,
    alternates: list[str],
    chromosome: str,
    position: int,
) -> list[float]:
    """Parse one selected AF field or fail with a record-specific error."""
    raw_value = values_by_tag.get(tag)
    location = f"{chromosome}:{position} INFO/{tag}"
    if raw_value in {None, "", "."}:
        raise ValueError(f"Population VCF has missing selected AF at {location}")
    values = raw_value.split(",")
    if number == "A" and len(values) == len(alternates):
        selected = values
    elif number == "R" and len(values) == len(alternates) + 1:
        selected = values[1:]
    elif number == "1" and len(alternates) == 1 and len(values) == 1:
        selected = values
    elif number == "1" and len(alternates) > 1:
        raise ValueError(f"Population VCF has scalar selected AF on a multiallelic record at {location}")
    else:
        raise ValueError(f"Population VCF has invalid selected AF cardinality at {location}")
    parsed: list[float] = []
    for value in selected:
        try:
            frequency = float(value)
        except ValueError as exc:
            raise ValueError(f"Population VCF has nonnumeric selected AF at {location}") from exc
        if not math.isfinite(frequency) or not 0 <= frequency <= 1:
            raise ValueError(f"Population VCF has invalid selected AF outside [0, 1] at {location}")
        parsed.append(frequency)
    return parsed


def population_rows(
    path: str | Path,
    index: str | Path | None,
    chromosome: str,
    start: int,
    end: int,
    af_tags: dict[str, str],
    x_af_policy: str = "overall_proxy",
    hemizygote: bool = False,
) -> pl.DataFrame:
    """Extract observed/rejected population states for one genomic window."""
    if x_af_policy not in {"sex_specific", "overall_proxy"}:
        raise ValueError(f"Unsupported X chromosome AF policy {x_af_policy!r}")
    if not af_tags.get("overall"):
        raise ValueError("Population AF configuration requires an overall INFO tag")
    configured_tags = {"overall": af_tags["overall"]}
    if x_af_policy == "sex_specific":
        if not af_tags.get("female"):
            raise ValueError("sex_specific X policy requires generation.af_source.female")
        configured_tags["female"] = af_tags["female"]
        if hemizygote:
            if not af_tags.get("male"):
                raise ValueError("HEM generation with sex_specific X policy requires generation.af_source.male")
            configured_tags["male"] = af_tags["male"]
    definitions = {role: population_info_definition(path, tag) for role, tag in configured_tags.items()}
    records: list[dict[str, object]] = []
    for raw in indexed_lines(path, index, chromosome, start, end):
        if raw.startswith("#"):
            continue
        fields = raw.rstrip("\n").split("\t")
        if len(fields) < 8:
            raise ValueError("Population VCF record has fewer than eight columns")
        record_chromosome = canonical_chromosome(fields[0])
        position = int(fields[1])
        if record_chromosome != chromosome or not start < position <= end:
            continue
        reference, alternates, filter_value, info = fields[3], fields[4].split(","), fields[6], fields[7]
        keys = [allele_key(record_chromosome, position, reference, alternate) for alternate in alternates]
        if filter_value != "PASS":
            records.extend(
                {
                    "CHROM": key[0], "POS": key[1], "REF": key[2], "ALT": key[3],
                    "allele_state": "rejected", "af": None, "af_overall": None,
                    "af_female": None, "af_male": None, "state_reason": "non_pass_filter",
                }
                for key in keys
            )
            continue
        values_by_tag = {
            item.split("=", 1)[0]: item.split("=", 1)[1]
            for item in info.split(";") if "=" in item
        }
        required_roles = ["overall"]
        if record_chromosome == "X" and x_af_policy == "sex_specific":
            required_roles.append("female")
            if hemizygote:
                required_roles.append("male")
        try:
            selected = {
                role: selected_af_values(values_by_tag, configured_tags[role], definitions[role][0], alternates, record_chromosome, position)
                for role in required_roles
            }
        except ValueError as exc:
            records.extend(
                {
                    "CHROM": key[0], "POS": key[1], "REF": key[2], "ALT": key[3],
                    "allele_state": "invalid", "af": None, "af_overall": None,
                    "af_female": None, "af_male": None, "state_reason": str(exc),
                }
                for key in keys
            )
            continue
        for alternate_index, key in enumerate(keys):
            records.append(
                {
                    "CHROM": key[0], "POS": key[1], "REF": key[2], "ALT": key[3],
                    "allele_state": "observed", "af": selected["overall"][alternate_index],
                    "af_overall": selected["overall"][alternate_index],
                    "af_female": selected.get("female", [None] * len(keys))[alternate_index],
                    "af_male": selected.get("male", [None] * len(keys))[alternate_index],
                    "state_reason": "valid_selected_af",
                }
            )
    return pl.DataFrame(records, schema=POPULATION_SCHEMA, strict=False)
