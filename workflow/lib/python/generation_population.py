"""Indexed population-panel parsing for partition-local generation work."""

from __future__ import annotations

import math
from bisect import bisect_right
from pathlib import Path
from typing import Iterable, Iterator
import warnings

import polars as pl

from genomic_contracts import allele_key, canonical_chromosome
from generation_resource_policy import check_unindexed_fixture
from io_contracts import open_text

try:
    import pysam
except ImportError:  # The fixture reader supports plain, unindexed text files.
    pysam = None


POPULATION_SCHEMA = {
    "CHROM": pl.String,
    "POS": pl.Int64,
    "REF": pl.String,
    "ALT": pl.String,
    "allele_state": pl.String,
    "af": pl.Float64,
    "af_overall": pl.Float64,
    "af_female": pl.Float64,
    "af_male": pl.Float64,
    "state_reason": pl.String,
    "af_zero_imputed_tags": pl.String,
}


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
    Tabix can return one long record for every interval its reference span
    overlaps; retain it only for the interval containing its one-based POS
    anchor so disjoint fetches cannot duplicate an allele state.
    """
    regions = [(start, end)] if intervals is None else intervals
    if not regions:
        return
    if index:
        if pysam is None:
            raise RuntimeError(
                "Indexed generation resources require pysam in the generation environment"
            )
        with pysam.TabixFile(str(path), index=str(index)) as reader:
            names = {canonical_chromosome(name): name for name in reader.contigs}
            source_chromosome = names.get(chromosome)
            if source_chromosome is not None:
                for region_start, region_end in regions:
                    for raw in reader.fetch(
                        source_chromosome, region_start, region_end
                    ):
                        fields = raw.split("\t", 2)
                        try:
                            position = int(fields[1])
                        except (IndexError, ValueError):
                            # Preserve the caller's existing malformed-row diagnostic.
                            yield raw
                            continue
                        if region_start < position <= region_end:
                            yield raw
        return

    check_unindexed_fixture(path)
    with open_text(path) as handle:
        yield from handle


def window_records(
    path: str | Path,
    index: str | Path | None,
    chromosome: str,
    start: int,
    end: int,
    *,
    min_columns: int,
    intervals: list[tuple[int, int]] | None = None,
) -> Iterator[tuple[int, list[str], tuple[str, int, str, str]]]:
    """Yield parsed allele rows owned by the one-based window ``(start, end]``.

    The indexed query remains zero-based and half-open, while resource alleles
    are VCF-style one-based positions.  ``intervals`` narrows indexed fetches;
    callers that need exact disjoint-interval membership retain that policy.
    """
    chromosome = canonical_chromosome(chromosome)
    for line_number, raw in enumerate(
        indexed_lines(path, index, chromosome, start, end, intervals=intervals),
        start=1,
    ):
        if not raw.strip() or raw.startswith("#"):
            continue
        fields = raw.rstrip("\n").split("\t")
        if len(fields) < min_columns:
            raise ValueError(f"{path}: line {line_number} needs {min_columns} columns")
        key = allele_key(fields[0], int(fields[1]), fields[2], fields[3])
        if key[0] == chromosome and start < key[1] <= end:
            yield line_number, fields, key


def population_info_definitions(
    path: str | Path, info_tags: Iterable[str]
) -> dict[str, tuple[str, str]]:
    """Return validated INFO definitions for the requested tags in one header pass."""
    remaining = set(info_tags)
    definitions: dict[str, tuple[str, str]] = {}
    with open_text(path) as handle:
        for raw in handle:
            if raw.startswith("##INFO=<ID="):
                fields = dict(
                    field.split("=", 1)
                    for field in raw.rstrip("\n")
                    .removeprefix("##INFO=<")
                    .removesuffix(">")
                    .split(",")
                    if "=" in field
                )
                info_tag = fields.get("ID")
                if info_tag not in remaining:
                    continue
                number, value_type = fields.get("Number"), fields.get("Type")
                if number not in {"A", "R", "1"}:
                    raise ValueError(
                        f"Population VCF INFO {info_tag!r} must declare Number=A, R, or 1; got {number!r}"
                    )
                if value_type not in {"Float", "Integer"}:
                    raise ValueError(
                        f"Population VCF INFO {info_tag!r} must have numeric Type=Float or Integer; got {value_type!r}"
                    )
                definitions[info_tag] = (number, value_type)
                remaining.remove(info_tag)
                if not remaining:
                    return definitions
            if raw.startswith("#CHROM"):
                break
    missing = sorted(remaining)[0]
    raise ValueError(
        f"Population VCF has no INFO header for configured tag {missing!r}"
    )


def population_info_definition(path: str | Path, info_tag: str) -> tuple[str, str]:
    """Return the declared cardinality and numeric type for a selected INFO tag."""
    return population_info_definitions(path, [info_tag])[info_tag]


def validate_population_af_headers(path: str | Path, af_tags: dict[str, str]) -> None:
    """Confirm every configured population AF tag is declared before partition work."""
    population_info_definitions(path, set(af_tags.values()))


def selected_info_values(info: str, tags: Iterable[str]) -> dict[str, str]:
    """Extract just the configured INFO fields without materializing every tag."""
    padded = ";" + info
    values: dict[str, str] = {}
    for tag in tags:
        marker = ";" + tag + "="
        hit = padded.find(marker)
        if hit < 0:
            continue
        value_start = hit + len(marker)
        value_end = padded.find(";", value_start)
        values[tag] = padded[value_start : value_end if value_end >= 0 else None]
    return values


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
        raise ValueError(
            f"Population VCF has scalar selected AF on a multiallelic record at {location}"
        )
    else:
        raise ValueError(
            f"Population VCF has invalid selected AF cardinality at {location}"
        )
    parsed: list[float] = []
    for value in selected:
        try:
            frequency = float(value)
        except ValueError as exc:
            raise ValueError(
                f"Population VCF has nonnumeric selected AF at {location}"
            ) from exc
        if not math.isfinite(frequency) or not 0 <= frequency <= 1:
            raise ValueError(
                f"Population VCF has invalid selected AF outside [0, 1] at {location}"
            )
        parsed.append(frequency)
    return parsed


def selected_af_values_or_zero(
    values_by_tag: dict[str, str],
    tag: str,
    number: str,
    alternates: list[str],
    chromosome: str,
    position: int,
) -> tuple[list[float], bool]:
    """Read an AF field, treating an absent INFO tag as zero for every ALT.

    gnomAD omits population-specific INFO tags when no allele was observed in
    that population.  A present-but-empty or malformed value remains invalid.
    """
    if tag not in values_by_tag:
        return [0.0] * len(alternates), True
    return selected_af_values(
        values_by_tag, tag, number, alternates, chromosome, position
    ), False


def population_rows(
    path: str | Path,
    index: str | Path | None,
    chromosome: str,
    start: int,
    end: int,
    af_tags: dict[str, str],
    x_af_policy: str = "overall_proxy",
    hemizygote: bool = False,
    *,
    intervals: list[tuple[int, int]] | None = None,
) -> pl.DataFrame:
    """Extract observed/rejected population states for one genomic window."""
    if x_af_policy not in {"sex_specific", "overall_proxy"}:
        raise ValueError(f"Unsupported X chromosome AF policy {x_af_policy!r}")
    if not af_tags.get("overall"):
        raise ValueError("Population AF configuration requires an overall INFO tag")
    configured_tags = {"overall": af_tags["overall"]}
    if x_af_policy == "sex_specific":
        if not af_tags.get("female"):
            raise ValueError(
                "sex_specific X policy requires generation.af_source.female"
            )
        configured_tags["female"] = af_tags["female"]
        if hemizygote:
            if not af_tags.get("male"):
                raise ValueError(
                    "HEM generation with sex_specific X policy requires generation.af_source.male"
                )
            configured_tags["male"] = af_tags["male"]
    tag_definitions = population_info_definitions(path, set(configured_tags.values()))
    definitions = {role: tag_definitions[tag] for role, tag in configured_tags.items()}
    columns: dict[str, list[object]] = {name: [] for name in POPULATION_SCHEMA}
    imputed_records: list[tuple[str, int, tuple[str, ...]]] = []
    interval_starts = (
        [region_start for region_start, _ in intervals] if intervals is not None else []
    )

    def in_requested_intervals(position: int) -> bool:
        if intervals is None:
            return True
        interval_index = bisect_right(interval_starts, position - 1) - 1
        return interval_index >= 0 and position <= intervals[interval_index][1]

    def append_rows(
        keys: list[tuple[str, int, str, str]],
        allele_state: str,
        state_reason: str,
        af_values: dict[str, list[float]] | None = None,
        imputed_tags: tuple[str, ...] = (),
    ) -> None:
        for alternate_index, key in enumerate(keys):
            columns["CHROM"].append(key[0])
            columns["POS"].append(key[1])
            columns["REF"].append(key[2])
            columns["ALT"].append(key[3])
            columns["allele_state"].append(allele_state)
            if af_values is None:
                columns["af"].append(None)
                columns["af_overall"].append(None)
                columns["af_female"].append(None)
                columns["af_male"].append(None)
                columns["af_zero_imputed_tags"].append(None)
            else:
                columns["af"].append(af_values["overall"][alternate_index])
                columns["af_overall"].append(af_values["overall"][alternate_index])
                columns["af_female"].append(
                    af_values.get("female", [None] * len(keys))[alternate_index]
                )
                columns["af_male"].append(
                    af_values.get("male", [None] * len(keys))[alternate_index]
                )
                columns["af_zero_imputed_tags"].append(
                    ",".join(imputed_tags) if imputed_tags else None
                )
            columns["state_reason"].append(state_reason)

    for raw in indexed_lines(path, index, chromosome, start, end, intervals=intervals):
        if raw.startswith("#"):
            continue
        fields = raw.rstrip("\n").split("\t")
        if len(fields) < 8:
            raise ValueError("Population VCF record has fewer than eight columns")
        record_chromosome = canonical_chromosome(fields[0])
        position = int(fields[1])
        if (
            record_chromosome != chromosome
            or not start < position <= end
            or not in_requested_intervals(position)
        ):
            continue
        reference, alternates, filter_value, info = (
            fields[3],
            fields[4].split(","),
            fields[6],
            fields[7],
        )
        keys = [
            allele_key(record_chromosome, position, reference, alternate)
            for alternate in alternates
        ]
        if filter_value != "PASS":
            append_rows(keys, "rejected", "non_pass_filter")
            continue
        required_roles = ["overall"]
        if record_chromosome == "X" and x_af_policy == "sex_specific":
            required_roles.append("female")
            if hemizygote:
                required_roles.append("male")
        values_by_tag = selected_info_values(
            info, (configured_tags[role] for role in required_roles)
        )
        try:
            selected_with_provenance = {
                role: selected_af_values_or_zero(
                    values_by_tag,
                    configured_tags[role],
                    definitions[role][0],
                    alternates,
                    record_chromosome,
                    position,
                )
                for role in required_roles
            }
        except ValueError as exc:
            append_rows(keys, "invalid", str(exc))
            continue
        selected = {
            role: values for role, (values, _) in selected_with_provenance.items()
        }
        invalid_sum_roles = [
            role for role, values in selected.items() if sum(values) > 1 + 1e-12
        ]
        if invalid_sum_roles:
            reason = (
                "Population VCF selected AF sum exceeds 1 for "
                + ", ".join(sorted(invalid_sum_roles))
                + f" at {record_chromosome}:{position}"
            )
            append_rows(keys, "invalid", reason)
            continue
        imputed_tags = tuple(
            configured_tags[role]
            for role, (_, imputed) in selected_with_provenance.items()
            if imputed
        )
        if imputed_tags:
            imputed_records.append((record_chromosome, position, imputed_tags))
        append_rows(
            keys,
            "observed",
            "selected_af_zero_imputed" if imputed_tags else "valid_selected_af",
            selected,
            imputed_tags,
        )
    if imputed_records:
        examples = ", ".join(
            f"{chromosome}:{position} INFO/{','.join(tags)}"
            for chromosome, position, tags in imputed_records[:10]
        )
        warnings.warn(
            "Population VCF omitted configured AF INFO tag(s) for "
            f"{len(imputed_records)} passing record(s); using AF=0 for their alternate allele(s). "
            f"The final bundle will record any affected calibration alleles. Examples: {examples}",
            RuntimeWarning,
            stacklevel=2,
        )
    return pl.DataFrame(columns, schema=POPULATION_SCHEMA, strict=False)
