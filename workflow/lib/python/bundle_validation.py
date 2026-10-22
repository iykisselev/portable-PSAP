"""Semantic validation and certification of completed generated bundles."""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path

import jsonschema
import yaml

from generation_bindings import BundleValidationInputs, BundleValidationOutputs
from bundle_certificate import BundleCertificateError, validate_bundle_component
from bundle_contract import fixed_width_bin_count
from genomic_contracts import canonical_chromosome
from model_catalog import enabled_models, models_for_chromosome, validate_status
from io_contracts import open_text, sha256_file, write_text_atomic


def validate_lookup(path: str | Path, bin_count: int) -> set[str]:
    """Validate one compressed lookup row at a time and retain only unit IDs."""
    identifiers: set[str] = set()
    with open_text(path, newline="", require_gzip=True) as handle:
        rows = csv.reader(handle, delimiter="\t")
        header = next(rows, None)
        if not header or header[0] != "testing_unit_id" or len(header) != bin_count + 1:
            raise ValueError(f"Lookup header does not match {bin_count} bins: {path}")
        for row in rows:
            if len(row) != bin_count + 1 or not row[0] or row[0] in identifiers:
                raise ValueError(f"Lookup contains an invalid or duplicate row: {path}")
            values = [float(value) for value in row[1:]]
            if any(
                not math.isfinite(value) or value < 0 or value > 1 for value in values
            ):
                raise ValueError(f"Lookup contains values outside [0, 1]: {path}")
            if any(
                left < right for left, right in zip(values, values[1:], strict=False)
            ):
                raise ValueError(f"Lookup is not nonincreasing: {path}")
            identifiers.add(row[0])
    return identifiers


def data_row_count(path: str | Path) -> int:
    """Count compressed data rows without retaining their values."""
    with open_text(path, newline="", require_gzip=True) as handle:
        rows = csv.reader(handle, delimiter="\t")
        next(rows, None)
        return sum(1 for _ in rows)


def _read_partition_bed(path: str | Path, label: str) -> list[tuple[str, int, int, str]]:
    """Read one certified four-column BED while retaining each interval."""
    rows: list[tuple[str, int, int, str]] = []
    with open_text(path, newline="", require_gzip=True) as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader, None)
        if header != ["#CHROM", "START", "END", "testing_unit_id"]:
            raise ValueError(f"{label} has an invalid BED header")
        for row in reader:
            if len(row) != 4 or not row[3]:
                raise ValueError(f"{label} contains an invalid BED row")
            try:
                start, end = int(row[1]), int(row[2])
            except ValueError as exc:
                raise ValueError(f"{label} contains non-integer coordinates") from exc
            if start < 0 or end <= start:
                raise ValueError(f"{label} contains invalid coordinates")
            rows.append((canonical_chromosome(row[0]), start, end, row[3]))
    return rows


def _validate_resolved_ownership(
    partitions: list[dict[str, object]], components: dict[str, str | Path]
) -> dict[str, set[tuple[str, int, int, str]]]:
    """Require resolved ownership to be non-overlapping across all partitions."""
    by_partition: dict[str, set[tuple[str, int, int, str]]] = {}
    all_rows: list[tuple[str, int, int, str, str]] = []
    for partition in partitions:
        identifier = str(partition["id"])
        rows = _read_partition_bed(
            components[str(partition["resolved_component"])],
            f"resolved partition {identifier}",
        )
        row_set = set(rows)
        if len(row_set) != len(rows):
            raise ValueError(f"Resolved partition {identifier} contains duplicate intervals")
        by_partition[identifier] = row_set
        all_rows.extend((*row, identifier) for row in rows)
    for (chromosome, unit), rows in _group_intervals(all_rows).items():
        previous_end = -1
        for start, end, partition in rows:
            if start < previous_end:
                raise ValueError(
                    "Resolved partition intervals overlap for "
                    f"{chromosome}/{unit}; ownership is not disjoint"
                )
            previous_end = end
    return by_partition


def _group_intervals(rows: list[tuple[str, int, int, str, str]]) -> dict[tuple[str, str], list[tuple[int, int, str]]]:
    grouped: dict[tuple[str, str], list[tuple[int, int, str]]] = {}
    for chromosome, start, end, unit, partition in rows:
        grouped.setdefault((chromosome, unit), []).append((start, end, partition))
    for values in grouped.values():
        values.sort()
    return grouped


def _validate_eligible_subset(
    partitions: list[dict[str, object]], components: dict[str, str | Path],
    resolved: dict[str, set[tuple[str, int, int, str]]],
) -> dict[str, set[str]]:
    """Require each eligible interval to be fully covered by its resolved owner."""
    eligible_ids: dict[str, set[str]] = {}
    for partition in partitions:
        identifier = str(partition["id"])
        eligible = _read_partition_bed(
            components[str(partition["eligible_component"])],
            f"eligible partition {identifier}",
        )
        resolved_by_key: dict[tuple[str, str], list[tuple[int, int]]] = {}
        for chromosome, start, end, unit in resolved[identifier]:
            resolved_by_key.setdefault((chromosome, unit), []).append((start, end))
        for intervals in resolved_by_key.values():
            intervals.sort()
        identifiers: set[str] = set()
        for chromosome, start, end, unit in eligible:
            intervals = resolved_by_key.get((chromosome, unit), [])
            if not any(left <= start and end <= right for left, right in intervals):
                raise ValueError(
                    f"Eligible partition {identifier} has an interval outside resolved ownership"
                )
            identifiers.add(unit)
        eligible_ids[identifier] = identifiers
    return eligible_ids


def _read_status_rows(path: str | Path, partition_ids: list[str]) -> list[dict[str, str]]:
    required = {
        "testing_unit_id", "partition_id", "CHROM", "raw_bases",
        "priority_removed_bases", "resolved_bases", "supported_bases",
        "retained_bases", "retained_fraction", "status", "failed_roles_json",
    }
    with open_text(path, newline="", require_gzip=True) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None or required - set(reader.fieldnames):
            raise ValueError("Partition status table has an invalid header")
        rows = list(reader)
    seen: set[tuple[str, str]] = set()
    for row in rows:
        key = (row["testing_unit_id"], row["partition_id"])
        if not all(key) or key in seen or key[1] not in partition_ids:
            raise ValueError("Partition status table has invalid or duplicate pair keys")
        seen.add(key)
        try:
            failed_roles = json.loads(row["failed_roles_json"])
        except json.JSONDecodeError as exc:
            raise ValueError("Partition status table has invalid failed_roles_json") from exc
        if not isinstance(failed_roles, list) or any(not isinstance(role, str) for role in failed_roles):
            raise ValueError("Partition status failed_roles_json must be an array of strings")
        validate_status(row["status"], "scope", failed_roles=tuple(failed_roles))
        numeric = ("raw_bases", "priority_removed_bases", "resolved_bases", "supported_bases", "retained_bases", "retained_fraction")
        try:
            values = {name: float(row[name]) for name in numeric}
        except ValueError as exc:
            raise ValueError("Partition status table has nonnumeric scope fields") from exc
        if any(not math.isfinite(value) or value < 0 for value in values.values()):
            raise ValueError("Partition status table has invalid scope values")
        if values["priority_removed_bases"] + values["resolved_bases"] != values["raw_bases"]:
            raise ValueError("Partition status row has inconsistent ownership bases")
        if values["retained_bases"] > values["supported_bases"] or values["supported_bases"] > values["resolved_bases"]:
            raise ValueError("Partition status row has inconsistent retained bases")
    return rows


def _validate_overlap_audit(path: str | Path, partition_ids: list[str]) -> None:
    expected = {
        "left_partition_id", "right_partition_id", "left_unique_id_count",
        "right_unique_id_count", "intersection_count", "union_count",
        "left_overlap_fraction", "right_overlap_fraction", "minimum_overlap_fraction",
        "empty_set_status", "low_overlap_warning",
    }
    with open_text(path, newline="", require_gzip=True) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None or set(reader.fieldnames) != expected:
            raise ValueError("Testing-unit ID-overlap audit has an invalid header")
        rows = list(reader)
    pairs = {(row["left_partition_id"], row["right_partition_id"]) for row in rows}
    expected_pairs = {
        (partition_ids[left], partition_ids[right])
        for left in range(len(partition_ids)) for right in range(left + 1, len(partition_ids))
    }
    if pairs != expected_pairs or len(rows) != len(pairs):
        raise ValueError("Testing-unit ID-overlap audit does not cover configured partition pairs")
    for row in rows:
        try:
            left, right, intersection, union = (
                int(row[name]) for name in (
                    "left_unique_id_count", "right_unique_id_count", "intersection_count", "union_count"
                )
            )
        except ValueError as exc:
            raise ValueError("Testing-unit ID-overlap audit has noninteger counts") from exc
        if min(left, right, intersection, union) < 0 or intersection > min(left, right) or union != left + right - intersection:
            raise ValueError("Testing-unit ID-overlap audit has inconsistent counts")


def validate_generated_bundle(
    inputs: BundleValidationInputs, outputs: BundleValidationOutputs
) -> None:
    """Verify bytes and semantic row contracts before writing the marker last."""
    if inputs.metadata is None or inputs.schema is None:
        raise ValueError("Partitioned bundle validation requires metadata and schema")
    metadata_path = Path(inputs.metadata)
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    schema = yaml.safe_load(Path(inputs.schema).read_text(encoding="utf-8"))
    jsonschema.Draft202012Validator(schema).validate(metadata)
    if metadata["contract"] != "partitioned_null_bundle_v1":
        raise ValueError("Only partitioned_null_bundle_v1 may be certified")
    binning = metadata["binning"]
    expected_count = fixed_width_bin_count(binning["score_max"], binning["bin_width"])
    if binning["bin_count"] != expected_count:
        raise ValueError(
            "Bundle bin_count is inconsistent with score_max and bin_width"
        )
    partitions = metadata["partitions"]
    partition_ids = [entry["id"] for entry in partitions]
    if partition_ids != list(dict.fromkeys(partition_ids)):
        raise ValueError("Bundle partition IDs must be unique and ordered")
    if [entry["order"] for entry in partitions] != list(range(len(partitions))):
        raise ValueError("Bundle partition orders must be contiguous from zero")
    hemizygote = "hem_chrX" in metadata["model_assumptions"]["roles"]
    expected_lookup_count = len(partitions) * len(enabled_models(hemizygote))
    if len(inputs.partition_lookups) != expected_lookup_count:
        raise ValueError("Partition/model lookup input count does not match metadata")
    if len(inputs.resolved_partitions) != len(partitions) or len(inputs.eligible_partitions) != len(partitions):
        raise ValueError("Partition BED input count does not match metadata")
    if inputs.partition_statuses is None or inputs.testing_unit_id_overlap is None:
        raise ValueError("Partitioned bundle validation requires scope audit components")
    components = {
        "calibration_regions": inputs.calibration_regions,
        "excluded_alleles": inputs.excluded_alleles,
        "candidate_audit": inputs.candidate_audit,
        "population_af_zero_imputations": inputs.population_af_zero_imputations,
        "unit_partition_status": inputs.partition_statuses,
        "testing_unit_id_overlap": inputs.testing_unit_id_overlap,
    }
    lookup_inputs = iter(inputs.partition_lookups)
    resolved_inputs = iter(inputs.resolved_partitions)
    eligible_inputs = iter(inputs.eligible_partitions)
    for partition in partitions:
        components[partition["resolved_component"]] = next(resolved_inputs)
        components[partition["eligible_component"]] = next(eligible_inputs)
        for model in enabled_models(hemizygote):
            components[partition["lookups"].get(model.identifier)] = next(lookup_inputs)
    for name, path in components.items():
        try:
            validate_bundle_component(name, path, metadata)
        except BundleCertificateError as exc:
            raise ValueError(str(exc)) from exc
    resolved_by_partition = _validate_resolved_ownership(partitions, components)
    eligible_ids_from_bed = _validate_eligible_subset(
        partitions, components, resolved_by_partition
    )
    status_rows = _read_status_rows(inputs.partition_statuses, partition_ids)
    _validate_overlap_audit(inputs.testing_unit_id_overlap, partition_ids)
    eligible_by_partition = {
        partition: {
            row["testing_unit_id"]: row["CHROM"]
            for row in status_rows
            if row["partition_id"] == partition and row["status"] == "eligible"
        }
        for partition in partition_ids
    }
    if {
        partition: set(rows)
        for partition, rows in eligible_by_partition.items()
    } != eligible_ids_from_bed:
        raise ValueError("Eligible partition BED membership differs from scope statuses")
    for partition in partitions:
        identifier = str(partition["id"])
        declared = sum(row["partition_id"] == identifier for row in status_rows)
        eligible = len(eligible_by_partition[identifier])
        if partition["pair_counts"] != {"declared": declared, "eligible": eligible}:
            raise ValueError(
                f"Partition pair counts differ from statuses for {identifier}"
            )
    for partition in partitions:
        ids_by_model = {
            model.identifier: validate_lookup(
                components[partition["lookups"][model.identifier]], binning["bin_count"]
            ) for model in enabled_models(hemizygote)
        }
        for model in enabled_models(hemizygote):
            expected = {
                unit for unit, chromosome in eligible_by_partition[partition["id"]].items()
                if any(item.identifier == model.identifier for item in models_for_chromosome(chromosome, hemizygote))
            }
            if ids_by_model[model.identifier] != expected:
                raise ValueError(
                    f"Lookup membership differs from eligible scope for {partition['id']}/{model.identifier}"
                )
    observed_counts = {
        "eligible_pairs": sum(len(value) for value in eligible_by_partition.values()),
        "declared_pairs": len(status_rows),
        "excluded_alleles": data_row_count(components["excluded_alleles"]),
        "calibration_intervals": data_row_count(components["calibration_regions"]),
        "population_af_zero_imputations": data_row_count(
            components["population_af_zero_imputations"]
        ),
    }
    if metadata["counts"] != observed_counts:
        raise ValueError(
            f"Bundle exclusion/scope counts differ from metadata: {observed_counts!r}"
        )
    metadata_sha256 = sha256_file(metadata_path)
    report = {
        "status": "ok",
        "contract": metadata["contract"],
        "schema_version": metadata["schema_version"],
        "metadata_sha256": metadata_sha256,
        "component_count": len(components),
        "counts": observed_counts,
    }
    write_text_atomic(
        outputs.validation, json.dumps(report, indent=2, sort_keys=True) + "\n"
    )
    validation_sha256 = sha256_file(outputs.validation)
    write_text_atomic(
        outputs.marker,
        f"{metadata['contract']}\t{metadata_sha256}\t{validation_sha256}\n",
    )
