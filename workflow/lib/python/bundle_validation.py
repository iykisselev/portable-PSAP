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
from model_catalog import enabled_models, models_for_chromosome
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


def eligible_units(path: str | Path) -> dict[str, set[str]]:
    """Read eligible IDs with their canonical chromosomes from the bundle BED."""
    identifiers: dict[str, set[str]] = {}
    with open_text(path, newline="", require_gzip=True) as handle:
        rows = csv.reader(handle, delimiter="\t")
        next(rows, None)
        for row in rows:
            if len(row) < 4 or not row[3]:
                raise ValueError("Eligible testing-unit BED contains an invalid row")
            identifiers.setdefault(row[3], set()).add(canonical_chromosome(row[0]))
    return identifiers


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
    status_rows = []
    with open_text(inputs.partition_statuses, newline="", require_gzip=True) as handle:
        status_rows = list(csv.DictReader(handle, delimiter="\t"))
    eligible_by_partition = {
        partition: {
            row["testing_unit_id"]: row["CHROM"]
            for row in status_rows
            if row["partition_id"] == partition and row["status"] == "eligible"
        }
        for partition in partition_ids
    }
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
