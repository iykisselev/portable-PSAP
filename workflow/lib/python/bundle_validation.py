"""Semantic validation and certification of completed generated bundles."""

from __future__ import annotations

import csv
from decimal import Decimal, ROUND_CEILING
import gzip
import json
import math
from pathlib import Path

import jsonschema
import yaml

from generation_bindings import BundleValidationInputs, BundleValidationOutputs
from bundle_certificate import BundleCertificateError, validate_bundle_component
from bundle_contract import HEM_COMPONENT, VALIDATION_ARTIFACT_FILENAMES
from indel_bundle_contract import validate_indel_bundle_metadata
from io_contracts import sha256_file, write_text_atomic


def validate_lookup(path: str | Path, bin_count: int) -> set[str]:
    """Validate one compressed lookup row at a time and retain only unit IDs."""
    identifiers: set[str] = set()
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        rows = csv.reader(handle, delimiter="\t")
        header = next(rows, None)
        if not header or header[0] != "testing_unit_id" or len(header) != bin_count + 1:
            raise ValueError(f"Lookup header does not match {bin_count} bins: {path}")
        for row in rows:
            if len(row) != bin_count + 1 or not row[0] or row[0] in identifiers:
                raise ValueError(f"Lookup contains an invalid or duplicate row: {path}")
            values = [float(value) for value in row[1:]]
            if any(not math.isfinite(value) or value < 0 or value > 1 for value in values):
                raise ValueError(f"Lookup contains values outside [0, 1]: {path}")
            if any(left < right for left, right in zip(values, values[1:])):
                raise ValueError(f"Lookup is not nonincreasing: {path}")
            identifiers.add(row[0])
    return identifiers


def data_row_count(path: str | Path) -> int:
    """Count compressed data rows without retaining their values."""
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        rows = csv.reader(handle, delimiter="\t")
        next(rows, None)
        return sum(1 for _ in rows)


def eligible_units(path: str | Path) -> dict[str, set[str]]:
    """Read eligible IDs with their canonical chromosomes from the bundle BED."""
    identifiers: dict[str, set[str]] = {}
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        rows = csv.reader(handle, delimiter="\t")
        next(rows, None)
        for row in rows:
            if len(row) < 4 or not row[3]:
                raise ValueError("Eligible testing-unit BED contains an invalid row")
            chromosome = row[0][3:] if row[0].lower().startswith("chr") else row[0]
            identifiers.setdefault(row[3], set()).add(chromosome.upper())
    return identifiers


def validate_generated_bundle(
    inputs: BundleValidationInputs, outputs: BundleValidationOutputs
) -> None:
    """Verify bytes and semantic row contracts before writing the marker last."""
    metadata_path = Path(inputs.metadata)
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    schema = yaml.safe_load(Path(inputs.schema).read_text(encoding="utf-8"))
    jsonschema.Draft202012Validator(schema).validate(metadata)
    if metadata["contract"] == "generated_bundle_v3":
        validate_indel_bundle_metadata(metadata)
    binning = metadata["binning"]
    expected_count = int((Decimal(str(binning["score_max"])) / Decimal(str(binning["bin_width"]))).to_integral_value(rounding=ROUND_CEILING))
    if binning["bin_count"] != expected_count:
        raise ValueError("Bundle bin_count is inconsistent with score_max and bin_width")
    components = {
        "lookup_het": inputs.lookup_het, "lookup_hom": inputs.lookup_hom,
        "eligible_testing_units": inputs.eligible_bed,
        "calibration_regions": inputs.calibration_regions,
        "excluded_alleles": inputs.excluded_alleles,
        "population_af_zero_imputations": inputs.population_af_zero_imputations,
        "excluded_units": inputs.excluded_units,
    }
    hemizygote = "hem_chrX" in metadata["model_assumptions"]["roles"]
    if hemizygote:
        components[HEM_COMPONENT] = inputs.lookup_hem
    for name, path in components.items():
        try:
            validate_bundle_component(name, path, metadata)
        except BundleCertificateError as exc:
            raise ValueError(str(exc)) from exc
    lookup_ids = {role: validate_lookup(components[f"lookup_{role}"], binning["bin_count"]) for role in ("het", "hom")}
    unit_chromosomes = eligible_units(components["eligible_testing_units"])
    expected_ids = set(unit_chromosomes)
    if not expected_ids:
        raise ValueError("Bundle contains no eligible testing units")
    if lookup_ids["het"] != expected_ids or lookup_ids["hom"] != expected_ids:
        raise ValueError("Lookup role IDs do not exactly match eligible testing units")
    if hemizygote:
        expected_x_ids = {unit for unit, chromosomes in unit_chromosomes.items() if chromosomes == {"X"}}
        hem_ids = validate_lookup(components[HEM_COMPONENT], binning["bin_count"])
        if not expected_x_ids or hem_ids != expected_x_ids:
            raise ValueError("HEM lookup IDs do not exactly match eligible non-PAR chromosome-X testing units")
    observed_counts = {
        "eligible_units": len(expected_ids), "excluded_units": data_row_count(components["excluded_units"]),
        "excluded_alleles": data_row_count(components["excluded_alleles"]),
        "calibration_intervals": data_row_count(components["calibration_regions"]),
        "population_af_zero_imputations": data_row_count(components["population_af_zero_imputations"]),
    }
    if metadata["counts"] != observed_counts:
        raise ValueError(f"Bundle exclusion/scope counts differ from metadata: {observed_counts!r}")
    metadata_sha256 = sha256_file(metadata_path)
    report = {
        "status": "ok", "contract": metadata["contract"], "schema_version": metadata["schema_version"],
        "metadata_sha256": metadata_sha256, "component_count": len(components), "counts": observed_counts,
    }
    write_text_atomic(outputs.validation, json.dumps(report, indent=2, sort_keys=True) + "\n")
    validation_sha256 = sha256_file(outputs.validation)
    write_text_atomic(outputs.marker, f"{metadata['contract']}\t{metadata_sha256}\t{validation_sha256}\n")
