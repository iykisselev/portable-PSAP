"""Genomic-resource validation for the apply barrier."""

from __future__ import annotations

from pathlib import Path
from typing import Sequence

from validation import (
    IndexedResource,
    read_vcf,
    validate_annotation_chromosome_names,
    validate_cadd_chromosome_names,
    validate_cadd_resources,
    validate_reference_tables,
    validate_region_chromosomes,
)
from validation_reporting import build_validation_report


def validate_resource_inputs(
    *,
    vcf: str | Path,
    calibration_regions: IndexedResource | None,
    exclusions: Sequence[IndexedResource],
    testing_units: IndexedResource,
    par_regions: Sequence[IndexedResource],
    cadd_snv: IndexedResource,
    cadd_indel: IndexedResource | None,
    unit: str,
    include_indels: bool,
) -> dict[str, object]:
    """Validate reference formats and chromosome compatibility in one category."""
    errors: list[str] = []
    warnings: list[str] = []
    reference_chromosomes = validate_reference_tables(
        calibration_regions=calibration_regions,
        exclusions=exclusions,
        testing_units=testing_units,
        par_regions=par_regions,
        unit=unit,
        errors=errors,
        warnings=warnings,
    )
    vcf_errors: list[str] = []
    _, vcf_chromosomes = read_vcf(vcf, vcf_errors)
    errors.extend(vcf_errors)
    validate_region_chromosomes(
        reference_chromosomes["exclusions"],
        vcf_chromosomes,
        "Exclusion interval files",
        warnings,
    )
    validate_annotation_chromosome_names(
        reference_chromosomes["testing_unit"],
        vcf_chromosomes,
        f"Configured {unit} testing-unit BED file",
        warnings,
    )
    cadd_chromosomes = validate_cadd_resources(
        cadd_snv.path,
        cadd_indel.path if cadd_indel else None,
        include_indels,
        errors,
        warnings,
        cadd_snv.index,
        cadd_indel.index if cadd_indel else None,
    )
    for filename, chromosomes in cadd_chromosomes.items():
        validate_cadd_chromosome_names(
            chromosomes,
            vcf_chromosomes,
            f"CADD resource {filename}",
            errors,
            warnings,
        )
    return build_validation_report(
        errors=errors,
        warnings=warnings,
        summary={
            "include_indels": include_indels,
            "resource_vcf_chromosomes": vcf_chromosomes,
            "cadd_chromosomes": {
                filename: sorted(chromosomes)
                for filename, chromosomes in cadd_chromosomes.items()
            },
        },
        sections=("resources",),
    )
