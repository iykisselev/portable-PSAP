"""Genomic-resource validation for the apply barrier."""

from __future__ import annotations

from types import SimpleNamespace

from validation import (
    parse_bool,
    read_vcf,
    validate_annotation_chromosome_names,
    validate_cadd_chromosome_names,
    validate_cadd_resources,
    validate_reference_tables,
    validate_region_chromosomes,
)
from validation_reporting import build_validation_report


def validate_resource_inputs(args: SimpleNamespace) -> dict[str, object]:
    """Validate reference formats and chromosome compatibility in one category."""
    errors: list[str] = []
    warnings: list[str] = []
    include_indels = parse_bool(args.include_indels, "include_indels", errors)
    reference_chromosomes = validate_reference_tables(args, errors, warnings)
    vcf_errors: list[str] = []
    _, vcf_chromosomes = read_vcf(args.vcf, vcf_errors)
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
        f"Configured {args.unit} testing-unit BED file",
        warnings,
    )
    cadd_chromosomes = validate_cadd_resources(
        args.cadd_snv,
        args.cadd_indel,
        include_indels,
        errors,
        warnings,
        args.cadd_snv_index,
        args.cadd_indel_index,
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
