#!/usr/bin/env python3
"""Validate interval, exclusion, and CADD resources independently."""

from __future__ import annotations

import argparse
from pathlib import Path

from validate_inputs import (
    parse_bool,
    read_vcf,
    validate_annotation_chromosome_names,
    validate_cadd_chromosome_names,
    validate_cadd_resources,
    validate_reference_tables,
    validate_region_chromosomes,
    write_report,
)


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--vcf", required=True)
    result.add_argument("--coverage", default="")
    result.add_argument("--variants-exclude", default="")
    result.add_argument("--testing-unit-bed", required=True)
    result.add_argument("--par-regions", default="")
    result.add_argument("--cadd-path", required=True)
    result.add_argument("--unit", required=True)
    result.add_argument("--include-indels", required=True)
    result.add_argument("--out-json", required=True)
    return result


def main() -> int:
    args = parser().parse_args()
    errors: list[str] = []
    warnings: list[str] = []
    include_indels = parse_bool(args.include_indels, "include_indels", errors)

    reference_chromosomes = validate_reference_tables(args, errors, warnings)
    vcf_read_errors: list[str] = []
    _, vcf_chromosomes = read_vcf(args.vcf, vcf_read_errors)
    errors.extend(vcf_read_errors)
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
        args.cadd_path,
        include_indels,
        errors,
        warnings,
    )
    for filename, chromosomes in cadd_chromosomes.items():
        validate_cadd_chromosome_names(
            chromosomes,
            vcf_chromosomes,
            f"CADD resource {filename}",
            errors,
            warnings,
        )

    report = {
        "status": "failed" if errors else "ok",
        "errors": errors,
        "warnings": warnings,
        "summary": {
            "include_indels": include_indels,
            "resource_vcf_chromosomes": vcf_chromosomes,
            "cadd_chromosomes": {
                filename: sorted(chromosomes)
                for filename, chromosomes in cadd_chromosomes.items()
            },
        },
    }
    write_report(report, args.out_json)
    print(
        f"Resource validation wrote {Path(args.out_json)} with "
        f"{len(errors)} error(s) and {len(warnings)} warning(s)."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
