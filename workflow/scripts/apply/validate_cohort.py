#!/usr/bin/env python3
"""Validate cohort VCF/PED/FASTA inputs and write the sample manifest."""

from __future__ import annotations

import argparse
from pathlib import Path

from validate_inputs import (
    parse_bool,
    read_ped,
    read_vcf,
    validate_fasta,
    validate_hemizygote_chromosomes,
    validate_normalized_split_vcf,
    validate_vcf_resource,
    validate_vcf_samples,
    write_manifest,
    write_report,
)


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--vcf", required=True)
    result.add_argument("--ped", required=True)
    result.add_argument("--fasta", required=True)
    result.add_argument("--hemizygote-model", required=True)
    result.add_argument("--out-json", required=True)
    result.add_argument("--out-manifest", required=True)
    return result


def main() -> int:
    args = parser().parse_args()
    errors: list[str] = []
    warnings: list[str] = []

    hemizygote = parse_bool(args.hemizygote_model, "hemizygote_model", errors)
    ped = read_ped(args.ped, errors, warnings)
    vcf_samples, vcf_chromosomes = read_vcf(args.vcf, errors)
    validate_hemizygote_chromosomes(hemizygote, vcf_chromosomes, errors)
    validate_vcf_samples(vcf_samples, ped["sample_ids"], errors)
    validate_vcf_resource(args.vcf, "Input VCF", errors, warnings)
    validate_fasta(args.fasta, errors)
    validate_normalized_split_vcf(args.vcf, args.fasta, errors, warnings)

    report = {
        "status": "failed" if errors else "ok",
        "errors": errors,
        "warnings": warnings,
        "summary": {
            "n_ped_samples": len(ped["all_sample_ids"]),
            "n_processed_samples": len(ped["sample_ids"]),
            "n_excluded_phenotype_samples": len(
                ped["excluded_phenotype_sample_ids"]
            ),
            "n_vcf_samples": len(vcf_samples),
            "n_affected_samples": ped["affected_count"],
            "n_male_samples": ped["male_count"],
            "sex_counts": ped["sex_counts"],
            "phenotype_counts": ped["phenotype_counts"],
            "vcf_chromosomes": vcf_chromosomes,
        },
    }
    write_report(report, args.out_json)
    write_manifest(ped["manifest"], args.out_manifest)
    print(
        f"Cohort validation wrote {Path(args.out_json)} with "
        f"{len(errors)} error(s) and {len(warnings)} warning(s)."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
