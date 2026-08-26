#!/usr/bin/env python3
"""Validate cohort VCF/PED/FASTA inputs and write the sample manifest."""

from __future__ import annotations

import argparse
from pathlib import Path

from validate_inputs import (
    parse_bool,
    read_ped,
    read_samples_tsv,
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
    result.add_argument("--sample-metadata", required=True)
    result.add_argument("--sample-metadata-kind", choices=["ped", "samples"], required=True)
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
    metadata = (
        read_ped(args.sample_metadata, errors, warnings)
        if args.sample_metadata_kind == "ped"
        else read_samples_tsv(args.sample_metadata, errors, warnings)
    )
    vcf_samples, vcf_chromosomes = read_vcf(args.vcf, errors)
    validate_hemizygote_chromosomes(hemizygote, vcf_chromosomes, errors)
    validate_vcf_samples(
        vcf_samples,
        metadata["sample_ids"],
        errors,
        metadata_label="PED" if args.sample_metadata_kind == "ped" else "Sample sheet",
    )
    validate_vcf_resource(args.vcf, "Input VCF", errors, warnings)
    validate_fasta(args.fasta, errors)
    validate_normalized_split_vcf(args.vcf, args.fasta, errors, warnings)

    report = {
        "status": "failed" if errors else "ok",
        "errors": errors,
        "warnings": warnings,
        "summary": {
            "n_ped_samples": len(metadata["all_sample_ids"]),
            "n_processed_samples": len(metadata["sample_ids"]),
            "n_excluded_phenotype_samples": len(
                metadata["excluded_phenotype_sample_ids"]
            ),
            "n_vcf_samples": len(vcf_samples),
            "n_affected_samples": metadata["affected_count"],
            "n_male_samples": metadata["male_count"],
            "sex_counts": metadata["sex_counts"],
            "phenotype_counts": metadata["phenotype_counts"],
            "vcf_chromosomes": vcf_chromosomes,
        },
    }
    write_report(report, args.out_json)
    write_manifest(metadata["manifest"], args.out_manifest)
    print(
        f"Cohort validation wrote {Path(args.out_json)} with "
        f"{len(errors)} error(s) and {len(warnings)} warning(s)."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
