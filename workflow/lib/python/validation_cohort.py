"""Cohort metadata, VCF, and reference validation for the apply barrier."""

from __future__ import annotations

from pathlib import Path

from sample_manifest import read_sample_manifest_index, write_sample_manifest
from validation import (
    parse_bool,
    read_ped,
    read_samples_tsv,
    read_vcf,
    validate_fasta,
    validate_hemizygote_chromosomes,
    validate_normalized_split_vcf,
    validate_vcf_resource,
    validate_vcf_samples,
)
from validation_reporting import build_validation_report


def validate_cohort_inputs(
    *,
    metadata_path: str | Path,
    metadata_kind: str,
    vcf: str | Path,
    vcf_index: str | Path | None,
    fasta: str | Path,
    fasta_index: str | Path,
    fasta_gzi: str | Path | None,
    hemizygote_model: str,
    manifest_path: str | Path,
) -> dict[str, object]:
    """Publish and validate one canonical sample manifest against the cohort VCF."""
    errors: list[str] = []
    warnings: list[str] = []
    hemizygote = parse_bool(hemizygote_model, "hemizygote_model", errors)
    if metadata_kind == "ped":
        metadata = read_ped(metadata_path, errors, warnings)
        metadata_label = "PED"
    elif metadata_kind == "samples":
        metadata = read_samples_tsv(metadata_path, errors, warnings)
        metadata_label = "Sample sheet"
    else:
        raise ValueError(f"Unsupported sample metadata kind: {metadata_kind!r}")

    # The physical manifest is the only serialized sample contract. Re-read it
    # once so VCF validation and all later workflow rules agree on its content.
    write_sample_manifest(metadata["manifest"], manifest_path, create_parents=False)
    manifest = read_sample_manifest_index(manifest_path)
    vcf_samples, vcf_chromosomes = read_vcf(vcf, errors)
    validate_hemizygote_chromosomes(hemizygote, vcf_chromosomes, errors)
    validate_vcf_samples(
        vcf_samples,
        [record.pid for record in manifest.records],
        errors,
        metadata_label=metadata_label,
    )
    validate_vcf_resource(vcf, "Input VCF", errors, warnings, vcf_index)
    validate_fasta(fasta, fasta_index, fasta_gzi, errors)
    validate_normalized_split_vcf(vcf, fasta, fasta_index, errors, warnings)
    return build_validation_report(
        errors=errors,
        warnings=warnings,
        summary={
            "n_ped_samples": len(metadata["all_sample_ids"]),
            "n_processed_samples": len(manifest.records),
            "n_excluded_phenotype_samples": len(
                metadata["excluded_phenotype_sample_ids"]
            ),
            "n_vcf_samples": len(vcf_samples),
            "n_affected_samples": sum(
                record.phenotype == "affected" for record in manifest.records
            ),
            "n_male_samples": sum(record.sex == "1" for record in manifest.records),
            "sex_counts": metadata["sex_counts"],
            "phenotype_counts": metadata["phenotype_counts"],
            "vcf_chromosomes": vcf_chromosomes,
        },
        sections=("cohort",),
    )
