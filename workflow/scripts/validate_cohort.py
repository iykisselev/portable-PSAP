"""Validate cohort inputs and write the manifest via injected Snakemake I/O."""

from __future__ import annotations

from apply.validate_inputs import (
    parse_bool, read_ped, read_samples_tsv, read_vcf, validate_fasta,
    validate_hemizygote_chromosomes, validate_normalized_split_vcf,
    validate_vcf_resource, validate_vcf_samples, write_manifest, write_report,
)


def main() -> None:
    errors: list[str] = []
    warnings: list[str] = []
    kind = str(snakemake.params.metadata_kind)
    hemizygote = parse_bool(str(snakemake.params.hemizygote_model), "hemizygote_model", errors)
    metadata = read_ped(snakemake.input.metadata, errors, warnings) if kind == "ped" else read_samples_tsv(snakemake.input.metadata, errors, warnings)
    vcf_samples, vcf_chromosomes = read_vcf(snakemake.input.vcf, errors)
    validate_hemizygote_chromosomes(hemizygote, vcf_chromosomes, errors)
    validate_vcf_samples(vcf_samples, metadata["sample_ids"], errors, metadata_label="PED" if kind == "ped" else "Sample sheet")
    validate_vcf_resource(snakemake.input.vcf, "Input VCF", errors, warnings)
    validate_fasta(snakemake.input.fasta, errors)
    validate_normalized_split_vcf(snakemake.input.vcf, snakemake.input.fasta, errors, warnings)
    write_report({"status": "failed" if errors else "ok", "errors": errors, "warnings": warnings, "summary": {"n_ped_samples": len(metadata["all_sample_ids"]), "n_processed_samples": len(metadata["sample_ids"]), "n_excluded_phenotype_samples": len(metadata["excluded_phenotype_sample_ids"]), "n_vcf_samples": len(vcf_samples), "n_affected_samples": metadata["affected_count"], "n_male_samples": metadata["male_count"], "sex_counts": metadata["sex_counts"], "phenotype_counts": metadata["phenotype_counts"], "vcf_chromosomes": vcf_chromosomes}}, snakemake.output.report)
    write_manifest(metadata["manifest"], snakemake.output.manifest)


if __name__ == "__main__":
    main()
