"""Validate resource contracts through injected Snakemake inputs and parameters."""

from types import SimpleNamespace

from _library import add_library_path

add_library_path(snakemake.scriptdir)

from validation import (
    parse_bool, read_vcf, validate_annotation_chromosome_names,
    validate_cadd_chromosome_names, validate_cadd_resources,
    validate_reference_tables, validate_region_chromosomes, write_report,
)


def main() -> None:
    """Validate reference-resource formats and chromosome compatibility.

    Converts Snakemake bindings into the shared validation contract, checks
    required indexes and annotation resources, and writes the resulting JSON
    report for the workflow validation gate.
    """
    args = SimpleNamespace(
        vcf=str(snakemake.input.vcf),
        coverage=str(snakemake.input.coverage[0]) if snakemake.input.coverage else "None",
        variants_exclude=",".join(map(str, snakemake.input.variants_exclude)) or "None",
        testing_unit_bed=str(snakemake.input.testing_unit_bed),
        par_regions=",".join(map(str, snakemake.input.par_regions)) or "None",
        cadd_snv=str(snakemake.input.cadd_snv),
        cadd_indel=str(snakemake.input.cadd_indel[0]) if snakemake.input.cadd_indel else "None",
        coverage_index=str(snakemake.params.coverage_index) or None,
        variant_exclude_indexes=list(snakemake.params.variant_exclude_indexes),
        testing_unit_index=str(snakemake.params.testing_unit_index) or None,
        par_region_indexes=list(snakemake.params.par_region_indexes),
        cadd_snv_index=str(snakemake.params.cadd_snv_index) or None,
        cadd_indel_index=str(snakemake.params.cadd_indel_index) or None,
        unit=str(snakemake.params.unit),
        include_indels=str(snakemake.params.include_indels),
    )
    errors: list[str] = []
    warnings: list[str] = []
    include_indels = parse_bool(args.include_indels, "include_indels", errors)
    reference_chromosomes = validate_reference_tables(args, errors, warnings)
    vcf_errors: list[str] = []
    _, vcf_chromosomes = read_vcf(args.vcf, vcf_errors)
    errors.extend(vcf_errors)
    validate_region_chromosomes(reference_chromosomes["exclusions"], vcf_chromosomes, "Exclusion interval files", warnings)
    validate_annotation_chromosome_names(reference_chromosomes["testing_unit"], vcf_chromosomes, f"Configured {args.unit} testing-unit BED file", warnings)
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
        validate_cadd_chromosome_names(chromosomes, vcf_chromosomes, f"CADD resource {filename}", errors, warnings)
    write_report({"status": "failed" if errors else "ok", "errors": errors, "warnings": warnings, "summary": {"include_indels": include_indels, "resource_vcf_chromosomes": vcf_chromosomes, "cadd_chromosomes": {filename: sorted(chromosomes) for filename, chromosomes in cadd_chromosomes.items()}}}, snakemake.output.report, create_parents=False)


if __name__ == "__main__":
    main()
