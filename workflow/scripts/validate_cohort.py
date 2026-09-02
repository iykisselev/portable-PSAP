"""Validate cohort inputs and write the manifest via injected Snakemake I/O."""

from _library import add_library_path

add_library_path(snakemake.scriptdir)

from validation_cohort import validate_cohort_inputs
from validation_reporting import write_validation_report


def main() -> None:
    """Validate cohort inputs and write the processed-sample manifest.

    Checks metadata, VCF, reference FASTA, indexes, genotype normalization, and
    hemizygote chromosome requirements supplied by Snakemake. Always writes the
    accumulated validation report and writes the manifest for downstream dynamic
    sample selection.
    """
    kind = str(snakemake.params.metadata_kind)
    vcf_index = str(snakemake.input.vcf_index[0]) if snakemake.input.vcf_index else None
    fasta_gzi = str(snakemake.input.fasta_gzi[0]) if snakemake.input.fasta_gzi else None
    report = validate_cohort_inputs(
        metadata_path=snakemake.input.metadata,
        metadata_kind=kind,
        vcf=snakemake.input.vcf,
        vcf_index=vcf_index,
        fasta=snakemake.input.fasta,
        fasta_index=snakemake.input.fasta_index,
        fasta_gzi=fasta_gzi,
        hemizygote_model=str(snakemake.params.hemizygote_model),
        manifest_path=snakemake.output.manifest,
    )
    write_validation_report(report, snakemake.output.report, create_parents=False)


if __name__ == "__main__":
    main()
