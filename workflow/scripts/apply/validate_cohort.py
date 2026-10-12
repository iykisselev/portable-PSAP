"""Validate cohort inputs and write the sample manifest."""

from pathlib import Path
import sys

sys.path.insert(0, str(Path(snakemake.scriptdir).parents[1] / "shared"))
from _bootstrap import bootstrap

bootstrap(snakemake)
from apply_bindings import (
    CohortValidationInputs,
    CohortValidationOutputs,
    CohortValidationParameters,
)
from _library import bind


from validation_cohort import validate_cohort_inputs
from validation_reporting import write_validation_report


def main() -> None:
    """Validate cohort inputs and write the processed-sample manifest.

    Checks metadata, VCF, reference FASTA, indexes, genotype normalization, and
    hemizygote chromosome requirements supplied by Snakemake. Always writes the
    accumulated validation report and writes the manifest for downstream dynamic
    sample selection.
    """
    inputs = bind(CohortValidationInputs, snakemake.input)
    outputs = bind(CohortValidationOutputs, snakemake.output)
    params = bind(CohortValidationParameters, snakemake.params)
    kind = str(params.metadata_kind)
    report = validate_cohort_inputs(
        metadata_path=inputs.metadata,
        metadata_kind=kind,
        vcf=inputs.vcf,
        vcf_index=inputs.vcf_index,
        fasta=inputs.fasta,
        fasta_index=inputs.fasta_index,
        fasta_gzi=inputs.fasta_gzi,
        hemizygote_model=params.hemizygote_model,
        manifest_path=outputs.manifest,
    )
    write_validation_report(report, outputs.report, create_parents=False)


if __name__ == "__main__":
    main()
