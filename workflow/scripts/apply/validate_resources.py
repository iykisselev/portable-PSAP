"""Validate resource contracts through injected Snakemake inputs and parameters."""

from pathlib import Path
import sys

SHARED_SCRIPT_DIR = Path(snakemake.scriptdir).parents[1] / "shared"
sys.path.insert(0, str(SHARED_SCRIPT_DIR))
from _library import add_library_path

add_library_path(snakemake.scriptdir)

from validation import IndexedResource
from validation_reporting import write_validation_report
from validation_resources import validate_resource_inputs


def indexed_resource(binding) -> IndexedResource:
    """Build one typed binding without inferring an index by position."""
    path, index = binding
    return IndexedResource(path=str(path), index=str(index) if index else None)


def main() -> None:
    """Validate reference-resource formats and chromosome compatibility.

    Converts Snakemake bindings into the shared validation contract, checks
    required indexes and annotation resources, and writes the resulting JSON
    report for the workflow validation gate.
    """
    calibration_binding = snakemake.params.calibration_regions
    cadd_indel_binding = snakemake.params.cadd_indel
    write_validation_report(
        validate_resource_inputs(
            vcf=str(snakemake.input.vcf),
            calibration_regions=(
                indexed_resource(calibration_binding)
                if calibration_binding
                else None
            ),
            exclusions=[
                indexed_resource(binding) for binding in snakemake.params.exclusions
            ],
            testing_units=indexed_resource(snakemake.params.testing_units),
            par_regions=[
                indexed_resource(binding) for binding in snakemake.params.par_regions
            ],
            cadd_snv=indexed_resource(snakemake.params.cadd_snv),
            cadd_indel=(
                indexed_resource(cadd_indel_binding) if cadd_indel_binding else None
            ),
            unit=str(snakemake.params.unit),
            include_indels=snakemake.params.include_indels,
        ),
        snakemake.output.report,
        create_parents=False,
    )


if __name__ == "__main__":
    main()
