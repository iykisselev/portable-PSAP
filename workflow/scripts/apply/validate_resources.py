"""Validate resource contracts through injected Snakemake inputs and parameters."""

from types import SimpleNamespace
from pathlib import Path
import sys

SHARED_SCRIPT_DIR = Path(snakemake.scriptdir).parents[1] / "shared"
sys.path.insert(0, str(SHARED_SCRIPT_DIR))
from _library import add_library_path

add_library_path(snakemake.scriptdir)

from validation_reporting import write_validation_report
from validation_resources import validate_resource_inputs


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
    write_validation_report(
        validate_resource_inputs(args),
        snakemake.output.report,
        create_parents=False,
    )


if __name__ == "__main__":
    main()
