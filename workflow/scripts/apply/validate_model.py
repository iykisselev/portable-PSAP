"""Validate model resources through injected Snakemake inputs and parameters."""

from pathlib import Path
import sys

sys.path.insert(0, str(Path(snakemake.scriptdir).parents[1] / "shared"))
from _bootstrap import bootstrap

bootstrap(snakemake)
from apply_bindings import (
    ModelValidationInputs,
    ReportOutput,
    ModelValidationParameters,
)
from _library import bind


from bundle_contract import is_supported_contract
from validation_model import (
    certificate_components_with_exclusions,
    validate_model_inputs,
)
from validation_reporting import write_validation_report


def main() -> None:
    """Validate model settings and lookup tables against testing-unit annotations.

    Converts Snakemake settings into the validation library's argument contract,
    checks every enabled lookup table, and writes a JSON report used by the
    shared validation gate.
    """
    inputs = bind(ModelValidationInputs, snakemake.input)
    outputs = bind(ReportOutput, snakemake.output)
    params = bind(ModelValidationParameters, snakemake.params)
    configured = {
        "het": str(inputs.lookup_het),
        "hom": str(inputs.lookup_hom),
        "hem_chrX": inputs.lookup_hem_chrX or "",
    }
    contract = str(params.lookup_contract)
    certificate = None
    components: dict[str, str] = {}
    if is_supported_contract(contract):
        certificate = {
            "metadata": str(inputs.lookup_metadata),
            "validation": str(inputs.lookup_validation),
            "marker": str(inputs.lookup_marker),
        }
        components = certificate_components_with_exclusions(
            {
                "lookup_het": configured["het"],
                "lookup_hom": configured["hom"],
                "eligible_testing_units": str(inputs.lookup_eligible_units),
                "calibration_regions": str(inputs.lookup_calibration_mask),
                "excluded_units": str(inputs.lookup_excluded_units),
            },
            inputs.lookup_excluded_alleles,
        )
        if params.hemizygote_model:
            components["lookup_hem_chrX"] = configured["hem_chrX"]
    report = validate_model_inputs(
        score_max=params.score_max,
        assembly=str(params.assembly),
        unit=str(params.unit),
        hemizygote_model=params.hemizygote_model,
        include_indels=params.include_indels,
        testing_unit_bed=inputs.testing_unit_bed,
        lookup_tables=configured,
        lookup_contract=contract,
        lookup_bin_width=params.lookup_bin_width,
        bundle_certificate=certificate,
        bundle_components=components,
    )
    write_validation_report(report, outputs.report, create_parents=False)


if __name__ == "__main__":
    main()
