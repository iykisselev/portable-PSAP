"""Validate model resources through injected Snakemake inputs and parameters."""

from types import SimpleNamespace

from _library import add_library_path

add_library_path(snakemake.scriptdir)

from bundle_contract import is_supported_contract
from validation_model import validate_model_inputs
from validation_reporting import write_validation_report


def main() -> None:
    """Validate model settings and lookup tables against testing-unit annotations.

    Converts Snakemake settings into the validation library's argument contract,
    checks every enabled lookup table, and writes a JSON report used by the
    shared validation gate.
    """
    args = SimpleNamespace(
        score_max=str(snakemake.params.score_max),
        assembly=str(snakemake.params.assembly),
        unit=str(snakemake.params.unit),
        hemizygote_model=str(snakemake.params.hemizygote_model),
        include_indels=str(snakemake.params.include_indels),
    )
    configured = {
        "het": str(snakemake.input.lookup_het),
        "hom": str(snakemake.input.lookup_hom),
        "hem_chrX": str(snakemake.input.lookup_hem_chrx[0]) if snakemake.input.lookup_hem_chrx else "",
    }
    contract = str(snakemake.params.lookup_contract)
    certificate = None
    components: dict[str, str] = {}
    if is_supported_contract(contract):
        certificate = {
            "metadata": str(snakemake.input.lookup_metadata[0]),
            "validation": str(snakemake.input.lookup_validation[0]),
            "marker": str(snakemake.input.lookup_marker[0]),
        }
        components = {
            "lookup_het": configured["het"],
            "lookup_hom": configured["hom"],
            "eligible_testing_units": str(snakemake.input.lookup_eligible_units[0]),
            "calibration_regions": str(snakemake.input.lookup_calibration_mask[0]),
            "excluded_alleles": str(snakemake.input.lookup_excluded_alleles[0]),
            "excluded_units": str(snakemake.input.lookup_excluded_units[0]),
        }
        if str(snakemake.params.hemizygote_model).lower() == "true":
            components["lookup_hem_chrX"] = configured["hem_chrX"]
    report = validate_model_inputs(
        args=args,
        testing_unit_bed=snakemake.input.testing_unit_bed,
        lookup_tables=configured,
        lookup_contract=contract,
        lookup_bin_width=str(snakemake.params.lookup_bin_width),
        bundle_certificate=certificate,
        bundle_components=components,
    )
    write_validation_report(report, snakemake.output.report, create_parents=False)


if __name__ == "__main__":
    main()
