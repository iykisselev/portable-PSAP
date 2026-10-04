"""Snakemake adapter for generated-bundle certification."""

from pathlib import Path
import sys

SHARED_SCRIPT_DIR = Path(snakemake.scriptdir).parents[1] / "shared"
sys.path.insert(0, str(SHARED_SCRIPT_DIR))
from _library import add_library_path

add_library_path(snakemake.scriptdir)

from bundle_validation import validate_generated_bundle
from generation_bindings import BundleValidationInputs, BundleValidationOutputs


def _optional_path(values) -> str | None:
    if not values:
        return None
    if isinstance(values, (str, Path)):
        return str(values)
    return str(values[0])


def main() -> None:
    """Decode Snakemake bindings and delegate semantic bundle validation."""
    validate_generated_bundle(
        BundleValidationInputs(
            lookup_het=str(snakemake.input.lookup_het),
            lookup_hom=str(snakemake.input.lookup_hom),
            lookup_hem=_optional_path(snakemake.input.lookup_hem),
            eligible_bed=str(snakemake.input.eligible_bed),
            calibration_regions=str(snakemake.input.calibration_regions),
            excluded_alleles=str(snakemake.input.excluded_alleles),
            excluded_units=str(snakemake.input.excluded_units),
            metadata=str(snakemake.input.metadata),
            schema=str(snakemake.input.schema),
        ),
        BundleValidationOutputs(
            validation=str(snakemake.output.validation),
            marker=str(snakemake.output.marker),
        ),
    )


if __name__ == "__main__":
    main()
