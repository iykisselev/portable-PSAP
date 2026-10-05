"""Snakemake adapter for generated-bundle finalization."""

import os
from pathlib import Path
import sys

if "snakemake" in globals():
    os.environ.setdefault("POLARS_MAX_THREADS", str(snakemake.threads))

SCRIPT_DIR = Path(snakemake.scriptdir) if "snakemake" in globals() else Path(__file__).resolve().parent
SHARED_SCRIPT_DIR = SCRIPT_DIR.parents[1] / "shared"
sys.path.insert(0, str(SHARED_SCRIPT_DIR))
from _library import add_library_path

add_library_path(SCRIPT_DIR)

from bundle_finalization import finalize_generated_bundle
from generation_bindings import (
    BundleFinalizationInputs,
    BundleFinalizationOutputs,
    BundleFinalizationParameters,
)


def _paths(values) -> list[str]:
    return [str(value) for value in values]


def _optional_path(values) -> str | None:
    if not values:
        return None
    if isinstance(values, (str, Path)):
        return str(values)
    return str(values[0])


def main() -> None:
    """Decode Snakemake bindings and publish the bundle through the library."""
    hem_output = getattr(snakemake.output, "lookup_hem", [])
    finalize_generated_bundle(
        BundleFinalizationInputs(
            lookup_het=tuple(_paths(snakemake.input.lookup_het)),
            lookup_hom=tuple(_paths(snakemake.input.lookup_hom)),
            lookup_hem=tuple(_paths(snakemake.input.lookup_hem)),
            diagnostics=tuple(_paths(snakemake.input.diagnostics)),
            units=str(snakemake.input.units),
            masked_unit_intervals=str(snakemake.input.masked_unit_intervals),
            allele_states=tuple(_paths(snakemake.input.allele_states)),
            effective_mask=str(snakemake.input.effective_mask),
            resource_identities=str(snakemake.input.resource_identities),
        ),
        BundleFinalizationOutputs(
            lookup_het=str(snakemake.output.lookup_het),
            lookup_hom=str(snakemake.output.lookup_hom),
            lookup_hem=_optional_path(hem_output),
            eligible_bed=str(snakemake.output.eligible_bed),
            calibration_regions=str(snakemake.output.calibration_regions),
            excluded_alleles=str(snakemake.output.excluded_alleles),
            population_af_zero_imputations=str(snakemake.output.population_af_zero_imputations),
            excluded_units=str(snakemake.output.excluded_units),
            metadata=str(snakemake.output.metadata),
        ),
        BundleFinalizationParameters(
            bin_count=int(snakemake.params.bin_count),
            minimum_retained_fraction=float(snakemake.params.minimum_retained_fraction),
            bin_width=float(snakemake.params.bin_width),
            score_max=float(snakemake.params.score_max),
            policy=str(snakemake.params.policy),
            conversion=str(snakemake.params.conversion),
            assembly=str(snakemake.params.assembly),
            testing_unit_type=str(snakemake.params.testing_unit_type),
            af_tags=dict(snakemake.params.af_tags),
            x_af_policy=str(snakemake.params.x_af_policy),
            x_mutation_rate_policy=str(snakemake.params.x_mutation_rate_policy),
            hemizygote=snakemake.params.hemizygote,
            variant_classes=tuple(snakemake.params.variant_classes),
            resource_ids=dict(snakemake.params.resource_ids),
        ),
    )


if __name__ == "__main__":
    main()
