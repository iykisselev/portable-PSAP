"""Snakemake adapter for generated-bundle certification."""

from pathlib import Path
import sys

sys.path.insert(0, str(Path(snakemake.scriptdir).parents[1] / "shared"))
from _bootstrap import bootstrap

bootstrap(snakemake)

from _library import bind

from bundle_validation import validate_generated_bundle
from generation_bindings import BundleValidationInputs, BundleValidationOutputs


def main() -> None:
    """Decode Snakemake bindings and delegate semantic bundle validation."""
    validate_generated_bundle(
        bind(BundleValidationInputs, snakemake.input),
        bind(BundleValidationOutputs, snakemake.output),
    )


if __name__ == "__main__":
    main()
