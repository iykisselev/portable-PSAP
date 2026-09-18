"""Bundle finalization with explicit zero-event exclusions."""

import math

GENERATION_BUNDLE_DIR = GENERATION_OUTDIR / "bundle"

rule finalize_null_distribution_bundle:
    input:
        lookup_het=rules.calculate_null_distribution.output.lookup_het,
        lookup_hom=rules.calculate_null_distribution.output.lookup_hom,
        diagnostics=rules.calculate_null_distribution.output.diagnostics,
        units=rules.prepare_generation_inputs.output.units,
        testing_units=config["generation"]["testing_unit_bed"]["data"],
        coverage=config["generation"]["coverage"]["data"],
        par_regions=config["generation"]["par_regions"]["data"]
    output:
        lookup_het=ensure(str(GENERATION_BUNDLE_DIR / "lookup_het.tsv.gz"), non_empty=True),
        lookup_hom=ensure(str(GENERATION_BUNDLE_DIR / "lookup_hom.tsv.gz"), non_empty=True),
        eligible_bed=ensure(str(GENERATION_BUNDLE_DIR / "eligible_testing_units.bed.gz"), non_empty=True),
        calibration_regions=ensure(str(GENERATION_BUNDLE_DIR / "calibration_regions.bed.gz"), non_empty=True),
        excluded_units=ensure(str(GENERATION_BUNDLE_DIR / "excluded_units.tsv.gz"), non_empty=True),
        metadata=ensure(str(GENERATION_BUNDLE_DIR / "metadata.json"), non_empty=True),
        validation=ensure(str(GENERATION_BUNDLE_DIR / "validation.json"), non_empty=True),
        marker=ensure(str(GENERATION_BUNDLE_DIR / "validation.ok"), non_empty=True)
    params:
        bin_count=lambda wc: int(math.ceil(config["model"]["score_max"] / config["generation_settings"]["score_bin_width"])),
        bin_width=config["generation_settings"]["score_bin_width"],
        score_max=config["model"]["score_max"],
        policy=config["generation"]["unobserved_policy"],
        conversion=config["generation_settings"]["mutation_conversion"]
    conda:
        str(CONDA_ENV_DIR / "validation.yaml")
    envmodules: *mod("validation")
    container: sif("validation")
    log:
        str(GENERATION_LOG_DIR / "finalize.log")
    script:
        str(Path(workflow.basedir) / "scripts" / "generate" / "finalize_bundle.py")
