"""Zero-policy event calculation; mutation conversion remains deliberately absent."""

GENERATION_CALCULATED_DIR = GENERATION_OUTDIR / "calculated"


rule calculate_zero_null_distribution:
    """Calculate HET/HOM distributions from observed AF plus explicit zero absence."""
    input:
        observed_af=rules.prepare_generation_inputs.output.observed_af,
        score_universe=rules.prepare_generation_inputs.output.score_universe,
        units=rules.prepare_generation_inputs.output.units
    output:
        lookup_het=ensure(str(GENERATION_CALCULATED_DIR / "lookup_het.long.tsv"), non_empty=True),
        lookup_hom=ensure(str(GENERATION_CALCULATED_DIR / "lookup_hom.long.tsv"), non_empty=True),
        diagnostics=ensure(str(GENERATION_CALCULATED_DIR / "distribution_diagnostics.tsv"), non_empty=True)
    params:
        score_max=config["model"]["score_max"],
        bin_count=config["generation_settings"]["bin_count"]
    threads: 1
    resources:
        mem_mb=2048,
        runtime=60
    conda:
        str(CONDA_ENV_DIR / "r.yaml")
    envmodules: *mod("r")
    container: sif("r")
    log:
        str(GENERATION_LOG_DIR / "calculate_zero.log")
    script:
        str(Path(workflow.basedir) / "scripts" / "generate" / "calculate_null_distribution.R")
