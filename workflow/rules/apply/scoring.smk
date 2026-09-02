LOOKUP_TABLES = {
    name: path
    for name, path in LOOKUP_SOURCE.tables.items()
}
HEMIZYGOUS_MODEL_ENABLED = WORKFLOW_PLAN.model.hemizygote


rule apply_psap:
    """
    Calculate per-sample PSAP results
    """
    input:
        annotated_gt=rules.join_sample_annotations.output.annotated_gt,
        lookup_het=LOOKUP_TABLES["het"],
        lookup_hom=LOOKUP_TABLES["hom"],
        lookup_hem_chrx=(
            [LOOKUP_TABLES["hem_chrX"]] if HEMIZYGOUS_MODEL_ENABLED else []
        ),
        lookup_metadata=[LOOKUP_SOURCE.metadata] if LOOKUP_SOURCE.metadata else [],
        lookup_validation=[LOOKUP_SOURCE.validation] if LOOKUP_SOURCE.validation else [],
        lookup_marker=[LOOKUP_SOURCE.marker] if LOOKUP_SOURCE.marker else [],
        eligible_units=[LOOKUP_SOURCE.eligible_units] if LOOKUP_SOURCE.eligible_units else [],
        calibration_mask=[LOOKUP_SOURCE.calibration_mask] if LOOKUP_SOURCE.calibration_mask else [],
        excluded_alleles=list(LOOKUP_SOURCE.excluded_alleles),
        excluded_units=[LOOKUP_SOURCE.excluded_units] if LOOKUP_SOURCE.excluded_units else []
    output:
        psap=ensure(
            str(SCORED_DIR / "{sample}.{testing_unit}.tsv.gz"),
            non_empty=True,
        )

    params:
        score_max=str(WORKFLOW_PLAN.model.score_max),
        hemizygote=WORKFLOW_PLAN.model.hemizygote,
        lookup_contract=LOOKUP_SOURCE.contract,
        lookup_bin_width=LOOKUP_SOURCE.bin_width or ""
    threads: 4
    resources:
        mem_mb=16000,
        runtime=360
    conda:
        str(CONDA_ENV_DIR / "r.yaml")
    envmodules: *mod("r")
    container: sif("r")
    benchmark:
        str(BENCHMARK_DIR / "apply_psap_{sample}_{testing_unit}.tsv")
    log:
        str(
            APPLY_LOG_DIR / "samples" /
            "apply_psap_{sample}_{testing_unit}.log"
        )
    script:
        str(APPLY_SCRIPT_DIR / "apply_psap.R")
