rule build_sample_candidates:
    """
    Construct patient-specific candidates for one affected individual.
    """
    input:
        psap=rules.apply_psap.output.psap,
        sample_manifest=validated_manifest_path
    output:
        candidates=ensure(
            str(CANDIDATE_DIR / "{sample}.{testing_unit}.tsv.gz"),
            non_empty=True,
        )
    wildcard_constraints:
        sample=SAMPLE_PATTERN,
        testing_unit=TESTING_UNIT_PATTERN
    params:
        pid=lambda wc: wc.sample,
        sex=affected_sample_sex,
        candidate_id_separator=WORKFLOW_PLAN.report.candidate_id_separator,
        candidate_id_ordinal_width=WORKFLOW_PLAN.report.candidate_id_ordinal_width
    threads: 2
    resources:
        mem_mb=4096,
        runtime=60
    conda:
        str(CONDA_ENV_DIR / "r.yaml")
    envmodules: *mod("r")
    container: sif("r")
    log:
        str(
            APPLY_LOG_DIR
            / "cohort"
            / "build_sample_candidates_{sample}_{testing_unit}.log"
        )
    script:
        str(APPLY_SCRIPT_DIR / "build_sample_candidates.R")


rule build_control_reference:
    """
    Retain model evidence for every carried variant in unaffected individuals.
    """
    input:
        validated=validation_marker,
        sample_manifest=validated_manifest_path,
        annotated=control_annotated_tables,
        psap=control_psap_tables
    output:
        control_variants=ensure(
            str(COHORT_DIR / "control_variants.tsv.gz"),
            non_empty=True,
        )
    params:
        pids=control_pids,
        sexes=control_sexes,
        hemizygote=WORKFLOW_PLAN.model.hemizygote,
        deduplicate=WORKFLOW_PLAN.report.deduplicate_controls
    threads: 4
    resources:
        mem_mb=lambda wc, input, attempt: execution_mem_mb(
            wc, input, attempt, "cohort_aggregate"
        ),
        runtime=lambda wc, input, attempt: execution_runtime_minutes(
            wc, input, attempt, "cohort_aggregate"
        )
    conda:
        str(CONDA_ENV_DIR / "r.yaml")
    envmodules: *mod("r")
    container: sif("r")
    benchmark:
        str(BENCHMARK_DIR / "build_control_reference.tsv")
    log:
        str(
            APPLY_LOG_DIR
            / "cohort"
            / "build_control_reference.log"
        )
    script:
        str(APPLY_SCRIPT_DIR / "build_control_reference.R")


rule build_cohort_report:
    """
    Validate affected candidates against controls and summarize the cohort.
    """
    input:
        candidates=affected_candidate_tables,
        controls=rules.build_control_reference.output.control_variants
    output:
        candidates=ensure(
            str(COHORT_DIR / "cohort_candidates.tsv.gz"),
            non_empty=True,
        ),
        summary=ensure(
            str(COHORT_DIR / "cohort_summary.tsv.gz"),
            non_empty=True,
        )
    threads: 4
    resources:
        mem_mb=lambda wc, input, attempt: execution_mem_mb(
            wc, input, attempt, "cohort_aggregate"
        ),
        runtime=lambda wc, input, attempt: execution_runtime_minutes(
            wc, input, attempt, "cohort_aggregate"
        )
    conda:
        str(CONDA_ENV_DIR / "r.yaml")
    envmodules: *mod("r")
    container: sif("r")
    benchmark:
        str(BENCHMARK_DIR / "build_cohort_report.tsv")
    log:
        str(
            APPLY_LOG_DIR
            / "cohort"
            / "build_cohort_report.log"
        )
    script:
        str(APPLY_SCRIPT_DIR / "build_cohort_report.R")
