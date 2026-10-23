rule build_sample_candidates:
    """
    Construct patient-specific candidates for one affected individual.
    """
    input:
        psap=rules.apply_psap.output.psap,
        sample_manifest=str(SAMPLE_MANIFEST_PATH),
    output:
        candidates=ensure(
            str(CANDIDATE_DIR / f"{{sample}}.{TESTING_UNIT}.tsv.gz"),
            non_empty=True,
        ),
    log:
        str(
            APPLY_LOG_DIR
            / "cohort"
            / f"build_sample_candidates_{{sample}}_{TESTING_UNIT}.log"
        ),
    wildcard_constraints:
        sample=SAMPLE_PATTERN,
    conda:
        env("r")
    container:
        sif("r")
    envmodules:
        *mod("r"),
    threads: 2
    resources:
        mem_mb=policy_resources("build_sample_candidates").mem_mb,
        runtime=policy_resources("build_sample_candidates").runtime,
    params:
        pid=lambda wc: wc.sample,
        sex=affected_sample_sex,
        candidate_id_separator=WORKFLOW_PLAN.report.candidate_id_separator,
        candidate_id_ordinal_width=WORKFLOW_PLAN.report.candidate_id_ordinal_width,
    script:
        str(APPLY_SCRIPT_DIR / "build_sample_candidates.R")


rule build_control_reference:
    """
    Retain model evidence for every carried variant in unaffected individuals.
    """
    input:
        validated=rules.require_valid_inputs.output.ok,
        sample_manifest=str(SAMPLE_MANIFEST_PATH),
        annotated=per_sample(ANNOTATED_DIR, "unaffected"),
        psap=per_sample(SCORED_DIR, "unaffected"),
    output:
        control_variants=ensure(
            COHORT_OUTPUTS["control_variants"],
            non_empty=True,
        ),
    log:
        str(APPLY_LOG_DIR / "cohort" / "build_control_reference.log"),
    benchmark:
        str(BENCHMARK_DIR / "build_control_reference.tsv")
    conda:
        env("r")
    container:
        sif("r")
    envmodules:
        *mod("r"),
    threads: 4
    resources:
        mem_mb=policy_resources("build_control_reference").mem_mb,
        runtime=policy_resources("build_control_reference").runtime,
    params:
        pids=control_pids,
        sexes=control_sexes,
        hemizygote=WORKFLOW_PLAN.model.hemizygote,
        deduplicate=WORKFLOW_PLAN.report.deduplicate_controls,
    script:
        str(APPLY_SCRIPT_DIR / "build_control_reference.R")


rule build_cohort_report:
    """
    Validate affected candidates against controls and summarize the cohort.
    """
    input:
        candidates=per_sample(CANDIDATE_DIR, "affected"),
        evidence=per_sample(ANNOTATED_DIR, "affected"),
        controls=rules.build_control_reference.output.control_variants,
        partition_statuses=LOOKUP_SOURCE.partition_statuses,
    output:
        candidates=ensure(
            COHORT_OUTPUTS["cohort_candidates"],
            non_empty=True,
        ),
        summary=ensure(
            COHORT_OUTPUTS["cohort_summary"],
            non_empty=True,
        ),
        report_metadata=ensure(REPORT_METADATA, non_empty=True),
        aggregation_diagnostics=ensure(AGGREGATION_DIAGNOSTICS, non_empty=True),
        partition_summaries=[
            ensure(path, non_empty=True) for path in PARTITION_SUMMARIES.values()
        ],
    log:
        str(APPLY_LOG_DIR / "cohort" / "build_cohort_report.log"),
    benchmark:
        str(BENCHMARK_DIR / "build_cohort_report.tsv")
    conda:
        env("r")
    container:
        sif("r")
    envmodules:
        *mod("r"),
    threads: 4
    resources:
        mem_mb=policy_resources("build_cohort_report").mem_mb,
        runtime=policy_resources("build_cohort_report").runtime,
    params:
        partition_ids=LOOKUP_SOURCE.partition_ids,
        affected_pids=affected_pids,
        affected_sexes=affected_sexes,
    script:
        str(APPLY_SCRIPT_DIR / "build_cohort_report.R")


rule build_endpoint_diagnostics:
    """Measure discrete lookup endpoints without loading all lookup tables."""
    input:
        validated=rules.require_valid_inputs.output.ok,
        lookup_tables=list(LOOKUP_SOURCE.tables.values()),
    output:
        diagnostics=ensure(ENDPOINT_DIAGNOSTICS, non_empty=True),
    log:
        str(APPLY_LOG_DIR / "cohort" / "build_endpoint_diagnostics.log"),
    conda:
        env("validation")
    container:
        sif("validation")
    envmodules:
        *mod("validation"),
    threads: 1
    resources:
        mem_mb=policy_resources("build_endpoint_diagnostics").mem_mb,
        runtime=policy_resources("build_endpoint_diagnostics").runtime,
    params:
        lookup_keys=tuple(
            f"{partition_id}\t{model_id}"
            for partition_id, model_id in LOOKUP_SOURCE.tables
        ),
        near_one_thresholds=(0.999, 0.999999, 0.999999999),
    script:
        str(APPLY_SCRIPT_DIR / "build_endpoint_diagnostics.py")
