rule build_sample_candidates:
    """
    Construct patient-specific candidates for one affected individual.
    """
    input:
        psap=rules.apply_psap.output.psap,
        sample_manifest=validated_manifest_path,
        script=str(APPLY_SCRIPT_DIR / "build_sample_candidates.R"),
        libraries=R_RUNTIME_FILES
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
        sex=affected_sample_sex
    threads: 2
    resources:
        mem_mb=4096,
        runtime=60
    conda:
        str(CONDA_ENV_DIR / "r.yaml")
    envmodules: mod("r")
    container: sif("r")
    log:
        str(
            LOG_DIR
            / "apply"
            / config["run"]["apply_prefix"]
            / "cohort"
            / "build_sample_candidates_{sample}_{testing_unit}.log"
        )
    shell:
        """
        Rscript {input.script:q} \
            --input {input.psap:q} \
            --pid {params.pid:q} \
            --sex {params.sex:q} \
            --output {output.candidates:q} \
            --threads {threads} \
            > {log:q} 2>&1
        """


rule build_control_reference:
    """
    Retain model evidence for every carried variant in unaffected individuals.
    """
    input:
        validated=validation_marker,
        sample_manifest=validated_manifest_path,
        annotated=control_annotated_tables,
        psap=control_psap_tables,
        script=str(APPLY_SCRIPT_DIR / "build_control_reference.R"),
        libraries=R_RUNTIME_FILES
    output:
        control_variants=ensure(
            str(COHORT_DIR / "control_variants.tsv.gz"),
            non_empty=True,
        )
    params:
        pids=control_pids,
        sexes=control_sexes,
        hemizygote=config["model"]["hemizygote"],
        deduplicate=(
            "--deduplicate"
            if config.get("report", {}).get("deduplicate_controls", True)
            else ""
        )
    threads: 4
    resources:
        mem_mb=lambda wc, input, attempt: scaled_mem_mb(wc, input, attempt, 16000, 256),
        runtime=lambda wc, input, attempt: scaled_runtime_minutes(wc, input, attempt, 360, 20)
    conda:
        str(CONDA_ENV_DIR / "r.yaml")
    envmodules: mod("r")
    container: sif("r")
    benchmark:
        str(BENCHMARK_DIR / "build_control_reference.tsv")
    log:
        str(
            LOG_DIR
            / "apply"
            / config["run"]["apply_prefix"]
            / "cohort"
            / "build_control_reference.log"
        )
    shell:
        """
        Rscript {input.script:q} \
            --input {input.annotated:q} \
            --psap-input {input.psap:q} \
            --pid {params.pids:q} \
            --sex {params.sexes:q} \
            --hemizygote {params.hemizygote} \
            {params.deduplicate} \
            --output {output.control_variants:q} \
            --threads {threads} \
            > {log:q} 2>&1
        """


rule build_cohort_report:
    """
    Validate affected candidates against controls and summarize the cohort.
    """
    input:
        candidates=affected_candidate_tables,
        controls=rules.build_control_reference.output.control_variants,
        script=str(APPLY_SCRIPT_DIR / "build_cohort_report.R"),
        libraries=R_RUNTIME_FILES
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
        mem_mb=lambda wc, input, attempt: scaled_mem_mb(wc, input, attempt, 16000, 256),
        runtime=lambda wc, input, attempt: scaled_runtime_minutes(wc, input, attempt, 360, 20)
    conda:
        str(CONDA_ENV_DIR / "r.yaml")
    envmodules: mod("r")
    container: sif("r")
    benchmark:
        str(BENCHMARK_DIR / "build_cohort_report.tsv")
    log:
        str(
            LOG_DIR
            / "apply"
            / config["run"]["apply_prefix"]
            / "cohort"
            / "build_cohort_report.log"
        )
    shell:
        """
        Rscript {input.script:q} \
            --candidates {input.candidates:q} \
            --controls {input.controls:q} \
            --output-candidates {output.candidates:q} \
            --output-summary {output.summary:q} \
            --threads {threads} \
            > {log:q} 2>&1
        """
