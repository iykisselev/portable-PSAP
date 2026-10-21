rule validate_cohort:
    """Validate cohort VCF/metadata/FASTA inputs and write the sample manifest."""
    input:
        vcf=apply_resource("cohort").source,
        vcf_index=as_inputs(apply_resource("cohort").source_index),
        metadata=INPUT_BINDINGS.samples or INPUT_BINDINGS.ped or "",
        fasta=apply_resource("fasta").source,
        fasta_index=apply_resource("fasta").source_index,
        fasta_gzi=as_inputs(
            apply_resource("fasta").companion_indexes[0]
            if apply_resource("fasta").companion_indexes
            else ""
        ),
        libraries=VALIDATION_LIBRARIES,
        schema=VALIDATION_REPORT_SCHEMA,
    output:
        report=ensure(str(COHORT_VALIDATION_REPORT), non_empty=True),
        manifest=ensure(str(SAMPLE_MANIFEST_PATH), non_empty=True),
    log:
        str(APPLY_LOG_DIR / "validation" / "cohort.log"),
    conda:
        env("validation")
    container:
        sif("validation")
    envmodules:
        *mod("validation"),
    threads: 1
    resources:
        mem_mb=policy_resources("validate_cohort").mem_mb,
        runtime=policy_resources("validate_cohort").runtime,
    params:
        hemizygote_model=WORKFLOW_PLAN.model.hemizygote,
        metadata_kind="samples" if INPUT_BINDINGS.samples else "ped",
    script:
        str(APPLY_SCRIPT_DIR / "validate_cohort.py")


rule validate_resources:
    """Validate interval, exclusion, testing-unit, and CADD resources."""
    input:
        **APPLY_RESOURCE_REGISTRY.validation_inputs(),
        vcf=apply_resource("cohort").source,
        libraries=VALIDATION_LIBRARIES,
        schema=VALIDATION_REPORT_SCHEMA,
    output:
        report=ensure(str(RESOURCE_VALIDATION_REPORT), non_empty=True),
    log:
        str(APPLY_LOG_DIR / "validation" / "resources.log"),
    conda:
        env("validation")
    container:
        sif("validation")
    envmodules:
        *mod("validation"),
    threads: 1
    resources:
        mem_mb=policy_resources("validate_resources").mem_mb,
        runtime=policy_resources("validate_resources").runtime,
    params:
        unit=WORKFLOW_PLAN.model.unit,
        include_indels=WORKFLOW_PLAN.model.include_indels,
        resource_fields=APPLY_RESOURCE_REGISTRY.validation_fields(),
    script:
        str(APPLY_SCRIPT_DIR / "validate_resources.py")


rule validate_model:
    """Validate model configuration and enabled lookup resources."""
    input:
        lookup_tables=list(LOOKUP_SOURCE.tables.values()),
        lookup_metadata=as_inputs(LOOKUP_SOURCE.metadata),
        lookup_validation=as_inputs(LOOKUP_SOURCE.validation),
        lookup_marker=as_inputs(LOOKUP_SOURCE.marker),
        lookup_calibration_mask=(
            as_inputs(LOOKUP_SOURCE.calibration_mask)
            if LOOKUP_SOURCE.is_bundle
            else []
        ),
        lookup_excluded_alleles=(
            list(LOOKUP_SOURCE.excluded_alleles) if LOOKUP_SOURCE.is_bundle else []
        ),
        resolved_partitions=list(LOOKUP_SOURCE.resolved_partitions.values()),
        eligible_partitions=list(LOOKUP_SOURCE.eligible_partitions.values()),
        partition_statuses=str(Path(LOOKUP_SOURCE.metadata).parent / "unit_partition_status.tsv.gz"),
        testing_unit_id_overlap=str(Path(LOOKUP_SOURCE.metadata).parent / "testing_unit_id_overlap.tsv.gz"),
        libraries=VALIDATION_LIBRARIES,
        schema=VALIDATION_REPORT_SCHEMA,
    output:
        report=ensure(str(MODEL_VALIDATION_REPORT), non_empty=True),
    log:
        str(APPLY_LOG_DIR / "validation" / "model.log"),
    conda:
        env("validation")
    container:
        sif("validation")
    envmodules:
        *mod("validation"),
    threads: 1
    resources:
        mem_mb=policy_resources("validate_model").mem_mb,
        runtime=policy_resources("validate_model").runtime,
    params:
        score_max=LOOKUP_SOURCE.score_max,
        assembly=WORKFLOW_PLAN.reference.assembly,
        unit=WORKFLOW_PLAN.model.unit,
        hemizygote_model=WORKFLOW_PLAN.model.hemizygote,
        include_indels=WORKFLOW_PLAN.model.include_indels,
        lookup_contract=LOOKUP_SOURCE.contract,
        lookup_bin_width=LOOKUP_SOURCE.bin_width,
        lookup_component_keys=tuple(f"lookup.{partition_id}.{model_id}" for partition_id, model_id in LOOKUP_SOURCE.tables),
    script:
        str(APPLY_SCRIPT_DIR / "validate_model.py")


rule merge_validation:
    """Merge independent validation fragments into validation.json."""
    input:
        cohort=rules.validate_cohort.output.report,
        resources=rules.validate_resources.output.report,
        model=rules.validate_model.output.report,
    output:
        report=protected(ensure(str(COMBINED_VALIDATION_REPORT), non_empty=True)),
    log:
        str(APPLY_LOG_DIR / "validation" / "merge.log"),
    conda:
        env("validation")
    container:
        sif("validation")
    envmodules:
        *mod("validation"),
    threads: 1
    resources:
        mem_mb=policy_resources("merge_validation").mem_mb,
        runtime=policy_resources("merge_validation").runtime,
    script:
        str(APPLY_SCRIPT_DIR / "merge_validation.py")


rule require_valid_inputs:
    """Gate expensive workflow rules on the merged validation report."""
    input:
        report=rules.merge_validation.output.report,
        manifest=rules.validate_cohort.output.manifest,
    output:
        ok=ensure(str(Path(APPLY_OUTPUT_DIR) / "validation.ok"), non_empty=True),
    log:
        str(APPLY_LOG_DIR / "validation" / "require_valid_inputs.log"),
    conda:
        env("validation")
    container:
        sif("validation")
    envmodules:
        *mod("validation"),
    threads: 1
    resources:
        mem_mb=policy_resources("require_valid_inputs").mem_mb,
        runtime=policy_resources("require_valid_inputs").runtime,
    params:
        failure_message="Input validation failed",
        include_errors=True,
        success_marker="validated\n",
    script:
        str(APPLY_SCRIPT_DIR / "require_valid_report.py")
