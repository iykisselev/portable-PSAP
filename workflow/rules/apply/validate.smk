rule validate_cohort:
    """Validate cohort VCF/metadata/FASTA inputs and write the sample manifest."""
    input:
        vcf=apply_resource("cohort").source,
        vcf_index=optional_input(apply_resource("cohort").source_index),
        metadata=INPUT_BINDINGS.samples or INPUT_BINDINGS.ped or "",
        fasta=apply_resource("fasta").source,
        fasta_index=apply_resource("fasta").source_index,
        fasta_gzi=optional_input(
            apply_resource("fasta").companion_indexes[0]
            if apply_resource("fasta").companion_indexes
            else ""
        ),
        libraries=VALIDATION_LIBRARIES,
        schema=VALIDATION_REPORT_SCHEMA
    output:
        report=ensure(str(COHORT_VALIDATION_REPORT), non_empty=True),
        manifest=protected(ensure(str(SAMPLE_MANIFEST_PATH), non_empty=True))
    params:
        hemizygote_model=WORKFLOW_PLAN.model.hemizygote,
        metadata_kind="samples" if INPUT_BINDINGS.samples else "ped"
    threads: 1
    resources:
        mem_mb=policy_resources("validate_cohort").mem_mb,
        runtime=policy_resources("validate_cohort").runtime
    conda:
        env("validation")
    envmodules: *mod("validation")
    container: sif("validation")
    log:
        str(APPLY_LOG_DIR / "validation" / "cohort.log")
    script:
        str(APPLY_SCRIPT_DIR / "validate_cohort.py")


rule validate_resources:
    """Validate interval, exclusion, testing-unit, and CADD resources."""
    input:
        vcf=apply_resource("cohort").source,
        calibration_regions=optional_input(
            optional_resource_field("calibration_regions", "source")
        ),
        calibration_regions_index=optional_input(
            optional_resource_field("calibration_regions", "source_index")
        ),
        variants_exclude=[
            resource_spec.source
            for resource_spec in resources_with_role("variants_exclude")
        ],
        variant_exclude_indexes=[
            resource_spec.source_index
            for resource_spec in resources_with_role("variants_exclude")
            if resource_spec.source_index
        ],
        testing_unit_bed=apply_resource("testing_units").source,
        testing_unit_index=optional_input(apply_resource("testing_units").source_index),
        par_regions=[
            resource_spec.source
            for resource_spec in resources_with_role("par_regions")
        ],
        par_region_indexes=[
            resource_spec.source_index
            for resource_spec in resources_with_role("par_regions")
            if resource_spec.source_index
        ],
        cadd_snv=apply_resource("cadd_snv").source,
        cadd_snv_index=optional_input(apply_resource("cadd_snv").source_index),
        cadd_indel=(
            optional_input(
                optional_resource_field("cadd_indel", "source")
            )
            if WORKFLOW_PLAN.model.include_indels
            else []
        ),
        cadd_indel_index=(
            optional_input(
                optional_resource_field("cadd_indel", "source_index")
            )
            if WORKFLOW_PLAN.model.include_indels
            else []
        ),
        libraries=VALIDATION_LIBRARIES,
        schema=VALIDATION_REPORT_SCHEMA
    output:
        report=ensure(str(RESOURCE_VALIDATION_REPORT), non_empty=True)
    params:
        unit=WORKFLOW_PLAN.model.unit,
        include_indels=WORKFLOW_PLAN.model.include_indels,
        calibration_regions=optional_apply_resource_binding("calibration_regions"),
        exclusions=resource_bindings_with_role("variants_exclude"),
        testing_units=apply_resource_binding("testing_units"),
        par_regions=resource_bindings_with_role("par_regions"),
        cadd_snv=apply_resource_binding("cadd_snv"),
        cadd_indel=optional_apply_resource_binding("cadd_indel")
    threads: 1
    resources:
        mem_mb=policy_resources("validate_resources").mem_mb,
        runtime=policy_resources("validate_resources").runtime
    conda:
        env("validation")
    envmodules: *mod("validation")
    container: sif("validation")
    log:
        str(APPLY_LOG_DIR / "validation" / "resources.log")
    script:
        str(APPLY_SCRIPT_DIR / "validate_resources.py")


rule validate_model:
    """Validate model configuration and enabled lookup resources."""
    input:
        testing_unit_bed=apply_resource("testing_units").source,
        lookup_het=LOOKUP_SOURCE.tables["het"],
        lookup_hom=LOOKUP_SOURCE.tables["hom"],
        lookup_hem_chrx=(
            [LOOKUP_SOURCE.tables["hem_chrX"]]
            if WORKFLOW_PLAN.model.hemizygote
            else []
        ),
        lookup_metadata=optional_input(LOOKUP_SOURCE.metadata),
        lookup_validation=optional_input(LOOKUP_SOURCE.validation),
        lookup_marker=optional_input(LOOKUP_SOURCE.marker),
        lookup_eligible_units=(
            optional_input(LOOKUP_SOURCE.eligible_units)
            if LOOKUP_SOURCE.is_bundle
            else []
        ),
        lookup_calibration_mask=(
            optional_input(LOOKUP_SOURCE.calibration_mask)
            if LOOKUP_SOURCE.is_bundle
            else []
        ),
        lookup_excluded_alleles=(
            list(LOOKUP_SOURCE.excluded_alleles) if LOOKUP_SOURCE.is_bundle else []
        ),
        lookup_excluded_units=(
            optional_input(LOOKUP_SOURCE.excluded_units)
            if LOOKUP_SOURCE.is_bundle
            else []
        ),
        libraries=VALIDATION_LIBRARIES,
        schema=VALIDATION_REPORT_SCHEMA
    output:
        report=ensure(str(MODEL_VALIDATION_REPORT), non_empty=True)
    params:
        score_max=WORKFLOW_PLAN.model.score_max,
        assembly=WORKFLOW_PLAN.reference.assembly,
        unit=WORKFLOW_PLAN.model.unit,
        hemizygote_model=WORKFLOW_PLAN.model.hemizygote,
        include_indels=WORKFLOW_PLAN.model.include_indels,
        lookup_contract=LOOKUP_SOURCE.contract,
        lookup_bin_width=LOOKUP_SOURCE.bin_width
    threads: 1
    resources:
        mem_mb=policy_resources("validate_model").mem_mb,
        runtime=policy_resources("validate_model").runtime
    conda:
        env("validation")
    envmodules: *mod("validation")
    container: sif("validation")
    log:
        str(APPLY_LOG_DIR / "validation" / "model.log")
    script:
        str(APPLY_SCRIPT_DIR / "validate_model.py")


rule merge_validation:
    """Merge independent validation fragments into validation.json."""
    input:
        cohort=rules.validate_cohort.output.report,
        resources=rules.validate_resources.output.report,
        model=rules.validate_model.output.report
    output:
        report=protected(ensure(str(COMBINED_VALIDATION_REPORT), non_empty=True))
    threads: 1
    resources:
        mem_mb=policy_resources("merge_validation").mem_mb,
        runtime=policy_resources("merge_validation").runtime
    conda:
        env("validation")
    envmodules: *mod("validation")
    container: sif("validation")
    log:
        str(APPLY_LOG_DIR / "validation" / "merge.log")
    script:
        str(APPLY_SCRIPT_DIR / "merge_validation.py")


rule require_valid_inputs:
    """Gate expensive workflow rules on the merged validation report."""
    input:
        report=rules.merge_validation.output.report,
        manifest=rules.validate_cohort.output.manifest
    output:
        ok=ensure(str(Path(APPLY_OUTPUT_DIR) / "validation.ok"), non_empty=True)
    params:
        failure_message="Input validation failed",
        include_errors=True,
        success_marker="validated\n"
    threads: 1
    resources:
        mem_mb=policy_resources("require_valid_inputs").mem_mb,
        runtime=policy_resources("require_valid_inputs").runtime
    conda:
        env("validation")
    envmodules: *mod("validation")
    container: sif("validation")
    log:
        str(APPLY_LOG_DIR / "validation" / "require_valid_inputs.log")
    script:
        str(APPLY_SCRIPT_DIR / "require_valid_report.py")
