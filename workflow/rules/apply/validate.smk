rule validate_cohort:
    """Validate cohort VCF/metadata/FASTA inputs and write the sample manifest."""
    input:
        vcf=VALIDATION_VCF,
        vcf_index=VALIDATION_VCF_INDEX_INPUT,
        metadata=VALIDATION_SAMPLES or VALIDATION_PED,
        fasta=VALIDATION_FASTA,
        fasta_index=VALIDATION_FASTA_INDEX,
        fasta_gzi=VALIDATION_FASTA_GZI_INPUT,
        library=str(VALIDATION_LIBRARY)
    output:
        report=ensure(str(COHORT_VALIDATION), non_empty=True),
        manifest=protected(ensure(str(SAMPLE_MANIFEST_PATH), non_empty=True))
    params:
        hemizygote_model=str(config["model"]["hemizygote"]),
        metadata_kind="samples" if VALIDATION_SAMPLES else "ped"
    threads: 1
    resources:
        mem_mb=4096,
        runtime=30
    conda:
        str(CONDA_ENV_DIR / "validation.yaml")
    envmodules: *mod("validation")
    container: sif("validation")
    log:
        str(LOG_DIR / "apply" / config["run"]["apply_prefix"] / "validation" / "cohort.log")
    script:
        str(Path(workflow.basedir) / "scripts" / "validate_cohort.py")


rule validate_resources:
    """Validate interval, exclusion, testing-unit, and CADD resources."""
    input:
        vcf=VALIDATION_VCF,
        coverage=VALIDATION_COVERAGE_INPUT,
        coverage_index=VALIDATION_COVERAGE_INDEX_INPUT,
        variants_exclude=VALIDATION_VARIANT_EXCLUSIONS,
        variant_exclude_indexes=[index for index in VALIDATION_VARIANT_EXCLUSION_INDEXES if index],
        testing_unit_bed=VALIDATION_TESTING_UNIT_BED,
        testing_unit_index=VALIDATION_TESTING_UNIT_INDEX_INPUT,
        par_regions=VALIDATION_PAR_REGIONS,
        par_region_indexes=[index for index in VALIDATION_PAR_REGION_INDEXES if index],
        cadd_snv=VALIDATION_CADD_SNV,
        cadd_snv_index=VALIDATION_CADD_SNV_INDEX_INPUT,
        cadd_indel=VALIDATION_CADD_INDEL_INPUT,
        cadd_indel_index=VALIDATION_CADD_INDEL_INDEX_INPUT,
        library=str(VALIDATION_LIBRARY)
    output:
        report=ensure(str(RESOURCE_VALIDATION), non_empty=True)
    params:
        unit=str(config["model"]["unit"]),
        include_indels=str(config["model"]["include_indels"]),
        coverage_index=VALIDATION_COVERAGE_INDEX,
        variant_exclude_indexes=VALIDATION_VARIANT_EXCLUSION_INDEXES,
        testing_unit_index=VALIDATION_TESTING_UNIT_INDEX,
        par_region_indexes=VALIDATION_PAR_REGION_INDEXES,
        cadd_snv_index=VALIDATION_CADD_SNV_INDEX,
        cadd_indel_index=VALIDATION_CADD_INDEL_INDEX
    threads: 1
    resources:
        mem_mb=4096,
        runtime=30
    conda:
        str(CONDA_ENV_DIR / "validation.yaml")
    envmodules: *mod("validation")
    container: sif("validation")
    log:
        str(LOG_DIR / "apply" / config["run"]["apply_prefix"] / "validation" / "resources.log")
    script:
        str(Path(workflow.basedir) / "scripts" / "validate_resources.py")


rule validate_model:
    """Validate model configuration and enabled lookup resources."""
    input:
        testing_unit_bed=VALIDATION_TESTING_UNIT_BED,
        lookup_het=VALIDATION_LOOKUP_HET,
        lookup_hom=VALIDATION_LOOKUP_HOM,
        lookup_hem_chrx=VALIDATION_LOOKUP_HEM_CHRX_INPUT,
        library=str(VALIDATION_LIBRARY)
    output:
        report=ensure(str(MODEL_VALIDATION), non_empty=True)
    params:
        score_max=str(config["model"]["score_max"]),
        assembly=str(config["reference"]["assembly"]),
        unit=str(config["model"]["unit"]),
        hemizygote_model=str(config["model"]["hemizygote"]),
        include_indels=str(config["model"]["include_indels"])
    threads: 1
    resources:
        mem_mb=512,
        runtime=10
    conda:
        str(CONDA_ENV_DIR / "validation.yaml")
    envmodules: *mod("validation")
    container: sif("validation")
    log:
        str(LOG_DIR / "apply" / config["run"]["apply_prefix"] / "validation" / "model.log")
    script:
        str(Path(workflow.basedir) / "scripts" / "validate_model.py")


rule merge_validation:
    """Merge independent validation fragments into validation.json."""
    input:
        cohort=rules.validate_cohort.output.report,
        resources=rules.validate_resources.output.report,
        model=rules.validate_model.output.report
    output:
        report=protected(ensure(str(VALIDATION_REPORT), non_empty=True))
    threads: 1
    resources:
        mem_mb=512,
        runtime=5
    conda:
        str(CONDA_ENV_DIR / "validation.yaml")
    envmodules: *mod("validation")
    container: sif("validation")
    log:
        str(LOG_DIR / "apply" / config["run"]["apply_prefix"] / "validation" / "merge.log")
    script:
        str(Path(workflow.basedir) / "scripts" / "merge_validation.py")


checkpoint require_valid_inputs:
    """Gate expensive workflow rules on the merged validation report."""
    input:
        report=rules.merge_validation.output.report,
        manifest=rules.validate_cohort.output.manifest
    output:
        ok=ensure(str(Path(APPLY_OUTDIR) / "validation.ok"), non_empty=True)
    params:
        failure_message="Input validation failed",
        include_errors=True,
        success_marker="validated\n"
    threads: 1
    resources:
        mem_mb=512,
        runtime=5
    conda:
        str(CONDA_ENV_DIR / "validation.yaml")
    envmodules: *mod("validation")
    container: sif("validation")
    log:
        str(LOG_DIR / "apply" / config["run"]["apply_prefix"] / "validation" / "require_valid_inputs.log")
    script:
        str(Path(workflow.basedir) / "scripts" / "require_valid_report.py")
