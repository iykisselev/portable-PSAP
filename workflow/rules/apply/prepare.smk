rule prepare_vcf_resource:
    """Stage either the cohort or projected exclusions at a stable VCF target."""
    input:
        source=vcf_resource_source,
        source_index=vcf_resource_source_index,
        validated=vcf_resource_validation_marker
    output:
        data=ensure(str(PREPARED_DIR / "{resource}.vcf.gz"), non_empty=True),
        index=ensure(str(PREPARED_DIR / "{resource}.vcf.gz.tbi"), non_empty=True)
    params:
        kind="vcf",
        tsv_layout="auto"
    wildcard_constraints:
        resource=VCF_RESOURCE_PATTERN
    threads: 1
    resources:
        mem_mb=lambda wc, input, attempt: resource_mem_mb(
            "stage_vcf", input, attempt
        ),
        runtime=lambda wc, input, attempt: resource_runtime_minutes(
            "stage_vcf", input, attempt
        )
    conda:
        str(CONDA_ENV_DIR / "bioinformatics.yaml")
    envmodules: *mod("bioinformatics")
    container: sif("bioinformatics")
    log:
        str(APPLY_LOG_DIR / "prepare" / "vcf_{resource}.log")
    script:
        str(Path(workflow.basedir) / "scripts" / "prepare_resource.py")


rule prepare_table_resource:
    """Stage one configured table resource using metadata-selected constants only."""
    input:
        source=table_resource_source,
        source_index=table_resource_source_index,
        cohort=table_resource_cohort,
        cohort_index=table_resource_cohort_index,
        validated=validation_marker
    output:
        data=ensure(str(PREPARED_TABLE_DIR / "{resource}.{extension}.gz"), non_empty=True),
        index=ensure(str(PREPARED_TABLE_DIR / "{resource}.{extension}.gz.tbi"), non_empty=True)
    params:
        kind=table_resource_kind,
        tsv_layout=table_resource_tsv_layout
    wildcard_constraints:
        resource=TABLE_RESOURCE_PATTERN,
        extension="(?:bed|tsv|vcf)"
    threads: 1
    resources:
        mem_mb=table_resource_mem_mb,
        runtime=table_resource_runtime_minutes
    conda:
        str(CONDA_ENV_DIR / "bioinformatics.yaml")
    envmodules: *mod("bioinformatics")
    container: sif("bioinformatics")
    benchmark:
        str(BENCHMARK_DIR / "prepare_{resource}_{extension}.tsv")
    log:
        str(APPLY_LOG_DIR / "prepare" / "table_{resource}_{extension}.log")
    script:
        str(Path(workflow.basedir) / "scripts" / "prepare_resource.py")


rule combine_exclusions:
    """Project full-span interval and exact-allele exclusions onto the cohort VCF."""
    input:
        cohort=WORKFLOW_INPUT_VCF,
        cohort_index=WORKFLOW_INPUT_VCF_INDEX,
        exclusions=WORKFLOW_EXCLUSIONS,
        validated=validation_marker
    output:
        vcf=ensure(COMBINED_EXCLUSIONS_VCF, non_empty=True)
    params:
        sources=[
            f"{EXCLUSION_KIND_BY_ID[item]}:{WORKFLOW_EXCLUSIONS[position]}"
            for position, item in enumerate(EXCLUSION_IDS)
        ]
    threads: 1
    resources:
        mem_mb=12000,
        runtime=180
    conda:
        str(CONDA_ENV_DIR / "validation.yaml")
    envmodules: *mod("validation")
    container: sif("validation")
    log:
        str(APPLY_LOG_DIR / "prepare" / "combine_exclusions.log")
    script:
        str(Path(workflow.basedir) / "scripts" / "combine_exclusions.py")
