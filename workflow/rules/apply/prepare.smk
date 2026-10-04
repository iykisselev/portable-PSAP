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
        mem_mb=policy_resources("prepare_vcf_resource").mem_mb,
        runtime=policy_resources("prepare_vcf_resource").runtime
    conda:
        env("bioinformatics")
    envmodules: *mod("bioinformatics")
    container: sif("bioinformatics")
    log:
        str(APPLY_LOG_DIR / "prepare" / "vcf_{resource}.log")
    script:
        str(APPLY_SCRIPT_DIR / "prepare_resource.py")


rule prepare_table_resource:
    """Stage one configured table resource using metadata-selected constants only."""
    input:
        source=table_resource_source,
        source_index=table_resource_source_index,
        cohort=table_resource_cohort,
        cohort_index=table_resource_cohort_index,
        validated=validation_marker
    output:
        data=ensure(
            str(APPLY_RESOURCE_REGISTRY.prepared_table_dir / "{resource}.{extension}.gz"),
            non_empty=True,
        ),
        index=ensure(
            str(APPLY_RESOURCE_REGISTRY.prepared_table_dir / "{resource}.{extension}.gz.tbi"),
            non_empty=True,
        )
    params:
        kind=table_resource_kind,
        tsv_layout=table_resource_tsv_layout
    wildcard_constraints:
        resource=TABLE_RESOURCE_PATTERN,
        extension="(?:bed|tsv|vcf)"
    threads: 1
    resources:
        mem_mb=TABLE_RESOURCE_CALLBACKS.mem_mb,
        runtime=TABLE_RESOURCE_CALLBACKS.runtime
    conda:
        env("bioinformatics")
    envmodules: *mod("bioinformatics")
    container: sif("bioinformatics")
    benchmark:
        str(BENCHMARK_DIR / "prepare_{resource}_{extension}.tsv")
    log:
        str(APPLY_LOG_DIR / "prepare" / "table_{resource}_{extension}.log")
    script:
        str(APPLY_SCRIPT_DIR / "prepare_resource.py")


rule combine_exclusions:
    """Project full-span interval and exact-allele exclusions onto the cohort VCF."""
    input:
        cohort=apply_resource("cohort").staged_data,
        cohort_index=apply_resource("cohort").staged_index,
        exclusions=list(APPLY_RESOURCE_REGISTRY.workflow_exclusions),
        validated=validation_marker
    output:
        vcf=ensure(APPLY_RESOURCE_REGISTRY.combined_exclusions_vcf, non_empty=True)
    params:
        sources=[
            f"{APPLY_RESOURCE_REGISTRY.exclusion_kind_by_id[item]}:"
            f"{APPLY_RESOURCE_REGISTRY.workflow_exclusions[position]}"
            for position, item in enumerate(APPLY_RESOURCE_REGISTRY.exclusion_ids)
        ]
    threads: 1
    resources:
        mem_mb=policy_resources("combine_exclusions").mem_mb,
        runtime=policy_resources("combine_exclusions").runtime
    conda:
        env("validation")
    envmodules: *mod("validation")
    container: sif("validation")
    log:
        str(APPLY_LOG_DIR / "prepare" / "combine_exclusions.log")
    script:
        str(APPLY_SCRIPT_DIR / "combine_exclusions.py")
