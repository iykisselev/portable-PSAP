rule prepare_vcf_resource:
    """Stage the cohort VCF at a stable BGZF target."""
    input:
        source=vcf_resource_source,
        source_index=vcf_resource_source_index,
        validated=vcf_resource_validation_marker,
    output:
        data=ensure(str(PREPARED_DIR / "{resource}.vcf.gz"), non_empty=True),
        index=ensure(str(PREPARED_DIR / "{resource}.vcf.gz.tbi"), non_empty=True),
    log:
        str(APPLY_LOG_DIR / "prepare" / "vcf_{resource}.log"),
    wildcard_constraints:
        resource=VCF_RESOURCE_PATTERN,
    conda:
        env("bioinformatics")
    container:
        sif("bioinformatics")
    envmodules:
        *mod("bioinformatics"),
    threads: 1
    resources:
        mem_mb=policy_resources("prepare_vcf_resource").mem_mb,
        runtime=policy_resources("prepare_vcf_resource").runtime,
    params:
        kind="vcf",
        tsv_layout="auto",
    script:
        str(APPLY_SCRIPT_DIR / "prepare_resource.py")


rule prepare_table_resource:
    """Stage one configured table resource using metadata-selected constants only."""
    input:
        source=table_resource_source,
        source_index=table_resource_source_index,
        cohort=table_resource_cohort,
        cohort_index=table_resource_cohort_index,
        validated=rules.require_valid_inputs.output.ok,
    output:
        data=ensure(
            str(
                APPLY_RESOURCE_REGISTRY.prepared_table_dir
                / "{resource}.{extension}.gz"
            ),
            non_empty=True,
        ),
        index=ensure(
            str(
                APPLY_RESOURCE_REGISTRY.prepared_table_dir
                / "{resource}.{extension}.gz.tbi"
            ),
            non_empty=True,
        ),
    log:
        str(APPLY_LOG_DIR / "prepare" / "table_{resource}_{extension}.log"),
    benchmark:
        str(BENCHMARK_DIR / "prepare_{resource}_{extension}.tsv")
    wildcard_constraints:
        resource=TABLE_RESOURCE_PATTERN,
        extension="(?:bed|tsv|vcf)",
    conda:
        env("bioinformatics")
    container:
        sif("bioinformatics")
    envmodules:
        *mod("bioinformatics"),
    threads: 1
    resources:
        mem_mb=TABLE_RESOURCE_CALLBACKS.mem_mb,
        runtime=TABLE_RESOURCE_CALLBACKS.runtime,
    params:
        kind=table_resource_kind,
        tsv_layout=table_resource_tsv_layout,
    script:
        str(APPLY_SCRIPT_DIR / "prepare_resource.py")
