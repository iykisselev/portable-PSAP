import shlex


PREPARATION_SCRIPT = str(APPLY_SCRIPT_DIR / "prepare_genomic_resource.py")
COMBINE_EXCLUSIONS_SCRIPT = str(APPLY_SCRIPT_DIR / "combine_exclusions.py")


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
        script=PREPARATION_SCRIPT
    wildcard_constraints:
        resource=VCF_RESOURCE_PATTERN
    threads: 1
    resources:
        mem_mb=4096,
        runtime=180
    conda:
        str(CONDA_ENV_DIR / "bioinformatics.yaml")
    envmodules: mod("bioinformatics")
    container: sif("bioinformatics")
    log:
        str(LOG_DIR / "apply" / config["run"]["apply_prefix"] / "prepare" / "vcf_{resource}.log")
    shell:
        """
        python {params.script:q} --input {input.source:q} --output {output.data:q} \
            --kind vcf > {log:q} 2>&1
        """


rule prepare_table_resource:
    """Stage one configured table resource using metadata-selected constants only."""
    input:
        source=table_resource_source,
        source_index=table_resource_source_index,
        cohort=table_resource_cohort,
        cohort_index=table_resource_cohort_index,
        validated=validation_marker
    output:
        data=ensure(str(PREPARED_TABLE_DIR / "{resource}.gz"), non_empty=True),
        index=ensure(str(PREPARED_TABLE_DIR / "{resource}.gz.tbi"), non_empty=True)
    params:
        script=PREPARATION_SCRIPT,
        kind=lambda wc: TABLE_RESOURCE_SPECS[wc.resource]["kind"],
        chromosome_map=table_resource_chromosome_map,
        tsv_layout=table_resource_tsv_layout
    wildcard_constraints:
        resource=TABLE_RESOURCE_PATTERN
    threads: 1
    resources:
        mem_mb=table_resource_mem_mb,
        runtime=table_resource_runtime_minutes
    conda:
        str(CONDA_ENV_DIR / "bioinformatics.yaml")
    envmodules: mod("bioinformatics")
    container: sif("bioinformatics")
    benchmark:
        str(BENCHMARK_DIR / "prepare_{resource}.tsv")
    log:
        str(LOG_DIR / "apply" / config["run"]["apply_prefix"] / "prepare" / "table_{resource}.log")
    shell:
        """
        python {params.script:q} --input {input.source:q} --output {output.data:q} \
            --kind {params.kind:q} {params.tsv_layout} {params.chromosome_map} \
            > {log:q} 2>&1
        """


rule combine_exclusions:
    """Project full-span interval and exact-allele exclusions onto the cohort VCF."""
    input:
        cohort=WORKFLOW_INPUT_VCF,
        cohort_index=WORKFLOW_INPUT_VCF_INDEX,
        exclusions=WORKFLOW_EXCLUSIONS,
        indexes=WORKFLOW_EXCLUSION_INDEXES
    output:
        vcf=ensure(COMBINED_EXCLUSIONS_VCF, non_empty=True)
    params:
        script=COMBINE_EXCLUSIONS_SCRIPT,
        source_args=" ".join(
            "--source "
            + shlex.quote(
                f"{EXCLUSION_KIND_BY_ID[item]}:{WORKFLOW_EXCLUSIONS[position]}"
            )
            for position, item in enumerate(EXCLUSION_IDS)
        )
    threads: 1
    resources:
        mem_mb=4096,
        runtime=180
    conda:
        str(CONDA_ENV_DIR / "validation.yaml")
    envmodules: mod("validation")
    container: sif("validation")
    log:
        str(LOG_DIR / "apply" / config["run"]["apply_prefix"] / "prepare" / "combine_exclusions.log")
    shell:
        """
        python {params.script:q} \
            --cohort-vcf {input.cohort:q} \
            {params.source_args} \
            --output {output.vcf:q} > {log:q} 2>&1
        """
