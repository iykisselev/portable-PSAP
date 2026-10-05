"""Rules that materialize generation preparation artifacts."""


checkpoint prepare_generation_scope:
    """Prepare the shared mask and testing units before partition work begins."""
    input:
        population_vcf=GENERATION_CONTEXT.plan.population_vcf.data,
        population_vcf_index=lambda wc: generation_index(GENERATION_CONTEXT.plan.population_vcf),
        score_snv=GENERATION_CONTEXT.plan.score_snv.data,
        score_snv_index=lambda wc: generation_index(GENERATION_CONTEXT.plan.score_snv),
        score_indel=(
            [GENERATION_CONTEXT.plan.score_indel.data]
            if WORKFLOW_PLAN.model.include_indels and GENERATION_CONTEXT.plan.score_indel
            else []
        ),
        score_indel_index=(
            generation_index(GENERATION_CONTEXT.plan.score_indel)
            if WORKFLOW_PLAN.model.include_indels and GENERATION_CONTEXT.plan.score_indel
            else []
        ),
        exclusion=(
            [GENERATION_CONTEXT.plan.exclusion.data]
            if GENERATION_CONTEXT.plan.exclusion.data
            else []
        ),
        testing_units=GENERATION_CONTEXT.plan.testing_unit_bed.data,
        par_regions=GENERATION_CONTEXT.plan.par_regions.data,
        exclusions=lambda wc: generation_exclusion_data(),
        fasta=GENERATION_CONTEXT.plan.fasta.data,
        fasta_fai=GENERATION_CONTEXT.plan.fasta.fai,
        mutation_rate=lambda wc: generation_mutation_rate_data(),
        mutation_rate_index=lambda wc: generation_mutation_rate_index()
    output:
        partitions=GENERATION_CONTEXT.scope_paths["partitions"],
        buckets=GENERATION_CONTEXT.scope_paths["buckets"],
        units=GENERATION_CONTEXT.scope_paths["units"],
        unit_intervals=GENERATION_CONTEXT.scope_paths["unit_intervals"],
        masked_unit_intervals=GENERATION_CONTEXT.scope_paths["masked_unit_intervals"],
        mask=GENERATION_CONTEXT.scope_paths["effective_mask"],
        resource_identities=GENERATION_CONTEXT.scope_paths["resource_identities"]
    params:
        minimum_retained_fraction=WORKFLOW_PLAN.generation_settings.minimum_retained_fraction,
        partition_size=GENERATION_CONTEXT.plan.partition_size_bp,
        unit_bucket_size=GENERATION_CONTEXT.plan.unit_bucket_size,
        af_tags=dict(GENERATION_CONTEXT.plan.af_source),
        compression=GENERATION_CONTEXT.plan.parquet_compression,
        include_indels=WORKFLOW_PLAN.model.include_indels,
        indel_max_reference_length=WORKFLOW_PLAN.generation_settings.indel_max_reference_length,
        exclusions=GENERATION_CONTEXT.exclusion_paths_by_id()
    threads: 1
    resources:
        mem_mb=policy_resources("prepare_generation_scope").mem_mb,
        runtime=policy_resources("prepare_generation_scope").runtime
    conda:
        env("validation")
    envmodules: *mod("validation")
    container: sif("validation")
    log:
        str(GENERATION_CONTEXT.log_dir / "prepare_scope.log")
    script:
        str(GENERATION_SCRIPT_DIR / "prepare_generation_scope.py")


rule prepare_generation_partition:
    """Create allele states, unit memberships, and scores for one genomic window."""
    input:
        partitions=lambda wc: generation_scope_artifact(wc, "partitions"),
        units=lambda wc: generation_scope_artifact(wc, "units"),
        unit_intervals=lambda wc: generation_scope_artifact(wc, "unit_intervals"),
        mask=lambda wc: generation_scope_artifact(wc, "effective_mask"),
        population_vcf=GENERATION_CONTEXT.plan.population_vcf.data,
        population_vcf_index=lambda wc: generation_index(GENERATION_CONTEXT.plan.population_vcf),
        score_snv=GENERATION_CONTEXT.plan.score_snv.data,
        score_snv_index=lambda wc: generation_index(GENERATION_CONTEXT.plan.score_snv),
        fasta=GENERATION_CONTEXT.plan.fasta.data,
        fasta_fai=GENERATION_CONTEXT.plan.fasta.fai
    output:
        allele_states=str(GENERATION_CONTEXT.partition_dir / "{partition}.allele_states.parquet"),
        allele_membership=str(GENERATION_CONTEXT.partition_dir / "{partition}.allele_membership.parquet"),
        allele_scores=str(GENERATION_CONTEXT.partition_dir / "{partition}.allele_scores.parquet")
    params:
        af_tags=dict(GENERATION_CONTEXT.plan.af_source),
        x_af_policy=WORKFLOW_PLAN.generation_settings.x_af_policy,
        hemizygote=WORKFLOW_PLAN.model.hemizygote,
        compression=GENERATION_CONTEXT.plan.parquet_compression
    threads: 1
    resources:
        mem_mb=policy_resources("prepare_generation_partition").mem_mb,
        runtime=policy_resources("prepare_generation_partition").runtime
    conda:
        env("validation")
    envmodules: *mod("validation")
    container: sif("validation")
    log:
        str(GENERATION_CONTEXT.log_dir / "partitions" / "{partition}.prepare.log")
    script:
        str(GENERATION_SCRIPT_DIR / "prepare_generation_partition.py")
