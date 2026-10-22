"""Rules that materialize generation preparation artifacts."""


checkpoint prepare_generation_scope:
    """Prepare the shared mask and testing units before window work begins."""
    input:
        **GENERATION_CONTEXT.raw_inputs("scope"),
    output:
        windows=GENERATION_CONTEXT.scope_paths["windows"],
        buckets=GENERATION_CONTEXT.scope_paths["buckets"],
        effective_mask=GENERATION_CONTEXT.scope_paths["effective_mask"],
        resolved_partition_intervals=GENERATION_CONTEXT.scope_paths["resolved_partition_intervals"],
        partition_statuses=GENERATION_CONTEXT.scope_paths["partition_statuses"],
        testing_unit_id_overlap=GENERATION_CONTEXT.scope_paths["testing_unit_id_overlap"],
    log:
        str(GENERATION_CONTEXT.plan.log_dir / "prepare_scope.log"),
    conda:
        env("validation")
    container:
        sif("validation")
    envmodules:
        *mod("validation"),
    threads: 1
    resources:
        mem_mb=policy_resources("prepare_generation_scope").mem_mb,
        runtime=policy_resources("prepare_generation_scope").runtime,
    params:
        minimum_retained_fraction=WORKFLOW_PLAN.generation_settings.minimum_retained_fraction,
        window_size=GENERATION_CONTEXT.plan.window_size_bp,
        unit_bucket_size=GENERATION_CONTEXT.plan.unit_bucket_size,
        af_tags=dict(GENERATION_CONTEXT.plan.af_source),
        compression=GENERATION_CONTEXT.plan.parquet_compression,
        include_indels=WORKFLOW_PLAN.model.include_indels,
        indel_max_reference_length=WORKFLOW_PLAN.generation_settings.indel_max_reference_length,
        partition_ids=tuple(partition.identifier for partition in GENERATION_CONTEXT.plan.testing_unit_partitions),
    script:
        str(GENERATION_SCRIPT_DIR / "prepare_generation_scope.py")


rule prepare_generation_window:
    """Create allele states, unit memberships, and scores for one genomic window."""
    input:
        **GENERATION_CONTEXT.raw_inputs("window"),
        windows=lambda wc: generation_scope_artifact(wc, "windows"),
        resolved_partition_intervals=lambda wc: generation_scope_artifact(wc, "resolved_partition_intervals"),
        partition_statuses=lambda wc: generation_scope_artifact(wc, "partition_statuses"),
        effective_mask=lambda wc: generation_scope_artifact(wc, "effective_mask"),
    output:
        allele_states=str(
            GENERATION_CONTEXT.window_dir / "{window}.allele_states.parquet"
        ),
        allele_membership=str(
            GENERATION_CONTEXT.window_dir / "{window}.allele_membership.parquet"
        ),
        allele_scores=str(
            GENERATION_CONTEXT.window_dir / "{window}.allele_scores.parquet"
        ),
        candidate_audit=str(
            GENERATION_CONTEXT.window_dir / "{window}.candidate_audit.parquet"
        ),
    log:
        str(GENERATION_CONTEXT.plan.log_dir / "windows" / "{window}.prepare.log"),
    conda:
        env("validation")
    container:
        sif("validation")
    envmodules:
        *mod("validation"),
    threads: 1
    resources:
        mem_mb=policy_resources("prepare_generation_window").mem_mb,
        runtime=policy_resources("prepare_generation_window").runtime,
    params:
        af_tags=dict(GENERATION_CONTEXT.plan.af_source),
        x_af_policy=WORKFLOW_PLAN.generation_settings.x_af_policy,
        hemizygote=WORKFLOW_PLAN.model.hemizygote,
        compression=GENERATION_CONTEXT.plan.parquet_compression,
        include_indels=WORKFLOW_PLAN.model.include_indels,
        indel_max_reference_length=WORKFLOW_PLAN.generation_settings.indel_max_reference_length,
        indel_max_alternate_length=WORKFLOW_PLAN.generation_settings.indel_max_alternate_length,
    script:
        str(GENERATION_SCRIPT_DIR / "prepare_generation_window.py")
