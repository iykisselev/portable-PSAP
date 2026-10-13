"""Streaming assembly and independent validation of generated bundles."""


rule finalize_null_distribution_bundle:
    """Render public bundle files from direct bucket lookup rows."""
    input:
        lookup_het=lambda wc: generation_bucket_artifacts(wc, "lookup_het"),
        lookup_hom=lambda wc: generation_bucket_artifacts(wc, "lookup_hom"),
        lookup_hem_chrX=lambda wc: generation_bucket_artifacts(wc, "lookup_hem_chrX"),
        diagnostics=lambda wc: generation_bucket_artifacts(wc, "diagnostics"),
        units=lambda wc: generation_scope_artifact(wc, "units"),
        masked_unit_intervals=lambda wc: generation_scope_artifact(
            wc, "masked_unit_intervals"
        ),
        allele_states=lambda wc: generation_partition_artifacts(wc, "allele_states"),
        candidate_audit=lambda wc: generation_partition_artifacts(wc, "candidate_audit"),
        effective_mask=lambda wc: generation_scope_artifact(wc, "effective_mask"),
    output:
        **{
            name: ensure(path, non_empty=True) if path else []
            for name, path in GENERATION_CONTEXT.bundle_outputs().items()
        },
    log:
        str(GENERATION_CONTEXT.plan.log_dir / "finalize.log"),
    conda:
        env("validation")
    container:
        sif("validation")
    envmodules:
        *mod("validation"),
    threads: 1
    resources:
        mem_mb=policy_resources("finalize_null_distribution_bundle").mem_mb,
        runtime=policy_resources("finalize_null_distribution_bundle").runtime,
    params:
        minimum_retained_fraction=WORKFLOW_PLAN.generation_settings.minimum_retained_fraction,
        bin_count=GENERATION_CONTEXT.execution_layout.bin_count,
        bin_width=WORKFLOW_PLAN.generation_settings.score_bin_width,
        score_max=WORKFLOW_PLAN.model.score_max,
        policy=GENERATION_CONTEXT.plan.unobserved_policy,
        conversion=WORKFLOW_PLAN.generation_settings.mutation_conversion,
        assembly=WORKFLOW_PLAN.reference.assembly,
        testing_unit_type=WORKFLOW_PLAN.model.unit,
        af_tags=dict(GENERATION_CONTEXT.plan.af_source),
        x_af_policy=WORKFLOW_PLAN.generation_settings.x_af_policy,
        x_mutation_rate_policy=WORKFLOW_PLAN.generation_settings.x_mutation_rate_policy,
        hemizygote=WORKFLOW_PLAN.model.hemizygote,
        variant_classes=list(WORKFLOW_PLAN.generation_settings.variant_classes),
        indel_candidate_universe=WORKFLOW_PLAN.generation_settings.indel_candidate_universe,
        indel_membership=WORKFLOW_PLAN.generation_settings.indel_membership,
        indel_mask_policy=WORKFLOW_PLAN.generation_settings.indel_mask_policy,
        indel_max_reference_length=WORKFLOW_PLAN.generation_settings.indel_max_reference_length,
        indel_max_alternate_length=WORKFLOW_PLAN.generation_settings.indel_max_alternate_length,
        resource_identities={
            name: dict(identity)
            for name, identity in GENERATION_CONTEXT.bundle_resource_identities.items()
        },
    script:
        str(GENERATION_SCRIPT_DIR / "finalize_bundle.py")


rule validate_null_distribution_bundle:
    """Certify completed bundle components and write validation.ok last."""
    input:
        **{
            name: getattr(rules.finalize_null_distribution_bundle.output, name)
            for name in GENERATION_CONTEXT.bundle_outputs()
        },
        schema=str(
            Path(workflow.basedir) / "schemas" / "null_distribution_bundle.schema.yaml"
        ),
    output:
        validation=ensure(GENERATION_CONTEXT.bundle_paths["validation"], non_empty=True),
        marker=ensure(GENERATION_CONTEXT.bundle_paths["marker"], non_empty=True),
    log:
        str(GENERATION_CONTEXT.plan.log_dir / "validate_bundle.log"),
    conda:
        env("validation")
    container:
        sif("validation")
    envmodules:
        *mod("validation"),
    threads: 1
    resources:
        mem_mb=policy_resources("validate_null_distribution_bundle").mem_mb,
        runtime=policy_resources("validate_null_distribution_bundle").runtime,
    script:
        str(GENERATION_SCRIPT_DIR / "validate_bundle.py")
