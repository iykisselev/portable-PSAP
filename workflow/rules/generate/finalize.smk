"""Streaming assembly and independent validation of generated bundles."""


rule finalize_null_distribution_bundle:
    """Render public bundle files from direct bucket lookup rows."""
    input:
        lookup_het=lambda wc: generation_bucket_artifacts(wc, "lookup_het"),
        lookup_hom=lambda wc: generation_bucket_artifacts(wc, "lookup_hom"),
        lookup_hem=(
            lambda wc: generation_bucket_artifacts(wc, "lookup_hem")
            if GENERATION_CONTEXT.hem_enabled
            else []
        ),
        diagnostics=lambda wc: generation_bucket_artifacts(wc, "diagnostics"),
        units=lambda wc: generation_scope_artifact(wc, "units"),
        masked_unit_intervals=lambda wc: generation_scope_artifact(wc, "masked_unit_intervals"),
        allele_states=lambda wc: generation_partition_artifacts(wc, "allele_states"),
        effective_mask=lambda wc: generation_scope_artifact(wc, "effective_mask"),
        resource_identities=lambda wc: generation_scope_artifact(wc, "resource_identities")
    output:
        lookup_het=ensure(GENERATION_CONTEXT.bundle_paths["lookup_het"], non_empty=True),
        lookup_hom=ensure(GENERATION_CONTEXT.bundle_paths["lookup_hom"], non_empty=True),
        lookup_hem=(
            ensure(GENERATION_CONTEXT.bundle_paths["lookup_hem"], non_empty=True)
            if GENERATION_CONTEXT.hem_enabled
            else []
        ),
        eligible_bed=ensure(GENERATION_CONTEXT.bundle_paths["eligible_bed"], non_empty=True),
        calibration_regions=ensure(
            GENERATION_CONTEXT.bundle_paths["calibration_regions"], non_empty=True
        ),
        excluded_alleles=ensure(
            GENERATION_CONTEXT.bundle_paths["excluded_alleles"], non_empty=True
        ),
        population_af_zero_imputations=ensure(
            GENERATION_CONTEXT.bundle_paths["population_af_zero_imputations"], non_empty=True
        ),
        excluded_units=ensure(GENERATION_CONTEXT.bundle_paths["excluded_units"], non_empty=True),
        metadata=ensure(GENERATION_CONTEXT.bundle_paths["metadata"], non_empty=True)
    params:
        minimum_retained_fraction=WORKFLOW_PLAN.generation_settings.minimum_retained_fraction,
        bin_count=GENERATION_CONTEXT.bin_count,
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
        resource_ids={
            "fasta": WORKFLOW_PLAN.reference.fasta_id,
            "exclusion": WORKFLOW_PLAN.reference.exclusion_id,
            "testing_units": WORKFLOW_PLAN.reference.testing_unit_id,
            "score_snv": WORKFLOW_PLAN.reference.cadd_snv_id,
            "population_vcf": f"population_vcf:{GENERATION_CONTEXT.plan.af_source['overall']}",
            "par_regions": "par_regions",
            "mutation_rate": "mutation_rate",
        }
    threads: 1
    resources:
        mem_mb=policy_resources("finalize_null_distribution_bundle").mem_mb,
        runtime=policy_resources("finalize_null_distribution_bundle").runtime
    conda:
        env("validation")
    envmodules: *mod("validation")
    container: sif("validation")
    log:
        str(GENERATION_CONTEXT.log_dir / "finalize.log")
    script:
        str(GENERATION_SCRIPT_DIR / "finalize_bundle.py")


rule validate_null_distribution_bundle:
    """Certify completed bundle components and write validation.ok last."""
    input:
        lookup_het=rules.finalize_null_distribution_bundle.output.lookup_het,
        lookup_hom=rules.finalize_null_distribution_bundle.output.lookup_hom,
        lookup_hem=(
            rules.finalize_null_distribution_bundle.output.lookup_hem
            if GENERATION_CONTEXT.hem_enabled
            else []
        ),
        eligible_bed=rules.finalize_null_distribution_bundle.output.eligible_bed,
        calibration_regions=rules.finalize_null_distribution_bundle.output.calibration_regions,
        excluded_alleles=rules.finalize_null_distribution_bundle.output.excluded_alleles,
        population_af_zero_imputations=(
            rules.finalize_null_distribution_bundle.output.population_af_zero_imputations
        ),
        excluded_units=rules.finalize_null_distribution_bundle.output.excluded_units,
        metadata=rules.finalize_null_distribution_bundle.output.metadata,
        schema=str(Path(workflow.basedir) / "schemas" / "null_distribution_bundle.schema.yaml")
    output:
        validation=ensure(GENERATION_CONTEXT.bundle_paths["validation"], non_empty=True),
        marker=ensure(GENERATION_CONTEXT.bundle_paths["marker"], non_empty=True)
    threads: 1
    resources:
        mem_mb=policy_resources("validate_null_distribution_bundle").mem_mb,
        runtime=policy_resources("validate_null_distribution_bundle").runtime
    conda:
        env("validation")
    envmodules: *mod("validation")
    container: sif("validation")
    log:
        str(GENERATION_CONTEXT.log_dir / "validate_bundle.log")
    script:
        str(GENERATION_SCRIPT_DIR / "validate_bundle.py")
