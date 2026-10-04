"""Partition-local event aggregation and bounded unit-bucket reduction."""


rule calculate_partition_null_distribution:
    """Aggregate events once within a genomic partition and write bucket partials."""
    input:
        partitions=lambda wc: generation_scope_artifact(wc, "partitions"),
        buckets=lambda wc: generation_scope_artifact(wc, "buckets"),
        allele_states=str(GENERATION_PARTITION_DIR / "{partition}.allele_states.parquet"),
        allele_membership=str(GENERATION_PARTITION_DIR / "{partition}.allele_membership.parquet"),
        allele_scores=str(GENERATION_PARTITION_DIR / "{partition}.allele_scores.parquet"),
        mutation_rate=lambda wc: generation_mutation_rate_data(),
        mutation_rate_index=lambda wc: generation_mutation_rate_index()
    output:
        statistics=str(
            GENERATION_STATISTICS_DIR / "{partition}.bucket-{bucket}.parquet"
        )
    params:
        score_max=WORKFLOW_PLAN.model.score_max,
        bin_width=WORKFLOW_PLAN.generation_settings.score_bin_width,
        unobserved_policy=GENERATION_PLAN.unobserved_policy,
        roles=GENERATION_ROLES,
        x_af_policy=WORKFLOW_PLAN.generation_settings.x_af_policy,
        x_mutation_rate_policy=WORKFLOW_PLAN.generation_settings.x_mutation_rate_policy,
        compression=GENERATION_PARQUET_COMPRESSION
    threads: 1
    resources:
        mem_mb=policy_resources("calculate_partition_null_distribution").mem_mb,
        runtime=policy_resources("calculate_partition_null_distribution").runtime
    conda:
        env("validation")
    envmodules: *mod("validation")
    container: sif("validation")
    log:
        str(
            GENERATION_LOG_DIR
            / "partitions"
            / "{partition}.{bucket}.calculate.log"
        )
    script:
        str(GENERATION_SCRIPT_DIR / "calculate_partition_statistics.py")


rule reduce_null_distribution_bucket:
    """Produce direct lookup rows for one bounded range of unit identifiers."""
    input:
        buckets=lambda wc: generation_scope_artifact(wc, "buckets"),
        partition_statistics=generation_statistics_shards
    output:
        lookup_het=str(GENERATION_LOOKUP_BUCKET_DIR / "{bucket}.het.tsv"),
        lookup_hom=str(GENERATION_LOOKUP_BUCKET_DIR / "{bucket}.hom.tsv"),
        lookup_hem=(
            str(GENERATION_LOOKUP_BUCKET_DIR / "{bucket}.hem_chrX.tsv")
            if GENERATION_HEM_ENABLED
            else []
        ),
        diagnostics=str(GENERATION_LOOKUP_BUCKET_DIR / "{bucket}.diagnostics.tsv")
    params:
        bin_count=GENERATION_BIN_COUNT,
        roles=GENERATION_ROLES
    threads: 1
    resources:
        mem_mb=policy_resources("reduce_null_distribution_bucket").mem_mb,
        runtime=policy_resources("reduce_null_distribution_bucket").runtime
    conda:
        env("validation")
    envmodules: *mod("validation")
    container: sif("validation")
    log:
        str(GENERATION_LOG_DIR / "buckets" / "{bucket}.reduce.log")
    script:
        str(GENERATION_SCRIPT_DIR / "reduce_partition_statistics.py")
