"""Rules that materialize generation preparation artifacts."""


checkpoint prepare_generation_scope:
    """Prepare the shared mask and testing units before partition work begins."""
    input:
        population_vcf=GENERATION_POPULATION_VCF_BINDING.data,
        score_snv=GENERATION_SCORE_SNV_BINDING.data,
        coverage=GENERATION_COVERAGE_BINDING.data,
        testing_units=GENERATION_TESTING_UNIT_BINDING.data,
        par_regions=GENERATION_PAR_REGIONS_BINDING.data,
        exclusions=lambda wc: generation_exclusion_data(),
        fasta=GENERATION_FASTA_BINDING.data,
        fasta_fai=GENERATION_FASTA_BINDING.fai,
        mutation_rate=lambda wc: generation_mutation_rate_data()
    output:
        partitions=GENERATION_SCOPE_PATHS["partitions"],
        buckets=GENERATION_SCOPE_PATHS["buckets"],
        units=GENERATION_SCOPE_PATHS["units"],
        unit_intervals=GENERATION_SCOPE_PATHS["unit_intervals"],
        masked_unit_intervals=GENERATION_SCOPE_PATHS["masked_unit_intervals"],
        mask=GENERATION_SCOPE_PATHS["effective_mask"],
        resource_identities=GENERATION_SCOPE_PATHS["resource_identities"]
    params:
        coverage_minimum=WORKFLOW_PLAN.generation_settings.coverage_minimum,
        partition_size=GENERATION_PARTITION_SIZE,
        unit_bucket_size=GENERATION_UNIT_BUCKET_SIZE,
        compression=GENERATION_PARQUET_COMPRESSION
    threads: 1
    resources:
        mem_mb=4096,
        runtime=120
    conda:
        str(Path(workflow.basedir) / "envs" / "validation.yaml")
    envmodules: *mod("validation")
    container: sif("validation")
    log:
        str(GENERATION_LOG_DIR / "prepare_scope.log")
    script:
        str(Path(workflow.basedir) / "scripts" / "generate" / "prepare_generation_scope.py")


rule prepare_generation_partition:
    """Create allele states, unit memberships, and scores for one genomic window."""
    input:
        partitions=lambda wc: generation_scope_artifact(wc, "partitions"),
        units=lambda wc: generation_scope_artifact(wc, "units"),
        unit_intervals=lambda wc: generation_scope_artifact(wc, "unit_intervals"),
        mask=lambda wc: generation_scope_artifact(wc, "effective_mask"),
        population_vcf=GENERATION_POPULATION_VCF_BINDING.data,
        population_vcf_index=lambda wc: generation_index(GENERATION_POPULATION_VCF_BINDING),
        score_snv=GENERATION_SCORE_SNV_BINDING.data,
        score_snv_index=lambda wc: generation_index(GENERATION_SCORE_SNV_BINDING),
        fasta=GENERATION_FASTA_BINDING.data,
        fasta_fai=GENERATION_FASTA_BINDING.fai
    output:
        allele_states=str(GENERATION_PARTITION_DIR / "{partition}.allele_states.parquet"),
        allele_membership=str(GENERATION_PARTITION_DIR / "{partition}.allele_membership.parquet"),
        allele_scores=str(GENERATION_PARTITION_DIR / "{partition}.allele_scores.parquet")
    params:
        af_tags=dict(GENERATION_PLAN.af_source),
        x_af_policy=WORKFLOW_PLAN.generation_settings.x_af_policy,
        hemizygote=WORKFLOW_PLAN.model.hemizygote,
        compression=GENERATION_PARQUET_COMPRESSION
    threads: 1
    resources:
        mem_mb=8192,
        runtime=240
    conda:
        str(Path(workflow.basedir) / "envs" / "validation.yaml")
    envmodules: *mod("validation")
    container: sif("validation")
    log:
        str(GENERATION_LOG_DIR / "partitions" / "{partition}.prepare.log")
    script:
        str(Path(workflow.basedir) / "scripts" / "generate" / "prepare_generation_partition.py")
