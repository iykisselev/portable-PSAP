"""Plan-derived generation paths and short Snakemake wildcard adapters."""

from generation_topology import (
    build_generation_execution_layout,
    build_generation_scope_artifacts,
)
from generation_manifests import bucket_ids, partition_ids


if WORKFLOW_PLAN.generation is None:
    raise ValueError("Generation rules require a generate or full WorkflowPlan")

GENERATION_PLAN = WORKFLOW_PLAN.generation
GENERATION_SCOPE = build_generation_scope_artifacts(GENERATION_PLAN)
GENERATION_EXECUTION = build_generation_execution_layout(
    GENERATION_PLAN,
    WORKFLOW_PLAN.model,
    WORKFLOW_PLAN.generation_settings,
)
GENERATION_PARTITION_DIR = GENERATION_PLAN.scalable_dir / "partitions"
GENERATION_LOG_DIR = GENERATION_PLAN.log_dir
GENERATION_PARTITION_SIZE = GENERATION_PLAN.partition_size_bp
GENERATION_UNIT_BUCKET_SIZE = GENERATION_PLAN.unit_bucket_size
GENERATION_PARQUET_COMPRESSION = GENERATION_PLAN.parquet_compression
GENERATION_POPULATION_VCF_BINDING = GENERATION_PLAN.population_vcf
GENERATION_SCORE_SNV_BINDING = GENERATION_PLAN.score_snv
GENERATION_COVERAGE_BINDING = GENERATION_PLAN.coverage
GENERATION_TESTING_UNIT_BINDING = GENERATION_PLAN.testing_unit_bed
GENERATION_FASTA_BINDING = GENERATION_PLAN.fasta
GENERATION_PAR_REGIONS_BINDING = GENERATION_PLAN.par_regions
GENERATION_EXCLUSION_BINDINGS = GENERATION_PLAN.exclusions
GENERATION_ARTIFACTS = GENERATION_SCOPE.artifacts
GENERATION_STATISTICS_DIR = GENERATION_EXECUTION.statistics_dir
GENERATION_LOOKUP_BUCKET_DIR = GENERATION_EXECUTION.lookup_bucket_dir
GENERATION_BIN_COUNT = GENERATION_EXECUTION.bin_count
GENERATION_ROLES = GENERATION_EXECUTION.roles
GENERATION_HEM_ENABLED = WORKFLOW_PLAN.model.hemizygote
GENERATION_BUNDLE_PATHS = GENERATION_EXECUTION.bundle_paths


def generation_index(binding):
    return [binding.index] if binding.index else []


def generation_mutation_rate_data():
    return [GENERATION_PLAN.mutation_rate.data] if GENERATION_PLAN.mutation_rate else []


def generation_mutation_rate_index():
    if GENERATION_PLAN.mutation_rate and GENERATION_PLAN.mutation_rate.index:
        return [GENERATION_PLAN.mutation_rate.index]
    return []


def generation_exclusion_data():
    return [binding.data for binding in GENERATION_EXCLUSION_BINDINGS]


def generation_scope_artifact(wildcards, artifact):
    """Resolve one scope output after its topology checkpoint completes."""
    scope = checkpoints.prepare_generation_scope.get()
    output_name = {"effective_mask": "mask"}.get(artifact, artifact)
    return str(getattr(scope.output, output_name))


def generation_partition_ids(wildcards):
    """Expand partition work from the validated scope manifest only."""
    return partition_ids(generation_scope_artifact(wildcards, "partitions"))


def generation_bucket_ids(wildcards):
    """Expand unit-bucket work from the validated scope manifest only."""
    return bucket_ids(generation_scope_artifact(wildcards, "buckets"))


def generation_partition_artifacts(wildcards, artifact):
    """Return one explicit output path for every discovered partition."""
    suffixes = {
        "allele_states": "allele_states.parquet",
        "allele_membership": "allele_membership.parquet",
        "allele_scores": "allele_scores.parquet",
    }
    return [
        str(GENERATION_PARTITION_DIR / f"{partition}.{suffixes[artifact]}")
        for partition in generation_partition_ids(wildcards)
    ]


def generation_bucket_artifacts(wildcards, artifact):
    """Return one explicit reduction output path for every discovered bucket."""
    suffixes = {
        "lookup_het": "het.tsv",
        "lookup_hom": "hom.tsv",
        "lookup_hem": "hem_chrX.tsv",
        "diagnostics": "diagnostics.tsv",
    }
    return [
        str(GENERATION_LOOKUP_BUCKET_DIR / f"{bucket}.{suffixes[artifact]}")
        for bucket in generation_bucket_ids(wildcards)
    ]


def generation_statistics_shards(wildcards):
    """Return all explicit partition shards needed for one unit bucket."""
    bucket = str(wildcards.bucket)
    if bucket not in generation_bucket_ids(wildcards):
        raise ValueError(f"Unknown generation bucket {bucket!r}")
    return [
        str(GENERATION_STATISTICS_DIR / f"{partition}.bucket-{bucket}.parquet")
        for partition in generation_partition_ids(wildcards)
    ]


def generation_statistics_artifacts(wildcards):
    """Return every explicit partition–bucket statistic shard for recovery."""
    return [
        str(GENERATION_STATISTICS_DIR / f"{partition}.bucket-{bucket}.parquet")
        for partition in generation_partition_ids(wildcards)
        for bucket in generation_bucket_ids(wildcards)
    ]
