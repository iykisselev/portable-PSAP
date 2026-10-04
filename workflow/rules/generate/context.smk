"""Shared paths and wildcard helpers for generation rules."""

from generation_topology import (
    build_generation_execution_layout,
    build_generation_scope_artifacts,
)
from generation_manifests import bucket_ids, partition_ids


if WORKFLOW_PLAN.generation is None:
    raise ValueError("Generation rules require a generate or full WorkflowPlan")

GENERATION_PLAN = WORKFLOW_PLAN.generation
GENERATION_SCOPE_LAYOUT = build_generation_scope_artifacts(GENERATION_PLAN)
GENERATION_EXECUTION_LAYOUT = build_generation_execution_layout(
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
GENERATION_EXCLUSION_BINDING = GENERATION_PLAN.exclusion
GENERATION_TESTING_UNIT_BINDING = GENERATION_PLAN.testing_unit_bed
GENERATION_FASTA_BINDING = GENERATION_PLAN.fasta
GENERATION_PAR_REGIONS_BINDING = GENERATION_PLAN.par_regions
GENERATION_EXCLUSION_BINDINGS = GENERATION_PLAN.exclusions
GENERATION_SCOPE_PATHS = GENERATION_SCOPE_LAYOUT.artifacts
GENERATION_STATISTICS_DIR = GENERATION_EXECUTION_LAYOUT.statistics_dir
GENERATION_LOOKUP_BUCKET_DIR = GENERATION_EXECUTION_LAYOUT.lookup_bucket_dir
GENERATION_BIN_COUNT = GENERATION_EXECUTION_LAYOUT.bin_count
GENERATION_ROLES = GENERATION_EXECUTION_LAYOUT.roles
GENERATION_HEM_ENABLED = WORKFLOW_PLAN.model.hemizygote
GENERATION_BUNDLE_PATHS = GENERATION_EXECUTION_LAYOUT.bundle_paths


def generation_index(resource_binding):
    return [resource_binding.index] if resource_binding.index else []


def generation_mutation_rate_data():
    return [GENERATION_PLAN.mutation_rate.data] if GENERATION_PLAN.mutation_rate else []


def generation_mutation_rate_index():
    if GENERATION_PLAN.mutation_rate and GENERATION_PLAN.mutation_rate.index:
        return [GENERATION_PLAN.mutation_rate.index]
    return []


def generation_exclusion_data():
    return [exclusion_binding.data for exclusion_binding in GENERATION_EXCLUSION_BINDINGS]


def generation_scope_artifact(wildcards, artifact_name):
    """Return one scope file after the scope checkpoint has completed."""
    scope_checkpoint = checkpoints.prepare_generation_scope.get()
    output_name = {"effective_mask": "mask"}.get(artifact_name, artifact_name)
    return str(getattr(scope_checkpoint.output, output_name))


def generation_partition_ids(wildcards):
    """Read partition IDs from the completed scope manifest."""
    return partition_ids(generation_scope_artifact(wildcards, "partitions"))


def generation_bucket_ids(wildcards):
    """Read unit-bucket IDs from the completed scope manifest."""
    return bucket_ids(generation_scope_artifact(wildcards, "buckets"))


def generation_partition_artifacts(wildcards, artifact_name):
    """Return the requested file for every discovered genomic partition."""
    artifact_suffixes = {
        "allele_states": "allele_states.parquet",
        "allele_membership": "allele_membership.parquet",
        "allele_scores": "allele_scores.parquet",
    }
    return [
        str(GENERATION_PARTITION_DIR / f"{partition_id}.{artifact_suffixes[artifact_name]}")
        for partition_id in generation_partition_ids(wildcards)
    ]


def generation_bucket_artifacts(wildcards, artifact_name):
    """Return the requested lookup file for every discovered unit bucket."""
    artifact_suffixes = {
        "lookup_het": "het.tsv",
        "lookup_hom": "hom.tsv",
        "lookup_hem": "hem_chrX.tsv",
        "diagnostics": "diagnostics.tsv",
    }
    return [
        str(GENERATION_LOOKUP_BUCKET_DIR / f"{bucket_id}.{artifact_suffixes[artifact_name]}")
        for bucket_id in generation_bucket_ids(wildcards)
    ]


def generation_statistics_shards(wildcards):
    """Return the partition statistics needed to reduce one unit bucket."""
    bucket_id = str(wildcards.bucket)
    if bucket_id not in generation_bucket_ids(wildcards):
        raise ValueError(f"Unknown generation bucket {bucket_id!r}")
    return [
        str(GENERATION_STATISTICS_DIR / f"{partition_id}.bucket-{bucket_id}.parquet")
        for partition_id in generation_partition_ids(wildcards)
    ]


def generation_statistics_artifacts(wildcards):
    """Return every partition-by-bucket statistic file needed for recovery."""
    return [
        str(GENERATION_STATISTICS_DIR / f"{partition_id}.bucket-{bucket_id}.parquet")
        for partition_id in generation_partition_ids(wildcards)
        for bucket_id in generation_bucket_ids(wildcards)
    ]
