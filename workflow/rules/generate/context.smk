"""Shared paths and wildcard helpers for generation rules."""

from generation_context import GenerationRuleContext
from generation_manifests import bucket_ids, partition_ids


GENERATION_CONTEXT = GenerationRuleContext.from_workflow_plan(WORKFLOW_PLAN)
GENERATION_SCRIPT_DIR = Path(workflow.basedir) / "scripts" / "generate"


def generation_index(resource_binding):
    return GENERATION_CONTEXT.index_input(resource_binding)


def generation_mutation_rate_data():
    return GENERATION_CONTEXT.mutation_rate_data()


def generation_mutation_rate_index():
    return GENERATION_CONTEXT.mutation_rate_index()


def generation_exclusion_data():
    return GENERATION_CONTEXT.exclusion_data()


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
    return [
        GENERATION_CONTEXT.partition_artifact(partition_id, artifact_name)
        for partition_id in generation_partition_ids(wildcards)
    ]


def generation_bucket_artifacts(wildcards, artifact_name):
    """Return the requested lookup file for every discovered unit bucket."""
    return [
        GENERATION_CONTEXT.bucket_artifact(bucket_id, artifact_name)
        for bucket_id in generation_bucket_ids(wildcards)
    ]


def generation_statistics_shards(wildcards):
    """Return the partition statistics needed to reduce one unit bucket."""
    bucket_id = str(wildcards.bucket)
    if bucket_id not in generation_bucket_ids(wildcards):
        raise ValueError(f"Unknown generation bucket {bucket_id!r}")
    return [
        GENERATION_CONTEXT.statistics_artifact(partition_id, bucket_id)
        for partition_id in generation_partition_ids(wildcards)
    ]


def generation_statistics_artifacts(wildcards):
    """Return every partition-by-bucket statistic file needed for recovery."""
    return [
        GENERATION_CONTEXT.statistics_artifact(partition_id, bucket_id)
        for partition_id in generation_partition_ids(wildcards)
        for bucket_id in generation_bucket_ids(wildcards)
    ]


def generation_outputs(wildcards):
    """Return recoverable outputs of one generation run after certification."""
    outputs = []
    for artifact_name in ("allele_states", "allele_membership", "allele_scores"):
        outputs.extend(generation_partition_artifacts(wildcards, artifact_name))
    outputs.extend(generation_statistics_artifacts(wildcards))
    for artifact_name in ("lookup_het", "lookup_hom", "diagnostics"):
        outputs.extend(generation_bucket_artifacts(wildcards, artifact_name))
    if GENERATION_CONTEXT.hem_enabled:
        outputs.extend(generation_bucket_artifacts(wildcards, "lookup_hem"))
    return outputs
