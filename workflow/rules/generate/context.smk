"""Shared paths and wildcard helpers for generation rules."""

from generation_context import GenerationRuleContext, as_inputs
from generation_manifests import (
    bucket_from_manifest,
    bucket_ids,
    partition_ids,
    read_bucket_manifest,
)


GENERATION_CONTEXT = GenerationRuleContext.from_workflow_plan(WORKFLOW_PLAN)
GENERATION_SCRIPT_DIR = Path(workflow.basedir) / "scripts" / "generate"


def generation_index(resource_binding):
    return as_inputs(resource_binding)


def optional_binding(binding, *, enabled=True):
    """Return optional data/index Snakemake bindings without placeholder paths."""
    if not enabled or binding is None or not binding.data:
        return [], []
    return [binding.data], ([binding.index] if binding.index else [])


SCORE_INDEL, SCORE_INDEL_INDEX = optional_binding(
    GENERATION_CONTEXT.plan.score_indel,
    enabled=WORKFLOW_PLAN.model.include_indels,
)


def hem_artifact(value):
    """Retain a named optional hem binding while disabling its file dependency."""
    return value if GENERATION_CONTEXT.hem_enabled else []


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
    if artifact_name == "lookup_hem" and not GENERATION_CONTEXT.hem_enabled:
        return []
    return [
        GENERATION_CONTEXT.bucket_artifact(bucket_id, artifact_name)
        for bucket_id in generation_bucket_ids(wildcards)
    ]


def generation_statistics_shards(wildcards):
    """Return the partition statistics needed to reduce one unit bucket."""
    bucket_id = str(wildcards.bucket)
    bucket = bucket_from_manifest(
        generation_scope_artifact(wildcards, "buckets"), bucket_id
    )
    return [
        GENERATION_CONTEXT.statistics_artifact(partition_id, bucket_id)
        for partition_id in bucket.partition_ids
    ]


def generation_statistics_artifacts(wildcards):
    """Return every partition-by-bucket statistic file needed for recovery."""
    return [
        GENERATION_CONTEXT.statistics_artifact(partition_id, bucket.identifier)
        for bucket in read_bucket_manifest(
            generation_scope_artifact(wildcards, "buckets")
        )
        for partition_id in bucket.partition_ids
    ]


def generation_outputs(wildcards):
    """Return recoverable outputs of one generation run after certification."""
    outputs = []
    for artifact_name in ("allele_states", "allele_membership", "allele_scores", "candidate_audit"):
        outputs.extend(generation_partition_artifacts(wildcards, artifact_name))
    outputs.extend(generation_statistics_artifacts(wildcards))
    for artifact_name in ("lookup_het", "lookup_hom", "lookup_hem", "diagnostics"):
        outputs.extend(generation_bucket_artifacts(wildcards, artifact_name))
    return outputs
