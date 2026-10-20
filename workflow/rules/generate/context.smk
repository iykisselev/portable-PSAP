"""Shared paths and wildcard helpers for generation rules."""

from generation_context import GenerationRuleContext
from generation_manifests import (
    bucket_from_manifest,
    bucket_ids,
    window_ids,
    read_bucket_manifest,
)

GENERATION_CONTEXT = GenerationRuleContext.from_workflow_plan(WORKFLOW_PLAN)
GENERATION_SCRIPT_DIR = Path(workflow.basedir) / "scripts" / "generate"


def hem_artifact(value):
    """Retain a named optional hem binding while disabling its file dependency."""
    return value if GENERATION_CONTEXT.model.hemizygote else []


def generation_scope_artifact(wildcards, artifact_name):
    """Return one scope file after the scope checkpoint has completed."""
    scope_checkpoint = checkpoints.prepare_generation_scope.get()
    return str(getattr(scope_checkpoint.output, artifact_name))


def generation_window_ids(wildcards):
    """Read window IDs from the completed scope manifest."""
    return window_ids(generation_scope_artifact(wildcards, "windows"))


def generation_bucket_ids(wildcards):
    """Read unit-bucket IDs from the completed scope manifest."""
    return bucket_ids(generation_scope_artifact(wildcards, "buckets"))


def generation_window_artifacts(wildcards, artifact_name):
    """Return the requested file for every discovered genomic window."""
    return [
        GENERATION_CONTEXT.window_artifact(window_id, artifact_name)
        for window_id in generation_window_ids(wildcards)
    ]


def generation_bucket_artifacts(wildcards, artifact_name):
    """Return the requested lookup file for every discovered unit bucket."""
    if artifact_name == "lookup_hem_chrX" and not GENERATION_CONTEXT.model.hemizygote:
        return []
    return [
        GENERATION_CONTEXT.bucket_artifact(bucket_id, artifact_name)
        for bucket_id in generation_bucket_ids(wildcards)
    ]


def generation_bucket_partition_model_artifacts(wildcards):
    return [
        GENERATION_CONTEXT.bucket_partition_model_artifact(bucket_id, partition_id, model_id)
        for bucket_id in generation_bucket_ids(wildcards)
        for partition_id, model_id in GENERATION_CONTEXT.execution_layout.partition_lookup_paths
    ]


def generation_statistics_shards(wildcards):
    """Return the window statistics needed to reduce one unit bucket."""
    bucket_id = str(wildcards.bucket)
    bucket = bucket_from_manifest(
        generation_scope_artifact(wildcards, "buckets"), bucket_id
    )
    return [
        GENERATION_CONTEXT.statistics_artifact(window_id, bucket_id)
        for window_id in bucket.window_ids
    ]


def generation_statistics_artifacts(wildcards):
    """Return every window-by-bucket statistic file needed for recovery."""
    return [
        GENERATION_CONTEXT.statistics_artifact(window_id, bucket.identifier)
        for bucket in read_bucket_manifest(
            generation_scope_artifact(wildcards, "buckets")
        )
        for window_id in bucket.window_ids
    ]


def generation_outputs(wildcards):
    """Return recoverable outputs of one generation run after certification."""
    outputs = []
    for artifact_name in (
        "allele_states",
        "allele_membership",
        "allele_scores",
        "candidate_audit",
    ):
        outputs.extend(generation_window_artifacts(wildcards, artifact_name))
    outputs.extend(generation_statistics_artifacts(wildcards))
    for artifact_name in (
        "lookup_het", "lookup_hom", "lookup_hem_chrX", "diagnostics",
        "partition_model_diagnostics",
    ):
        outputs.extend(generation_bucket_artifacts(wildcards, artifact_name))
    return outputs
