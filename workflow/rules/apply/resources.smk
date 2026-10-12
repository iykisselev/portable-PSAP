"""Expose the run's shared resource catalog to apply rules.

The catalog is built once from ``WORKFLOW_PLAN``. Validation and staging then
use the same paths instead of resolving a resource independently.
"""

from apply_resources import build_apply_resource_registry
from execution_policy import table_resource_policy
from workflow_bindings import as_inputs
from functools import partial

APPLY_RESOURCE_REGISTRY = build_apply_resource_registry(WORKFLOW_PLAN, PREPARED_DIR)
TABLE_RESOURCE_PATTERN = APPLY_RESOURCE_REGISTRY.table_resource_pattern
TABLE_RESOURCE_POLICY_BY_ID = {
    resource_id: table_resource_policy(resource_spec.scientific_role)
    for resource_id, resource_spec in APPLY_RESOURCE_REGISTRY.table_resources.items()
}
VCF_RESOURCE_PATTERN = "cohort"


def apply_resource(name):
    """Return one required resource without flattening its fields into globals."""
    return APPLY_RESOURCE_REGISTRY.resource(name)


def optional_apply_resource(name):
    """Return one optional resource, if this run enables it."""
    return APPLY_RESOURCE_REGISTRY.resources.get(name)


def table_resource(wildcards):
    return APPLY_RESOURCE_REGISTRY.table_resources[wildcards.resource]


def table_resource_source(wildcards):
    return table_resource(wildcards).source


def table_resource_source_index(wildcards):
    return as_inputs(table_resource(wildcards).source_index)


def table_resource_cohort(wildcards, *, field="staged_data"):
    return (
        APPLY_RESOURCE_REGISTRY.inputs("cohort", field)
        if table_resource(wildcards).requires_chromosome_normalization
        else []
    )


table_resource_cohort_index = partial(table_resource_cohort, field="staged_index")


def table_resource_policy_name(wildcards):
    return TABLE_RESOURCE_POLICY_BY_ID[wildcards.resource]


TABLE_RESOURCE_CALLBACKS = policy_resources(table_resource_policy_name)


def table_resource_kind(wildcards):
    return table_resource(wildcards).kind


def table_resource_tsv_layout(wildcards):
    return table_resource(wildcards).tsv_layout or "auto"


def vcf_resource_source(wildcards):
    return APPLY_RESOURCE_REGISTRY.resource(wildcards.resource).source


def vcf_resource_source_index(wildcards):
    index = APPLY_RESOURCE_REGISTRY.resource(wildcards.resource).source_index
    return as_inputs(index)


def vcf_resource_validation_marker(wildcards):
    return (
        [rules.require_valid_inputs.output.ok] if wildcards.resource == "cohort" else []
    )
