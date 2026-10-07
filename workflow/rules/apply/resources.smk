"""Expose the run's shared resource catalog to apply rules.

The catalog is built once from ``WORKFLOW_PLAN``. Validation and staging then
use the same paths instead of resolving a resource independently.
"""

from apply_resources import build_apply_resource_registry
from execution_policy import table_resource_policy


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


def optional_resource_field(name, field):
    """Read an explicit field from an optional resource without a mirror global."""
    resource_spec = optional_apply_resource(name)
    return getattr(resource_spec, field) if resource_spec else ""


def optional_input(value):
    """Represent one optional scalar as a Snakemake input list."""
    return [value] if value else []


def resources_with_role(role):
    """Return resources in their configured order for one scientific role."""
    return [
        resource_spec
        for resource_spec in APPLY_RESOURCE_REGISTRY.resources.values()
        if resource_spec.scientific_role == role
    ]


def resource_binding(resource_spec):
    """Keep a resource path and optional index associated for script adapters."""
    return (resource_spec.source, resource_spec.source_index or None)


def apply_resource_binding(name):
    """Return one required resource as a path/index pair."""
    return resource_binding(apply_resource(name))


def optional_apply_resource_binding(name):
    """Return an enabled optional resource as a path/index pair."""
    resource_spec = optional_apply_resource(name)
    return resource_binding(resource_spec) if resource_spec else None


def resource_bindings_with_role(role):
    """Return ordered bindings without reconstructing indexes by position."""
    return [resource_binding(resource_spec) for resource_spec in resources_with_role(role)]


def table_resource(wildcards):
    return APPLY_RESOURCE_REGISTRY.table_resources[wildcards.resource]


def table_resource_source(wildcards):
    return table_resource(wildcards).source


def table_resource_source_index(wildcards):
    return optional_input(table_resource(wildcards).source_index)


def table_resource_cohort(wildcards):
    return (
        optional_input(apply_resource("cohort").staged_data)
        if table_resource(wildcards).requires_chromosome_normalization
        else []
    )


def table_resource_cohort_index(wildcards):
    return (
        optional_input(apply_resource("cohort").staged_index)
        if table_resource(wildcards).requires_chromosome_normalization
        else []
    )


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
    return [index] if index else []


def vcf_resource_validation_marker(wildcards):
    return [validation_marker(wildcards)] if wildcards.resource == "cohort" else []
