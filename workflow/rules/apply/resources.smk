"""Expose the run's shared resource catalog to apply rules.

The catalog is built once from ``WORKFLOW_PLAN``. Validation and staging then
use the same paths instead of resolving a resource independently.
"""

from apply_resources import build_apply_resource_registry
from execution_policy import resource_mem_mb, resource_runtime_minutes


APPLY_RESOURCE_REGISTRY = build_apply_resource_registry(WORKFLOW_PLAN, PREPARED_DIR)
RESOURCE_SPECS = APPLY_RESOURCE_REGISTRY.resources
INPUT_VCF_SOURCE = APPLY_RESOURCE_REGISTRY.input_vcf_source
PREPARED_INPUT_VCF = APPLY_RESOURCE_REGISTRY.prepared_input_vcf
WORKFLOW_INPUT_VCF = PREPARED_INPUT_VCF
WORKFLOW_INPUT_VCF_INDEX = f"{PREPARED_INPUT_VCF}.tbi"
PREPARED_TABLE_DIR = APPLY_RESOURCE_REGISTRY.prepared_table_dir
TABLE_RESOURCE_SPECS = APPLY_RESOURCE_REGISTRY.table_resources
TABLE_RESOURCE_IDS = APPLY_RESOURCE_REGISTRY.table_resource_ids
TABLE_RESOURCE_PATTERN = APPLY_RESOURCE_REGISTRY.table_resource_pattern
TABLE_RESOURCE_DATA = APPLY_RESOURCE_REGISTRY.table_resource_data
TABLE_RESOURCE_INDEX = APPLY_RESOURCE_REGISTRY.table_resource_index
EXCLUSION_IDS = APPLY_RESOURCE_REGISTRY.exclusion_ids
EXCLUSION_SOURCE_BY_ID = APPLY_RESOURCE_REGISTRY.exclusion_source_by_id
EXCLUSION_KIND_BY_ID = APPLY_RESOURCE_REGISTRY.exclusion_kind_by_id

WORKFLOW_COVERAGE = TABLE_RESOURCE_DATA.get("coverage", "")
WORKFLOW_COVERAGE_INDEX = TABLE_RESOURCE_INDEX.get("coverage", "")
WORKFLOW_TESTING_UNIT = TABLE_RESOURCE_DATA["testing_units"]
WORKFLOW_TESTING_UNIT_INDEX = TABLE_RESOURCE_INDEX["testing_units"]
WORKFLOW_CADD_SNV = TABLE_RESOURCE_DATA["cadd_snv"]
WORKFLOW_CADD_SNV_INDEX = TABLE_RESOURCE_INDEX["cadd_snv"]
WORKFLOW_CADD_INDEL = TABLE_RESOURCE_DATA.get("cadd_indel", "")
WORKFLOW_CADD_INDEL_INDEX = TABLE_RESOURCE_INDEX.get("cadd_indel", "")
WORKFLOW_EXCLUSIONS = list(APPLY_RESOURCE_REGISTRY.workflow_exclusions)
WORKFLOW_EXCLUSION_INDEXES = list(APPLY_RESOURCE_REGISTRY.workflow_exclusion_indexes)
EXCLUSION_SOURCES = list(APPLY_RESOURCE_REGISTRY.exclusion_sources)

COMBINED_EXCLUSIONS_VCF = APPLY_RESOURCE_REGISTRY.combined_exclusions_vcf
COMBINED_EXCLUSIONS = APPLY_RESOURCE_REGISTRY.combined_exclusions
COMBINED_EXCLUSIONS_INDEX = APPLY_RESOURCE_REGISTRY.combined_exclusions_index
VCF_RESOURCE_IDS = ("cohort", "combined_exclusions")
VCF_RESOURCE_PATTERN = "(?:cohort|combined_exclusions)"


def table_resource_source(wildcards):
    return TABLE_RESOURCE_SPECS[wildcards.resource].source


def table_resource_source_index(wildcards):
    index = TABLE_RESOURCE_SPECS[wildcards.resource].source_index
    return [index] if index else []


def table_resource_cohort(wildcards):
    resource_spec = TABLE_RESOURCE_SPECS[wildcards.resource]
    return [WORKFLOW_INPUT_VCF] if resource_spec.requires_chromosome_normalization else []


def table_resource_cohort_index(wildcards):
    resource_spec = TABLE_RESOURCE_SPECS[wildcards.resource]
    return [WORKFLOW_INPUT_VCF_INDEX] if resource_spec.requires_chromosome_normalization else []


def table_resource_mem_mb(wildcards, input, attempt):
    return resource_mem_mb("stage_" + TABLE_RESOURCE_SPECS[wildcards.resource].scientific_role, input, attempt)


def table_resource_runtime_minutes(wildcards, input, attempt):
    return resource_runtime_minutes("stage_" + TABLE_RESOURCE_SPECS[wildcards.resource].scientific_role, input, attempt)


def table_resource_kind(wildcards):
    return TABLE_RESOURCE_SPECS[wildcards.resource].kind


def table_resource_tsv_layout(wildcards):
    return TABLE_RESOURCE_SPECS[wildcards.resource].tsv_layout or "auto"


def vcf_resource_source(wildcards):
    return APPLY_RESOURCE_REGISTRY.resource(wildcards.resource).source


def vcf_resource_source_index(wildcards):
    index = APPLY_RESOURCE_REGISTRY.resource(wildcards.resource).source_index
    return [index] if index else []


def vcf_resource_validation_marker(wildcards):
    return [validation_marker(wildcards)] if wildcards.resource == "cohort" else []
