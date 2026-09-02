"""Static sample-index adapters for apply rules.

The configured PED/sample sheet determines safe wildcard values before DAG
construction. The ordinary validation marker still blocks every expensive rule
until the physical canonical manifest and cohort metadata have been validated.
"""

from sample_manifest import (
    configured_sample_index,
)


SAMPLE_PATTERN = r"[^/]+"
SAMPLE_METADATA_KIND = "samples" if INPUT_BINDINGS.samples else "ped"
SAMPLE_METADATA_PATH = INPUT_BINDINGS.samples or INPUT_BINDINGS.ped
SAMPLE_INDEX = configured_sample_index(SAMPLE_METADATA_PATH, SAMPLE_METADATA_KIND)


def validation_marker(wildcards):
    return rules.require_valid_inputs.output.ok


def validated_manifest_path(wildcards):
    return str(SAMPLE_MANIFEST_PATH)


def sample_sex(wildcards):
    return SAMPLE_INDEX.sample(wildcards.sample).sex


def affected_sample_sex(wildcards):
    return SAMPLE_INDEX.sample(
        wildcards.sample, expected_phenotype="affected"
    ).sex


def sample_ids_for(wildcards, phenotype):
    return SAMPLE_INDEX.sample_ids(phenotype)


def control_annotated_tables(wildcards):
    return expand(
        str(ANNOTATED_DIR / "{sample}.{testing_unit}.tsv.gz"),
        sample=sample_ids_for(wildcards, "unaffected"),
        testing_unit=[TESTING_UNIT],
    )


def control_psap_tables(wildcards):
    return expand(
        str(SCORED_DIR / "{sample}.{testing_unit}.tsv.gz"),
        sample=sample_ids_for(wildcards, "unaffected"),
        testing_unit=[TESTING_UNIT],
    )


def affected_candidate_tables(wildcards):
    return expand(
        str(CANDIDATE_DIR / "{sample}.{testing_unit}.tsv.gz"),
        sample=sample_ids_for(wildcards, "affected"),
        testing_unit=[TESTING_UNIT],
    )


def control_pids(wildcards):
    return sample_ids_for(wildcards, "unaffected")


def control_sexes(wildcards):
    return SAMPLE_INDEX.sample_sexes("unaffected")
