"""Read the configured sample list early enough to build a safe DAG.

The validation rules still publish the canonical sample manifest and block
expensive work until the cohort metadata has been checked.
"""

from sample_manifest import (
    configured_sample_index,
)

SAMPLE_PATTERN = r"[^/]+"
SAMPLE_METADATA_KIND = "samples" if INPUT_BINDINGS.samples else "ped"
SAMPLE_METADATA_PATH = INPUT_BINDINGS.samples or INPUT_BINDINGS.ped
SAMPLE_INDEX = configured_sample_index(SAMPLE_METADATA_PATH, SAMPLE_METADATA_KIND)


def sample_sex(wildcards):
    return SAMPLE_INDEX.sample(wildcards.sample).sex


def affected_sample_sex(wildcards):
    return SAMPLE_INDEX.sample(wildcards.sample, expected_phenotype="affected").sex


def per_sample(directory, phenotype):
    return [
        str(directory / f"{sample}.{TESTING_UNIT}.tsv.gz")
        for sample in SAMPLE_INDEX.sample_ids(phenotype)
    ]


def control_pids(wildcards):
    return SAMPLE_INDEX.sample_ids("unaffected")


def control_sexes(wildcards):
    return SAMPLE_INDEX.sample_sexes("unaffected")
