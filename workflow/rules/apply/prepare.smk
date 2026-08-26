import shlex


PREPARATION_SCRIPT = str(APPLY_SCRIPT_DIR / "prepare_genomic_resource.py")
COMBINE_EXCLUSIONS_SCRIPT = str(APPLY_SCRIPT_DIR / "combine_exclusions.py")

INPUT_VCF_SOURCE = configured_path(config["inputs"]["vcf"])
PREPARED_INPUT_VCF = str(PREPARED_DIR / "cohort.vcf.gz")
WORKFLOW_INPUT_VCF = (
    INPUT_VCF_SOURCE
    if genomic_resource_is_ready(INPUT_VCF_SOURCE)
    else PREPARED_INPUT_VCF
)
WORKFLOW_INPUT_VCF_INDEX = (
    tabix_index(INPUT_VCF_SOURCE)
    if WORKFLOW_INPUT_VCF == INPUT_VCF_SOURCE
    else f"{PREPARED_INPUT_VCF}.tbi"
)

COVERAGE_SOURCE = optional_reference_path("coverage")
COVERAGE_KIND = genomic_resource_kind(COVERAGE_SOURCE) if COVERAGE_SOURCE else ""
PREPARED_COVERAGE = str(
    PREPARED_DIR / f"coverage.{COVERAGE_KIND or 'unknown'}.gz"
)
WORKFLOW_COVERAGE = PREPARED_COVERAGE if COVERAGE_SOURCE else ""
WORKFLOW_COVERAGE_INDEX = f"{PREPARED_COVERAGE}.tbi" if COVERAGE_SOURCE else ""

TESTING_UNIT_SOURCE = configured_reference_path("testing_unit_bed")
PREPARED_TESTING_UNIT = str(PREPARED_DIR / "testing_units.bed.gz")
WORKFLOW_TESTING_UNIT = PREPARED_TESTING_UNIT
WORKFLOW_TESTING_UNIT_INDEX = f"{PREPARED_TESTING_UNIT}.tbi"

CADD_SNV_SOURCE = configured_path(REFERENCE_BINDINGS["cadd"]["snv"])
CADD_INDEL_SOURCE = optional_config_path(REFERENCE_BINDINGS["cadd"]["indel"])
PREPARED_CADD_SNV = str(PREPARED_DIR / "whole_genome_SNVs.tsv.gz")
PREPARED_CADD_INDEL = str(PREPARED_DIR / "InDels.tsv.gz")
# CADD is always staged because an otherwise valid/indexed source can still use
# a different contig spelling (for example, 1 versus chr1) than the cohort VCF.
WORKFLOW_CADD_SNV = PREPARED_CADD_SNV
WORKFLOW_CADD_INDEL = PREPARED_CADD_INDEL
WORKFLOW_CADD_SNV_INDEX = f"{PREPARED_CADD_SNV}.tbi"
WORKFLOW_CADD_INDEL_INDEX = f"{PREPARED_CADD_INDEL}.tbi"

# Ordered list bindings replace comma-delimited exclusions.  The legacy PAR
# resource retains its own role in metadata while sharing the same projection path.
EXCLUSION_SOURCES = [
    configured_path(path) for path in REFERENCE_BINDINGS["exclusions"]
]
EXCLUSION_IDS = [f"{position:03d}" for position in range(len(EXCLUSION_SOURCES))]
EXCLUSION_SOURCE_BY_ID = dict(zip(EXCLUSION_IDS, EXCLUSION_SOURCES))
EXCLUSION_KIND_BY_ID = {
    item: genomic_resource_kind(path)
    for item, path in EXCLUSION_SOURCE_BY_ID.items()
}
PREPARED_EXCLUSION_BY_ID = {
    item: str(PREPARED_DIR / "exclusions" / f"{item}.gz")
    for item in EXCLUSION_IDS
}
WORKFLOW_EXCLUSION_BY_ID = {
    item: (
        source
        if genomic_resource_is_ready(source)
        else PREPARED_EXCLUSION_BY_ID[item]
    )
    for item, source in EXCLUSION_SOURCE_BY_ID.items()
}
WORKFLOW_EXCLUSIONS = [WORKFLOW_EXCLUSION_BY_ID[item] for item in EXCLUSION_IDS]
WORKFLOW_EXCLUSION_INDEXES = [
    (
        tabix_index(EXCLUSION_SOURCE_BY_ID[item])
        if WORKFLOW_EXCLUSION_BY_ID[item] == EXCLUSION_SOURCE_BY_ID[item]
        else f"{PREPARED_EXCLUSION_BY_ID[item]}.tbi"
    )
    for item in EXCLUSION_IDS
]
COMBINED_EXCLUSIONS_VCF = str(PREPARED_DIR / "combined_exclusions.vcf")
COMBINED_EXCLUSIONS = str(PREPARED_DIR / "combined_exclusions.vcf.gz")
COMBINED_EXCLUSIONS_INDEX = f"{COMBINED_EXCLUSIONS}.tbi"


rule prepare_input_vcf:
    input:
        source=INPUT_VCF_SOURCE,
        validated=validation_marker
    output:
        data=ensure(PREPARED_INPUT_VCF, non_empty=True),
        index=ensure(f"{PREPARED_INPUT_VCF}.tbi", non_empty=True)
    params:
        script=PREPARATION_SCRIPT
    threads: 1
    conda:
        str(CONDA_ENV_DIR / "bioinformatics.yaml")
    envmodules:
        "samtools/1.20",
        "python/3.11.8"
    log:
        str(LOG_DIR / "apply" / config["run"]["apply_prefix"] / "prepare" / "input_vcf.log")
    shell:
        """
        python {params.script:q} --input {input.source:q} --output {output.data:q} \
            --kind vcf > {log:q} 2>&1
        """


rule prepare_coverage:
    input:
        source=lambda wc: COVERAGE_SOURCE,
        cohort=WORKFLOW_INPUT_VCF,
        cohort_index=WORKFLOW_INPUT_VCF_INDEX,
        validated=validation_marker
    output:
        data=ensure(PREPARED_COVERAGE, non_empty=True),
        index=ensure(f"{PREPARED_COVERAGE}.tbi", non_empty=True)
    params:
        script=PREPARATION_SCRIPT,
        kind=COVERAGE_KIND
    threads: 1
    conda:
        str(CONDA_ENV_DIR / "bioinformatics.yaml")
    envmodules:
        "samtools/1.20",
        "python/3.11.8"
    log:
        str(LOG_DIR / "apply" / config["run"]["apply_prefix"] / "prepare" / "coverage.log")
    shell:
        """
        python {params.script:q} --input {input.source:q} --output {output.data:q} \
            --kind {params.kind:q} --chromosome-map-vcf {input.cohort:q} \
            > {log:q} 2>&1
        """


rule prepare_testing_unit_bed:
    input:
        source=TESTING_UNIT_SOURCE,
        cohort=WORKFLOW_INPUT_VCF,
        cohort_index=WORKFLOW_INPUT_VCF_INDEX,
        validated=validation_marker
    output:
        data=ensure(PREPARED_TESTING_UNIT, non_empty=True),
        index=ensure(f"{PREPARED_TESTING_UNIT}.tbi", non_empty=True)
    params:
        script=PREPARATION_SCRIPT
    threads: 1
    conda:
        str(CONDA_ENV_DIR / "bioinformatics.yaml")
    envmodules:
        "samtools/1.20",
        "python/3.11.8"
    log:
        str(LOG_DIR / "apply" / config["run"]["apply_prefix"] / "prepare" / "testing_units.log")
    shell:
        """
        python {params.script:q} --input {input.source:q} --output {output.data:q} \
            --kind bed --chromosome-map-vcf {input.cohort:q} > {log:q} 2>&1
        """


rule prepare_cadd_snv:
    input:
        source=CADD_SNV_SOURCE,
        source_index=existing_genomic_indexes(CADD_SNV_SOURCE),
        cohort=WORKFLOW_INPUT_VCF,
        cohort_index=WORKFLOW_INPUT_VCF_INDEX,
        validated=validation_marker
    output:
        data=ensure(PREPARED_CADD_SNV, non_empty=True),
        index=ensure(f"{PREPARED_CADD_SNV}.tbi", non_empty=True)
    params:
        script=PREPARATION_SCRIPT
    threads: 1
    conda:
        str(CONDA_ENV_DIR / "bioinformatics.yaml")
    envmodules:
        "samtools/1.20",
        "python/3.11.8"
    log:
        str(LOG_DIR / "apply" / config["run"]["apply_prefix"] / "prepare" / "cadd_snv.log")
    shell:
        """
        python {params.script:q} --input {input.source:q} --output {output.data:q} \
            --kind tsv --tsv-layout position \
            --chromosome-map-vcf {input.cohort:q} > {log:q} 2>&1
        """


rule prepare_cadd_indel:
    input:
        source=CADD_INDEL_SOURCE,
        source_index=existing_genomic_indexes(CADD_INDEL_SOURCE),
        cohort=WORKFLOW_INPUT_VCF,
        cohort_index=WORKFLOW_INPUT_VCF_INDEX,
        validated=validation_marker
    output:
        data=ensure(PREPARED_CADD_INDEL, non_empty=True),
        index=ensure(f"{PREPARED_CADD_INDEL}.tbi", non_empty=True)
    params:
        script=PREPARATION_SCRIPT
    threads: 1
    conda:
        str(CONDA_ENV_DIR / "bioinformatics.yaml")
    envmodules:
        "samtools/1.20",
        "python/3.11.8"
    log:
        str(LOG_DIR / "apply" / config["run"]["apply_prefix"] / "prepare" / "cadd_indel.log")
    shell:
        """
        python {params.script:q} --input {input.source:q} --output {output.data:q} \
            --kind tsv --tsv-layout position \
            --chromosome-map-vcf {input.cohort:q} > {log:q} 2>&1
        """


rule prepare_exclusion_resource:
    input:
        source=lambda wc: EXCLUSION_SOURCE_BY_ID[wc.item],
        validated=validation_marker
    output:
        data=ensure(str(PREPARED_DIR / "exclusions" / "{item}.gz"), non_empty=True),
        index=ensure(str(PREPARED_DIR / "exclusions" / "{item}.gz.tbi"), non_empty=True)
    params:
        script=PREPARATION_SCRIPT,
        kind=lambda wc: EXCLUSION_KIND_BY_ID[wc.item]
    wildcard_constraints:
        item=r"\d{3}"
    threads: 1
    conda:
        str(CONDA_ENV_DIR / "bioinformatics.yaml")
    envmodules:
        "samtools/1.20",
        "python/3.11.8"
    log:
        str(LOG_DIR / "apply" / config["run"]["apply_prefix"] / "prepare" / "exclusion_{item}.log")
    shell:
        """
        python {params.script:q} --input {input.source:q} --output {output.data:q} \
            --kind {params.kind:q} > {log:q} 2>&1
        """


rule combine_exclusions:
    input:
        cohort=WORKFLOW_INPUT_VCF,
        cohort_index=WORKFLOW_INPUT_VCF_INDEX,
        exclusions=WORKFLOW_EXCLUSIONS,
        indexes=WORKFLOW_EXCLUSION_INDEXES
    output:
        vcf=ensure(COMBINED_EXCLUSIONS_VCF, non_empty=True)
    params:
        script=COMBINE_EXCLUSIONS_SCRIPT,
        source_args=" ".join(
            "--source "
            + shlex.quote(
                f"{EXCLUSION_KIND_BY_ID[item]}:{WORKFLOW_EXCLUSION_BY_ID[item]}"
            )
            for item in EXCLUSION_IDS
        )
    threads: 1
    conda:
        str(CONDA_ENV_DIR / "validation.yaml")
    envmodules:
        "python/3.11.8"
    log:
        str(LOG_DIR / "apply" / config["run"]["apply_prefix"] / "prepare" / "combine_exclusions.log")
    shell:
        """
        python {params.script:q} \
            --cohort-vcf {input.cohort:q} \
            {params.source_args} \
            --output {output.vcf:q} > {log:q} 2>&1
        """


rule prepare_combined_exclusions:
    input:
        vcf=rules.combine_exclusions.output.vcf
    output:
        vcf=ensure(COMBINED_EXCLUSIONS, non_empty=True),
        index=ensure(COMBINED_EXCLUSIONS_INDEX, non_empty=True)
    params:
        script=PREPARATION_SCRIPT
    threads: 1
    conda:
        str(CONDA_ENV_DIR / "bioinformatics.yaml")
    envmodules:
        "samtools/1.20",
        "python/3.11.8"
    log:
        str(LOG_DIR / "apply" / config["run"]["apply_prefix"] / "prepare" / "combined_exclusions.log")
    shell:
        """
        python {params.script:q} --input {input.vcf:q} --output {output.vcf:q} \
            --kind vcf > {log:q} 2>&1
        """
