from pathlib import Path


VALIDATION_DIR = Path(APPLY_OUTDIR) / "validation"
COHORT_VALIDATION = VALIDATION_DIR / "cohort.json"
RESOURCE_VALIDATION = VALIDATION_DIR / "resources.json"
MODEL_VALIDATION = VALIDATION_DIR / "model.json"
VALIDATION_REPORT = Path(APPLY_OUTDIR) / "validation.json"
VALIDATION_LIBRARY = Path(APPLY_SCRIPT_DIR) / "validate_inputs.py"

VALIDATION_VCF = configured_path(config["inputs"]["vcf"])
VALIDATION_PED = configured_path(config["inputs"]["ped"])
VALIDATION_FASTA = configured_path(config["reference"]["fasta"])
VALIDATION_FASTA_INDEX = f"{VALIDATION_FASTA}.fai"
VALIDATION_FASTA_GZI = (
    [f"{VALIDATION_FASTA}.gzi"] if VALIDATION_FASTA.endswith(".gz") else []
)

VALIDATION_COVERAGE = optional_config_path(config["reference"].get("coverage"))
VALIDATION_COVERAGE_INPUT = [VALIDATION_COVERAGE] if VALIDATION_COVERAGE else []
VALIDATION_TESTING_UNIT_BED = configured_path(
    config["reference"]["testing_unit_bed"]
)
VALIDATION_VARIANT_EXCLUSIONS = comma_separated_config_paths(
    config["reference"].get("variants_exclude")
)
VALIDATION_PAR_REGIONS = comma_separated_config_paths(
    config["reference"].get("par_regions")
)

VALIDATION_CADD_DIR = Path(configured_path(config["reference"]["cadd_dir"]))
VALIDATION_CADD_SNV = str(VALIDATION_CADD_DIR / "whole_genome_SNVs.tsv.gz")
VALIDATION_CADD_INDEL_INPUT = (
    [str(VALIDATION_CADD_DIR / "InDels.tsv.gz")]
    if config["model"]["include_indels"]
    else []
)

VALIDATION_LOOKUP_HET = configured_path(
    config["reference"]["lookup_tables"]["het"]
)
VALIDATION_LOOKUP_HOM = configured_path(
    config["reference"]["lookup_tables"]["hom"]
)
VALIDATION_LOOKUP_HEM_CHRX_INPUT = (
    [configured_path(config["reference"]["lookup_tables"]["hem_chrX"])]
    if config["model"]["hemizygote"]
    else []
)

VALIDATION_INDEXED_RESOURCES = (
    [VALIDATION_VCF, VALIDATION_TESTING_UNIT_BED, VALIDATION_CADD_SNV]
    + VALIDATION_COVERAGE_INPUT
    + VALIDATION_VARIANT_EXCLUSIONS
    + VALIDATION_PAR_REGIONS
    + VALIDATION_CADD_INDEL_INPUT
)


rule validate_cohort:
    """Validate cohort VCF/PED/FASTA inputs and write the sample manifest."""
    input:
        vcf=VALIDATION_VCF,
        vcf_indexes=existing_genomic_indexes(VALIDATION_VCF),
        ped=VALIDATION_PED,
        fasta=VALIDATION_FASTA,
        fasta_index=VALIDATION_FASTA_INDEX,
        fasta_gzi=VALIDATION_FASTA_GZI,
        script=str(Path(APPLY_SCRIPT_DIR) / "validate_cohort.py"),
        library=str(VALIDATION_LIBRARY)
    output:
        report=ensure(str(COHORT_VALIDATION), non_empty=True),
        manifest=protected(ensure(str(SAMPLE_MANIFEST_PATH), non_empty=True))
    params:
        hemizygote_model=str(config["model"]["hemizygote"])
    threads: 1
    conda:
        str(CONDA_ENV_DIR / "validation.yaml")
    envmodules:
        "samtools/1.20",
        "python/3.11.8"
    log:
        str(LOG_DIR / "apply" / config["run"]["apply_prefix"] / "validation" / "cohort.log")
    shell:
        """
        python {input.script:q} \
            --vcf {input.vcf:q} \
            --ped {input.ped:q} \
            --fasta {input.fasta:q} \
            --hemizygote-model {params.hemizygote_model:q} \
            --out-json {output.report:q} \
            --out-manifest {output.manifest:q} \
            > {log:q} 2>&1
        """


rule validate_resources:
    """Validate interval, exclusion, testing-unit, and CADD resources."""
    input:
        vcf=VALIDATION_VCF,
        coverage=VALIDATION_COVERAGE_INPUT,
        variants_exclude=VALIDATION_VARIANT_EXCLUSIONS,
        testing_unit_bed=VALIDATION_TESTING_UNIT_BED,
        par_regions=VALIDATION_PAR_REGIONS,
        cadd_snv=VALIDATION_CADD_SNV,
        cadd_indel=VALIDATION_CADD_INDEL_INPUT,
        indexes=existing_genomic_indexes(VALIDATION_INDEXED_RESOURCES),
        script=str(Path(APPLY_SCRIPT_DIR) / "validate_resources.py"),
        library=str(VALIDATION_LIBRARY)
    output:
        report=ensure(str(RESOURCE_VALIDATION), non_empty=True)
    params:
        coverage=lambda wc, input: (
            str(input.coverage[0]) if input.coverage else "None"
        ),
        variants_exclude=lambda wc, input: (
            ",".join(map(str, input.variants_exclude))
            if input.variants_exclude
            else "None"
        ),
        par_regions=lambda wc, input: (
            ",".join(map(str, input.par_regions))
            if input.par_regions
            else "None"
        ),
        cadd_path=lambda wc, input: str(Path(str(input.cadd_snv)).parent),
        unit=str(config["model"]["unit"]),
        include_indels=str(config["model"]["include_indels"])
    threads: 1
    conda:
        str(CONDA_ENV_DIR / "validation.yaml")
    envmodules:
        "samtools/1.20",
        "python/3.11.8"
    log:
        str(LOG_DIR / "apply" / config["run"]["apply_prefix"] / "validation" / "resources.log")
    shell:
        """
        python {input.script:q} \
            --vcf {input.vcf:q} \
            --coverage {params.coverage:q} \
            --variants-exclude {params.variants_exclude:q} \
            --testing-unit-bed {input.testing_unit_bed:q} \
            --par-regions {params.par_regions:q} \
            --cadd-path {params.cadd_path:q} \
            --unit {params.unit:q} \
            --include-indels {params.include_indels:q} \
            --out-json {output.report:q} \
            > {log:q} 2>&1
        """


rule validate_model:
    """Validate model configuration and enabled lookup resources."""
    input:
        testing_unit_bed=VALIDATION_TESTING_UNIT_BED,
        lookup_het=VALIDATION_LOOKUP_HET,
        lookup_hom=VALIDATION_LOOKUP_HOM,
        lookup_hem_chrx=VALIDATION_LOOKUP_HEM_CHRX_INPUT,
        script=str(Path(APPLY_SCRIPT_DIR) / "validate_model.py"),
        library=str(VALIDATION_LIBRARY)
    output:
        report=ensure(str(MODEL_VALIDATION), non_empty=True)
    params:
        lookup_hem_chrx=lambda wc, input: (
            str(input.lookup_hem_chrx[0]) if input.lookup_hem_chrx else "None"
        ),
        score_max=str(config["model"]["score_max"]),
        assembly=str(config["reference"]["assembly"]),
        unit=str(config["model"]["unit"]),
        hemizygote_model=str(config["model"]["hemizygote"]),
        include_indels=str(config["model"]["include_indels"])
    threads: 1
    conda:
        str(CONDA_ENV_DIR / "validation.yaml")
    envmodules:
        "python/3.11.8"
    log:
        str(LOG_DIR / "apply" / config["run"]["apply_prefix"] / "validation" / "model.log")
    shell:
        """
        python {input.script:q} \
            --testing-unit-bed {input.testing_unit_bed:q} \
            --lookup-het {input.lookup_het:q} \
            --lookup-hom {input.lookup_hom:q} \
            --lookup-hem-chrx {params.lookup_hem_chrx:q} \
            --score-max {params.score_max:q} \
            --assembly {params.assembly:q} \
            --unit {params.unit:q} \
            --hemizygote-model {params.hemizygote_model:q} \
            --include-indels {params.include_indels:q} \
            --out-json {output.report:q} \
            > {log:q} 2>&1
        """


rule merge_validation:
    """Merge independent validation fragments into validation.json."""
    input:
        cohort=rules.validate_cohort.output.report,
        resources=rules.validate_resources.output.report,
        model=rules.validate_model.output.report,
        script=str(Path(APPLY_SCRIPT_DIR) / "merge_validation.py")
    output:
        report=protected(ensure(str(VALIDATION_REPORT), non_empty=True))
    threads: 1
    conda:
        str(CONDA_ENV_DIR / "validation.yaml")
    envmodules:
        "python/3.11.8"
    log:
        str(LOG_DIR / "apply" / config["run"]["apply_prefix"] / "validation" / "merge.log")
    shell:
        """
        python {input.script:q} \
            --fragment {input.cohort:q} \
            --fragment {input.resources:q} \
            --fragment {input.model:q} \
            --output {output.report:q} \
            > {log:q} 2>&1
        """


checkpoint require_valid_inputs:
    """Gate expensive workflow rules on the merged validation report."""
    input:
        report=rules.merge_validation.output.report,
        manifest=rules.validate_cohort.output.manifest,
        script=str(Path(APPLY_SCRIPT_DIR) / "require_valid_inputs.py")
    output:
        ok=ensure(str(Path(APPLY_OUTDIR) / "validation.ok"), non_empty=True)
    threads: 1
    conda:
        str(CONDA_ENV_DIR / "validation.yaml")
    envmodules:
        "python/3.11.8"
    log:
        str(LOG_DIR / "apply" / config["run"]["apply_prefix"] / "validation" / "require_valid_inputs.log")
    shell:
        """
        python {input.script:q} \
            --report {input.report:q} \
            --output {output.ok:q} \
            > {log:q} 2>&1
        """
