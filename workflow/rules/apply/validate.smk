from pathlib import Path


VALIDATION_DIR = Path(APPLY_OUTDIR) / "validation"
COHORT_VALIDATION = VALIDATION_DIR / "cohort.json"
RESOURCE_VALIDATION = VALIDATION_DIR / "resources.json"
MODEL_VALIDATION = VALIDATION_DIR / "model.json"
VALIDATION_REPORT = Path(APPLY_OUTDIR) / "validation.json"
VALIDATION_LIBRARY = APPLY_LIBRARY_DIR / "validate_inputs.py"

VALIDATION_VCF = INPUT_VCF_SOURCE
VALIDATION_VCF_INDEX = [INPUT_VCF_SOURCE_INDEX] if INPUT_VCF_SOURCE_INDEX else []
VALIDATION_PED = optional_config_path(config["inputs"].get("ped"))
VALIDATION_SAMPLES = optional_config_path(config["inputs"].get("samples"))
VALIDATION_FASTA = configured_reference_path("fasta")
VALIDATION_FASTA_INDEX = configured_path(REFERENCE_BINDINGS["fasta"]["fai"])
VALIDATION_FASTA_GZI = [
    configured_path(REFERENCE_BINDINGS["fasta"]["gzi"])
] if REFERENCE_BINDINGS["fasta"]["gzi"] is not None else []

VALIDATION_COVERAGE = optional_reference_path("coverage")
VALIDATION_COVERAGE_INPUT = [VALIDATION_COVERAGE] if VALIDATION_COVERAGE else []
VALIDATION_TESTING_UNIT_BED = configured_reference_path("testing_unit_bed")
VALIDATION_VARIANT_EXCLUSIONS = [
    binding_data(path)
    for resource, path in zip(
        config["reference"]["exclusions"], REFERENCE_BINDINGS["exclusions"]
    )
    if resource["role"] == "variants_exclude"
]
VALIDATION_PAR_REGIONS = [
    binding_data(path)
    for resource, path in zip(
        config["reference"]["exclusions"], REFERENCE_BINDINGS["exclusions"]
    )
    if resource["role"] == "par_regions"
]

VALIDATION_CADD_SNV = binding_data(REFERENCE_BINDINGS["cadd"]["snv"])
VALIDATION_CADD_INDEL_INPUT = (
    [binding_data(REFERENCE_BINDINGS["cadd"]["indel"])]
    if config["model"]["include_indels"]
    else []
)

VALIDATION_LOOKUP_HET = configured_path(
    REFERENCE_BINDINGS["lookup_tables"]["het"]
)
VALIDATION_LOOKUP_HOM = configured_path(
    REFERENCE_BINDINGS["lookup_tables"]["hom"]
)
VALIDATION_LOOKUP_HEM_CHRX_INPUT = (
    [configured_path(REFERENCE_BINDINGS["lookup_tables"]["hem_chrX"])]
    if config["model"]["hemizygote"]
    else []
)

VALIDATION_CONFIGURED_INDEXES = [
    index
    for index in [
        INPUT_VCF_SOURCE_INDEX,
        binding_index(REFERENCE_BINDINGS["coverage"]),
        binding_index(REFERENCE_BINDINGS["testing_unit_bed"]),
        *(binding_index(binding) for binding in REFERENCE_BINDINGS["exclusions"]),
        binding_index(REFERENCE_BINDINGS["cadd"]["snv"]),
        binding_index(REFERENCE_BINDINGS["cadd"]["indel"]),
    ]
    if index
]


rule validate_cohort:
    """Validate cohort VCF/metadata/FASTA inputs and write the sample manifest."""
    input:
        vcf=VALIDATION_VCF,
        vcf_indexes=VALIDATION_VCF_INDEX,
        metadata=VALIDATION_SAMPLES or VALIDATION_PED,
        fasta=VALIDATION_FASTA,
        fasta_index=VALIDATION_FASTA_INDEX,
        fasta_gzi=VALIDATION_FASTA_GZI,
        library=str(VALIDATION_LIBRARY)
    output:
        report=ensure(str(COHORT_VALIDATION), non_empty=True),
        manifest=protected(ensure(str(SAMPLE_MANIFEST_PATH), non_empty=True))
    params:
        hemizygote_model=str(config["model"]["hemizygote"]),
        metadata_kind="samples" if VALIDATION_SAMPLES else "ped"
    threads: 1
    resources:
        mem_mb=4096,
        runtime=30
    conda:
        str(CONDA_ENV_DIR / "validation.yaml")
    envmodules: mod("validation")
    container: sif("validation")
    log:
        str(LOG_DIR / "apply" / config["run"]["apply_prefix"] / "validation" / "cohort.log")
    script:
        str(Path(workflow.basedir) / "scripts" / "validate_cohort.py")


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
        indexes=VALIDATION_CONFIGURED_INDEXES,
        library=str(VALIDATION_LIBRARY)
    output:
        report=ensure(str(RESOURCE_VALIDATION), non_empty=True)
    params:
        unit=str(config["model"]["unit"]),
        include_indels=str(config["model"]["include_indels"])
    threads: 1
    resources:
        mem_mb=4096,
        runtime=30
    conda:
        str(CONDA_ENV_DIR / "validation.yaml")
    envmodules: mod("validation")
    container: sif("validation")
    log:
        str(LOG_DIR / "apply" / config["run"]["apply_prefix"] / "validation" / "resources.log")
    script:
        str(Path(workflow.basedir) / "scripts" / "validate_resources.py")


rule validate_model:
    """Validate model configuration and enabled lookup resources."""
    input:
        testing_unit_bed=VALIDATION_TESTING_UNIT_BED,
        lookup_het=VALIDATION_LOOKUP_HET,
        lookup_hom=VALIDATION_LOOKUP_HOM,
        lookup_hem_chrx=VALIDATION_LOOKUP_HEM_CHRX_INPUT,
        library=str(VALIDATION_LIBRARY)
    output:
        report=ensure(str(MODEL_VALIDATION), non_empty=True)
    params:
        score_max=str(config["model"]["score_max"]),
        assembly=str(config["reference"]["assembly"]),
        unit=str(config["model"]["unit"]),
        hemizygote_model=str(config["model"]["hemizygote"]),
        include_indels=str(config["model"]["include_indels"])
    threads: 1
    resources:
        mem_mb=512,
        runtime=10
    conda:
        str(CONDA_ENV_DIR / "validation.yaml")
    envmodules: mod("validation")
    container: sif("validation")
    log:
        str(LOG_DIR / "apply" / config["run"]["apply_prefix"] / "validation" / "model.log")
    script:
        str(Path(workflow.basedir) / "scripts" / "validate_model.py")


rule merge_validation:
    """Merge independent validation fragments into validation.json."""
    input:
        cohort=rules.validate_cohort.output.report,
        resources=rules.validate_resources.output.report,
        model=rules.validate_model.output.report
    output:
        report=protected(ensure(str(VALIDATION_REPORT), non_empty=True))
    threads: 1
    resources:
        mem_mb=512,
        runtime=5
    conda:
        str(CONDA_ENV_DIR / "validation.yaml")
    envmodules: mod("validation")
    container: sif("validation")
    log:
        str(LOG_DIR / "apply" / config["run"]["apply_prefix"] / "validation" / "merge.log")
    script:
        str(Path(workflow.basedir) / "scripts" / "merge_validation.py")


checkpoint require_valid_inputs:
    """Gate expensive workflow rules on the merged validation report."""
    input:
        report=rules.merge_validation.output.report,
        manifest=rules.validate_cohort.output.manifest
    output:
        ok=ensure(str(Path(APPLY_OUTDIR) / "validation.ok"), non_empty=True)
    params:
        failure_message="Input validation failed",
        include_errors=True,
        success_marker="validated\n"
    threads: 1
    resources:
        mem_mb=512,
        runtime=5
    conda:
        str(CONDA_ENV_DIR / "validation.yaml")
    envmodules: mod("validation")
    container: sif("validation")
    log:
        str(LOG_DIR / "apply" / config["run"]["apply_prefix"] / "validation" / "require_valid_inputs.log")
    script:
        str(Path(workflow.basedir) / "scripts" / "require_valid_report.py")
