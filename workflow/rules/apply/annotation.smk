CADD_SNV = WORKFLOW_CADD_SNV
CADD_INDEL = WORKFLOW_CADD_INDEL
INCLUDE_INDELS = WORKFLOW_PLAN.model.include_indels
TESTING_UNIT_BED = WORKFLOW_TESTING_UNIT
ANNOTATED_VCF = ANNOTATION_DIR / f"{TESTING_UNIT}.vcf.gz"
ANNOTATED_VCF_CSI = Path(f"{ANNOTATED_VCF}.csi")
RAW_UNIT_ANNOTATIONS = ANNOTATION_DIR / f"{TESTING_UNIT}.annotated.tsv.gz"
ANNOTATION_QC_REPORT = ANNOTATION_DIR / f"{TESTING_UNIT}.annotation_qc.json"
ANNOTATION_QC_OK = ANNOTATION_DIR / f"{TESTING_UNIT}.annotation_qc.ok"


rule annotate_variants:
    """
    Attach CADD and the configured testing-unit annotation to the filtered VCF.
    """
    input:
        vcf=rules.filter_excluded_records.output.vcf,
        vcf_csi=rules.filter_excluded_records.output.csi,
        cadd_snv=str(CADD_SNV),
        cadd_snv_index=WORKFLOW_CADD_SNV_INDEX,
        cadd_indel=[str(CADD_INDEL)] if INCLUDE_INDELS else [],
        cadd_indel_index=(
            [WORKFLOW_CADD_INDEL_INDEX] if INCLUDE_INDELS else []
        ),
        unit_bed=TESTING_UNIT_BED,
        unit_bed_index=WORKFLOW_TESTING_UNIT_INDEX
    output:
        vcf=ensure(str(ANNOTATED_VCF), non_empty=True),
        csi=ensure(str(ANNOTATED_VCF_CSI), non_empty=True),
        annotations=ensure(str(RAW_UNIT_ANNOTATIONS), non_empty=True)
    params:
        include_indels=INCLUDE_INDELS,
        testing_unit=TESTING_UNIT
    threads: 8
    resources:
        mem_mb=lambda wc, input, attempt: execution_mem_mb(
            wc, input, attempt, "annotate_variants"
        ),
        runtime=lambda wc, input, attempt: execution_runtime_minutes(
            wc, input, attempt, "annotate_variants"
        )
    conda:
        str(CONDA_ENV_DIR / "bioinformatics.yaml")
    envmodules: *mod("bioinformatics")
    container: sif("bioinformatics")
    benchmark:
        str(BENCHMARK_DIR / "annotate_variants.tsv")
    log:
        str(
            APPLY_LOG_DIR /
            "annotation" / "annotate_variants.log"
        )
    script:
        str(APPLY_SCRIPT_DIR / "annotate_variants.py")


rule summarize_annotations:
    """Report complete annotation misses before any per-sample scoring."""
    input:
        annotations=rules.annotate_variants.output.annotations
    output:
        report=ensure(str(ANNOTATION_QC_REPORT), non_empty=True)
    params:
        include_indels=INCLUDE_INDELS
    threads: 1
    resources:
        mem_mb=1024,
        runtime=30
    conda:
        str(CONDA_ENV_DIR / "validation.yaml")
    envmodules: *mod("validation")
    container: sif("validation")
    log:
        str(
            APPLY_LOG_DIR /
            "annotation" / "summarize_annotations.log"
        )
    script:
        str(APPLY_SCRIPT_DIR / "summarize_annotations.py")


rule require_valid_annotations:
    """Stop the workflow when an enabled annotation source matched nothing."""
    input:
        report=rules.summarize_annotations.output.report
    output:
        ok=ensure(str(ANNOTATION_QC_OK), non_empty=True)
    params:
        failure_message="Annotation validation failed",
        include_errors=False,
        success_marker="annotation validation passed\n"
    threads: 1
    resources:
        mem_mb=512,
        runtime=5
    conda:
        str(CONDA_ENV_DIR / "validation.yaml")
    envmodules: *mod("validation")
    container: sif("validation")
    log:
        str(
            APPLY_LOG_DIR /
            "annotation" / "require_valid_annotations.log"
        )
    script:
        str(APPLY_SCRIPT_DIR / "require_valid_report.py")


rule prepare_testing_units:
    """
    Normalize the configured testing-unit annotation table.
    """
    input:
        annotations=rules.annotate_variants.output.annotations,
        annotation_qc=rules.require_valid_annotations.output.ok
    output:
        unit_table=ensure(str(TESTING_UNIT_TABLE), non_empty=True)
    threads: 4
    resources:
        mem_mb=16000,
        runtime=180
    conda:
        str(CONDA_ENV_DIR / "r.yaml")
    envmodules: *mod("r")
    container: sif("r")
    log:
        str(
            APPLY_LOG_DIR /
            "annotation" / "prepare_testing_units.log"
        )
    script:
        str(APPLY_SCRIPT_DIR / "prepare_testing_units.R")
