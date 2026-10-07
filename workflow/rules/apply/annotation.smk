import shlex


ANNOTATED_VCF = ANNOTATION_DIR / f"{TESTING_UNIT}.vcf.gz"
ANNOTATED_VCF_CSI = Path(f"{ANNOTATED_VCF}.csi")
RAW_UNIT_ANNOTATIONS = ANNOTATION_DIR / f"{TESTING_UNIT}.annotated.tsv.gz"
ANNOTATION_QC_REPORT = ANNOTATION_DIR / f"{TESTING_UNIT}.annotation_qc.json"
ANNOTATION_QC_OK = ANNOTATION_DIR / f"{TESTING_UNIT}.annotation_qc.ok"


rule annotate_variants:
    """
    Attach CADD and the configured testing-unit annotation to retained variants.
    """
    input:
        vcf=FILTERED_VCF,
        vcf_csi=FILTERED_VCF_INDEX,
        cadd_snv=apply_resource("cadd_snv").staged_data,
        cadd_snv_index=apply_resource("cadd_snv").staged_index,
        cadd_indel=(
            optional_input(
                optional_resource_field("cadd_indel", "staged_data")
            )
            if WORKFLOW_PLAN.model.include_indels
            else []
        ),
        cadd_indel_index=(
            optional_input(
                optional_resource_field("cadd_indel", "staged_index")
            )
            if WORKFLOW_PLAN.model.include_indels
            else []
        ),
        unit_bed=apply_resource("testing_units").staged_data,
        unit_bed_index=apply_resource("testing_units").staged_index
    output:
        vcf=ensure(str(ANNOTATED_VCF), non_empty=True),
        csi=ensure(str(ANNOTATED_VCF_CSI), non_empty=True),
        annotations=ensure(str(RAW_UNIT_ANNOTATIONS), non_empty=True)
    params:
        cadd_header=(
            "##INFO=<ID=PSAP_CADD_PHRED,Number=1,Type=Float,"
            "Description=\"CADD PHRED score\">"
        ),
        unit_header=(
            "##INFO=<ID=TESTING_UNIT,Number=.,Type=String,"
            f"Description=\"PSAP {TESTING_UNIT} unit IDs\">"
        ),
        indel_stage=lambda wc, input: (
            "bcftools annotate "
            f"--annotations {shlex.quote(str(input.cadd_indel[0]))} "
            "--columns CHROM,POS,REF,ALT,-,+INFO/PSAP_CADD_PHRED "
            "--pair-logic exact -Ou - |"
            if input.cadd_indel
            else ""
        )
    threads: 8
    resources:
        mem_mb=policy_resources("annotate_variants").mem_mb,
        runtime=policy_resources("annotate_variants").runtime
    conda:
        env("bioinformatics")
    envmodules: *mod("bioinformatics")
    container: sif("bioinformatics")
    benchmark:
        str(BENCHMARK_DIR / "annotate_variants.tsv")
    log:
        str(APPLY_LOG_DIR / "annotation" / "annotate_variants.log")
    shell:
        """
        {{
            bcftools annotate \
                --annotations {input.cadd_snv:q} \
                --columns CHROM,POS,REF,ALT,-,INFO/PSAP_CADD_PHRED \
                --header-line {params.cadd_header:q} \
                --pair-logic exact \
                -Ou \
                {input.vcf:q} |
            {params.indel_stage}
            bcftools annotate \
                --annotations {input.unit_bed:q} \
                --columns CHROM,FROM,TO,INFO/TESTING_UNIT \
                --header-line {params.unit_header:q} \
                --merge-logic TESTING_UNIT:unique \
                --threads {threads} \
                -Oz \
                -o {output.vcf:q} \
                -

            bcftools index --threads {threads} --force {output.vcf:q}

            {{
                printf 'CHROM\\tPOS\\tREF\\tALT\\tCADD_PHRED\\tunit\\n'
                bcftools query \
                    --format '%CHROM\\t%POS\\t%REF\\t%ALT\\t%INFO/PSAP_CADD_PHRED\\t%INFO/TESTING_UNIT\\n' \
                    {output.vcf:q}
            }} | bgzip --threads {threads} --stdout > {output.annotations:q}
        }} > {log:q} 2>&1
        """


rule summarize_annotations:
    """Report complete annotation misses before any per-sample scoring."""
    input:
        annotations=rules.annotate_variants.output.annotations
    output:
        report=ensure(str(ANNOTATION_QC_REPORT), non_empty=True)
    params:
        include_indels=WORKFLOW_PLAN.model.include_indels
    threads: 1
    resources:
        mem_mb=policy_resources("summarize_annotations").mem_mb,
        runtime=policy_resources("summarize_annotations").runtime
    conda:
        env("validation")
    envmodules: *mod("validation")
    container: sif("validation")
    log:
        str(APPLY_LOG_DIR / "annotation" / "summarize_annotations.log")
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
        mem_mb=policy_resources("require_valid_annotations").mem_mb,
        runtime=policy_resources("require_valid_annotations").runtime
    conda:
        env("validation")
    envmodules: *mod("validation")
    container: sif("validation")
    log:
        str(APPLY_LOG_DIR / "annotation" / "require_valid_annotations.log")
    script:
        str(APPLY_SCRIPT_DIR / "require_valid_report.py")
