import shlex


CADD_SNV = WORKFLOW_CADD_SNV
CADD_INDEL = WORKFLOW_CADD_INDEL
INCLUDE_INDELS = config["model"]["include_indels"]
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
        info_header=(
            "##INFO=<ID=TESTING_UNIT,Number=.,Type=String,"
            f'Description="PSAP {TESTING_UNIT} unit IDs">'
        ),
        query_format=(
            "%CHROM\\t%POS\\t%REF\\t%ALT\\t%INFO/PSAP_CADD_PHRED"
            "\\t%INFO/TESTING_UNIT\\n"
        ),
        indel_annotation=(
            "bcftools annotate "
            f"--annotations {shlex.quote(str(CADD_INDEL))} "
            "--columns CHROM,POS,REF,ALT,-,+INFO/PSAP_CADD_PHRED "
            "--pair-logic exact -Ou -"
            if INCLUDE_INDELS
            else "bcftools view -Ou -"
        )
    threads: 8
    resources:
        mem_mb=lambda wc, input, attempt: scaled_mem_mb(wc, input, attempt, 32000, 512),
        runtime=lambda wc, input, attempt: scaled_runtime_minutes(wc, input, attempt, 720, 30)
    conda:
        str(CONDA_ENV_DIR / "bioinformatics.yaml")
    envmodules: mod("bioinformatics")
    container: sif("bioinformatics")
    benchmark:
        str(BENCHMARK_DIR / "annotate_variants.tsv")
    log:
        str(
            LOG_DIR / "apply" / config["run"]["apply_prefix"] /
            "annotation" / "annotate_variants.log"
        )
    shell:
        r"""
        {{
            bcftools annotate \
                --annotations {input.cadd_snv:q} \
                --columns CHROM,POS,REF,ALT,-,INFO/PSAP_CADD_PHRED \
                --header-line '##INFO=<ID=PSAP_CADD_PHRED,Number=1,Type=Float,Description="CADD PHRED score">' \
                --pair-logic exact \
                -Ou \
                {input.vcf:q} | \
            {params.indel_annotation} | \
            bcftools annotate \
                --annotations {input.unit_bed:q} \
                --columns CHROM,FROM,TO,INFO/TESTING_UNIT \
                --header-line {params.info_header:q} \
                --merge-logic TESTING_UNIT:unique \
                --threads {threads} \
                -Oz \
                -o {output.vcf:q} \
                -

            bcftools index \
                --threads {threads} \
                --force \
                {output.vcf:q}

            {{
                printf 'CHROM\tPOS\tREF\tALT\tCADD_PHRED\tunit\n'
                bcftools query \
                    --format {params.query_format:q} \
                    {output.vcf:q}
            }} | bgzip --threads {threads} --stdout > {output.annotations:q}
        }} > {log:q} 2>&1
        """


rule summarize_annotations:
    """Report complete annotation misses before any per-sample scoring."""
    input:
        annotations=rules.annotate_variants.output.annotations,
        script=str(APPLY_SCRIPT_DIR / "summarize_annotations.py")
    output:
        report=ensure(str(ANNOTATION_QC_REPORT), non_empty=True)
    params:
        include_indels=str(INCLUDE_INDELS)
    threads: 1
    resources:
        mem_mb=1024,
        runtime=30
    conda:
        str(CONDA_ENV_DIR / "validation.yaml")
    envmodules: mod("validation")
    container: sif("validation")
    log:
        str(
            LOG_DIR / "apply" / config["run"]["apply_prefix"] /
            "annotation" / "summarize_annotations.log"
        )
    shell:
        """
        python {input.script:q} \
            --annotations {input.annotations:q} \
            --include-indels {params.include_indels:q} \
            --output {output.report:q} \
            > {log:q} 2>&1
        """


rule require_valid_annotations:
    """Stop the workflow when an enabled annotation source matched nothing."""
    input:
        report=rules.summarize_annotations.output.report,
        script=str(APPLY_SCRIPT_DIR / "require_valid_annotations.py")
    output:
        ok=ensure(str(ANNOTATION_QC_OK), non_empty=True)
    threads: 1
    resources:
        mem_mb=512,
        runtime=5
    conda:
        str(CONDA_ENV_DIR / "validation.yaml")
    envmodules: mod("validation")
    container: sif("validation")
    log:
        str(
            LOG_DIR / "apply" / config["run"]["apply_prefix"] /
            "annotation" / "require_valid_annotations.log"
        )
    shell:
        """
        python {input.script:q} \
            --report {input.report:q} \
            --output {output.ok:q} \
            > {log:q} 2>&1
        """


rule prepare_testing_units:
    """
    Normalize the configured testing-unit annotation table.
    """
    input:
        annotations=rules.annotate_variants.output.annotations,
        annotation_qc=rules.require_valid_annotations.output.ok,
        script=str(APPLY_SCRIPT_DIR / "prepare_testing_units.R"),
        libraries=R_RUNTIME_FILES
    output:
        unit_table=ensure(str(TESTING_UNIT_TABLE), non_empty=True)
    threads: 4
    resources:
        mem_mb=16000,
        runtime=180
    conda:
        str(CONDA_ENV_DIR / "r.yaml")
    envmodules: mod("r")
    container: sif("r")
    log:
        str(
            LOG_DIR / "apply" / config["run"]["apply_prefix"] /
            "annotation" / "prepare_testing_units.log"
        )
    shell:
        """
        Rscript {input.script:q} \
            --input {input.annotations:q} \
            --output {output.unit_table:q} \
            --threads {threads} \
            > {log:q} 2>&1
        """
