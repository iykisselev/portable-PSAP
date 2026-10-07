SPLIT_GENOTYPE_DIR = RAW_DIR / "split_genotypes"
PROCESSING_SAMPLE_IDS = [record.pid for record in SAMPLE_INDEX.records]
SPLIT_GENOTYPE_VCFS = [
    temp(ensure(str(SPLIT_GENOTYPE_DIR / f"{sample}.vcf.gz"), non_empty=True))
    for sample in PROCESSING_SAMPLE_IDS
]
SPLIT_GENOTYPE_INDICES = [
    temp(ensure(str(SPLIT_GENOTYPE_DIR / f"{sample}.vcf.gz.csi"), non_empty=True))
    for sample in PROCESSING_SAMPLE_IDS
]


rule split_genotypes:
    """Split all non-reference sample genotypes in one cohort pass."""
    input:
        vcf=rules.annotate_variants.output.vcf,
        vcf_index=rules.annotate_variants.output.csi,
        validated=validation_marker,
        annotation_qc=rules.require_valid_annotations.output.ok
    output:
        vcf=SPLIT_GENOTYPE_VCFS,
        csi=SPLIT_GENOTYPE_INDICES
    params:
        samples=PROCESSING_SAMPLE_IDS
    threads: 1
    resources:
        mem_mb=policy_resources("split_genotypes").mem_mb,
        runtime=policy_resources("split_genotypes").runtime
    conda:
        env("bioinformatics")
    envmodules: *mod("bioinformatics")
    container: sif("bioinformatics")
    log:
        str(APPLY_LOG_DIR / "samples" / "split_genotypes.log")
    script:
        str(APPLY_SCRIPT_DIR / "split_genotypes.py")


rule extract_gt:
    """
    Extract genotype information for every individual in the PED file
    """
    input:
        vcf=str(SPLIT_GENOTYPE_DIR / "{sample}.vcf.gz"),
        vcf_index=str(SPLIT_GENOTYPE_DIR / "{sample}.vcf.gz.csi")
    output:
        gt=temp(ensure(str(RAW_DIR / "{sample}.tsv.gz"), non_empty=True))
    wildcard_constraints:
        sample=SAMPLE_PATTERN
    threads: 2
    resources:
        mem_mb=policy_resources("extract_gt").mem_mb,
        runtime=policy_resources("extract_gt").runtime
    conda:
        env("bioinformatics")
    envmodules: *mod("bioinformatics")
    container: sif("bioinformatics")
    log:
        str(APPLY_LOG_DIR / "samples" / "extract_gt_{sample}.log")
    shell:
        """
        bcftools query \
            -f "%CHROM\t%POS\t%REF\t%ALT\t%INFO/PSAP_CADD_PHRED\t%INFO/TESTING_UNIT[\t%GT]\n" \
            {input.vcf:q} | \
            bgzip --threads {threads} --stdout > {output.gt:q}
        """

rule prepare_sample_genotypes:
    """Standardize, sex-classify, and annotate one sample in memory."""
    input:
        gt=rules.extract_gt.output.gt
    output:
        annotated_gt=ensure(
            str(ANNOTATED_DIR / "{sample}.{testing_unit}.tsv.gz"),
            non_empty=True,
        )
    wildcard_constraints:
        sample=SAMPLE_PATTERN,
        testing_unit=TESTING_UNIT_PATTERN
    params:
        sex=sample_sex,
        hemizygote_model=WORKFLOW_PLAN.model.hemizygote
    threads: 4
    resources:
        mem_mb=policy_resources("prepare_sample_genotypes").mem_mb,
        runtime=policy_resources("prepare_sample_genotypes").runtime
    conda:
        env("r")
    envmodules: *mod("r")
    container: sif("r")
    log:
        str(APPLY_LOG_DIR / "samples" / "prepare_sample_genotypes_{sample}_{testing_unit}.log")
    script:
        str(APPLY_SCRIPT_DIR / "prepare_sample_genotypes.R")
