rule extract_gt:
    """
    Extract genotype information for every individual in the PED file
    """
    input:
        vcf=rules.filter_excluded_records.output.vcf
    output:
        gt=ensure(str(RAW_DIR / "{sample}.tsv.gz"), non_empty=True)
    wildcard_constraints:
        sample=SAMPLE_PATTERN
    params:
        sample=lambda wc: wc.sample
    threads: 2
    resources:
        mem_mb=4096,
        runtime=120
    conda:
        str(CONDA_ENV_DIR / "bioinformatics.yaml")
    envmodules: mod("bioinformatics")
    container: sif("bioinformatics")
    log:
        str(
            LOG_DIR / "apply" / config["run"]["apply_prefix"] / "samples" /
            "extract_gt_{sample}.log"
        )
    shell:
        """
        bcftools query \
            -s {params.sample:q} \
            -i 'GT!="RR" && GT!="mis"' \
            -f "%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n" \
            {input.vcf:q} | \
            bgzip --threads {threads} --stdout > {output.gt:q}
        """

rule standardize_gt:
    """
    Standardize diploid calls and retain supported haploid X/Y calls.
    """
    input:
        gt=rules.extract_gt.output.gt
    output:
        gt=ensure(str(STANDARDIZED_DIR / "{sample}.tsv.gz"), non_empty=True)
    params:
        hemizygote_model=config["model"]["hemizygote"]
    wildcard_constraints:
        sample=SAMPLE_PATTERN
    threads: 2
    resources:
        mem_mb=8192,
        runtime=120
    conda:
        str(CONDA_ENV_DIR / "r.yaml")
    envmodules: mod("r")
    container: sif("r")
    log:
        str(
            LOG_DIR / "apply" / config["run"]["apply_prefix"] / "samples" /
            "standardize_gt_{sample}.log"
        )
    script:
        str(APPLY_SCRIPT_DIR / "standardize_gt.R")



rule apply_sex_model:
    """
    Apply sex-specific classification, including male haploid X/Y calls.
    """
    input:
        gt=rules.standardize_gt.output.gt,
        sample_manifest=validated_manifest_path
    output:
        gt=ensure(str(SEX_CLASSIFIED_DIR / "{sample}.tsv.gz"), non_empty=True)
    wildcard_constraints:
        sample=SAMPLE_PATTERN
    params:
        sex=sample_sex,
        hemizygote_model=config["model"]["hemizygote"]
    threads: 2
    resources:
        mem_mb=8192,
        runtime=120
    conda:
        str(CONDA_ENV_DIR / "r.yaml")
    envmodules: mod("r")
    container: sif("r")
    log:
        str(
            LOG_DIR / "apply" / config["run"]["apply_prefix"] / "samples" /
            "apply_sex_model_{sample}.log"
        )
    script:
        str(APPLY_SCRIPT_DIR / "classify_gt.R")



rule join_sample_annotations:
    """
    Join genotypes and testing unit annotations
    """
    input:
        gt=rules.apply_sex_model.output.gt,
        unit_tbl=str(TESTING_UNIT_TABLE)
    output:
        annotated_gt=ensure(
            str(ANNOTATED_DIR / "{sample}.{testing_unit}.tsv.gz"),
            non_empty=True,
        )
    wildcard_constraints:
        sample=SAMPLE_PATTERN,
        testing_unit=TESTING_UNIT_PATTERN
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
            LOG_DIR / "apply" / config["run"]["apply_prefix"] / "samples" /
            "join_sample_annotations_{sample}_{testing_unit}.log"
        )
    script:
        str(APPLY_SCRIPT_DIR / "join_sample_annotations.R")
