VARIANT_DIR = APPLY_OUTDIR / "variants"
FILTERED_BCF = VARIANT_DIR / f"{config['run']['apply_prefix']}.covered.bcf"
FILTERED_VCF = VARIANT_DIR / f"{config['run']['apply_prefix']}.filtered_pass.vcf.gz"


rule filter_low_coverage:
    """
    Retain variants in configured coverage regions, or pass through when absent.
    """
    input:
        vcf=WORKFLOW_INPUT_VCF,
        vcf_index=WORKFLOW_INPUT_VCF_INDEX,
        validated=validation_marker,
        covered_regions=[WORKFLOW_COVERAGE] if WORKFLOW_COVERAGE else [],
        covered_regions_index=(
            [WORKFLOW_COVERAGE_INDEX] if WORKFLOW_COVERAGE else []
        )
    output:
        bcf=ensure(str(FILTERED_BCF), non_empty=True),
        csi=ensure(str(FILTERED_BCF) + ".csi", non_empty=True)
    threads: 4
    resources:
        mem_mb=4096,
        runtime=180
    conda:
        str(CONDA_ENV_DIR / "bioinformatics.yaml")
    envmodules: *mod("bioinformatics")
    container: sif("bioinformatics")
    log:
        str(LOG_DIR / "apply" / config["run"]["apply_prefix"] / "variants" / "filter_low_coverage.log")
    params:
        coverage=(
            lambda wc, input: str(input.covered_regions[0])
            if input.covered_regions
            else ""
        ),
        coverage_enabled="true" if WORKFLOW_COVERAGE else "false"
    shell:
        """
        if {params.coverage_enabled}; then
            echo "Applying configured coverage-region filter" > {log:q}
            bcftools view \
                --threads {threads} \
                --regions-file {params.coverage:q} \
                -Ob \
                -o {output.bcf:q} \
                {input.vcf:q} \
                >> {log:q} 2>&1
        else
            echo "Coverage resource not configured; retaining all variants" > {log:q}
            bcftools view \
                --threads {threads} \
                -Ob \
                -o {output.bcf:q} \
                {input.vcf:q} \
                >> {log:q} 2>&1
        fi

        bcftools index --threads {threads} --force {output.bcf:q} >> {log:q} 2>&1
        """


rule filter_excluded_records:
    """
    Apply all interval and allele exclusions from one combined variant set once.
    """
    input:
        bcf=rules.filter_low_coverage.output.bcf,
        bcf_csi=rules.filter_low_coverage.output.csi,
        exclusions=[COMBINED_EXCLUSIONS] if EXCLUSION_SOURCES else [],
        exclusions_index=[COMBINED_EXCLUSIONS_INDEX] if EXCLUSION_SOURCES else []
    output:
        vcf=ensure(str(FILTERED_VCF), non_empty=True),
        csi=ensure(str(FILTERED_VCF) + ".csi", non_empty=True)
    threads: 4
    resources:
        mem_mb=4096,
        runtime=180
    conda:
        str(CONDA_ENV_DIR / "bioinformatics.yaml")
    envmodules: *mod("bioinformatics")
    container: sif("bioinformatics")
    log:
        str(LOG_DIR / "apply" / config["run"]["apply_prefix"] / "variants" / "filter_excluded_records.log")
    params:
        exclusions=(
            lambda wc, input: str(input.exclusions[0])
            if input.exclusions
            else ""
        ),
        exclusion_enabled="true" if EXCLUSION_SOURCES else "false"
    shell:
        """
        if {params.exclusion_enabled}; then
            echo "Applying the combined interval/variant exclusion set" > {log:q}
            bcftools isec \
                --threads {threads} \
                -C \
                -w1 \
                -Oz \
                -o {output.vcf:q} \
                {input.bcf:q} \
                {params.exclusions:q} \
                >> {log:q} 2>&1
        else
            echo "No exclusions configured; retaining all covered variants" > {log:q}
            bcftools view \
                --threads {threads} \
                -Oz \
                -o {output.vcf:q} \
                {input.bcf:q} \
                >> {log:q} 2>&1
        fi

        bcftools index --threads {threads} --force {output.vcf:q} >> {log:q} 2>&1
        """
