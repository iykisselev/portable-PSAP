VARIANT_DIR = APPLY_OUTPUT_DIR / "variants"
FILTERED_BCF = VARIANT_DIR / f"{APPLY_OUTPUT_DIR.name}.in_scope.bcf"
FILTERED_VCF = VARIANT_DIR / f"{APPLY_OUTPUT_DIR.name}.filtered_pass.vcf.gz"


rule filter_calibration_regions:
    """
    Retain variants in the generated bundle's calibration regions, if present.
    """
    input:
        vcf=WORKFLOW_INPUT_VCF,
        vcf_index=WORKFLOW_INPUT_VCF_INDEX,
        validated=validation_marker,
        calibration_regions=[WORKFLOW_CALIBRATION_REGIONS] if WORKFLOW_CALIBRATION_REGIONS else [],
        calibration_regions_index=(
            [WORKFLOW_CALIBRATION_REGIONS_INDEX] if WORKFLOW_CALIBRATION_REGIONS else []
        )
    output:
        bcf=ensure(str(FILTERED_BCF), non_empty=True),
        csi=ensure(str(FILTERED_BCF) + ".csi", non_empty=True)
    threads: 4
    resources:
        mem_mb=policy_resources("filter_calibration_regions").mem_mb,
        runtime=policy_resources("filter_calibration_regions").runtime
    conda:
        env("bioinformatics")
    envmodules: *mod("bioinformatics")
    container: sif("bioinformatics")
    log:
        str(APPLY_LOG_DIR / "variants" / "filter_calibration_regions.log")
    params:
        calibration_regions=(
            lambda wc, input: str(input.calibration_regions[0])
            if input.calibration_regions
            else ""
        ),
        calibration_regions_enabled="true" if WORKFLOW_CALIBRATION_REGIONS else "false"
    shell:
        """
        if {params.calibration_regions_enabled}; then
            echo "Restricting variants to generated calibration regions" > {log:q}
            bcftools view \
                --threads {threads} \
                --regions-file {params.calibration_regions:q} \
                -Ob \
                -o {output.bcf:q} \
                {input.vcf:q} \
                >> {log:q} 2>&1
        else
            echo "Calibration regions resource not configured; retaining all variants" > {log:q}
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
        bcf=rules.filter_calibration_regions.output.bcf,
        bcf_csi=rules.filter_calibration_regions.output.csi,
        exclusions=[COMBINED_EXCLUSIONS] if EXCLUSION_SOURCES else [],
        exclusions_index=[COMBINED_EXCLUSIONS_INDEX] if EXCLUSION_SOURCES else []
    output:
        vcf=ensure(str(FILTERED_VCF), non_empty=True),
        csi=ensure(str(FILTERED_VCF) + ".csi", non_empty=True)
    threads: 4
    resources:
        mem_mb=policy_resources("filter_excluded_records").mem_mb,
        runtime=policy_resources("filter_excluded_records").runtime
    conda:
        env("bioinformatics")
    envmodules: *mod("bioinformatics")
    container: sif("bioinformatics")
    log:
        str(APPLY_LOG_DIR / "variants" / "filter_excluded_records.log")
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
            echo "No exclusions configured; retaining all in-scope variants" > {log:q}
            bcftools view \
                --threads {threads} \
                -Oz \
                -o {output.vcf:q} \
                {input.bcf:q} \
                >> {log:q} 2>&1
        fi

        bcftools index --threads {threads} --force {output.vcf:q} >> {log:q} 2>&1
        """
