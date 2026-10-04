from shell_arguments import optional_regions_argument


VARIANT_DIR = APPLY_OUTPUT_DIR / "variants"
FILTERED_BCF = VARIANT_DIR / f"{APPLY_OUTPUT_DIR.name}.in_scope.bcf"
FILTERED_VCF = VARIANT_DIR / f"{APPLY_OUTPUT_DIR.name}.filtered_pass.vcf.gz"


rule filter_calibration_regions:
    """
    Retain variants in the generated bundle's calibration regions, if present.
    """
    input:
        vcf=apply_resource("cohort").staged_data,
        vcf_index=apply_resource("cohort").staged_index,
        validated=validation_marker,
        calibration_regions=optional_input(
            optional_resource_field("calibration_regions", "staged_data")
        ),
        calibration_regions_index=(
            optional_input(
                optional_resource_field("calibration_regions", "staged_index")
            )
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
        region_args=(
            lambda wc, input: optional_regions_argument(
                input.calibration_regions[0] if input.calibration_regions else None
            )
        ),
        calibration_message=(
            "Restricting variants to generated calibration regions"
            if optional_apply_resource("calibration_regions")
            else "Calibration regions resource not configured; retaining all variants"
        )
    shell:
        """
        echo "{params.calibration_message}" > {log:q}
        bcftools view \
            --threads {threads} \
            {params.region_args} \
            -Ob \
            -o {output.bcf:q} \
            {input.vcf:q} \
            >> {log:q} 2>&1

        bcftools index --threads {threads} --force {output.bcf:q} >> {log:q} 2>&1
        """


rule filter_excluded_records:
    """
    Apply all interval and allele exclusions from one combined variant set once.
    """
    input:
        bcf=rules.filter_calibration_regions.output.bcf,
        bcf_csi=rules.filter_calibration_regions.output.csi,
        exclusions=(
            optional_input(APPLY_RESOURCE_REGISTRY.combined_exclusions)
            if APPLY_RESOURCE_REGISTRY.exclusion_sources
            else []
        ),
        exclusions_index=(
            optional_input(APPLY_RESOURCE_REGISTRY.combined_exclusions_index)
            if APPLY_RESOURCE_REGISTRY.exclusion_sources
            else []
        )
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
        exclusion_enabled=(
            "true" if APPLY_RESOURCE_REGISTRY.exclusion_sources else "false"
        )
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
