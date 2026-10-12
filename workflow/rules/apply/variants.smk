from shell_arguments import optional_regions_argument

VARIANT_DIR = APPLY_OUTPUT_DIR / "variants"
FILTERED_BCF = VARIANT_DIR / f"{APPLY_OUTPUT_DIR.name}.in_scope.bcf"
EXCLUSION_FILTERED_VCF = VARIANT_DIR / f"{APPLY_OUTPUT_DIR.name}.filtered_pass.vcf.gz"
CALIBRATION_FILTER_ENABLED = optional_apply_resource("calibration_regions") is not None
EXCLUSION_FILTER_ENABLED = bool(APPLY_RESOURCE_REGISTRY.exclusion_sources)
IN_SCOPE_VCF = apply_resource("cohort").staged_data
IN_SCOPE_VCF_INDEX = apply_resource("cohort").staged_index


if CALIBRATION_FILTER_ENABLED:

    rule filter_calibration_regions:
        """Retain variants in the generated bundle's calibration regions."""
        input:
            vcf=apply_resource("cohort").staged_data,
            vcf_index=apply_resource("cohort").staged_index,
            validated=rules.require_valid_inputs.output.ok,
            calibration_regions=APPLY_RESOURCE_REGISTRY.inputs(
                "calibration_regions", "staged_data"
            ),
            calibration_regions_index=APPLY_RESOURCE_REGISTRY.inputs(
                "calibration_regions", "staged_index"
            ),
        output:
            bcf=ensure(str(FILTERED_BCF), non_empty=True),
            csi=ensure(str(FILTERED_BCF) + ".csi", non_empty=True),
        log:
            str(APPLY_LOG_DIR / "variants" / "filter_calibration_regions.log"),
        conda:
            env("bioinformatics")
        container:
            sif("bioinformatics")
        envmodules:
            *mod("bioinformatics"),
        threads: 4
        resources:
            mem_mb=policy_resources("filter_calibration_regions").mem_mb,
            runtime=policy_resources("filter_calibration_regions").runtime,
        params:
            region_args=lambda wc, input: optional_regions_argument(
                input.calibration_regions[0] if input.calibration_regions else None
            ),
        shell:
            """
            echo "Restricting variants to generated calibration regions" > {log:q}
            bcftools view \
                --threads {threads} \
                {params.region_args} \
                -Ob \
                -o {output.bcf:q} \
                {input.vcf:q} \
                >> {log:q} 2>&1

            bcftools index --threads {threads} --force {output.bcf:q} >> {log:q} 2>&1
            """

    IN_SCOPE_VCF = rules.filter_calibration_regions.output.bcf
    IN_SCOPE_VCF_INDEX = rules.filter_calibration_regions.output.csi


FILTERED_VCF = IN_SCOPE_VCF
FILTERED_VCF_INDEX = IN_SCOPE_VCF_INDEX

if EXCLUSION_FILTER_ENABLED:

    rule combine_exclusions:
        """Project exclusions from sites-only, calibration-filtered variants."""
        input:
            variants=IN_SCOPE_VCF,
            variants_index=IN_SCOPE_VCF_INDEX,
            exclusions=list(APPLY_RESOURCE_REGISTRY.workflow_exclusions),
        output:
            vcf=ensure(APPLY_RESOURCE_REGISTRY.combined_exclusions, non_empty=True),
            index=ensure(
                APPLY_RESOURCE_REGISTRY.combined_exclusions_index, non_empty=True
            ),
        log:
            str(APPLY_LOG_DIR / "prepare" / "combine_exclusions.log"),
        conda:
            env("validation")
        container:
            sif("validation")
        envmodules:
            *mod("validation"),
        threads: 1
        resources:
            mem_mb=policy_resources("combine_exclusions").mem_mb,
            runtime=policy_resources("combine_exclusions").runtime,
        params:
            sources=APPLY_RESOURCE_REGISTRY.exclusion_specs(),
            indexed_tsv_sources=[
                resource.staged_data
                for resource in APPLY_RESOURCE_REGISTRY.resources.values()
                if resource.origin == "generated_bundle"
                and resource.scientific_role == "variants_exclude"
                and resource.kind == "tsv"
            ],
        script:
            str(APPLY_SCRIPT_DIR / "combine_exclusions.py")

    rule filter_excluded_records:
        """Apply interval and exact-allele exclusions once."""
        input:
            variants=IN_SCOPE_VCF,
            variants_index=IN_SCOPE_VCF_INDEX,
            exclusions=APPLY_RESOURCE_REGISTRY.combined_exclusions,
            exclusions_index=APPLY_RESOURCE_REGISTRY.combined_exclusions_index,
        output:
            vcf=ensure(str(EXCLUSION_FILTERED_VCF), non_empty=True),
            csi=ensure(str(EXCLUSION_FILTERED_VCF) + ".csi", non_empty=True),
        log:
            str(APPLY_LOG_DIR / "variants" / "filter_excluded_records.log"),
        conda:
            env("bioinformatics")
        container:
            sif("bioinformatics")
        envmodules:
            *mod("bioinformatics"),
        threads: 4
        resources:
            mem_mb=policy_resources("filter_excluded_records").mem_mb,
            runtime=policy_resources("filter_excluded_records").runtime,
        shell:
            """
            echo "Applying the combined interval/variant exclusion set" > {log:q}
            bcftools isec \
                --threads {threads} \
                -C \
                -w1 \
                -Oz \
                -o {output.vcf:q} \
                {input.variants:q} \
                {input.exclusions:q} \
                >> {log:q} 2>&1

            bcftools index --threads {threads} --force {output.vcf:q} >> {log:q} 2>&1
            """

    FILTERED_VCF = rules.filter_excluded_records.output.vcf
    FILTERED_VCF_INDEX = rules.filter_excluded_records.output.csi
