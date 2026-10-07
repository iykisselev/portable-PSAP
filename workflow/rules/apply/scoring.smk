HEMIZYGOUS_MODEL_ENABLED = WORKFLOW_PLAN.model.hemizygote


rule prepare_lookups:
    """Parse and validate each shared PSAP lookup once into an RDS artifact."""
    input:
        lookup_het=LOOKUP_SOURCE.tables["het"],
        lookup_hom=LOOKUP_SOURCE.tables["hom"],
        lookup_hem_chrx=(
            [LOOKUP_SOURCE.tables["hem_chrX"]] if HEMIZYGOUS_MODEL_ENABLED else []
        ),
        validated=validation_marker
    output:
        lookup_het=ensure(str(PREPARED_LOOKUP_DIR / "het.rds"), non_empty=True),
        lookup_hom=ensure(str(PREPARED_LOOKUP_DIR / "hom.rds"), non_empty=True),
        lookup_hem_chrx=(
            ensure(str(PREPARED_LOOKUP_DIR / "hem_chrX.rds"), non_empty=True)
            if HEMIZYGOUS_MODEL_ENABLED
            else []
        )
    params:
        score_max=str(WORKFLOW_PLAN.model.score_max),
        hemizygote=WORKFLOW_PLAN.model.hemizygote,
        lookup_contract=LOOKUP_SOURCE.contract,
        lookup_bin_width=LOOKUP_SOURCE.bin_width or ""
    threads: 1
    resources:
        mem_mb=policy_resources("prepare_lookups").mem_mb,
        runtime=policy_resources("prepare_lookups").runtime
    conda:
        env("r")
    envmodules: *mod("r")
    container: sif("r")
    log:
        str(APPLY_LOG_DIR / "prepare" / "lookups.log")
    script:
        str(APPLY_SCRIPT_DIR / "prepare_lookups.R")


rule apply_psap:
    """
    Calculate per-sample PSAP results
    """
    input:
        annotated_gt=rules.prepare_sample_genotypes.output.annotated_gt,
        lookup_het=rules.prepare_lookups.output.lookup_het,
        lookup_hom=rules.prepare_lookups.output.lookup_hom,
        lookup_hem_chrx=(
            [rules.prepare_lookups.output.lookup_hem_chrx]
            if HEMIZYGOUS_MODEL_ENABLED
            else []
        )
    output:
        psap=ensure(
            str(SCORED_DIR / "{sample}.{testing_unit}.tsv.gz"),
            non_empty=True,
        )

    params:
        score_max=str(WORKFLOW_PLAN.model.score_max),
        hemizygote=WORKFLOW_PLAN.model.hemizygote,
        lookup_contract=LOOKUP_SOURCE.contract,
        lookup_bin_width=LOOKUP_SOURCE.bin_width or ""
    threads: 4
    resources:
        mem_mb=policy_resources("apply_psap").mem_mb,
        runtime=policy_resources("apply_psap").runtime
    conda:
        env("r")
    envmodules: *mod("r")
    container: sif("r")
    benchmark:
        str(BENCHMARK_DIR / "apply_psap_{sample}_{testing_unit}.tsv")
    log:
        str(APPLY_LOG_DIR / "samples" / "apply_psap_{sample}_{testing_unit}.log")
    script:
        str(APPLY_SCRIPT_DIR / "apply_psap.R")
