HEMIZYGOUS_MODEL_ENABLED = WORKFLOW_PLAN.model.hemizygote


rule apply_psap:
    """
    Calculate per-sample PSAP results
    """
    input:
        annotated_gt=rules.prepare_sample_genotypes.output.annotated_gt,
        lookup_tables=list(LOOKUP_SOURCE.tables.values()),
    output:
        psap=ensure(
            str(SCORED_DIR / f"{{sample}}.{TESTING_UNIT}.tsv.gz"),
            non_empty=True,
        ),
    log:
        str(APPLY_LOG_DIR / "samples" / f"apply_psap_{{sample}}_{TESTING_UNIT}.log"),
    benchmark:
        str(BENCHMARK_DIR / f"apply_psap_{{sample}}_{TESTING_UNIT}.tsv")
    conda:
        env("r")
    container:
        sif("r")
    envmodules:
        *mod("r"),
    threads: 4
    resources:
        mem_mb=policy_resources("apply_psap").mem_mb,
        runtime=policy_resources("apply_psap").runtime,
    params:
        score_max=str(LOOKUP_SOURCE.score_max),
        hemizygote=WORKFLOW_PLAN.model.hemizygote,
        lookup_contract=LOOKUP_SOURCE.contract,
        lookup_bin_width=LOOKUP_SOURCE.bin_width or "",
        partition_ids=LOOKUP_SOURCE.partition_ids,
        lookup_keys=tuple(
            f"{partition_id}\t{model_id}"
            for partition_id, model_id in LOOKUP_SOURCE.tables
        ),
    script:
        str(APPLY_SCRIPT_DIR / "apply_psap.R")
