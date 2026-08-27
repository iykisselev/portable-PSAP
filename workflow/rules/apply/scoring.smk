LOOKUP_TABLES = {
    name: configured_path(path)
    for name, path in REFERENCE_BINDINGS["lookup_tables"].items()
}
HEMIZYGOUS_MODEL_ENABLED = config["model"]["hemizygote"]


rule apply_psap:
    """
    Calculate per-sample PSAP results
    """
    input:
        annotated_gt=rules.join_sample_annotations.output.annotated_gt,
        lookup_het=LOOKUP_TABLES["het"],
        lookup_hom=LOOKUP_TABLES["hom"],
        lookup_hem_chrx=(
            [LOOKUP_TABLES["hem_chrX"]] if HEMIZYGOUS_MODEL_ENABLED else []
        )
    output:
        psap=ensure(
            str(SCORED_DIR / "{sample}.{testing_unit}.tsv.gz"),
            non_empty=True,
        )

    params:
        score_max=str(config["model"]["score_max"]),
        hemizygote=config["model"]["hemizygote"]
    threads: 4
    resources:
        mem_mb=16000,
        runtime=360
    conda:
        str(CONDA_ENV_DIR / "r.yaml")
    envmodules: mod("r")
    container: sif("r")
    benchmark:
        str(BENCHMARK_DIR / "apply_psap_{sample}_{testing_unit}.tsv")
    log:
        str(
            LOG_DIR / "apply" / f'{config["run"]["apply_prefix"]}' / "samples" /
            "apply_psap_{sample}_{testing_unit}.log"
        )
    script:
        str(APPLY_SCRIPT_DIR / "apply_psap.R")
