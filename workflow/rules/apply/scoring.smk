import shlex


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
        script=str(APPLY_SCRIPT_DIR / "apply_psap.R"),
        libraries=R_RUNTIME_FILES,
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
        lookup_hem_chrx_arg=lambda wc, input: (
            f"--lookup-hem-chrx {shlex.quote(str(input.lookup_hem_chrx[0]))}"
            if input.lookup_hem_chrx else ""
        ),
        score_max=str(config["model"]["score_max"]),
        hemizygote=config["model"]["hemizygote"]
    threads: 1
    conda:
        str(CONDA_ENV_DIR / "r.yaml")
    envmodules: "r/4.3.3"
    log:
        str(
            LOG_DIR / "apply" / f'{config["run"]["apply_prefix"]}' / "samples" /
            "apply_psap_{sample}_{testing_unit}.log"
        )
    shell:
        """
        Rscript {input.script:q} \
            --input {input.annotated_gt} \
            --lookup-het {input.lookup_het:q} \
            --lookup-hom {input.lookup_hom:q} \
            {params.lookup_hem_chrx_arg} \
            --score-max {params.score_max} \
            --hemizygote {params.hemizygote} \
            --output {output.psap} \
            --threads {threads} \
            > {log} 2>&1
        """
