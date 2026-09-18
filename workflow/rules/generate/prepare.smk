"""Cohort-independent preparation for null-distribution generation."""

GENERATION_OUTDIR = Path(config["outdir"]) / "generate" / config["run"]["generation_prefix"]
GENERATION_PREPARED_DIR = GENERATION_OUTDIR / "prepared"
GENERATION_LOG_DIR = Path(config["outdir"]) / "logs" / "generate" / config["run"]["generation_prefix"]
GENERATION_ARTIFACTS = {
    "report": str(GENERATION_OUTDIR / "validation.json"),
    "observed_af": str(GENERATION_PREPARED_DIR / "observed_af.tsv"),
    "excluded_alleles": str(GENERATION_PREPARED_DIR / "excluded_alleles.tsv"),
    "units": str(GENERATION_PREPARED_DIR / "provisional_units.tsv"),
    "score_universe": str(GENERATION_PREPARED_DIR / "score_universe.tsv"),
    "score_gaps": str(GENERATION_PREPARED_DIR / "score_gaps.tsv"),
}


def generation_index(name):
    index = config["generation"][name].get("index")
    return [index] if index else []


def generation_optional_data(name):
    binding = config["generation"].get(name)
    return [binding["data"]] if binding else []


rule prepare_generation_inputs:
    """Validate and normalize generation resources without cohort dependencies."""
    input:
        population_vcf=config["generation"]["population_vcf"]["data"],
        population_vcf_index=lambda wc: generation_index("population_vcf"),
        score_snv=config["generation"]["score_snv"]["data"],
        score_snv_index=lambda wc: generation_index("score_snv"),
        coverage=config["generation"]["coverage"]["data"],
        coverage_index=lambda wc: generation_index("coverage"),
        testing_units=config["generation"]["testing_unit_bed"]["data"],
        testing_units_index=lambda wc: generation_index("testing_unit_bed"),
        mutation_rate=lambda wc: generation_optional_data("mutation_rate"),
        fasta=config["generation"]["fasta"]["data"],
        fasta_fai=config["generation"]["fasta"]["fai"],
        fasta_gzi=lambda wc: [config["generation"]["fasta"]["gzi"]] if config["generation"]["fasta"]["gzi"] else []
    output:
        report=ensure(GENERATION_ARTIFACTS["report"], non_empty=True),
        observed_af=ensure(GENERATION_ARTIFACTS["observed_af"], non_empty=True),
        excluded_alleles=ensure(GENERATION_ARTIFACTS["excluded_alleles"], non_empty=True),
        units=ensure(GENERATION_ARTIFACTS["units"], non_empty=True),
        score_universe=ensure(GENERATION_ARTIFACTS["score_universe"], non_empty=True),
        score_gaps=ensure(GENERATION_ARTIFACTS["score_gaps"], non_empty=True)
    params:
        af_tag=config["generation"]["af_source"]["overall"],
        coverage_minimum=config["generation_settings"]["coverage_minimum"]
    threads: 1
    resources:
        mem_mb=2048,
        runtime=60
    conda:
        str(Path(workflow.basedir) / "envs" / "validation.yaml")
    envmodules: *mod("validation")
    container: sif("validation")
    log:
        str(GENERATION_LOG_DIR / "prepare.log")
    script:
        str(Path(workflow.basedir) / "scripts" / "generate" / "prepare_generation.py")
