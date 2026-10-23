# Portable PSAP

Portable PSAP is a Snakemake workflow [1] for prioritizing rare variants in a
cohort. It uses the PSAP (Population Sampling Probability) framework [3] and
its genomic-regions extension [7]. For each affected individual, it estimates
how unusually it is to observe a variant this deleterious (by CADD-PHRED
[2,5]) in a gene under the selected inheritance model.

The workflow has two stages:

1. Generate a null-distribution *bundle*: per-gene lookup tables that map a
   CADD score to a PSAP value. These are computed from a population VCF
   (e.g. gnomAD [4,6]), CADD scores, gene/region BEDs and exclusion regions.
   Build one bundle per resource set and reuse it.
2. Apply a bundle to a cohort VCF and get ranked candidates for affected
   samples, matched against carriers among unaffected controls.


## Main concepts

| Term | Meaning |
| --- | --- |
| Mode | Every run sets `workflow.mode` to `generate` (build a bundle), `apply` (score a cohort using an existing bundle) or `full` (both in one DAG). |
| Testing unit | The region a variant is tested within. This is usually a gene, identified by the 4th column of a BED file (e.g. an Ensembl gene ID). |
| Partition | A named BED of testing units, such as `cds`, `intronic`, `cadd-genomic-region`. A run can have several partitions, listed in order. If they overlap, the partition listed first owns the shared bases. Each partition gets its own PSAP. The per-partition PSAPs are combined per gene with ACAT [8]. |
| Disease model | The inheritance models: `DOM-het`, `REC-hom`, `X-linked-het`, `X-linked-hom` and `X-linked-hem`. The hemizygous model only runs when `model.hemizygote: true`. They are defined in `workflow/model_catalog.json`. |
| Bundle | The directory that `generate` writes, in the `partitioned_null_bundle_v1` format. It contains one lookup table per partition×model, the BEDs that define its scope, `metadata.json` with checksums, and a `validation.ok` certificate. `apply` refuses bundles that are uncertified, modified or incompatible. |
| Run manifest | A YAML file is created for each run (under `runs/`). It holds the mode, the input/output paths and the run name. It is passed with `--configfile`. |
| Scientific config | `config/config.yaml`: the assembly, resource release IDs, score grid and generation policy. It is loaded automatically and versioned with the code. Edit it only to deliberately change the model. |

## Repository layout

```text
config/config.yaml      scientific defaults (loaded automatically)
config/samples.tsv      sample-sheet template
runs/*.example.yaml     run-manifest templates and runnable demos
profiles/default/       local execution profile (auto-loaded by Snakemake)
profiles/slurm/         Slurm executor profile
profiles/grex/          site overlay for GREX (Univ. of Manitoba): account, partition, modules
workflow/Snakefile      entry point; rules in workflow/rules/{generate,apply}/
workflow/envs/          conda environments (controller + per-rule)
workflow/envmodules.yaml  Lmod module names used with --software-deployment-method env-modules
workflow/schemas/       JSON schemas the config is validated against
scripts/                standalone post-analysis utilities
tests/                  test suite and small fixture data
docs/                   detailed guides
```

## Quick start (local demos, ~minutes)

Run every command from the repository root. Relative paths in manifests are
resolved against the directory where you start `snakemake`, not against the
YAML file.

```bash
# 1. Controller environment: Snakemake 9.23.1 + Slurm executor plugin
mamba env create --file workflow/envs/snakemake.yaml   # or conda / micromamba
conda activate portable-psap-snakemake

# 2. Build a demo bundle from tiny fixture data
snakemake --configfile runs/generate.example.yaml --cores 2 --use-conda
#    -> results/generate/generate-demo/bundle/

# 3. Apply that bundle to a demo cohort
snakemake --configfile runs/apply.example.yaml --cores 2 --use-conda
#    -> results/apply/apply-demo/cohort/

# Or do both in one DAG:
snakemake --configfile runs/full.example.yaml --cores 2 --use-conda
```

Notes:

- On the first run, Snakemake builds the per-rule conda environments under
  `.snakemake/conda/`. These contain Python, R, bcftools and htslib, so that run
  needs network access or a channel mirror. You don't need to install those
  tools yourself.
- Snakemake automatically loads `profiles/default/` as the workflow profile.
  Its defaults are 8 cores, 32 GB, conda, `--rerun-incomplete` and printed shell
  commands. Values you give on the command line take precedence.
- The `apply` demo reads the bundle that the `generate` demo writes. Run
  `generate` first.
- Add `-n` (dry run) to any command to check the config and DAG without running
  anything. The DAG grows after the first checkpoint, so a dry run shows only
  the first few generation jobs.

## Configuration model

Two layers are merged and then validated against `workflow/schemas/`. A bad
key or a missing path fails at startup, before any job runs.

- `config/config.yaml` holds the scientific identity: `GRCh38`, CADD v1.7
  release IDs, `score_max: 70`, `score_bin_width: 0.05`, the minimum retained
  fraction of each unit, X-chromosome policies, and so on. It records *which*
  resources these are, but not *where* they are.
- A **run manifest** binds those identities to physical files and sets
  `outdir`, `workflow.mode` and the output prefix. Use absolute paths for
  production and cluster runs.

`outdir` plus `run.generation_prefix` and/or `run.apply_prefix` determine
every output path. Give each run a unique prefix.

## Applying an existing bundle to your cohort

For a production apply run, start from
`runs/apply.bundle.example.yaml` (minimal) or `runs/run.example.yaml` (with
exclusion resources):

```bash
cp runs/run.example.yaml runs/my_cohort.yaml
$EDITOR runs/my_cohort.yaml
snakemake --configfile runs/my_cohort.yaml -n      # dry run
snakemake --configfile runs/my_cohort.yaml --cores 16
```

What the manifest needs:

| Key | Content |
| --- | --- |
| `inputs.vcf` | `{data, index}`. A cohort VCF/VCF.GZ, **normalized and split** against the reference (`bcftools norm -m-any -f ref.fa`). Only `CHROM POS REF ALT` and `FORMAT/GT` are used. `index` may be `null`; the workflow then sorts, compresses and indexes a working copy. |
| `inputs.samples` **or** `inputs.ped` | Set exactly one and set the other to `null`. The sample sheet is a headered TSV with `family_id pid sex phenotype`. A PED file uses the standard first 6 columns. Sex is `1`=male, `2`=female. Phenotype is `1`=unaffected, `2`=affected; other values are excluded. Every processed `pid` must appear exactly once in both the VCF and the metadata. |
| `bindings.reference.fasta` | `{data, fai, gzi}`. The `.fai` index is required. `.gzi` is required for a bgzipped FASTA and can be `null` for a plain one. |
| `bindings.reference.cadd.snv` / `.indel` | CADD tables whose first six columns are `Chrom Pos Ref Alt RawScore PHRED`. The indel table is needed when `model.include_indels: true`. |
| `bindings.reference.exclusions` | Optional exclusion sets, keyed by the IDs in `config/config.yaml` (`encode-hg38-blacklist-v2`, `par-regions`). Each one is a BED, a TSV (intervals or exact alleles) or a VCF. A variant is dropped if any part of its REF span overlaps an excluded interval. |
| `bindings.reference.lookup_source` | `{type: bundle, directory: <path to bundle/>}` |
| `model.include_indels`, `model.hemizygote` | These must match the bundle. The bundle also fixes the assembly, the partitions and the score grid, and apply cannot override them. |

Apply mode can create source indexes (`.tbi`/`.csi`) for the VCF and tabular
resources. FASTA indexes are required. Treat inputs as read-only and keep them
outside `outdir`.

### What happens during apply

1. **Validation.** The workflow checks the VCF, metadata, FASTA, resources,
   bundle certificate and checksums, and model compatibility. The results go to
   `validation.json`. Errors stop the run before any heavy work. Warnings are
   recorded and the run continues.
2. **Staging and filtering.** Inputs are copied into `prepared/` and indexed
   there. Variants outside the bundle's calibration regions and variants that
   match an exclusion are dropped.
3. **Annotation.** Each variant gets a CADD-PHRED score and its gene/partition
   memberships. The run stops if CADD annotates no SNVs, if indel CADD
   annotates no indels when indels are enabled, or if no cohort variant falls
   in a partition. Each of these usually indicates a build or `chr`-prefix
   mismatch.
4. **Genotypes.** `GT` calls are classified as `het`, `hom` or `hem`, taking
   sex into account. Male X heterozygous calls are dropped when the hemizygote
   model is on. Haploid calls outside X/Y are rejected.
5. **Scoring.** For each sample, gene, partition and model, the workflow takes
   the highest-CADD compatible variant or variants (all ties are kept). It
   looks up the PSAP value for that score in the bundle.
6. **Reporting.** Each affected candidate is matched to controls carrying the
   identical allele in the same gene and partition. The partition PSAPs are
   then combined per individual, gene and model with ACAT.

## Generating a bundle

Use `runs/generate.example.yaml` as the template. The full annotated example is
in [docs/pipeline.md](docs/pipeline.md#generate-mode). Inputs:

| Key | Content |
| --- | --- |
| `generation.population_vcf` | Population sites VCF with AF INFO fields. Name the fields in `generation.af_source` (`overall`, plus `female` and `male` for chrX). |
| `generation.score_snv` / `score_indel` | CADD SNV and indel tables. The indel table defines the candidate indel universe. Set `generation_settings.variant_classes: [snv, indel]` together with `model.include_indels: true`. |
| `generation.testing_unit_partitions` | An ordered list of `{id, bed}`. Each BED has the columns `chrom start end unit_id` (0-based, half-open). Each `id` must also be declared under `reference.testing_unit_partitions` with a release identity. |
| `generation.exclusion`, `generation.par_regions` | BEDs of regions to remove (e.g. low-coverage or blacklisted regions) and the pseudoautosomal regions. |
| `generation.fasta` | The same structure as in apply. |
| `generation.unobserved_policy` | Controls how alleles missing from the population VCF are treated. `zero` gives them AF 0. `mutation_rate` uses a per-allele rate from `generation.mutation_rate` (a TSV with `CHROM POS REF ALT mu`) and optionally `generation.mutation_rate_indel`. |

Production generation **requires explicit `.tbi`/`.csi` indexes** for the
population VCF, the CADD tables and the mutation-rate tables. An index next to
the data file is not picked up automatically. The demos bypass this with
`generation.execution.allow_unindexed_fixture_inputs: true`. Never set that
option for real data.

Only autosomes and chrX are calibrated. A unit/partition is dropped if less
than `minimum_retained_fraction` (default 0.5) of its bases survive the
exclusions, and the reason is recorded in `unit_partition_status.tsv.gz`.

Generation is split into genomic windows × gene buckets, so it parallelizes
well on a cluster. You can tune the split under `generation.execution`
(`partition_size_bp`, `unit_bucket_size`).

### Full mode

`runs/full.example.yaml` generates a bundle and applies it in one DAG. The
FASTA, CADD tables and exclusions are declared once under
`bindings.reference.*`. Population and mutation-rate inputs go under
`generation.*`. The apply step automatically uses the new bundle's scope.

## Running on Slurm

The cluster route uses the `snakemake-executor-plugin-slurm` from the
controller environment. Run the controller on a login node, or inside a
long-running job if your site prefers that. It submits one Slurm job per rule
instance. Disable the auto-loaded local profile with `--workflow-profile none`.

**Conda on compute nodes** (generic site): create a small site overlay
containing your account and partition, then stack it on `profiles/slurm`.

```bash
mkdir -p profiles/mysite
cat > profiles/mysite/config.yaml <<'EOF'
default-resources:
  slurm_account: <acct>
  slurm_partition: <partition>
  mem_mb: 2048
  runtime: 60
EOF

snakemake --profile profiles/slurm --profile profiles/mysite --workflow-profile none \
  --configfile runs/my_cohort.yaml -n
```

`profiles/slurm/config.yaml` sets `jobs: 50`, conda deployment and
`rerun-incomplete`. Restate the default `mem_mb` and `runtime` in the overlay
so they remain in effect after the overlay is merged.
Conda environments are created on the first run under `.snakemake/conda/`, so
that directory must be on a shared filesystem.

**Environment modules + virtualenv** (GREX, or any Lmod site without conda): add
the site overlay.

```bash
snakemake --profile profiles/slurm --profile profiles/grex --workflow-profile none \
  --configfile runs/my_cohort.yaml -n
```

For this route:

- Edit `profiles/grex/config.yaml` to set `slurm_account`, `slurm_partition`
  and the node-local `tmpdir`.
- Check that the modules in `workflow/envmodules.yaml` exist (`module spider
  samtools`).
- Python rule dependencies (polars, pyarrow, pysam, jsonschema, …) come from a
  shared virtualenv. Compute jobs activate it via `precommand:` in the profile.

The exact `pip install` list and the checks for interpreter mix-ups are in
[docs/pipeline.md](docs/pipeline.md#running-on-grex-with-modules-and-a-virtualenv).
For your own cluster, copy `profiles/grex/` to a new directory and adapt it.

Remove `-n` to submit. Rules that only read or check small reports run locally
on the controller.

**Resuming:** if a job fails, fix the cause and rerun the same command.
Snakemake reruns incomplete jobs and reuses completed generation shards. Do not
create `validation.ok` by hand or edit generated manifests to get past a
failed check.

## Outputs

### Apply / full: `<outdir>/apply/<apply_prefix>/cohort/`

| File | One row per | Key columns |
| --- | --- | --- |
| `cohort_summary.tsv.gz` | affected individual × gene × model | `PSAP` (ACAT across partitions), `aggregation_status`, configured/eligible/contributing partition counts, per-partition PSAPs as JSON |
| `cohort_candidates.tsv.gz` | affected candidate variant | `candidate_id`, `pid`, `unit`, `partition_id`, `model`, `variant_key`, `PSAP`, `control_count`, `control_pids` |
| `control_variants.tsv.gz` | control × carried variant × compatible model | `pid`, `unit`, `partition_id`, `model`, `variant_key`, `PSAP` |
| `partitions/<id>.summary.tsv.gz` | per-partition evidence for affected individuals | |
| `endpoint_diagnostics.tsv.gz`, `aggregation_diagnostics.tsv.gz` | | Raw and mid-p diagnostics for calibration work. The production result uses raw survival probabilities. |
| `report_metadata.json` | | serialization details |

Start with `cohort_summary` to rank genes per individual. Then use
`cohort_candidates` to find the variant behind each result and to see whether
any controls carry it.

The `<apply_prefix>/` directory also holds `validation.json`,
`sample_manifest.tsv`, and intermediates: `prepared/`, `variants/`,
`annotations/` (including the annotation QC JSON), per-sample `raw/`,
`annotated/` and `scored/` tables, and `benchmarks/`. Logs are written to
`<outdir>/logs/apply/<apply_prefix>/`.

### Generate / full: `<outdir>/generate/<generation_prefix>/bundle/`

```text
lookup.<partition>.<model>.tsv.gz    one table per partition × model
resolved.<partition>.bed.gz          unit intervals after overlap resolution
eligible.<partition>.bed.gz          units that passed the retention threshold
calibration_regions.bed.gz           exact scope; apply uses it as an inclusion mask
excluded_alleles.tsv.gz              allele-level exclusions
candidate_audit.tsv.gz               indels rejected for structural reasons (e.g. REF crosses mask edge)
unit_partition_status.tsv.gz         which unit/partitions were kept or dropped, and why
metadata.json, validation.json, validation.ok
```

Copy or archive the **whole directory**. Apply verifies every component
against `metadata.json`. Generation intermediates under `scalable/` are kept
as restart points and can be deleted once the bundle is certified.

## Troubleshooting

| Symptom | Likely cause |
| --- | --- |
| Error before the DAG is built | A schema or plan error. Typical causes: both `ped` and `samples` are set, a prefix is missing, a production generation input has no index, or the bundle directory doesn't exist. The message names the offending key. |
| `No executor named slurm` | The Slurm plugin is missing from the Python that runs `snakemake`. |
| Annotation stops with "annotates no …" | Different builds or contig naming (`chr1` vs `1`) between the VCF, CADD and the bundle. |
| A gene is missing from the results | Look it up in the bundle's `unit_partition_status.tsv.gz`. It may have too little retained sequence or no compatible event mass. |
| Imports fail on compute nodes (modules route) | The virtualenv isn't activated by `precommand`, or a module swapped the Python interpreter. See [docs/pipeline.md](docs/pipeline.md#deployment-verification). |

## Tests

```bash
python tests/run_tier.py lint        # Ruff, snakefmt, and R lint; run before committing
python tests/run_tier.py unit        # fast, no Snakemake
python tests/run_tier.py workflow    # DAG construction and lint
python tests/run_tier.py e2e         # end-to-end runs of all modes on fixtures
python tests/run_tier.py all
```

Use the environment in `workflow/envs/ci.yaml` for these commands. The `all`
command runs formatting and lint checks before the full test suite.
GitHub Actions runs the functional test tiers; formatting and lint checks run
locally through the `lint` or `all` command.

## Further documentation

- [docs/project_concepts.md](docs/project_concepts.md): the statistical model
  (null probabilities, calibration scope, ACAT, interpretation limits).
- [docs/pipeline.md](docs/pipeline.md): the full operational reference for
  every input contract, stage, output and deployment route.
- [docs/generate_workflow_architecture.md](docs/generate_workflow_architecture.md)
  and [docs/apply_workflow_architecture.md](docs/apply_workflow_architecture.md):
  code structure for contributors.

## References

1. Köster, J. & Rahmann, S. Snakemake - a scalable bioinformatics workflow
   engine. *Bioinformatics* **28**, 2520–2522 (2012).
2. Kircher, M. *et al.* A general framework for estimating the relative
   pathogenicity of human genetic variants. *Nature Genetics* **46**,
   310–315 (2014).
3. Wilfert, A. B. *et al.* Genomewide significance testing of variation from
   single case exomes. *Nature Genetics* **48**, 1455–1461 (2016).
4. Karczewski, K. J. *et al.* The mutational constraint spectrum quantified
   from variation in 141,456 humans. *Nature* **581**, 434–443 (2020).
5. Rentzsch, P., Schubach, M., Shendure, J. & Kircher, M. CADD-Splice - improving
   genome-wide variant effect prediction using deep learning-derived splice
   scores. *Genome Medicine* **13**, 31 (2021).
6. Chen, S. *et al.* A genomic mutational constraint map using variation in
   76,156 human genomes. *Nature* **625**, 92–100 (2024).
7. Ogloblinsky, M.-S.C., *et al.* PSAP-Genomic-Regions: A method
   leveraging population data to prioritize coding and non-coding variants in
   whole-genome sequencing for rare disease diagnosis. *Genetic Epidemiology*
   **49**, e22593 (2025).
8. Liu, Y. *et al.* ACAT: A fast and powerful p value combination method for
   rare-variant analysis in sequencing studies. *The American Journal of Human
   Genetics* **104**, 410–421 (2019). https://doi.org/10.1016/j.ajhg.2019.01.002
