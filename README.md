# Portable PSAP

Portable PSAP is a Snakemake workflow [1] for rare-variant prioritization with
the PSAP (Population Sampling Probability) method [3] and its genomic-regions
extension [7]. Given a joint-called cohort VCF and a pedigree, it scores each
affected sample's most damaging variant per gene (or other testing unit) under
dominant, recessive and X-linked models, and reports whether unaffected
controls carry the same allele.

A PSAP value answers: *under this inheritance model, how likely is an event at
least this deleterious (by CADD-PHRED [2,5]) in this gene, given population
allele frequencies?* Smaller is more surprising. It is a ranking statistic, not
a pathogenicity call.

The workflow can also build the null-distribution lookup tables that PSAP
needs, from a population VCF (for example gnomAD [4,6]) and CADD scores, and
package them as a checksummed, self-describing "bundle".

- [What you need to know first](#what-you-need-to-know-first)
- [Quick start with the bundled demo](#quick-start-with-the-bundled-demo)
- [Scoring your own cohort (`apply`)](#scoring-your-own-cohort-apply)
- [Building lookup tables (`generate` and `full`)](#building-lookup-tables-generate-and-full)
- [Running on a Slurm cluster](#running-on-a-slurm-cluster)
- [Outputs](#outputs)
- [Troubleshooting](#troubleshooting)
- [Further documentation](#further-documentation)

## What you need to know first

### Three modes

Every run manifest sets `workflow.mode`. The target is always the default
`all`; there are no other supported targets.

| Mode | You supply | You get |
| --- | --- | --- |
| `apply` | Cohort VCF + sample metadata, reference resources, and lookup tables (third-party tables or a bundle from an earlier `generate`) | Cohort reports under `<outdir>/apply/<apply_prefix>/cohort/` |
| `generate` | Population VCF, CADD tables, testing-unit BED, exclusion BEDs, FASTA | A certified bundle under `<outdir>/generate/<generation_prefix>/bundle/` |
| `full` | Everything above in one manifest | The bundle, then the cohort scored against it, in one DAG |

Most users start with `apply`. Use `generate` once per reference/model
combination and reuse the bundle across cohorts.

### Two configuration layers

| File | Contains | Edit when |
| --- | --- | --- |
| `config/config.yaml` | The scientific model: assembly, testing-unit type, score grid, which models are on, resource *identities* (for example `CADD-v1.7-SNV`) | Rarely; only to change the science or record a new resource release |
| `runs/<name>.yaml` (passed with `--configfile`) | Mode, output location, file *paths*, and run-specific overrides | Every run |

Snakemake loads `config/config.yaml` first and merges the run manifest over
it. The merged result is schema-validated before any job is scheduled, so most
mistakes fail immediately with a message naming the offending key.

All relative paths in a manifest are resolved against the directory where you
run `snakemake`, not the location of the YAML file. Run from the repository
root, and prefer absolute paths for production and cluster runs.

### Key terms

- **Testing unit**: the region a PSAP value is computed for. The default is a
  gene, identified by column 4 of the testing-unit BED. The same IDs must
  appear in the lookup tables.
- **Lookup table**: one row per testing unit, one column per CADD-PHRED bin,
  holding the null probability for that role (`het`, `hom`, or `hem_chrX`).
- **Bundle**: the lookup tables plus the eligible units, calibration mask,
  excluded alleles, metadata and checksums that were used to build them. A
  bundle is usable only after its `validation.ok` marker exists. When applying
  a bundle, the workflow restricts cohort variants to the bundle's calibrated
  regions, so scores are never read outside the scope they were built for.
- **Disease models**: `DOM-het`, `REC-hom`, `X-linked-het`, `X-linked-hom` and
  `X-linked-hem`. Hemizygous scoring is controlled by `model.hemizygote`; indel
  scoring by `model.include_indels`.

### Repository layout

```text
config/config.yaml      scientific defaults and resource identities
runs/*.example.yaml     run manifests: small demos and production templates
profiles/               Snakemake profiles: default (local), slurm, grex (site overlay)
workflow/Snakefile      entry point; selects rules by workflow.mode
workflow/rules/         apply/ and generate/ rule files
workflow/scripts/       thin per-rule Python/R adapters
workflow/lib/           Python and R implementation (internal, not an installable package)
workflow/envs/          pinned Conda environments; envmodules.yaml and containers.yaml for other routes
workflow/schemas/       JSON schemas for config, manifests and bundles
scripts/                standalone post-analysis utilities
tests/                  fixtures and the tiered test runner
docs/                   pipeline guide, concepts, architecture notes
```

## Quick start with the bundled demo

The demo manifests use tiny fixtures from `tests/` and finish in a minute or
two on a laptop.

### 1. Install the controller environment

You need a Conda-compatible manager (`conda`, `mamba` or `micromamba`). The
controller environment provides Snakemake 9.23.1 and the Slurm executor plugin:

```bash
git clone https://github.com/iykisselev/portable-PSAP.git
cd portable-PSAP
mamba env create --file workflow/envs/snakemake.yaml
conda activate portable-psap-snakemake
snakemake --version   # 9.23.1
```

You do not need to install bcftools, htslib, R or the Python libraries
yourself. On the first run Snakemake builds the per-rule environments from
`workflow/envs/` into `.snakemake/conda/`. That needs internet access (or a
channel mirror) and a few GB of disk space in the repository directory.

### 2. Dry-run, then run

Snakemake picks up `workflow/Snakefile` and `profiles/default` automatically
when invoked from the repository root.

```bash
snakemake --configfile runs/apply.example.yaml -n                # show planned jobs
snakemake --configfile runs/apply.example.yaml --cores 2 --use-conda
```

The other demos work the same way:

```bash
snakemake --configfile runs/generate.example.yaml --cores 2 --use-conda
snakemake --configfile runs/full.example.yaml --cores 2 --use-conda
```

Results appear under `results/`, for example
`results/apply/apply-demo/cohort/cohort_candidates.tsv.gz` and
`results/generate/generate-demo/bundle/`.

`generate` and `full` start with a checkpoint, so a dry-run shows only the
first few jobs. The rest of the DAG (one job per genomic partition and unit
bucket) is added once scope discovery finishes.

## Scoring your own cohort (`apply`)

### 1. Make a manifest

```bash
cp runs/run.example.yaml runs/my_cohort.yaml           # third-party lookup tables
# or
cp runs/apply.bundle.example.yaml runs/my_cohort.yaml  # a bundle you generated
```

Edit every path, and give the run a unique `run.apply_prefix`. Outputs go to
`<outdir>/apply/<apply_prefix>/`, so reusing a prefix reuses (and may
overwrite) earlier results.

### 2. Prepare the inputs

Every resource must use the assembly in `config/config.yaml` (GRCh38 by
default). `chr`-prefixed and bare contig names are reconciled where possible,
but mixed builds are never valid.

| Manifest key | Requirement |
| --- | --- |
| `inputs.vcf` | VCF or bgzipped VCF with `GT`. It **must already be normalized and multiallelic-split** against the same FASTA (for example `bcftools norm -m -any -f ref.fa`). The index may be `null`; the workflow will sort and index a staged copy. |
| `inputs.samples` *or* `inputs.ped` | Exactly one; set the other to `null`. Sample sheet: headered TSV with `family_id pid sex phenotype`. PED: standard first six columns. Sex `1`=male, `2`=female; phenotype `1`=unaffected, `2`=affected. Each sample ID must appear exactly once in both the VCF and metadata, and IDs may not contain `/` or `\`. Samples with other phenotype codes are skipped. |
| `bindings.reference.fasta` | `data` plus `fai`. A bgzipped FASTA also needs `gzi`; for plain FASTA set `gzi: null`. These indexes are not created for you. |
| `bindings.reference.cadd.snv` | CADD SNV table. The first six columns must be `Chrom Pos Ref Alt RawScore PHRED`, as in the official CADD downloads. |
| `bindings.reference.cadd.indel` | Same layout for indels. Required when `model.include_indels: true` (the default); otherwise `null`. |
| `bindings.reference.lookup_source` | Either third-party tables (`type: external`, see `runs/run.example.yaml`) or a bundle: `{type: bundle, directory: <outdir>/generate/<prefix>/bundle}`. |
| `bindings.reference.testing_unit_bed` | BED with a unit ID in column 4. Needed for external tables; not needed for a bundle, which carries its own. |
| `bindings.reference.exclusion`, `bindings.reference.exclusions` | Regions or alleles to drop (blacklists, PAR, low-mappability regions). BED, coordinate or allele TSV, or VCF. A variant is dropped if any part of its REF span overlaps an excluded interval. Needed for external tables; a bundle carries its own. |

External lookup tables are tab-delimited: a unit ID followed by one numeric
column per CADD bin, with an optional header row. You need `het` and `hom`, and
`hem_chrX` when `model.hemizygote: true`.

When you apply a bundle, its model settings (`include_indels`, `hemizygote`,
unit type) and assembly must match the run. The bundle's own score grid is used
instead of the one in `config/config.yaml`. The certificate and component
checksums are verified before any annotation work starts.

### 3. Run

```bash
snakemake --profile profiles/default --configfile runs/my_cohort.yaml -n
snakemake --profile profiles/default --configfile runs/my_cohort.yaml --cores 16
```

`profiles/default` runs locally with Conda, a 32 GB total memory budget,
`--rerun-incomplete` and failed-log printing. `--cores` on the command line
overrides the profile's 8.

### What happens

1. **Validate.** The workflow checks the cohort, resources and model
   (sample IDs, normalization, FASTA agreement, CADD columns, unit IDs
   matching lookup rows). Errors stop the run here and are written to
   `validation.json`; warnings are recorded and the run continues.
2. **Stage.** Inputs are sorted, bgzipped and indexed into `prepared/` where
   needed. Source files are never modified.
3. **Filter.** Variants are restricted to bundle calibration regions (bundle
   mode only), then excluded records are removed.
4. **Annotate** with CADD-PHRED and every overlapping testing unit. The run
   stops if CADD or the BED annotates nothing, since that almost always means a
   build or contig-name mismatch.
5. **Classify genotypes per sample.** Calls are classified as `het`, `hom`,
   `hem` (male non-PAR X/Y, when hemizygote scoring is on), or female X
   `het`/`hom`. Hom-ref and missing calls are dropped.
6. **Score.** For each sample, unit and model, the workflow keeps the
   highest-CADD variant (all of them if tied), clips the score to the grid, and
   looks up its PSAP value.
7. **Report.** Affected samples produce candidates and unaffected samples
   produce control evidence. Each candidate is compared with controls that
   carry the identical allele in the same unit.

## Building lookup tables (`generate` and `full`)

Start from `runs/generate.example.yaml` (or `runs/full.example.yaml`) and the
production-shaped example in [docs/pipeline.md](docs/pipeline.md#generate-mode).
The main inputs are:

| Key (`generate` mode) | Content |
| --- | --- |
| `generation.population_vcf` | Sites VCF with AF INFO fields, for example gnomAD genomes. Name the fields in `generation.af_source` (`overall`, plus `female`/`male` for X). |
| `generation.score_snv` | CADD SNV table covering the calibration space. |
| `generation.score_indel` | CADD indel table; defines the indel candidate set when `generation_settings.variant_classes: [snv, indel]`. |
| `generation.testing_unit_bed` | Unit intervals (BED, ID in column 4). |
| `generation.exclusion`, `generation.par_regions`, `generation.exclusions` | Intervals to remove: low coverage or mappability, blacklists, PAR. Units keeping less than `minimum_retained_fraction` (default 0.5) of their bases are dropped. |
| `generation.fasta` | Reference FASTA and indexes. |
| `generation.unobserved_policy` | `zero` (alleles missing from the population VCF get AF 0) or `mutation_rate` (substitute a per-allele μ from `generation.mutation_rate`, and optionally `generation.mutation_rate_indel`). |

In `full` mode, the FASTA, CADD, BED and exclusion resources are given once
under `bindings.reference` and shared by both halves of the DAG. The generated
bundle is bound to the apply step automatically.

Practical notes for real-scale generation:

- **Indexes are mandatory.** The population VCF, CADD tables and mutation-rate
  tables must be bgzipped and have explicit `.tbi`/`.csi` `index:` bindings.
  An index sitting next to the file is not picked up automatically. The demos
  get around this with `allow_unindexed_fixture_inputs: true`; never set it
  for real data.
- **Generation covers autosomes 1–22 and chrX.** PAR is removed, and chrY is
  not calibrated.
- **Work is split into genomic partitions × unit buckets.** Tune
  `generation.execution.partition_size_bp` (default 5 Mb) and
  `unit_bucket_size` (default 250) to trade memory per job against job count.
  Neither changes the result.
- SNV-only runs write a `generated_bundle_v2` bundle; combined SNV+indel runs
  write `generated_bundle_v3`. Readers accept v1–v3, but bundles are never
  upgraded in place. Regenerate under a new prefix if you change the model or
  resources.
- Intermediate Parquet shards under `<outdir>/generate/<prefix>/scalable/` are
  kept on purpose as restart points. Delete them only after the bundle is
  certified and archived.

The scientific assumptions (`unobserved_policy`, X-chromosome AF and
mutation-rate policies, `score_bin_width`) live in `generation_settings` and
`generation` and are recorded in the bundle's `metadata.json`.

## Running on a Slurm cluster

Three profiles are provided and are meant to be layered:

| Profile | Purpose |
| --- | --- |
| `profiles/default` | Local executor, Conda environments. Snakemake loads it automatically as the workflow profile. |
| `profiles/slurm` | Slurm executor, up to 50 concurrent jobs, Conda environments, 2 GB / 60 min default resources. |
| `profiles/grex` | Example site overlay (University of Manitoba GREX): account, partition, node-local `tmpdir`, and environment modules instead of Conda. Copy and adapt it for your site. |

Pass `--workflow-profile none` with the Slurm profiles so that
`profiles/default` is not merged in as well.

### Option A: Slurm with Conda environments

This works wherever compute nodes can see the repository's `.snakemake/conda`
directory. Create a small site overlay with your scheduler identity, for
example `profiles/mysite/config.yaml`:

```yaml
default-resources:
  slurm_account: <account>
  slurm_partition: <partition>
  mem_mb: 2048
  runtime: 60
```

Activate the controller environment on the login node and run:

```bash
snakemake \
  --profile profiles/slurm \
  --profile profiles/mysite \
  --workflow-profile none \
  --configfile runs/my_cohort.yaml \
  -n                        # remove -n to submit
```

Build the Conda environments once on a node with internet access, for example
with `--conda-create-envs-only`, before submitting.

### Option B: Slurm with environment modules and a virtualenv (GREX-style)

Use this when your site provides modules and discourages Conda. External tools
come from `module load` (mapped in `workflow/envmodules.yaml`). Python
dependencies come from a shared virtualenv:

```bash
module purge
module load python/3.11.8 samtools/1.20

python -m venv .venv-grex
. .venv-grex/bin/activate
python -m pip install --upgrade pip wheel
python -m pip install "snakemake==9.23.1" snakemake-executor-plugin-slurm \
  "pandas>=2.0" "polars==1.44.1" pyarrow pysam "jsonschema==4.22.0" "PyYAML==6.0.1"
```

Then, in your copy of `profiles/grex/config.yaml`:

```yaml
software-deployment-method:
  - env-modules
precommand: ". /absolute/path/to/portable-PSAP/.venv-grex/bin/activate"
default-resources:
  slurm_account: <account>
  slurm_partition: <partition>
  tmpdir: /scratch            # node-local scratch
```

Check `module spider` for each name in `workflow/envmodules.yaml` and edit the
mapping to match your site. Then run:

```bash
snakemake --profile profiles/slurm --profile profiles/grex --workflow-profile none \
  --configfile runs/my_cohort.yaml -n          # remove -n to submit
```

Module and virtualenv setups often fail because a different `python` is
picked up on the compute node. Export `PYTHONNOUSERSITE=1`, clear stray
`PYTHONPATH` entries, and in a short interactive allocation confirm that
`python -c 'import snakemake, polars, pysam'` works after the modules load.

### Containers

`workflow/containers.yaml` maps each tool family to an Apptainer image path,
but no images are shipped and this route is not tested in CI. Conda is the
tested route.

### Operating notes

- Run the controller in `tmux` or `screen`, or as a small long-running batch
  job; it must stay alive for the whole run.
- Small validation and QC rules run on the controller host as `localrules`.
  Everything else goes to Slurm.
- Per-rule memory and runtime come from the rules and the profile's
  `default-resources`. Raise `mem_mb`/`runtime` in your overlay if your site
  needs it. Measured usage is written under `benchmarks/`.
- **Recovery:** fix the cause and rerun the same command. `rerun-incomplete` is
  on, and finished shards are reused. Never create `validation.ok` by hand or
  edit generated manifests to get past a failed check.

## Outputs

### Cohort reports: `<outdir>/apply/<apply_prefix>/cohort/`

All three files are gzip-compressed TSVs. `variant_key` is
`CHROM:POS-POS_REF>ALT`.

**`cohort_candidates.tsv.gz`**: one row per affected sample × unit × model ×
top-scoring variant. This is the main result.

| Column | Meaning |
| --- | --- |
| `candidate_id` | Stable ID such as `case_001:DOM:001` |
| `pid`, `phenotype`, `sex` | Sample metadata |
| `unit`, `model`, `variant_key` | Testing unit, disease model, variant |
| `PSAP` | Lookup value; lower is more surprising |
| `validation` | `control_match`: at least one unaffected sample carries the allele. `ok`: controls exist and none carry it. `no_controls`: the cohort has no control evidence. |
| `control_count`, `control_pids` | Matching unaffected carriers |

A `DOM-het` candidate matches control carriers under any model. Other models
match only controls carrying the allele under the same model.

**`cohort_summary.tsv.gz`**: candidates grouped by unit × model × variant
across affected samples (`affected_count`, `affected_pids`, `minimum_PSAP`, and
aggregated control fields). Use it for spotting recurrent hits.

**`control_variants.tsv.gz`**: the control evidence used above, one row per
unaffected sample × unit × model × variant, with the lowest available `PSAP`.

### Supporting files: `<outdir>/apply/<apply_prefix>/`

| Path | Content |
| --- | --- |
| `validation.json` | All validation errors and warnings; check this first |
| `sample_manifest.tsv` | Samples actually processed, with sex and phenotype |
| `prepared/`, `variants/`, `annotations/` | Staged inputs, filtered variants, annotated VCF and annotation QC |
| `raw/`, `standardized/`, `sex_classified/`, `annotated/`, `scored/` | Per-sample intermediates |
| `benchmarks/` | Per-rule runtime and memory |

Rule logs are in `<outdir>/logs/apply/<apply_prefix>/` and
`<outdir>/logs/generate/<generation_prefix>/`.

### Bundle: `<outdir>/generate/<generation_prefix>/bundle/`

```text
lookup_het.tsv.gz, lookup_hom.tsv.gz, lookup_hem_chrX.tsv.gz   lookup tables
eligible_testing_units.bed.gz, calibration_regions.bed.gz      calibrated scope
excluded_alleles.tsv.gz, excluded_units.tsv.gz                 what was removed, and why
candidate_audit.tsv.gz                                         v3 only: rejected indel candidates
metadata.json                                                  contract, policies, inputs, checksums
validation.json, validation.ok                                 certificate; ok is written last
```

Copy or archive the whole directory, since the checksums in `metadata.json`
cover every component.

### Optional review annotation

`scripts/annotate_cohort_summary.py` writes a copy of `cohort_summary.tsv.gz`
with added `gene_name` (from a GTF, treating `unit` as an Ensembl gene ID),
`cadd_phred` and `population_af` columns. It runs outside Snakemake and does
not touch workflow outputs. Its inputs must be bgzipped and indexed.

```bash
python scripts/annotate_cohort_summary.py \
  results/apply/my_cohort/cohort/cohort_summary.tsv.gz \
  results/apply/my_cohort/cohort/cohort_summary.review.tsv.gz \
  --cadd-snv cadd_snv.tsv.gz --cadd-indel cadd_indel.tsv.gz \
  --population-vcf gnomad.sites.vcf.gz --af-field AF \
  --gtf gencode.gtf.gz
```

## Troubleshooting

| Symptom | Likely cause |
| --- | --- |
| Error before any job runs | Schema or plan check failed. Typical causes: both or neither of `ped`/`samples` set, missing run prefix, or a generation resource without an `index:`. |
| Stops at `require_valid_inputs` | Read `validation.json`. Typical causes: sample IDs differing between VCF and metadata, an unnormalized VCF, or a FASTA/VCF contig mismatch. |
| Stops at `require_valid_annotations` | CADD or the testing-unit BED annotated nothing. Check the build, contig naming, and CADD column order. |
| `No executor named slurm` | `snakemake-executor-plugin-slurm` is not installed in the active controller environment. |
| Import or module errors on compute nodes | Wrong interpreter after `module load`; see the Option B notes above. |
| A gene is missing from the lookup | Look it up in the bundle's `excluded_units.tsv.gz`. It had low retained fraction or zero event mass. |
| Mutation-rate fallback fails on chrX | `x_chromosome.mutation_rate_policy` defaults to `reject`. Set `sex_averaged_proxy` if that assumption is acceptable. |

More cases are covered in
[docs/pipeline.md](docs/pipeline.md#validation-and-troubleshooting).

## Further documentation

- [docs/project_concepts.md](docs/project_concepts.md): the statistics and
  their interpretation limits, in plain language.
- [docs/pipeline.md](docs/pipeline.md): full operational reference, covering
  every input contract, configuration key, stage and output.
- [docs/apply_workflow_architecture.md](docs/apply_workflow_architecture.md),
  [docs/generate_workflow_architecture.md](docs/generate_workflow_architecture.md):
  code layout and rules for changing each half of the workflow.
- [tests/null_distribution/README.md](tests/null_distribution/README.md):
  generation test fixtures.

### Tests and contributing

```bash
python tests/run_tier.py unit        # fast library tests
python tests/run_tier.py workflow    # DAG construction and example manifests
python tests/run_tier.py all         # unit, schema, workflow, e2e, resilience
```

CI runs these tiers in the environment from `workflow/envs/ci.yaml`, and also
runs `ruff`, `snakefmt` and the R linter. `workflow/lib/python` is an internal
flat-module library that rule scripts import through
`workflow/shared/_library.py`; it is not a public package.

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

## License

MIT; see [LICENSE](LICENSE).
