# Portable PSAP pipeline guide

This document is the main technical and operational guide to Portable PSAP. It
explains what the pipeline computes, how the three workflow modes fit together,
what each stage expects and produces, and how to run or extend the workflow.

For a shorter scientific introduction, read
[project concepts](project_concepts.md). Developers working on one branch of
the DAG can go directly to the
[generate architecture](generate_workflow_architecture.md) or
[apply architecture](apply_workflow_architecture.md).

## What the pipeline does

Portable PSAP has two main jobs:

1. Build a null-distribution lookup bundle from population allele frequencies,
   variant scores, genomic exclusions, and testing-unit definitions.
2. Apply lookup values to a cohort, then organize the resulting evidence into
   affected-sample candidates and unaffected-control observations.

The workflow supports three modes:

| Mode | What runs | Default target | Main result |
| --- | --- | --- | --- |
| `generate` | Null-distribution generation only | `null_distributions` | A certified lookup bundle |
| `apply` | Cohort processing and scoring only | `all` | Cohort candidate, summary, and control reports |
| `full` | Generation followed by application in one DAG | `full_results` | Both the certified bundle and cohort reports |

At a high level, data move through the pipeline like this:

```text
Population AF + CADD + testing units + exclusions
                              |
                              v
              Generate conditional null distributions
                              |
                              v
                    Certified lookup bundle
                              |
                              v
Cohort VCF + sample metadata -> filter -> annotate -> classify -> score
                                                              |
                                                              v
                                         candidates + controls + summary
```

In `apply` mode, the lookup bundle may be replaced by external legacy lookup
tables. In `full` mode, application always consumes the bundle produced by the
generation branch of the same DAG.

## Scientific model

### Testing units and scores

A testing unit is the region over which an event is evaluated. The default is
a gene, although `model.unit` also supports `cadd_region` and
`coding_cadd_region`. The fourth BED column is the stable unit identifier. The
same identifier must appear in the lookup tables.

CADD-PHRED is the severity score. Scores are clipped to `[0, score_max]` and
placed on a fixed-width grid. During application, if several compatible
variants tie for the highest score in a sample/unit/model combination, the
workflow keeps all tied variants.

### Genotype roles

Generation converts allele frequency `f` into a role-specific event
probability `q`:

| Lookup role | Event probability |
| --- | --- |
| Heterozygous (`het`) | `q = 2f(1-f)` |
| Homozygous alternate (`hom`) | `q = f²` |
| Male non-PAR chromosome-X hemizygous (`hem_chrX`) | `q = f` |

For each unit and role, the workflow assumes independent allele events. For a
score threshold `t`, it calculates a conditional tail probability:

```text
P(event score >= t | at least one event in the unit)
```

Equivalently, the numerator is the probability of at least one event among
alleles at or above `t`, and the denominator is the probability of at least
one compatible event anywhere in the unit. This score-bin survival row becomes
the lookup table used by the apply workflow.

### Alleles absent from the population panel

`generation.unobserved_policy` controls alleles that are part of the
calibration space but absent from the population VCF:

- `zero`: use an effective frequency of zero;
- `mutation_rate`: substitute a per-allele mutation rate `mu` before applying
  the role probability formula.

The supported mutation conversion is `de_novo_proxy`. For chromosome X,
`generation_settings.x_chromosome.mutation_rate_policy` must be either:

- `reject`, which stops if mutation fallback is needed on chromosome X; or
- `sex_averaged_proxy`, which permits the supplied `mu` to serve as the X
  fallback.

These settings are scientific assumptions and are written into bundle
metadata.

## Configuration model

The workflow merges two configuration layers:

1. `config/config.yaml` contains versioned scientific defaults and resource
   identities.
2. A run manifest supplies the mode, physical input paths, output root, and
   run-specific overrides.

The merged configuration is validated against `workflow/schemas/config.schema.yaml`
before rules are loaded. `workflow/lib/python/workflow_plan.py` then resolves
it into an immutable `WorkflowPlan`. Rules consume this plan rather than
reading configuration keys independently.

Paths in manifests are interpreted relative to the directory where Snakemake
is invoked. For production and cluster runs, absolute paths are safer.

### Common settings

| Setting | Meaning |
| --- | --- |
| `outdir` | Root directory for all results and logs |
| `workflow.mode` | `generate`, `apply`, or `full` |
| `run.generation_prefix` | Names one generation output tree |
| `run.apply_prefix` | Names one apply output tree |
| `model.unit` | Testing-unit identifier class |
| `model.hemizygote` | Enables the chromosome-X HEM model |
| `model.include_indels` | Enables indel annotation in apply mode; generation currently requires `false` |
| `model.score_max` | Maximum CADD-PHRED value represented by the lookup grid |

### Resource identities versus paths

`config/config.yaml` records the declared identity of each scientific resource,
such as the assembly and CADD release. A run manifest binds those identities to
physical files. Changing a path does not change the declared scientific
identity; when replacing a release, update both intentionally.

## Input contracts

All genomic resources must use the configured assembly and compatible contig
names. The workflow canonicalizes common `chr` prefixes internally where the
relevant reader supports it, but mixed builds are never valid.

### Reference FASTA

The binding contains:

```yaml
fasta:
  data: /path/reference.fa.gz
  fai: /path/reference.fa.gz.fai
  gzi: /path/reference.fa.gz.gzi
```

An uncompressed FASTA requires `.fai` and may set `gzi: null`. A BGZF FASTA
requires both `.fai` and `.gzi`.

### Testing-unit BED

The testing-unit file needs at least four columns:

```text
CHROM  START  END  testing_unit_id
```

BED coordinates are zero-based and half-open. A unit may have several rows;
the workflow merges its intervals when calculating retained fractions.

### Exclusion resources and retained scope

Migration: rename `generation.coverage` to `generation.exclusion` in generate
mode, or `bindings.reference.coverage` to `bindings.reference.exclusion` in
apply/full mode. Rename scientific metadata `reference.coverage` to
`reference.exclusion`, and `generation_settings.coverage_minimum` to
`generation_settings.minimum_retained_fraction`. An existing exclusion BED
needs no content conversion. A former inclusion BED must not simply be renamed:
it needs the opposite intervals within the supported reference scope.
Regenerate the bundle under a new generation prefix after migrating; old bundles
retain their old scope. The new provisional-unit columns are `supported_bases`,
`retained_bases`, and `retained_fraction`.

The `exclusion` binding contains regions to remove, not regions to include.
It may describe low coverage, mapping problems, or another source of unreliable
sequence. Use `data: null, index: null` when no primary exclusion is needed.
Generation interprets interval inputs as zero-based, half-open BED coordinates.
It combines the primary binding with the other configured interval exclusions.

For each unit, the denominator is its union of supported, non-PAR bases before
other exclusions. The numerator is the part remaining after all other interval
exclusions. `generation_settings.minimum_retained_fraction: 0.5` requires at
least half those bases to remain. Units with no remaining bases are excluded.

Apply-mode exclusions may be BED, TSV, or VCF:

- BED: zero-based, half-open intervals;
- positional TSV: `CHROM, start, end` using the configured TSV layout;
- allele TSV: `CHROM, POS, REF, ALT`;
- VCF: exact allele matches.

Apply mode projects intervals onto the full span of each cohort variant,
including REF length and `INFO/END`, before constructing a combined exclusion
VCF.

### CADD score table

The SNV and optional indel tables use the first six CADD columns:

```text
Chrom  Pos  Ref  Alt  RawScore  PHRED
```

Generation currently calibrates SNVs only. Apply mode can annotate indels when
`model.include_indels: true` and an indel score resource is configured.

### Population VCF

The population VCF supplies the AF INFO fields named in `generation.af_source`:

```yaml
af_source:
  overall: AF
  female: AF_female
  male: AF_male
```

`overall` is always required. With `x_chromosome.af_policy: sex_specific`,
female AF is required; male AF is also required when the HEM model is enabled.
Selected frequencies must be finite and in `[0, 1]`.

### Mutation-rate resource

Mutation fallback uses a headerless or `#`-headered TSV:

```text
CHROM  POS  REF  ALT  mu
```

`mu` must be finite and in `[0, 1]`. Every absent calibration allele that is
encountered must have a matching mutation-rate row. Large resources should be
BGZF-compressed and tabix-indexed by columns 1 and 2.

### Cohort VCF

The cohort VCF must be normalized and split against the configured FASTA. The
workflow uses `CHROM`, `POS`, `REF`, `ALT`, and per-sample `GT`. Each selected
sample must appear exactly once in the VCF header and metadata.

### Sample metadata

Configure exactly one of `inputs.ped` and `inputs.samples`.

A sample sheet is a tab-delimited file with:

```text
family_id  pid  sex  phenotype
```

PED input uses the conventional first six columns. Sex is `1` for male and `2`
for female. Phenotype is `1` for unaffected and `2` for affected. Other values
are treated as unknown or excluded according to the validation logic.

## Generate mode

### Minimal zero-policy manifest

```yaml
outdir: /shared/results

workflow:
  mode: generate

run:
  generation_prefix: grch38-gene-v1

model:
  include_indels: false

generation:
  unobserved_policy: zero
  population_vcf: {data: /shared/resources/population.vcf.gz, index: /shared/resources/population.vcf.gz.tbi}
  score_snv: {data: /shared/resources/cadd.tsv.gz, index: /shared/resources/cadd.tsv.gz.tbi}
  exclusion: {data: /shared/resources/exclusion.bed.gz, index: /shared/resources/exclusion.bed.gz.tbi}
  testing_unit_bed: {data: /shared/resources/units.bed.gz, index: /shared/resources/units.bed.gz.tbi}
  par_regions: {data: /shared/resources/par.bed.gz, index: /shared/resources/par.bed.gz.tbi}
  fasta:
    data: /shared/resources/reference.fa
    fai: /shared/resources/reference.fa.fai
    gzi: null
  af_source:
    overall: AF
    female: AF_female
    male: AF_male
```

### Mutation-rate (`mu`) proxy

To use mutation rates for unobserved alleles, add or override:

```yaml
generation_settings:
  mutation_conversion: de_novo_proxy
  x_chromosome:
    mutation_rate_policy: sex_averaged_proxy

generation:
  unobserved_policy: mutation_rate
  mutation_rate:
    data: /shared/resources/mutation_rates.tsv.gz
    index: /shared/resources/mutation_rates.tsv.gz.tbi
```

There is no `mu_proxy` enum value. The executable configuration is the
combination of `unobserved_policy: mutation_rate` and
`mutation_conversion: de_novo_proxy`.

### Generation stages

#### 1. Scope discovery

`prepare_generation_scope` constructs the supported calibration space:

1. Read testing-unit intervals and the FASTA index.
2. Restrict the supported contigs to autosomes 1–22 and chromosome X.
3. Remove PAR to define the supported denominator for each unit.
4. Subtract the union of primary and additional interval exclusions.
5. Calculate retained fraction per unit and mark units below
   `minimum_retained_fraction` as ineligible.
6. Create deterministic genomic partitions and eligible-unit buckets.
7. Record source paths and SHA-256 identities.

Because this rule is a Snakemake checkpoint, the downstream partition and
bucket wildcards do not exist until its manifests have been written.

#### 2. Partition preparation

`prepare_generation_partition` works on one genomic window. It enumerates the
SNV alleles in callable testing-unit sequence, joins each allele to unit
membership and CADD score, and classifies its population state as present,
absent, or rejected. Duplicate population states, duplicate scores, invalid AF
values, and missing calibration scores stop the job.

The rule writes three Parquet shards:

- `allele_states`: population state and AF values;
- `allele_membership`: allele-to-unit relationships;
- `allele_scores`: CADD score per allele.

#### 3. Partition-by-bucket statistics

`calculate_partition_null_distribution` joins the three partition shards for
one bounded unit bucket. It applies the zero or mutation-rate fallback, derives
role probabilities, maps scores to fixed-width bins, and writes additive
no-event statistics. Storing counts and sums of `log(1-q)` allows the workflow
to combine partitions without creating a genome-wide event table.

#### 4. Bucket reduction

`reduce_null_distribution_bucket` combines all partition statistics for one
unit bucket. It uses compensated floating-point summation, calculates event
mass and the conditional survival row, and emits compact diagnostics.

Possible diagnostic states are:

- `ok`: finite positive event mass and a valid lookup row;
- `zero_event_mass`: no compatible event has positive mass;
- `numerical_failure`: invalid aggregate probability; finalization aborts.

#### 5. Bundle finalization and certification

`finalize_null_distribution_bundle` keeps a unit only if it meets retention and
has positive mass in every applicable role. For an eligible chromosome-X unit,
HEM must be valid when the model is enabled. It then streams the public lookup,
scope, exclusion, and metadata files.

`validate_null_distribution_bundle` independently checks the component
contract and writes `validation.ok` last. Downstream consumers should depend on
that marker rather than only checking for lookup files.

## Apply mode

### Lookup sources

External lookup tables use:

```yaml
lookup_source:
  type: external
  contract: external_legacy_v0
  tables:
    het: /shared/lookups/lookup_het.tsv.gz
    hom: /shared/lookups/lookup_hom.tsv.gz
    hem_chrX: /shared/lookups/lookup_hem_chrX.tsv.gz
```

A certified bundle uses `type: bundle` and supplies its tables, metadata,
validation report, marker, eligible units, calibration mask, excluded alleles,
and excluded units. Contracts `generated_bundle_v1` and `generated_bundle_v2`
are readable, although this version writes only `generated_bundle_v2`. A v1
bundle cannot be silently upgraded to v2 and should be regenerated when v2
semantics are required.

### Apply stages

#### 1. Validation

Three rules run independently:

- `validate_cohort` checks the VCF, metadata, FASTA, sample agreement,
  normalization, and sex-chromosome requirements, and writes the canonical
  sample manifest.
- `validate_resources` checks calibration regions, exclusions, testing units, CADD
  structure, contigs, and indexes.
- `validate_model` checks lookup roles, score-grid compatibility, testing-unit
  overlap, and bundle certification when applicable.

`merge_validation` combines the reports. `require_valid_inputs` writes
`validation.ok` only if no blocking error remains. Expensive stages depend on
this marker.

#### 2. Resource staging

`ApplyResourceRegistry` provides one source and one prepared target for every
resource. `prepare_vcf_resource` and `prepare_table_resource` create sorted,
BGZF-compressed, indexed working copies as needed. Validation and staging read
the same registry, preventing path drift between stages.

#### 3. Coverage and exclusion filtering

`filter_calibration_regions` keeps records overlapping a generated bundle's
calibration mask. External apply mode passes all variants through this stage;
the following exclusion filter removes the primary and additional exclusions.

`combine_exclusions` projects interval and exact-allele exclusions onto the
cohort. `filter_excluded_records` subtracts the combined set once.

#### 4. Annotation

`annotate_variants` attaches CADD-PHRED and every overlapping testing-unit ID.
SNVs use the SNV table; indels use the indel table only when enabled.

`summarize_annotations` and `require_valid_annotations` stop the pipeline when
an enabled source unexpectedly annotates no relevant records.
`prepare_testing_units` normalizes annotation output to one row for each
allele/unit pair.

#### 5. Per-sample genotype processing

For each sample:

1. `extract_gt` removes reference-homozygous and missing calls while extracting
   allele and GT fields.
2. `standardize_gt` normalizes phased and unphased diploid notation and retains
   supported haploid X/Y calls.
3. `apply_sex_model` assigns sex-aware classes. Male X/Y alternate calls can
   become `hem`; female X calls become `het_female` or `hom_female`; invalid
   ploidy/sex combinations are removed or rejected.
4. `join_sample_annotations` joins genotype rows to CADD and unit annotations.

#### 6. PSAP scoring

`apply_psap` selects genotypes compatible with each enabled disease model. For
each sample, testing unit, and model it keeps all variants tied at the highest
CADD-PHRED score, maps that score to the lookup grid, and records the selected
lookup value as `PSAP`.

#### 7. Cohort reporting

- `build_sample_candidates` creates deterministic candidate IDs for affected
  participants.
- `build_control_reference` retains carried-variant evidence from unaffected
  participants.
- `build_cohort_report` compares candidates with controls carrying the exact
  same `unit` and `variant_key`, then writes candidate and cohort summaries.

Dominant heterozygous candidates may use evidence from any configured control
model. Other candidate models require a matching control model.

## Full mode integration

Full mode is one connected DAG, not two independent commands. It uses one
configured reference binding set and changes the effective apply resources
only after generation has produced a certified bundle.

| Resource | Generation uses | Apply uses |
| --- | --- | --- |
| FASTA | `bindings.reference.fasta` | The same configured FASTA |
| CADD SNV | `bindings.reference.cadd.snv` | The same configured CADD table |
| Exclusion / retained scope | Configured exclusion intervals are subtracted | Generated `calibration_regions.bed.gz` is an inclusion mask |
| Testing units | Configured testing-unit BED | Generated `eligible_testing_units.bed.gz` |
| Variant exclusions | Configured `variants_exclude` resources | Generated `excluded_alleles.tsv.gz` |
| PAR | Exactly one configured `par_regions` resource | Already reflected in generated scope |
| Lookup tables | Generated from population inputs | Generated HET/HOM/HEM tables |

This handoff ensures that lookup values are applied only within the scope for
which they were generated.

## Output layout

### Generated bundle

```text
<outdir>/generate/<generation_prefix>/bundle/
├── lookup_het.tsv.gz
├── lookup_hom.tsv.gz
├── lookup_hem_chrX.tsv.gz        # when HEM is enabled
├── eligible_testing_units.bed.gz
├── calibration_regions.bed.gz
├── excluded_alleles.tsv.gz
├── excluded_units.tsv.gz
├── metadata.json
├── validation.json
└── validation.ok
```

`metadata.json` records the bundle contract, assembly, contig convention,
binning, roles, policies, resource identities and checksums, component
checksums, and row counts.

Generation intermediates live under `scalable/prepared`,
`scalable/partitions`, and `scalable/calculated`. They are intentionally kept
as recovery points for large runs.

### Apply results

```text
<outdir>/apply/<apply_prefix>/
├── validation.json
├── validation.ok
├── sample_manifest.tsv
├── prepared/
├── variants/
├── annotations/
├── raw/
├── standardized/
├── sex_classified/
├── annotated/
├── scored/
├── candidates/
├── benchmarks/
└── cohort/
    ├── cohort_candidates.tsv.gz
    ├── cohort_summary.tsv.gz
    └── control_variants.tsv.gz
```

The candidate report contains affected-sample evidence, PSAP values, exact
variant keys, control-match status, and matching control IDs. The summary
groups recurring unit/model/candidate signatures across affected samples. The
control report contains the evidence records used for comparison.

Logs are written below `<outdir>/logs/generate/<generation_prefix>/` and
`<outdir>/logs/apply/<apply_prefix>/`.

## Running locally

Use a mode-specific manifest from `runs/` as the starting point. A dry-run is
the fastest way to catch schema, path, and DAG errors:

```bash
snakemake \
  --profile profiles/default \
  --configfile runs/generate.example.yaml \
  --dry-run
```

Remove `--dry-run` to execute. Replace the manifest with
`runs/apply.example.yaml` or `runs/full.example.yaml` for the other modes.

## Running on GREX with modules and a virtualenv

The Slurm controller and Python rule dependencies can live in a shared
virtualenv; external programs are selected through environment modules.

Create the virtualenv once from the repository root:

```bash
module purge
module load python/3.11.8
module load samtools/1.20

python -m venv .venv-grex
. .venv-grex/bin/activate
python -m pip install --upgrade pip wheel
python -m pip install \
  "snakemake==9.23.1" \
  snakemake-executor-plugin-slurm \
  "pandas>=2.0" \
  "polars==1.44.1" \
  pyarrow pysam \
  "jsonschema==4.22.0" \
  "PyYAML==6.0.1"
```

Set the real account, partition, and scratch path in
`profiles/grex/config.yaml`. Add an absolute shared virtualenv path so compute
jobs activate it in their non-interactive shell:

```yaml
precommand: ". /absolute/path/to/portable-PSAP/.venv-grex/bin/activate"
```

Confirm that module names in `workflow/envmodules.yaml` exist on GREX and
adjust them if necessary.

Dry-run:

```bash
snakemake \
  --profile profiles/slurm \
  --profile profiles/grex \
  --workflow-profile none \
  --configfile runs/grex-generate.yaml \
  --software-deployment-method env-modules \
  --dry-run \
  null_distributions
```

Submit by removing `--dry-run`. Use `all` for apply mode or `full_results` for
full mode.

## Scaling, retries, and recovery

Generation requires explicit `.tbi` or `.csi` index bindings for the population
VCF, CADD SNV table, and any configured mutation-rate table. This is checked
while building the workflow plan, before jobs are submitted. An index next to
the data file is not automatically discovered. Full mode checks the shared
`bindings.reference.cadd.snv.index` binding for CADD.

The small demonstration manifests explicitly set
`generation.execution.allow_unindexed_fixture_inputs: true`. This exception
only permits uncompressed files up to 1 MiB; scope preparation checks that limit
before hashing inputs or scheduling partition jobs. Leave the option unset
(the default is `false`) for production runs. Coverage, testing-unit and exclusion
BEDs are read once during scope preparation and do not need indexes for generation.

Generation exposes three performance-tuning settings:

```yaml
generation:
  execution:
    partition_size_bp: 5000000
    unit_bucket_size: 250
    parquet_compression: zstd
```

- Smaller genomic partitions reduce memory per preparation job but create more
  scheduler work.
- Smaller unit buckets reduce reducer memory but increase the number of
  partition-by-bucket jobs.
- Keep the manifest order deterministic; changing partition or bucket sizing
  changes work layout, not the scientific calculation.

Both supplied execution profiles enable `rerun-incomplete`. If a cluster job
fails, correct the underlying problem and rerun the same Snakemake command.
Completed valid shards are reused. Do not manually create `validation.ok` or
edit generated manifests to bypass a failed check.

## Validation and troubleshooting

### Schema failure before DAG construction

Check the mode-specific shape of the manifest. Common causes are a missing
prefix, both `ped` and `samples` being set, a missing mutation-rate binding, or
generation being run with `include_indels: true`.

### No Slurm executor named `slurm`

Install `snakemake-executor-plugin-slurm` in the same virtualenv as Snakemake
and verify that the virtualenv is active on the submission host.

### Python import failure on a compute node

Verify that `precommand` activates the absolute shared virtualenv path and that
the Python module used to create the virtualenv is available on compute nodes.
Test the imports inside a short interactive allocation if necessary.

### Module not found

Run `module spider <tool>` on GREX, then update `workflow/envmodules.yaml` with
the site-provided module name. Generate mode uses the `validation` route;
apply mode also uses `bioinformatics` and `r`.

### Mutation-rate fallback fails

The mutation-rate table must contain every absent allele encountered in the
calibration space. If chromosome-X fallback is intended, set
`mutation_rate_policy: sex_averaged_proxy`; the default `reject` policy stops
on absent chromosome-X alleles.

### A unit is missing from the final lookup

Inspect `excluded_units.tsv.gz` and the bucket diagnostics. Units are removed
when retained fraction is below the configured minimum or when an applicable role has
zero event mass.

### Apply validation passes but annotation stops

Read the annotation QC JSON and rule log. The workflow deliberately stops if
CADD or testing-unit annotation matches no relevant variants, since continuing
would create empty or misleading scores.

## Testing changes

The test runner exposes increasingly broad tiers:

```bash
python tests/run_tier.py unit
python tests/run_tier.py schema
python tests/run_tier.py workflow
python tests/run_tier.py e2e
python tests/run_tier.py resilience
```

Use `python tests/run_tier.py all` for the complete suite. At minimum:

- library changes should pass `unit`;
- configuration and planning changes should pass `schema` and `workflow`;
- rule, script, or cross-mode changes should pass `e2e`;
- checkpoint, certification, or recovery changes should also pass
  `resilience`.

## Developer map

| Area | Main files |
| --- | --- |
| Entry point and mode selection | `workflow/Snakefile` |
| Immutable configuration plan | `workflow/lib/python/workflow_plan.py` |
| Apply resource catalog | `workflow/lib/python/apply_resources.py` |
| Generation rules | `workflow/rules/generate/` |
| Apply rules | `workflow/rules/apply/` |
| Generation adapters | `workflow/scripts/generate/` |
| Apply adapters | `workflow/scripts/apply/` |
| Shared adapter bootstraps | `workflow/shared/` |
| Python implementation | `workflow/lib/python/` |
| R implementation | `workflow/lib/r/` |
| Configuration schemas | `workflow/schemas/` |
| Software definitions | `workflow/envs/`, `workflow/envmodules.yaml`, `workflow/containers.yaml` |
| Execution profiles | `profiles/` |
| Fixtures and regression tests | `tests/` |

Keep Snakemake adapters thin. Put reusable scientific and file-format logic in
the Python or R libraries, resolve physical resources once through
`WorkflowPlan` and `ApplyResourceRegistry`, and treat the generated bundle
contract as a versioned public interface.
