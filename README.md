# Portable PSAP

Portable PSAP is a Snakemake workflow [1] that generates certified partitioned
PSAP null-distribution bundles, applies a certified bundle to a cohort VCF, or
performs both steps in one full workflow. It prioritizes variants under selected
disease models using the single-case PSAP framework and its genomic-regions
extension [3,7].

## Workflow modes

Every manifest explicitly selects one mode with `workflow.mode`:

| Mode | Inputs | Default target | Result |
| --- | --- | --- | --- |
| `apply` | Cohort, reference/annotation resources, certified partitioned bundle | `all` | Cohort candidate, summary, and control reports |
| `generate` | Population AF, score, scope, and reference resources | `all` | Certified generated lookup bundle |
| `full` | Cohort plus one shared generation/application binding set | `all` | Generated bundle followed by cohort reports |

These are the only supported explicit Snakemake targets. The workflow does not
maintain compatibility aliases for result targets.

`workflow/lib/python` is an internal flat-module library, not a public Python
package. Workflow scripts add that directory through `workflow/shared/_library.py`
and import modules such as `generation_scope`; package-qualified imports are
unsupported.

## Partitioned SNV/indel null generation

Generated null distributions can include SNVs alone or the combined SNV/indel
candidate universe. A combined run sets `model.include_indels: true`, declares
`generation_settings.variant_classes: [snv, indel]`, and supplies an indexed
indel CADD resource (`generation.score_indel` in standalone `generate` mode,
or `bindings.reference.cadd.indel` in `full` mode). Candidate indels are
normalized, left-aligned, reference-validated alleles from that score resource;
their complete REF spans must remain inside the calibration mask.

Every new generate or full run writes `partitioned_null_bundle_v1`, whether it
contains SNVs alone or SNVs and indels. The bundle records its score resources,
candidate-universe and span conventions, partition order, overlap ownership,
and one lookup table for every enabled partition/model pair. Apply accepts only
this certified contract; regenerate legacy bundles before use.

For `unobserved_policy: mutation_rate`, `generation.mutation_rate_indel` is an
optional indexed allele-specific indel rate source. If it is not supplied,
unobserved indels contribute AF 0; SNV mutation-rate requirements are
unchanged. Per-bucket generation diagnostics report `snv_event_count` and
`indel_event_count` for each unit/role, making the two contributions auditable
without changing lookup semantics.

Use the executable small demonstrations in
[`runs/apply.example.yaml`](runs/apply.example.yaml),
[`runs/generate.example.yaml`](runs/generate.example.yaml), and
[`runs/full.example.yaml`](runs/full.example.yaml) to exercise each mode. For a
production apply run, copy [`runs/run.example.yaml`](runs/run.example.yaml) and
replace its physical resource paths.

To apply an existing generated bundle, copy
[`runs/apply.bundle.example.yaml`](runs/apply.bundle.example.yaml). Supply its
directory under `bindings.reference.lookup_source`:

```yaml
lookup_source:
  type: bundle
  directory: results/generate/my-model/bundle
```

The workflow discovers the partition/model lookup tables, eligible partitions,
calibration mask, exclusions, contract, and score grid from the bundle
declarations. No source partition BED or repeated component paths are accepted
in apply mode. Paths are relative to the directory where you invoke Snakemake.
The bundle must already exist for this apply form; `full` generates and binds
its own bundle automatically. Model choices (`include_indels`, `hemizygote`,
testing-unit type) and assembly must match the bundle. Certificate and checksum
checks run before annotation. The bundle's score maximum and bin width govern
scoring; configuration defaults do not change them. See
[bundle configuration and deployment](docs/pipeline.md#applying-an-existing-bundle)
for migration and verification details.

## Workflow overview

### 1. Validation

Before processing variants, the workflow runs three validation stages:

- Cohort validation checks the VCF, sample metadata, reference FASTA and
  indexes, sample-ID agreement, VCF normalization/splitting, and sex-chromosome
  requirements.
- Resource validation checks resource formats, chromosome compatibility, CADD
  columns, certified partition artifacts, and configured exclusion files.
- Model validation checks model settings and verifies the bundle certificate,
  partition/model lookup tables, and testing-unit identities.

Each stage writes a JSON report; the reports are merged when validation
finishes. Errors stop the run before the computationally expensive rules begin.
Warnings remain in `validation.json` and do not stop the run.

### 2. Resource staging and variant filtering

The workflow creates indexed working copies of inputs below the run output 
directory. When necessary, it sorts and BGZF-compresses the VCF and
tabular genomic resources, then creates indexes for the staged files. It then:

1. retains variants in the generated bundle's calibration regions, when a
   bundle is provided;
2. removes variants matching the combined exclusion set.

Exclusions may be intervals or exact alleles. For interval exclusions, the
workflow considers the complete VCF span of a variant, including the `REF`
length and any `INFO/END` value. A deletion is therefore removed if any part of
its span overlaps an excluded region.

### 3. VCF annotation with CADD and owned partitions

The filtered VCF receives CADD-PHRED scores and owned testing-unit/partition
memberships. SNVs always use the configured SNV CADD table. Indels use the
configured indel table only when `model.include_indels` is `true`. CADD scores
and their splice-focused extension are described in [2,5].

Annotation QC stops the run if an expected resource annotates no relevant
records:

- CADD annotates no SNVs;
- enabled CADD indel annotation annotates no indels; or
- the certified partition scope annotates no cohort variants.

The annotation table is normalized to one row for each `variant_key`,
testing-unit, and partition combination. A variant may belong to several
testing units, but overlapping bases within a testing unit are owned by the
earliest configured partition.

### 4. Genotype extraction and classification

For each processed sample, the workflow extracts `GT` calls, normalizes phased
and unphased diploid notation to `a/b`, then assigns a genotype class:

- `het` for heterozygous diploid calls;
- `hom` for homozygous alternate diploid calls; and
- `haploid` for haploid X/Y calls when the hemizygote model is enabled.

Reference homozygous, missing, and duplicate calls are excluded or rejected. Haploid
calls outside X/Y are rejected. The sex model uses PED-coded sex values (`1` for
male, `2` for female):

- male X/Y homozygous or haploid alternate calls become `hem` when
  hemizygote scoring is enabled;
- male X/Y heterozygous calls are removed when the hemizygote model is
  enabled;
- female X heterozygous and homozygous calls become `het_female` and
  `hom_female`, respectively; and
- haploid calls are not accepted for samples with unknown sex.

### 5. Variant scoring with PSAP

Each classified genotype is inner-joined with its CADD and testing-unit/
partition annotations. For every enabled disease model, testing unit, and
partition, the workflow:

1. retains variants in the model's compatible genotype class;
2. keeps every variant tied for the highest CADD-PHRED value in that
   unit/partition;
3. clips CADD-PHRED values to the resolved lookup interval `[0, score_max]`;
4. maps each value to an evenly spaced column in the corresponding PSAP lookup
   grid; and
5. records the lookup value as `PSAP`.

Ties are retained and sorted by variant key. Every co-best variant proceeds to
candidate generation and comparison with controls.

### 6. Variant prioritization and aggregation

For every affected sample, the workflow creates one candidate for each scored
variant. Candidate IDs are deterministic and incorporate the partition-aware
candidate ordering.

For every unaffected sample, the workflow retains evidence for each carried
variant under compatible disease models. It attaches the lowest available PSAP
value for every exact `unit`/`partition`/`model`/`variant_key` combination.

The final cohort report compares each affected candidate with controls that
carry the same exact variant key in the same testing unit and partition. A
dominant-heterozygous candidate may use evidence from any configured control
model; all other models require a matching model. The main cohort summary
combines available partition PSAPs once per affected individual, testing unit,
and model with ACAT [8].

## Inputs

Create a production apply manifest by copying
[`runs/run.example.yaml`](runs/run.example.yaml):

```bash
cp runs/run.example.yaml runs/my_analysis.yaml
```

Paths in the run manifest are resolved relative to the directory from which
you invoke Snakemake. Together, `outdir` and `run.apply_prefix` determine all
output paths.

### Cohort VCF

`inputs.vcf.data` is a standard VCF or BGZF-compressed VCF containing the
cohort samples. It needs a valid header and sample columns whose IDs match the
selected metadata file. Variant records must be normalized and split against
the reference FASTA. The workflow uses `CHROM`, `POS`, `REF`, `ALT`, and `GT`.

`inputs.vcf.index` may point to a `.tbi` or `.csi` index, or be `null`. If the
source index is missing, the workflow issues a warning and creates an indexed
staged VCF internally.

### Cohort metadata

Set one of `inputs.ped` and `inputs.samples`.

#### PED

The PED file is tab-delimited and must provide at least the first six
conventional columns. The workflow reads:

| PED column | Required meaning |
| --- | --- |
| 1 | Family ID |
| 2 | Sample ID (`pid`) matching a VCF sample |
| 5 | Sex: `1` male, `2` female, another value unknown |
| 6 | Phenotype: `1` unaffected, `2` affected, another value excluded |

#### Sample sheet

The sample sheet is a headered, tab-delimited TSV with these required columns:

```tsv
family_id	pid	sex	phenotype
FAM001	case_001	1	2
FAM001	control_001	2	1
```

`pid` cannot contain a slash or backslash. Each processed `pid` must appear
exactly once in both the VCF and the metadata.

### Reference FASTA

`bindings.reference.fasta` requires:

```yaml
fasta:
  data: resources/reference.fa.gz
  fai: resources/reference.fa.gz.fai
  gzi: resources/reference.fa.gz.gzi
```

An uncompressed FASTA still requires `fai`, but may set `gzi: null`. A
compressed FASTA must be BGZF-compressed and has both `.fai` and `.gzi`.

### Genomic interval and annotation resources

The workflow can stage these resources without source indexes and create the
needed `.tbi` indexes itself. Providing an existing `.tbi` or `.csi` index can
avoid some staging work.

| Binding | Required structure |
| --- | --- |
| `exclusion` (optional) | Genomic intervals to remove, independent of source. Generation uses zero-based, half-open BED intervals. Set `data: null` and `index: null` when absent. |
| `testing_unit_partitions` (generate/full) | Ordered bindings, each with an ID and BED containing `CHROM`, zero-based start, end, and a non-empty testing-unit ID. Earlier bindings own overlapping bases. Apply reads these artifacts from the certified bundle and rejects source-partition overrides. |
| `exclusions` | An ordered list of BED, TSV, or VCF resources. BED rows are zero-based, half-open intervals. TSV rows are either `CHROM, start, end` with 1-based inclusive coordinates or `CHROM, POS, REF, ALT` exact variants. VCF exclusions match exact alleles. |
| `cadd.snv` | Tab-delimited CADD SNV table with `Chrom`, `Pos`, `Ref`, `Alt`, `RawScore`, and `PHRED` as the first six data columns. |
| `cadd.indel` | Same required CADD layout for indels when `model.include_indels: true`; otherwise it may be `null`. |

The exclusion identity is source-neutral; lookup identities name their
resource releases. When substituting a resource, reflect assembly and release 
changes in the run documentation; gnomAD constraint resources are described
in [4,6].

### Certified lookup bundle

Apply mode requires `bindings.reference.lookup_source` to name a certified
`partitioned_null_bundle_v1` directory. The bundle contains a
`testing_unit_id`-headed lookup table for every enabled partition/model pair,
along with resolved and eligible partition BEDs, the calibration mask,
exclusions, metadata, and certificate. Its partition order and score grid are
authoritative; apply manifests cannot override them.

## Outputs

The main reports are gzip-compressed TSVs located in:

```text
<outdir>/apply/<apply_prefix>/cohort/
```

### `cohort_candidates.tsv.gz`

One row per affected candidate variant after partition-aware control comparison.
It includes `partition_id` in addition to candidate, sample, unit, model,
variant, PSAP, and matching-control fields.

### `cohort_summary.tsv.gz`

One row per affected individual, testing unit, and model. `PSAP` is the ACAT
omnibus result across available partitions, not a minimum across candidates.
The report records the aggregation status, configured/eligible/contributing
partition counts, and ordered JSON descriptions of contributing IDs, raw
partition PSAPs, inputs, and statuses.

### Partition and diagnostic reports

`cohort/partitions/<partition>.summary.tsv.gz` contains the corresponding
per-partition affected evidence. `endpoint_diagnostics.tsv.gz` and
`aggregation_diagnostics.tsv.gz` retain raw and mid-p diagnostics. The
production aggregation policy is raw survival probabilities; the diagnostics
exist to support calibration of any future alternative policy against real
bundles and cohorts.

### Optional candidate-review annotations

`scripts/annotate_cohort_summary.py` is a standalone post-analysis utility. It
does not modify workflow outputs: it appends CADD PHRED, a selected population
AF INFO field, and `gene_name` from a GTF to a copy of `cohort_candidates.tsv.gz`.
It treats `unit` as an Ensembl gene ID (and ignores an optional version suffix).
All CADD tables and the population VCF must be BGZF-compressed and indexed.

```bash
python scripts/annotate_cohort_summary.py \
  results/apply/my-run/cohort/cohort_candidates.tsv.gz \
  results/apply/my-run/cohort/cohort_candidates.review.tsv.gz \
  --cadd-snv reference/cadd_snv.tsv.gz \
  --cadd-indel reference/cadd_indel.tsv.gz \
  --population-vcf reference/population.vcf.gz \
  --af-field AF \
  --gtf reference/gencode.gtf.gz
```

The output preserves every input row and adds `gene_name`, `cadd_phred`, and
`population_af` for the candidate's `variant_key`.

### `control_variants.tsv.gz`

One row for each carried variant/model evidence record from an unaffected
control.

| Column | Meaning |
| --- | --- |
| `pid`, `phenotype`, `sex` | Control sample metadata. `phenotype` is `unaffected`. |
| `unit`, `partition_id`, `model`, `variant_key` | The testing unit, owned partition, compatible disease model, and carried variant. |
| `PSAP` | Lowest available PSAP value for the exact partition-aware evidence record; may be missing when no score is available. |

### Supporting outputs

The same `<outdir>/apply/<apply_prefix>/` directory contains supporting files:

- `validation.json` and `validation.ok`: input/resource/model validation;
- `sample_manifest.tsv`: processed sample metadata;
- `prepared/`: staged VCFs, tables, and indexes;
- `variants/`: calibration-scope-filtered and exclusion-filtered variants;
- `annotations/`: annotated VCF, flat annotation table, and annotation QC;
- `raw/`, `standardized/`, `sex_classified/`, and `annotated/`: per-sample
  genotype intermediates;
- `scored/`: per-sample PSAP tables; and
- `benchmarks/` and `<outdir>/logs/`: rule runtime/memory reports and logs.

## Setup and run

Run the commands below from the repository root. Paths in a run manifest are
interpreted relative to the directory where Snakemake is started, so running
from another directory changes the meaning of relative paths.

### 1. Install the workflow runner

Install a Conda-compatible package manager (`mamba`, `micromamba`, or Conda)
and use the supplied environment to install Snakemake and its compatible 
executor plugin. `mamba` is used below; substitute your preferred manager 
if needed.

```bash
mamba env create --file workflow/envs/snakemake.yaml
conda activate portable-psap-snakemake
snakemake --version
```

The default profile creates the per-rule software environments automatically
from `workflow/envs/` on the first run. It therefore needs internet or a local
Conda channel mirror on that first run, as well as write access to
`.snakemake/conda` in the repository. The supplied environments contain the
required Python, R, bcftools, bgzip, and tabix tools. No manual installation
of those individual programs is needed for the local Conda route.

### 2. Create and configure a run manifest

Copy the template and give the invocation a unique output prefix:

```bash
cp runs/run.example.yaml runs/my_analysis.yaml
```

Edit `runs/my_analysis.yaml` and replace every example path. At a minimum,
configure the following:

| Manifest setting | What to provide |
| --- | --- |
| `outdir` and `run.apply_prefix` | A writable results directory and a unique name for this cohort/run. The pipeline writes to `<outdir>/apply/<apply_prefix>/`. |
| `inputs.vcf` | The cohort VCF and, when available, its `.tbi` or `.csi` index. |
| `inputs.ped` or `inputs.samples` | Exactly one metadata source. The selected sample IDs must match the VCF header exactly. |
| `bindings.reference.fasta` | The reference FASTA plus its required `.fai`; a BGZF-compressed FASTA also needs its `.gzi`. |
| `exclusion`, testing-unit, CADD, and lookup bindings | Resources compatible with the same assembly and testing-unit identifiers as the scientific configuration. |

Source indexes for the VCF and tabular genomic resources are optional: when an
index is omitted, the workflow can and index a working copy below the run
output directory. The source FASTA indexes are not optional. Keep input
resources outside the output directory and treat them as read-only.

The scientific settings are located in `config/config.yaml`, separately from the
run manifest. Change them only when deliberately changing the PSAP model or the
referenced resource release. In particular, if `model.include_indels` is enabled,
provide the indel CADD table; if `model.hemizygote` is enabled, provide the chromosome-X
hemizygous lookup table.

### 3. Local run

The default profile is intended for a workstation and enables Conda deployment,
8 cores, 32 GB of total memory, and automatic reruns of incomplete jobs. Change
`--cores` to the number of cores that may be used on the machine; command-line
values override the profile default.

```bash
snakemake \
  --profile profiles/default \
  --configfile runs/my_analysis.yaml \
  --cores 8
```

To execute a bundled demonstration, select its manifest directly. `generate`
publishes `<outdir>/generate/<generation_prefix>/bundle/`; `full` additionally
publishes `<outdir>/apply/<apply_prefix>/cohort/` after certifying and consuming
that bundle.

```bash
snakemake --configfile runs/generate.example.yaml --cores 2 --use-conda
snakemake --configfile runs/full.example.yaml --cores 2 --use-conda
```

### 4. Run on Slurm with environment modules

For the GREX Slurm configuration, edit
`profiles/grex/config.yaml` to set your account, partition, and node-local
temporary directory. Confirm that the module names in `workflow/envmodules.yaml`
are available at the site; adjust that mapping when the site provides different
module names. The Snakemake controller environment from step 1 is still needed
on the submission host.

Dry-run the scheduler configuration first:

```bash
snakemake \
  --profile profiles/slurm \
  --profile profiles/grex \
  --workflow-profile none \
  --configfile runs/my_analysis.yaml \
  -n
```

Submit the workflow by removing `-n`:

```bash
snakemake \
  --profile profiles/slurm \
  --profile profiles/grex \
  --workflow-profile none \
  --configfile runs/my_analysis.yaml
```

Tune `jobs`, `default-resources.mem_mb`, and `default-resources.runtime` in
the profile if they do not match the site or project allocation.

## Generated bundles and supporting outputs

`generate` and `full` publish a `partitioned_null_bundle_v1` under
`<outdir>/generate/<generation_prefix>/bundle/`. The bundle contains one
partition/model lookup table for every enabled pair, resolved and eligible
partition BEDs, calibration and exclusion artifacts, `metadata.json`,
`validation.json`, and `validation.ok`. The marker is written only after
component size/checksum and semantic validation succeeds. `full` automatically
binds this bundle to application; an `apply` manifest may name another
previously certified partitioned bundle explicitly.

<<<<<<< HEAD
=======
For a combined v3 bundle, retain `candidate_audit.tsv.gz` with the certificate
when copying or archiving the bundle. It records structural indel rejections
(for example, a candidate whose REF span crosses a calibration-mask boundary)
that are intentionally distinct from population FILTER exclusions.

>>>>>>> main
See the [pipeline guide](docs/pipeline.md#generated-bundle) for the generated
bundle interface, the [generate architecture](docs/generate_workflow_architecture.md)
for generation ownership, and the [apply architecture](docs/apply_workflow_architecture.md)
for application ownership.

## Testing and contributing

Run a documented test tier from the repository root with, for example,
`python tests/run_tier.py workflow`. The null-distribution fixtures are
described in [tests/null_distribution/README.md](tests/null_distribution/README.md).
Module boundaries and contract ownership are documented in the
[generate architecture](docs/generate_workflow_architecture.md) and
[apply architecture](docs/apply_workflow_architecture.md).

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
