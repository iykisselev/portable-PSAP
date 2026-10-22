# Apply workflow architecture

This is the developer guide to `workflow.mode: apply` and to the cohort-scoring
half of `workflow.mode: full`. Read [project concepts](project_concepts.md)
first if you want the scientific interpretation of PSAP.

## Entry point and resource resolution

The `workflow/Snakefile` creates one immutable `WorkflowPlan` before it loads
the apply rules. `rules/apply/paths.smk` establishes the output layout and the
single `APPLY_SCRIPT_DIR`; `resources.smk` builds an immutable
`ApplyResourceRegistry` from the plan. Validation, staging, and downstream
rules all use that registry, so each resource has one resolved source path.

| Layer | Location | Responsibility |
| --- | --- | --- |
| Configuration plan | `workflow/lib/python/workflow_plan.py` | Resolve apply/full inputs, configured versus effective references, lookup source, and output paths. |
| Resource catalog | `workflow/lib/python/apply_resources.py` and `rules/apply/resources.smk` | Describe sources, indexes, staging requirements, origins, and stable prepared paths. |
| Paths/sample discovery | `rules/apply/paths.smk`, `sample_discovery.smk` | Define output paths and create safe DAG callbacks from configured metadata. |
| Rule stages | `rules/apply/*.smk` | Validate, stage, filter, annotate, classify, score, and aggregate a cohort. |
| Script adapters | `workflow/scripts/apply/` | Thin Python/R adapters invoked by Snakemake rules. |
| Shared adapter support | `workflow/shared/_library.py`, `workflow/shared/bootstrap.R` | Locate the shared Python/R workflow libraries. |

Put new apply adapters in `workflow/scripts/apply/`. Put reusable
implementation in `workflow/lib/python` or `workflow/lib/r`, and reserve
`workflow/shared` for support code that is genuinely shared by adapters.

## Lookup source and full mode

Apply mode consumes a certified `partitioned_null_bundle_v1` directory. The
bundle supplies partition/model lookup tables, metadata, validation and marker
files, resolved and eligible partition BEDs, calibration regions, and excluded
alleles.

When the source is a bundle, `workflow_plan._effective_bundle_reference()`
uses the bundle's calibration mask and partition ownership artifacts. The
configured exclusion is already reflected in that mask and is not reapplied.
FASTA and CADD remain configured physical resources. In `full` mode, the source
is always the bundle produced by the same DAG.
The generated scope therefore becomes the apply scope, and the bundle's
`validation.ok` gates scoring.

## Inputs

| Input | Required | Contract used by apply |
| --- | --- | --- |
| Cohort VCF and optional index | Yes | Normalized/split VCF with `GT` and sample IDs matching metadata. |
| PED or sample sheet | Exactly one | Sample ID, phenotype, and sex used to build the DAG and classify genotypes. |
| FASTA, `.fai`, and optional `.gzi` | Yes | Assembly validation and VCF/resource compatibility. |
| Calibration-region resource | Bundle-derived | Callable interval filter. |
| Partition BEDs | Bundle-derived | Resolved ownership and eligible unit/partition membership. |
| CADD SNV and optional indel resource | SNV yes; indel when `include_indels` | Allele-level scores used for lookup-bin selection. |
| Exclusions | Optional | BED/TSV/VCF intervals or exact alleles removed from the cohort. |
| Lookup source | Yes | A certified partitioned generated bundle. |

Indexes are optional for stageable VCF, BED, and TSV resources. The workflow
creates indexed working copies under `prepared/`; FASTA companion indexes are
provided by the user and are not regenerated.

## DAG and transformations

The default `all` target writes the candidate, omnibus, and control reports,
one report for each configured partition, and endpoint/aggregation diagnostics.
The stages below are listed in dependency order, although the independent
validation rules can run in parallel.

| Stage / rules | Inputs | Transformation | Primary outputs |
| --- | --- | --- | --- |
| Validation: `validate_cohort`, `validate_resources`, `validate_model`, `merge_validation`, `require_valid_inputs` | Cohort/metadata/FASTA, resources, lookup source, schemas | Verify formats, contigs, sample agreement, model, lookup/bundle contract, and required role coverage; merge reports and gate expensive rules. | `validation/*.json`, `validation.json`, `validation.ok`, `sample_manifest.tsv` |
| Staging: `prepare_vcf_resource`, `prepare_table_resource` | Sources and optional indexes from `ApplyResourceRegistry` | Normalize chromosome names when needed, sort/compress/index stageable VCF/BED/TSV resources. | `prepared/cohort.vcf.gz`, staged tables and indexes |
| Exclusion projection: `combine_exclusions` | Staged cohort and every configured exclusion | Project full-span intervals and exact alleles to a cohort-aligned exclusion VCF. | `prepared/combined_exclusions.vcf` |
| Variant filtering: `filter_calibration_regions`, `filter_excluded_records` | Staged cohort, optional bundle calibration regions, combined exclusions | Restrict to generated calibration regions when present, then subtract excluded records. | In-scope BCF/CSI and filtered VCF/CSI |
| Annotation: `annotate_variants`, `summarize_annotations`, `require_valid_annotations` | Filtered VCF, CADD, resolved partition BEDs | Attach CADD and unit/partition memberships, then reject complete expected annotation misses. | Annotated VCF and annotation QC |
| Per-sample processing: `split_genotypes`, `extract_gt`, `prepare_sample_genotypes` | Annotated VCF | Split non-reference calls, retain the two INFO fields, normalize ploidy/alleles, assign sex-aware genotype classes, and expand unit/partition annotations in memory. | Temporary `raw/` calls and `annotated/` tables |
| Scoring: `apply_psap` | Per-sample annotations and lookup source/certificate | Keep compatible best-scoring variants and map score bins to partition/model PSAP values. | Partition-aware files under `scored/` |
| Cohort aggregation: `build_sample_candidates`, `build_control_reference`, `build_cohort_report` | Affected scores, unaffected annotations/scores, manifest | Build candidates, preserve control evidence, compare exact matching variants, and aggregate reports. | Candidate intermediates plus three final cohort reports |

## Data and join keys

The canonical allele key contains chromosome, position, reference, and
alternate allele. The normalized membership table adds the unit and owned
partition IDs, so a variant overlapping several units produces several rows.
Per-sample tables carry that allele/unit/partition tuple through classification
and scoring.

The final control comparison is intentionally strict. It uses the exact
`unit`/`partition_id`/`variant_key` tuple, not merely the same gene or position. The
dominant-heterozygous exception allows evidence from any configured control
model; other candidate models require a matching model. Keep this behavior if
you refactor report aggregation.

## Output layout

For `outdir` and `run.apply_prefix`, results are rooted at
`<outdir>/apply/<apply_prefix>/`; logs live below
`<outdir>/logs/apply/<apply_prefix>/`.

| Path | Contents |
| --- | --- |
| `validation.json`, `validation.ok` | Merged preflight result and success gate. |
| `sample_manifest.tsv` | Canonical validated sample metadata. |
| `prepared/` | Staged/indexed resources and cohort-aligned exclusions. |
| `variants/`, `annotations/` | Filtered variants, annotated VCF/table, and annotation QC. |
| `raw/`, `standardized/`, `sex_classified/`, `annotated/`, `scored/` | Per-sample intermediates. |
| `candidates/` | Affected-only candidate intermediates. |
| `cohort/cohort_candidates.tsv.gz` | Candidate-level affected evidence and control comparison. |
| `cohort/cohort_summary.tsv.gz` | Per-affected-individual, testing-unit/model ACAT omnibus results. |
| `cohort/partitions/<partition>.summary.tsv.gz` | Per-partition affected evidence summaries. |
| `cohort/{endpoint,aggregation}_diagnostics.tsv.gz` | Endpoint and raw-versus-mid-p ACAT diagnostics. |
| `cohort/control_variants.tsv.gz` | Unaffected participant evidence. |
| `benchmarks/` | Benchmark output for selected high-cost rules. |

## Rules worth protecting when changing apply

1. **Validate before cost.** Rules that read or transform full cohort data must
   depend, directly or through staging, on `validation.ok`.
2. **Use the registry.** Add resources through `ApplyResourceRegistry`; do not
   add a rule-local source path or a second path-resolution mechanism.
3. **Use the effective reference.** Bundle and full mode must stage/filter the
   bundle's calibration scope, eligible units, and excluded alleles—not the
   stale configured equivalents.
4. **Keep sample wildcards safe.** `sample_discovery.smk` builds wildcard lists
   from configured metadata, while the validated manifest gates processing.
   Preserve this split so malformed metadata fails clearly without weakening
   DAG construction.
5. **Preserve one-row-per-variant/unit/partition semantics.** Do not collapse
   owned memberships before scoring or control comparison.
6. **Keep adapters thin.** Put cross-rule logic in libraries and use
   `APPLY_SCRIPT_DIR` for every apply script path.

## Rule authoring conventions

Use `env`, `mod`, and `sif` from `rules/deployment.smk` together when adding a
rule, selecting one logical tool family. Keep rule-specific threads and
resource-policy declarations visible in the rule. Construct optional shell
arguments with a tested helper such as `optional_regions_argument()`; do not
hand-build unquoted optional fragments in a shell body.

## Where common changes belong

| Change | Main implementation points | Required follow-through |
| --- | --- | --- |
| Add a stageable resource | `workflow_plan.py`, `apply_resources.py`, schema, `resources.smk` | Define source/index/kind/origin; validate, stage, and add fixture tests. |
| Add an annotation column | Annotation adapter/library and downstream R/Python consumers | Document its missing-value policy and preserve the canonical allele/unit key. |
| Add a disease role | Model/lookup validation, genotype classification, scorer, and generation support | Define compatible genotypes/chromosomes and bundle contract compatibility. |
| Change reporting | `cohort.smk` adapters and report configuration | Keep candidate IDs stable or version the report contract; update expected fixture digests. |

After an apply change, run the end-to-end fixture:

```bash
python -m unittest discover -s tests -p 'test_apply_workflow.py'
```
