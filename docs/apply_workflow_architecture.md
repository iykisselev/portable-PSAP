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

Apply mode can use external lookup tables or a certified bundle. The source
always supplies HET and HOM tables, plus the chromosome-X HEM table when
`model.hemizygote` is enabled. A bundle also supplies metadata, validation and
marker files, eligible units, calibration regions, excluded alleles, and
excluded units.

When the source is a bundle, `workflow_plan._effective_bundle_reference()`
uses the bundle's inclusion mask as `calibration_regions`, and replaces the
testing-unit BED and variant exclusions with bundle artifacts. The configured
`exclusion` is already reflected in that mask and is not reapplied. FASTA and
CADD remain the configured physical resources. For external lookups, the
source-neutral `exclusion` binding joins the ordinary exclusion sources; it is
never interpreted as an inclusion mask.
In `full` mode, the source is always the bundle produced by the same DAG. The
generated scope therefore becomes the apply scope, and the bundle's
`validation.ok` gates scoring.

## Inputs

| Input | Required | Contract used by apply |
| --- | --- | --- |
| Cohort VCF and optional index | Yes | Normalized/split VCF with `GT` and sample IDs matching metadata. |
| PED or sample sheet | Exactly one | Sample ID, phenotype, and sex used to build the DAG and classify genotypes. |
| FASTA, `.fai`, and optional `.gzi` | Yes | Assembly validation and VCF/resource compatibility. |
| Coverage resource | Optional in external apply; bundle-derived in full/bundle apply | Callable interval filter. |
| Testing-unit BED | Yes in configured bindings; bundle-derived when applicable | Fourth column is the unit ID joined to lookup rows. |
| CADD SNV and optional indel resource | SNV yes; indel when `include_indels` | Allele-level scores used for lookup-bin selection. |
| Exclusions | Optional | BED/TSV/VCF intervals or exact alleles removed from the cohort. |
| Lookup source | Yes | External HET/HOM(/HEM) tables or a certified generated bundle. |

Indexes are optional for stageable VCF, BED, and TSV resources. The workflow
creates indexed working copies under `prepared/`; FASTA companion indexes are
provided by the user and are not regenerated.

## DAG and transformations

The default `all` target writes the three cohort reports. The stages below are
listed in dependency order, although the independent validation rules can run
in parallel.

| Stage / rules | Inputs | Transformation | Primary outputs |
| --- | --- | --- | --- |
| Validation: `validate_cohort`, `validate_resources`, `validate_model`, `merge_validation`, `require_valid_inputs` | Cohort/metadata/FASTA, resources, lookup source, schemas | Verify formats, contigs, sample agreement, model, lookup/bundle contract, and required role coverage; merge reports and gate expensive rules. | `validation/*.json`, `validation.json`, `validation.ok`, `sample_manifest.tsv` |
| Staging: `prepare_vcf_resource`, `prepare_table_resource` | Sources and optional indexes from `ApplyResourceRegistry` | Normalize chromosome names when needed, sort/compress/index stageable VCF/BED/TSV resources. | `prepared/cohort.vcf.gz`, staged tables and indexes |
| Exclusion projection: `combine_exclusions` | Staged cohort and every configured exclusion | Project full-span intervals and exact alleles to a cohort-aligned exclusion VCF. | `prepared/combined_exclusions.vcf` |
| Variant filtering: `filter_calibration_regions`, `filter_excluded_records` | Staged cohort, optional bundle calibration regions, combined exclusions | Restrict to generated calibration regions when present, then subtract excluded records. | In-scope BCF/CSI and filtered VCF/CSI |
| Annotation: `annotate_variants`, `summarize_annotations`, `require_valid_annotations`, `prepare_testing_units` | Filtered VCF, CADD, testing BED | Annotate VCF/flat table, reject complete expected annotation misses, normalize to variant/unit rows. | Annotated VCF, annotation QC, `<unit>.tsv.gz` |
| Per-sample processing: `extract_gt`, `standardize_gt`, `apply_sex_model`, `join_sample_annotations` | Filtered VCF, sample manifest, unit table | Extract calls, standardize ploidy/alleles, assign sex-aware genotype class, and join annotations. | `raw/`, `standardized/`, `sex_classified/`, `annotated/` tables |
| Scoring: `apply_psap` | Per-sample annotations and lookup source/certificate | Keep compatible best-scoring variants and map score bins to role-specific PSAP values. | `scored/<sample>.<unit>.tsv.gz` |
| Cohort aggregation: `build_sample_candidates`, `build_control_reference`, `build_cohort_report` | Affected scores, unaffected annotations/scores, manifest | Build candidates, preserve control evidence, compare exact matching variants, and aggregate reports. | Candidate intermediates plus three final cohort reports |

## Data and join keys

The canonical allele key contains chromosome, position, reference, and
alternate allele. The normalized testing-unit table adds the unit ID, so a
variant overlapping several units produces several rows. Per-sample tables
carry that allele/unit pair through classification and scoring.

The final control comparison is intentionally strict. It uses the exact
`unit`/`variant_key` pair, not merely the same gene or position. The
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
| `cohort/cohort_summary.tsv.gz` | Aggregated affected candidate signatures. |
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
5. **Preserve one-row-per-variant/unit semantics.** Do not collapse unit
   memberships before scoring or control comparison.
6. **Keep adapters thin.** Put cross-rule logic in libraries and use
   `APPLY_SCRIPT_DIR` for every apply script path.

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
