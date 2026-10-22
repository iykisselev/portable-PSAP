# Generate workflow architecture

This is the developer guide to `workflow.mode: generate` and to the generation
half of `workflow.mode: full`. Start with [project concepts](project_concepts.md)
if you want the scientific reasoning before the implementation details.

## Entry point and module ownership

The `workflow/Snakefile` first normalizes and validates the merged
configuration. It then builds one immutable `WorkflowPlan` in
`workflow/lib/python/workflow_plan.py`. Generation rules are loaded only in
`generate` and `full` modes.

| Layer | Location | Responsibility |
| --- | --- | --- |
| Configuration plan | `workflow_plan.py` | Resolve mode, bindings, policies, output roots, and full-mode resource ownership. |
| Rule context | `rules/generate/context.smk` | Turn the plan into paths, manifests, wildcards, roles, and callbacks. |
| DAG stages | `rules/generate/{prepare,calculate,finalize}.smk` | Declare the checkpoint, window/bucket jobs, bundle rendering, and certification. |
| Script adapters | `workflow/scripts/generate/` | Decode Snakemake values and call library code. |
| Domain library | `workflow/lib/python/generation_*.py` | Scope construction, population parsing, probability reduction, lookup rendering, and bundle finalization. |
| Shared adapter support | `workflow/shared/_library.py` | Add the workflow Python library to adapter import paths. |

Keep scientific and file-format logic in `workflow/lib/python`. The scripts
should stay thin: they decode Snakemake inputs and call library functions. That
keeps the useful logic testable without constructing a live Snakemake object.

## Inputs

In standalone `generate` mode, the following bindings live under `generation`:

| Input | Required | Used for |
| --- | --- | --- |
| `population_vcf` and optional index | Yes | Population allele states and configured AF INFO tags. |
| `score_snv` and optional index | Yes | Allele-level CADD score for each calibrated SNV. |
| `score_indel` and optional index | When `model.include_indels: true` | Normalized, finite indel candidate universe and CADD score. |
| `exclusion` | Binding required; data may be null | Intervals to subtract, independent of source. Combined with other interval exclusions. |
| `testing_unit_partitions` | Yes | Ordered biological-partition BEDs with stable unit IDs. Earlier entries own overlap bases. |
| `fasta` and `.fai` | Yes | Contig lengths, windowing, and reference allele enumeration. |
| `par_regions` | Yes | Pseudoautosomal intervals excluded before X-role assignment. |
| `exclusions` | Optional | Additional intervals excluded from calibration. |
| `mutation_rate` and optional index | Conditional | Required for an unobserved `mutation_rate` policy. |
| `mutation_rate_indel` and optional index | Optional | Allele-specific indel mutation rate; absent indels use AF 0 if it is omitted. |
| `af_source` | Yes | Names the overall, and when configured female/male, AF INFO tags. |

The current scope supports SNVs and, when enabled, normalized indels on contigs
`1`–`22` and `X`. Every calibration allele must have exactly one score row.
Indel candidates are owned by their anchor-position window, while a bounded
right halo ensures a REF span crossing that window boundary is evaluated
once against the complete mask and testing-unit span. Invalid selected AF
values, duplicate population states, duplicate scores, and missing SNV mutation
rates when they are needed are all hard failures.

### Full-mode binding behavior

Full mode deliberately does not define a second reference scope under
`generation`. It reuses FASTA, exclusion, testing units, CADD SNVs (and CADD
indels when enabled), and ordinary variant exclusions from
`bindings.reference`. Exactly one configured exclusion must have the
`par_regions` role; every `variants_exclude` binding becomes a generation
exclusion. `WorkflowPlan` enforces this while the DAG is being built.

## DAG and transformations

The scope checkpoint discovers the window and bucket wildcards. Every later
rule reads those JSON manifests instead of trying to rediscover work from the
filesystem.

| Rule | Inputs | Transformation | Outputs |
| --- | --- | --- | --- |
| `prepare_generation_scope` (checkpoint) | FASTA index, ordered partition BEDs, exclusion/PAR/other exclusions, population/score bindings | Resolve per-unit partition ownership in configured order; subtract exclusions; determine pair eligibility; divide the eligible owned scope into windows; assign one bucket per testing unit. | Window and bucket manifests; resolved ownership, pair-status, ID-overlap, and calibration-mask Parquet files. |
| `prepare_generation_window` | One window plus scope artifacts, FASTA, population VCF, and scores | Enumerate calibrated SNV alleles and, when enabled, anchor-owned indel catalog entries; validate spans, join unit membership, score, and population state. | One each of `allele_states`, `allele_membership`, `allele_scores`, and indel candidate-audit Parquet shards. |
| `calculate_window_null_distribution` | One window shard set and one unit bucket | Assign the absent-allele policy and role-specific effective AF; bin scores; aggregate additive no-event statistics and class-specific event counts. | One window-by-bucket statistics Parquet shard. |
| `reduce_null_distribution_bucket` | All statistics shards for one bucket | Compensated reduction across windows; calculate conditional score-bin survival rows and diagnostics by biological partition and public disease model. | Partition/model lookup rows plus diagnostics TSV with total, SNV, and indel event counts. |
| `finalize_null_distribution_bundle` | All bucket rows/diagnostics plus scope and state artifacts | Retain only units with positive mass for every applicable role; stream public bundle components and metadata. | Gzipped lookup, scope, exclusion, and metadata files. |
| `validate_null_distribution_bundle` | Finalized bundle and schema | Validate component checksums, row coverage, contracts, and semantics. | `validation.json` and `validation.ok`. |

The public `all` target exposes the certificate and recoverable generation
intermediates. In full mode, it also exposes the downstream apply reports.

## Intermediate-data contracts

These checkpoint artifacts are part of the internal DAG contract:

- `resolved_partition_intervals.parquet` has one row per owned biological
  partition interval and includes the configured `partition_order`.
- `unit_partition_status.parquet` has one row per testing-unit/partition pair
  with raw, priority-removed, resolved, supported, and retained bases plus its
  eligibility status. Retention is calculated after precedence subtraction.
- `testing_unit_id_overlap.parquet` records cross-partition raw-ID overlap
  diagnostics; it does not alter identifiers or eligibility.
- `effective_calibration_mask.parquet` records the exact calibration intervals.
- `windows.json` and `unit_buckets.json` are versioned, sorted manifests;
  they are the only source of window and bucket wildcards.
  The window manifest uses schema version 2, a `windows` array, `window_id`, and
  `window_size_bp`; the bucket manifest uses schema version 3 and `window_ids`.
  Older execution-partition manifests must be regenerated. The public config
  key remains `generation.execution.partition_size_bp`, resolved once into the
  internal `window_size_bp` plan field.
- Window `allele_states` holds `present`, `absent`, `rejected`, or invalid
  population state information. `allele_membership` maps alleles to units, and
  `allele_scores` supplies the matching score.
- Statistics shards are additive per `testing_unit_id`, `partition_id`, role,
  and score bin.
  They store `certain_count`, `log_none_sum`, total event count, and SNV/indel
  event counts so global event tables never need to be materialized. Bucket
  diagnostics expose those class counts as `snv_event_count` and
  `indel_event_count` for every unit/partition/model combination.

Downstream rules should not glob these files or infer wildcard names from
directory contents. Use the callbacks in `rules/generate/context.smk`; they
wait for the checkpoint and validate the manifests before returning paths.

## Rule authoring conventions

Use `env`, `mod`, and `sif` from `rules/deployment.smk` together when adding a
rule, selecting one logical tool family. Keep rule-specific threads and
resource-policy declarations at the call site. Add recoverable generation
target artifacts through `generation_outputs()` rather than copying
checkpoint-aware callbacks into `Snakefile`.

## Bundle contract

The public bundle lives at
`<outdir>/generate/<generation_prefix>/bundle/`. Its filenames are defined in
one place, `bundle_contract.py`. The public components are:

- One `lookup.<partition>.<model>.tsv.gz` table for every enabled partition/model pair;
- one resolved and one eligible BED for every partition, plus the calibration-region BED;
- excluded alleles, candidate audit, scope-status and testing-unit-overlap audit tables;
- `metadata.json`, `validation.json`, and `validation.ok`.

`metadata.json` records contract and schema versions, assembly, score-grid
parameters, model catalog, partition order, overlap policy, absent-allele and
X policies, source and component checksums, and counts. Every fresh generation
run writes `partitioned_null_bundle_v1`; apply accepts only this certified
contract. Downstream rules use `validation.ok` as the completion boundary.

## Rules worth protecting when changing generation

1. **One binding set per full run.** In full mode, do not add a second
   generation-only reference source. Change the plan resolver if a new shared
   resource is needed.
2. **Determinism.** Unit IDs are sorted before bucket assignment; windows
   and reductions are manifest ordered. Preserve that ordering when adding
   parallel work.
3. **Bounded memory.** Window jobs only read one genomic window; bucket
   reduction reads one bounded unit range. Avoid collecting genome-wide allele
   or event tables in a new rule.
4. **Explicit policies.** A new absence, chromosome, score-binning, or
   genotype policy must be represented in configuration and emitted in bundle
   metadata.
5. **Certificate last.** Never create `validation.ok` before all public
   components and semantic checks are complete.
6. **Role symmetry.** A unit is publishable only if every role that applies to
   it has valid positive event mass. HEM applies only to surviving chromosome-X
   units.

## Where common changes belong

| Change | Main implementation points | Required follow-through |
| --- | --- | --- |
| Add a generation input | `workflow_plan.py`, config schema, `context.smk`, scope/window adapter and library | Include it in resource identities and bundle metadata; add fixture coverage. |
| Add a role | Generation settings/model plan, `generation_probability.py`, context roles, lookup/finalization validation | Define compatible chromosomes/genotypes, table filename, certificate semantics, and apply support. |
| Change binning | `generation_probability.py` and contract constants | Version the contract if compatibility changes; update apply lookup validation. |
| Change a public artifact | `bundle_contract.py`, finalization, validation schema/library | Preserve backward compatibility or deliberately version the contract. |

After a generation change, run the small fixture:

```bash
python -m unittest discover -s tests -p 'test_generate_workflow.py'
```
