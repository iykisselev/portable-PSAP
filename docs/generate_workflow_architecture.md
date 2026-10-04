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
| DAG stages | `rules/generate/{prepare,calculate,finalize}.smk` | Declare the checkpoint, partition/bucket jobs, bundle rendering, and certification. |
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
| `exclusion` | Binding required; data may be null | Intervals to subtract, independent of source. Combined with other interval exclusions. |
| `testing_unit_bed` | Yes | Unit intervals and stable unit IDs. |
| `fasta` and `.fai` | Yes | Contig lengths, partitioning, and reference allele enumeration. |
| `par_regions` | Yes | Pseudoautosomal intervals excluded before X-role assignment. |
| `exclusions` | Optional | Additional intervals excluded from calibration. |
| `mutation_rate` and optional index | Conditional | Required for an unobserved `mutation_rate` policy. |
| `af_source` | Yes | Names the overall, and when configured female/male, AF INFO tags. |

The current scope supports SNVs on contigs `1`–`22` and `X`. Every calibration
allele must have exactly one score row. Invalid selected AF values, duplicate
population states, duplicate scores, and missing mutation rates when they are
needed are all hard failures.

### Full-mode binding behavior

Full mode deliberately does not define a second reference scope under
`generation`. It reuses FASTA, exclusion, testing units, CADD SNVs, and ordinary
variant exclusions from `bindings.reference`. Exactly one configured exclusion
must have the `par_regions` role; every `variants_exclude` binding becomes a
generation exclusion. `WorkflowPlan` enforces this while the DAG is being
built.

## DAG and transformations

The scope checkpoint discovers the partition and bucket wildcards. Every later
rule reads those JSON manifests instead of trying to rediscover work from the
filesystem.

| Rule | Inputs | Transformation | Outputs |
| --- | --- | --- | --- |
| `prepare_generation_scope` (checkpoint) | FASTA index, units, exclusion/PAR/other exclusions, population/score bindings | Build supported non-PAR scope; subtract all other interval exclusions; measure retained fractions; partition the genome; assign eligible-unit buckets; capture input checksums. | Partition and bucket manifests; provisional unit/interval/mask Parquet files; resource-identity JSON. |
| `prepare_generation_partition` | One partition plus scope artifacts, FASTA, population VCF, and scores | Enumerate calibrated SNV alleles in the masked unit intervals; join unit membership, CADD score, and population allele state. | One each of `allele_states`, `allele_membership`, and `allele_scores` Parquet shards. |
| `calculate_partition_null_distribution` | One partition shard set and one unit bucket | Assign the absent-allele policy and role-specific effective AF; bin scores; aggregate additive no-event statistics. | One partition-by-bucket statistics Parquet shard. |
| `reduce_null_distribution_bucket` | All statistics shards for one bucket | Compensated reduction across partitions; calculate conditional score-bin survival rows and diagnostics. | HET/HOM and optional HEM lookup rows plus diagnostics TSV for that bucket. |
| `finalize_null_distribution_bundle` | All bucket rows/diagnostics plus scope and state artifacts | Retain only units with positive mass for every applicable role; stream public bundle components and metadata. | Gzipped lookup, scope, exclusion, and metadata files. |
| `validate_null_distribution_bundle` | Finalized bundle and schema | Validate component checksums, row coverage, contracts, and semantics. | `validation.json` and `validation.ok`. |

The default `null_distributions` target exposes the certificate and the
recoverable generation intermediates. `generate_results` is the equivalent
explicit target. In full mode, `full_results` exposes those generation outputs
and the downstream apply reports together.

## Intermediate-data contracts

These checkpoint artifacts are part of the internal DAG contract:

- `provisional_units.parquet` has one row per unit with `supported_bases`,
  `retained_bases`, `retained_fraction`, chromosome, and eligibility status.
  Supported bases exclude PAR but not other interval exclusions; retained bases
  also exclude the union of all other exclusion sources. The default minimum
  retained fraction is 0.5. Statuses distinguish `no_calibration_scope`,
  `fully_excluded`, `retained_fraction_below_minimum`, and `eligible`.
- `unit_intervals.parquet` preserves merged raw unit intervals; the masked
  counterpart preserves only callable, supported, non-excluded segments.
- `effective_calibration_mask.parquet` records the exact calibration intervals.
- `partitions.json` and `unit_buckets.json` are versioned, sorted manifests;
  they are the only source of partition and bucket wildcards.
- Partition `allele_states` holds `present`, `absent`, `rejected`, or invalid
  population state information. `allele_membership` maps alleles to units, and
  `allele_scores` supplies the matching score.
- Statistics shards are additive per `testing_unit_id`, role, and score bin.
  They store `certain_count`, `log_none_sum`, and event count so global event
  tables never need to be materialized.

Downstream rules should not glob these files or infer wildcard names from
directory contents. Use the callbacks in `rules/generate/context.smk`; they
wait for the checkpoint and validate the manifests before returning paths.

## Bundle contract

The public bundle lives at
`<outdir>/generate/<generation_prefix>/bundle/`. Its filenames are defined in
one place, `bundle_contract.py`. The public components are:

- HET and HOM lookup tables; HEM chromosome-X lookup table when enabled;
- eligible testing-unit BED and calibration-region BED;
- excluded alleles and excluded units;
- `metadata.json`, `validation.json`, and `validation.ok`.

`metadata.json` records contract and schema versions, assembly, score-grid
parameters, roles, absent-allele and X policies, source and component
checksums, and counts. Downstream apply rules should use `validation.ok` as the
completion boundary.

## Rules worth protecting when changing generation

1. **One binding set per full run.** In full mode, do not add a second
   generation-only reference source. Change the plan resolver if a new shared
   resource is needed.
2. **Determinism.** Unit IDs are sorted before bucket assignment; partitions
   and reductions are manifest ordered. Preserve that ordering when adding
   parallel work.
3. **Bounded memory.** Partition jobs only read one genomic window; bucket
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
| Add a generation input | `workflow_plan.py`, config schema, `context.smk`, scope/partition adapter and library | Include it in resource identities and bundle metadata; add fixture coverage. |
| Add a role | Generation settings/model plan, `generation_probability.py`, context roles, lookup/finalization validation | Define compatible chromosomes/genotypes, table filename, certificate semantics, and apply support. |
| Change binning | `generation_probability.py` and contract constants | Version the contract if compatibility changes; update apply lookup validation. |
| Change a public artifact | `bundle_contract.py`, finalization, validation schema/library | Preserve backward compatibility or deliberately version the contract. |

After a generation change, run the small fixture:

```bash
python -m unittest discover -s tests -p 'test_generate_workflow.py'
```
