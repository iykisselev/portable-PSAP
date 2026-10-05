# Indel null-generation implementation plan

## Objective

Extend generated PSAP null distributions from the current SNV-only allele
universe to normalized small insertions and deletions while preserving:

- exact agreement between generation and apply-time allele identity;
- bounded, partition-local processing;
- sparse partition/bucket scheduling;
- explicit population-frequency and mutation-rate policies;
- deterministic, restartable artifacts; and
- certified bundle metadata that states the actual variant universe.

The first supported combined universe will be `variant_classes: [snv, indel]`.
Equal-length substitutions longer than one base, symbolic alleles, breakends,
star alleles, and structural variants remain out of scope.

## Central design decision: use a finite indel catalog

Do not attempt to enumerate indels directly from the reference sequence. Unlike
SNVs, there is no finite three-alternates-per-base construction: sequence and
length choices make the universe effectively unbounded.

Use the configured indel score resource as the authoritative candidate catalog.
An indel is eligible for null generation only when it:

1. appears as a normalized allele in the configured CADD indel resource;
2. satisfies the configured small-indel length policy;
3. has a reference allele matching the configured FASTA;
4. has a complete reference span inside the effective calibration mask; and
5. overlaps at least one eligible testing unit under the same membership rule
   used during apply.

This defines a finite, scoreable universe and prevents the generated null from
including alleles that the apply workflow cannot score. Record the catalog
identity and length policy in bundle metadata; they are scientific assumptions,
not merely execution settings.

## Required scientific contracts

### Allele identity and normalization

All generation, population, score, mutation-rate, exclusion, and apply inputs
must use the same biallelic VCF representation:

- uppercase `REF` and `ALT`;
- one-based `POS`;
- minimal alleles;
- left alignment against the configured FASTA; and
- one alternate allele per logical row.

Classify an allele after normalization:

| Class | Definition |
| --- | --- |
| SNV | `len(REF) == 1`, `len(ALT) == 1`, both canonical bases, and `REF != ALT` |
| indel | `len(REF) != len(ALT)`, both nonempty canonical DNA strings |
| unsupported | symbolic/breakend/star alleles, equal-length multi-base substitutions, or invalid bases |

Add shared helpers to `genomic_contracts.py` for classification, normalized
identity validation, and reference span. Do not independently reimplement
these rules in readers.

For an allele `(POS, REF, ALT)`, define its zero-based reference span as:

```text
[POS - 1, POS - 1 + len(REF))
```

Insertions therefore use their VCF anchor-base span. Deletions cover every
deleted reference base plus the anchor represented by `REF`.

### Mask and testing-unit membership

Use two distinct span rules:

- **Mask rule:** the complete reference span must be contained in one effective
  calibration-mask interval. An indel crossing an exclusion, PAR boundary, or
  unsupported region is rejected and audited.
- **Membership rule:** the reference span must overlap the testing-unit BED
  interval. One indel may belong to multiple testing units, producing one
  allele/unit membership for each.

Add conformance tests against apply-time `bcftools annotate` behavior. If its
interval matching does not exactly implement this contract, change apply-time
unit assignment to use the shared explicit span rule rather than weakening
generation semantics.

### Partition ownership and boundary handling

Assign every indel to exactly one partition using its normalized anchor
position (`POS`). This prevents duplicate score, state, and probability rows.

A deletion can overlap a unit while its anchor lies before that unit or in the
preceding partition. Introduce an indel candidate halo based on the maximum
supported `REF` length:

```text
left halo = max_reference_length - 1 bases
```

Use the halo only for partition discovery and sparse partition/bucket
incidence. The actual membership and mask decisions still use the exact allele
span. Extend partition preparation so the owner partition can associate a
boundary-crossing deletion with units in the adjacent partition.

The following invariant must hold:

> Repartitioning the same inputs must not change the allele set, unit
> memberships, event counts, or final lookup values.

### Population states and event probabilities

Reuse exact normalized `CHROM/POS/REF/ALT` joins and the existing population
states:

- `observed`: use configured overall or sex-specific population AF;
- `absent`: apply the configured unobserved policy;
- `rejected`: exclude before probability aggregation; and
- invalid AF/cardinality/value: fail the partition.

An omitted configured population-specific AF field continues to mean AF zero
and must remain present in `population_af_zero_imputations.tsv.gz`.

The existing role probabilities apply to a biallelic indel allele:

```text
HET: 2p(1-p)
HOM: p^2
HEM: p
```

The current generated-bundle model assumes allele events are independent.
Retain that approximation in the first indel implementation and state it
prominently in metadata. Add QC for multiallelic loci (`sum(AF) <= 1`) and
measure how often multiple retained indel alleles share a locus. A future
locus-level multinomial model would require a separate scientific review and
bundle-contract version.

### Unobserved indels

Never apply an SNV mutation rate to an indel.

- Under `unobserved_policy: zero`, a catalog indel absent from the population
  panel receives zero event probability, as SNVs do now.
- Under `unobserved_policy: mutation_rate`, require a separate normalized,
  allele-specific indel mutation-rate resource containing
  `CHROM, POS, REF, ALT, mu` for every absent retained indel.

For backward compatibility, retain the existing `generation.mutation_rate` as
the SNV rate binding and add `generation.mutation_rate_indel`. Reject an
indel-enabled mutation-rate run unless both required class-specific rate
resources are present. Record their identities separately.

## Configuration and plan changes

### User-facing configuration

Support this shape in standalone generation:

```yaml
model:
  include_indels: true

generation_settings:
  variant_classes: [snv, indel]
  indels:
    candidate_universe: score_resource
    membership: reference_span_overlap
    mask_policy: reference_span_contained
    max_reference_length: 50
    max_alternate_length: 50

generation:
  score_snv: {data: /path/to/whole_genome_SNVs.tsv.gz, index: /path/to/file.tbi}
  score_indel: {data: /path/to/InDels.tsv.gz, index: /path/to/file.tbi}
  mutation_rate: {data: /path/to/snv_rates.tsv.gz, index: /path/to/file.tbi}
  mutation_rate_indel: {data: /path/to/indel_rates.tsv.gz, index: /path/to/file.tbi}
```

`mutation_rate*` remains conditional on `unobserved_policy: mutation_rate`.
In full mode, resolve `score_indel` from `bindings.reference.cadd.indel`, just
as generation already reuses the SNV CADD binding.

Enforce these relationships in both JSON Schema and `WorkflowPlan`:

- `model.include_indels: false` requires `[snv]`;
- `model.include_indels: true` requires `[snv, indel]`;
- indel generation requires the indel score/catalog binding;
- indel plus mutation fallback requires the indel mutation-rate binding; and
- length limits must be positive and emitted in metadata.

Remove the current generate/full schema clauses that force
`model.include_indels: false` only after all downstream contracts are ready.

### Internal bindings

Extend the immutable plan and adapter dataclasses with:

- `score_indel` and index;
- `mutation_rate_indel` and index;
- enabled variant classes; and
- the resolved indel policy/length limits.

Thread these through `GenerationRuleContext`, preparation rules, script
adapters, resource inventories, and execution-policy sizing. Optional Snakemake
inputs must remain empty lists when indels are disabled.

## Data-flow implementation

### Phase 1: shared allele and resource validation

1. Add `variant_class()`, `reference_span()`, and canonical allele validation
   to `workflow/lib/python/genomic_contracts.py`.
2. Add a partition-local indel score reader in `generation_scores.py`. It must:
   - query only the owner partition;
   - reject unsupported/noncanonical alleles;
   - enforce length limits;
   - validate finite CADD-PHRED scores;
   - preserve genuine duplicate rows for explicit duplicate detection; and
   - verify `REF` against FASTA before accepting the candidate.
3. Extend resource preflight validation to inspect indel class composition,
   normalization, contig naming, numeric score values, and index compatibility.
4. Fail early when the population, CADD, mutation-rate, and FASTA assemblies or
   normalization conventions cannot produce exact allele joins.

If the raw CADD indel index is point-based, keep anchor ownership but perform
post-read span filtering. If production measurements show boundary queries are
insufficient, stage a derived BGZF table with explicit zero-based `START` and
`END` columns and a span-aware Tabix index; do not whole-file scan per
partition.

### Phase 2: scope and sparse topology

1. Build ordinary callable unit intervals as today.
2. When indels are enabled, derive candidate-owner intervals by expanding each
   eligible unit to the left by `max_reference_length - 1`, clipping to the
   chromosome and supported scope.
3. Build generation partitions from the union of SNV callable intervals and
   indel candidate-owner intervals.
4. Build sparse bucket/partition incidence from those class-aware owner
   intervals so a boundary deletion schedules its owner partition and bucket.
5. Version the partition/scope manifest if it gains variant-class or halo
   fields. Validate that every published sparse pair refers to a known
   partition and bucket.

Do not revert to a partition-by-bucket Cartesian product.

### Phase 3: partition preparation

Refactor `prepare_generation_partition()` into explicit candidate producers:

```text
SNV candidates   = FASTA enumeration in callable unit bases
indel candidates = normalized score/catalog rows owned by this partition
all candidates   = validated union keyed by CHROM/POS/REF/ALT
```

For indels:

1. Compute and retain `SPAN_START` and `SPAN_END` during validation.
2. Require complete span containment in the effective mask.
3. Interval-join the span to eligible testing-unit intervals.
4. Deduplicate `(testing_unit_id, allele key)` memberships.
5. Join the candidate to its catalog score exactly once.
6. Join exact population states; classify unmatched candidates as absent.
7. Write rejected candidates with a stable reason code.

Continue publishing the existing three Parquet artifact families, but add an
explicit `variant_class` column:

- `allele_states`;
- `allele_membership`; and
- `allele_scores`.

Add `SPAN_START/SPAN_END` to membership or to a separate candidate artifact if
needed for audits. Keep probability joins on the canonical allele columns so
SNV behavior remains unchanged.

Suggested stable rejection reasons include:

- `unsupported_allele_class`;
- `not_minimal_or_left_aligned`;
- `reference_mismatch`;
- `length_limit_exceeded`;
- `span_not_contained_in_calibration_mask`;
- `no_eligible_testing_unit_overlap`;
- `non_pass_population_filter`; and
- `missing_or_duplicate_score`.

Structural candidate failures should be reported separately from population
FILTER rejection; missing scores for an accepted candidate remain fatal.

### Phase 4: frequency assignment and statistics

1. Make `mutation_rate_rows()` class-aware or provide separate SNV and indel
   readers with the same returned schema.
2. Join mutation rates by exact normalized allele key and `variant_class`.
3. Require complete rates only for absent alleles with mutation fallback.
4. Keep HET/HOM/HEM probability formulas and additive statistics unchanged.
5. Add per-class event counts to diagnostics so SNV-only, indel-only, and
   combined contributions can be audited without changing lookup semantics.
6. Verify that a testing unit split across partitions still reduces to one row
   and that an indel crossing a partition boundary contributes exactly once.

### Phase 5: bundle and apply contracts

Adding indels changes the scientific meaning of generated lookup values.
Introduce `generated_bundle_v3` and schema version 3 rather than silently
reinterpreting v2.

Bundle v3 metadata must include:

- `variant_classes: [snv, indel]` in deterministic order;
- indel candidate-universe identity;
- normalization and reference-span conventions;
- membership and mask policies;
- maximum `REF`/`ALT` lengths;
- SNV and indel score-resource identities;
- separate mutation-rate identities when applicable;
- retained/rejected/observed/absent counts by variant class; and
- the existing allele-independence declaration.

Update bundle finalization and validation to verify these fields and component
checksums. Extend `excluded_alleles.tsv.gz` (or add a certified candidate-audit
component) so filtered indel candidates and reasons are preserved.

Update apply-time model validation so:

- `model.include_indels: true` rejects a bundle lacking `indel` in
  `variant_classes`;
- `model.include_indels: false` rejects or explicitly ignores an indel-enabled
  bundle according to one documented policy (prefer rejection); and
- the configured CADD indel identity matches the identity certified by the
  generated bundle in full mode.

Keep v1/v2 readers for existing SNV bundles. Prefer emitting v3 for all new
generation runs once implemented, including SNV-only runs, so there is one
current writer and explicit compatibility readers.

### Phase 6: documentation and operational controls

Update:

- `README.md` input and model descriptions;
- `docs/generate_workflow_architecture.md` candidate and partition semantics;
- example generate/full manifests;
- resource preparation documentation; and
- bundle migration notes.

Add preparation summaries containing candidate counts by class and rejection
reason. Add resource estimates for the indel catalog because partition I/O,
not probability reduction, will likely dominate runtime.

## Test plan

### Unit tests

- SNV/indel/unsupported classification.
- Minimal representation and FASTA-reference validation.
- Insertion and deletion span calculation.
- Complete-mask containment, including a deletion crossing a mask gap.
- Testing-unit overlap for an anchor inside, before, and at a boundary.
- Deterministic owner partition for a boundary-crossing deletion.
- Sparse bucket incidence includes the owner partition but no unrelated pair.
- Exact population AF matching for normalized biallelic and multiallelic rows.
- Missing population-specific AF becomes zero and remains audited.
- Duplicate/missing/nonfinite indel scores fail.
- Length-limit and symbolic-allele rejection reasons.
- Separate SNV and indel mutation-rate coverage.
- Multiallelic AF-sum validation.

### Numerical and invariance tests

- Compare Python aggregation with a small brute-force implementation containing
  SNVs, insertions, and deletions.
- Confirm SNV-only outputs remain byte-identical or numerically identical under
  the selected bundle-version policy.
- Confirm combined results are invariant to partition size and unit bucket
  size.
- Confirm a deletion spanning two partitions is counted once.
- Confirm a variant belonging to two genes contributes once to each gene.
- Confirm reduction is deterministic across restart and shard order.

### Workflow fixtures

Create a compact normalized fixture containing:

- one observed insertion;
- one observed deletion;
- one population-absent catalog indel;
- one multiallelic population record;
- one deletion crossing a partition boundary;
- one indel crossing an excluded interval;
- one indel overlapping two testing units; and
- chromosome-X indels for overall-proxy and sex-specific AF policies.

Run generate, apply, and full modes with both `zero` and `mutation_rate`
policies. Validate certificate corruption and recovery behavior for the new
resources and artifacts.

### Performance acceptance

On a representative whole-genome CADD indel resource:

- no partition performs a whole-file scan;
- peak memory remains within the partition resource policy;
- scheduled statistics jobs remain the sparse bucket/partition incidence;
- candidate counts and wall time scale with queried catalog rows;
- re-running a missing partition or bucket recomputes only its dependencies;
  and
- finalization streams audit/state artifacts rather than collecting them
  genome-wide.

## Recommended delivery sequence

1. **Contracts first:** allele helpers, schemas, plan fields, bindings, and
   early validation while generation still rejects enabled indels.
2. **Candidate preparation:** indel score reader, normalization/reference/span
   checks, halo-aware sparse topology, and partition artifacts.
3. **Population and rates:** exact state joins, zero policy, then explicit
   allele-specific indel mutation fallback.
4. **Aggregation:** combined per-class events, diagnostics, and partition
   invariance tests.
5. **Bundle v3:** metadata, audit output, semantic validation, and backward
   readers.
6. **Apply/full integration:** enforce bundle/model class compatibility and run
   combined end-to-end fixtures.
7. **Production benchmark:** tune partition memory/runtime only after the
   scientific and restart contracts pass.

Do not enable `model.include_indels` in generate or full mode until phases 1–6
are complete. Partial enablement could produce a bundle that scores cohort
indels against an SNV-only null distribution, which must remain a hard error.

## Definition of done

Indel null generation is complete when:

- one normalized candidate definition is used from resource validation through
  apply;
- every retained indel has exactly one score, one canonical population state,
  at least one unit membership, and at most one owner partition;
- every absent indel has the policy-required mutation assignment;
- boundary-crossing alleles are neither lost nor duplicated;
- bundle metadata completely identifies the indel universe and policies;
- apply rejects incompatible model/bundle combinations before cohort work;
- SNV-only compatibility tests pass; and
- unit, schema, workflow, end-to-end, resilience, numerical, and production
  performance tests all pass.
