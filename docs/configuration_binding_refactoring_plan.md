# Configuration schemas and exclusion bindings: implementation plan

## Review

The suggestion is valid, with qualifications about schema composition,
validation entry points, and configuration merging. Implement it in separate
stages; the schema cleanup can preserve accepted manifests, while keyed
exclusions are an intentional configuration-interface change.

Verified in the current tree:

- `workflow/schemas/run.schema.yaml` repeats `af_source` and generation
  `execution` three times. Apply/full reference definitions repeat the same
  properties except for `lookup_source`; standalone/full generation repeat
  their shared population, policy, mutation-rate, and execution properties.
- Top-level `generation` repeats shape checks supplied by the mode branches.
  Apply forbids it, and generate/full select different closed definitions.
- `_generation_plan()` checks exclusion lengths and zips metadata with paths
  for full mode. `build_apply_resource_registry()` independently checks and
  zips the effective lists. External apply can therefore build a plan whose
  exclusion mismatch is discovered only when the registry is constructed.
- Equal lengths cannot detect reordered bindings. IDs are read into
  `ExclusionMetadata` but are not used to join physical resources.
- `_effective_bundle_reference()` also produces parallel metadata/binding
  tuples, and registry names are positional (`exclusion_000`, etc.).
- Generation provenance uses another positional name,
  `interval_exclusion_000`, in `generation_scope.py`.

Three corrections to the proposal:

1. Keep shared schema bases closed only after they enumerate every
   mode-specific property. The installed JSON Schema validator does not safely
   track `unevaluatedProperties` across local YAML references, so concrete
   mode branches use `required` and `false` property constraints to select or
   forbid those enumerated fields. Verify this through both the JSON Schema
   validator and Snakemake's validator.
2. “Schema owns shape; plan owns every cross-field check” is too broad.
   The schema already expresses useful conditionals such as required HEM
   lookups, sex-specific AF tags, and mutation-rate bindings. Give each
   invariant one authority; keep these declarative conditions in the schema
   and reserve resolved resource relationships for the plan.
3. This does not require changing every run manifest. The breaking migration
   affects apply/full reference exclusions and science metadata. Standalone
   `generation.exclusions` has no parallel identity list today and can remain
   a list in this work.

Baseline: `python tests/run_tier.py schema` passes (28 tests). The repository
already contains other uncommitted refactoring; preserve those changes.

## Stage 1 — Factor schemas without changing the manifest contract

Add reusable `$defs` for `af_source`, `generation_execution`, and the shared
reference properties (including CADD bindings). Add a shared generation base
for population VCF, AF source, execution, absence policy, and optional
mutation-rate resource.

Compose concrete definitions as follows:

- Apply reference: reference base plus required `lookup_source`.
- Full reference: reference base, with `lookup_source` still forbidden.
- Standalone generation: generation base plus its required scope, score,
  FASTA, PAR, and optional extra-exclusion bindings.
- Full generation: generation base only; duplicate shared resources stay
  forbidden because full mode derives them from `bindings.reference`.

Keep the common bases closed after listing the full union of concrete fields.
Use mode definitions to require their applicable fields and set forbidden
shared fields to `false`; keep nested binding objects closed as they are now.
Replace top-level `generation` with a minimal declaration, leaving its complete
constraints to the selected mode definition. Preserve required fields,
nullability, bounds, enums, defaults, and mode conditionals exactly.

Verify acceptance and rejection for all current fixtures, plus misspelled
properties, full-mode `lookup_source`, duplicate shared generation resources,
missing `af_source.overall`, and malformed execution settings. Exercise the
merged `config.schema.yaml`, not just `run.schema.yaml`. Verify that default
insertion by Snakemake still behaves the same when `execution` is absent or
partially supplied. Do not change defaults during deduplication.

Acceptance: schema tests, all-mode lint and fresh DAG tests pass without
changing job counts or manifests. Treat line savings as incidental.

## Stage 2 — Establish one configuration-validation boundary

Recommended ownership:

| Authority | Responsibility |
| --- | --- |
| Schema | Object shape, required properties, primitive types, allowed values, mode restrictions, and declarative conditional requirements |
| Workflow plan | Binding correspondence, role selection, effective bundle substitutions, production index policy, and immutable resolved configuration |
| Runtime validators | Actual resource content, contigs, allele/sample compatibility, file integrity, and bundle certificates |
| Resource registry | Construct staging paths/specifications from an already valid plan |

Introduce a shared configuration-validation entry point in the Python library.
Have `build_workflow_plan()` normalize a copy and invoke it before resolving
records. Remove the separate schema-validation invocation from `Snakefile`
once this route is in place. Load the composed schema and its local references
from a repository-relative location, with no network resolution or file-content
validation. Preserve caller input immutability and current default behavior.

This protects direct callers in tests and other internal modules as well as
Snakemake. Wrap schema failures in `WorkflowPlanError`, including the useful
configuration path and underlying cause. Keep malformed boolean rejection at
plan construction, as required by the existing interface; never replace strict
booleans with truthiness conversion.

Inventory duplicated checks before deleting them. Remove schema-covered
presence/type checks only after every plan entry path reaches the shared
validator. Retain or separate transformations that are not validation:
`_required()` currently trims text, rejects sentinel strings, and expands user
paths, even for some non-path values. Removing it wholesale would change
behavior. Characterize those transformations and keep unrelated path/sentinel
cleanup out of this stage.

While the old list contract still exists, validate exclusion correspondence
once in plan construction, before generation selection or bundle substitution.
Remove the registry's independent length check. Report one consistent error
for apply/full mismatches; checking configured exclusions even for bundle
apply is an intentional earlier rejection and needs a regression test.

Acceptance: direct plan callers and Snakemake reject malformed configuration
before jobs run. Existing science conditionals still fail, with tests assigned
to their authoritative layer. Run schema, unit, and workflow tiers. The unit
tier must remain free of Snakemake subprocesses.

## Stage 3 — Replace positional exclusions with ID-keyed bindings

Use the mapping key as the ID on both sides; do not repeat an `id` property
inside each value:

```yaml
# Scientific configuration
reference:
  exclusions:
    gnomad-r2-notpass:
      kind: bed
      role: variants_exclude
    par-regions:
      kind: bed
      role: par_regions

# Run manifest
bindings:
  reference:
    exclusions:
      par-regions: {data: /reference/par.bed, index: null}
      gnomad-r2-notpass: {data: /reference/notpass.bed, index: null}
```

Require safe, nonempty IDs, for example
`^[A-Za-z0-9][A-Za-z0-9._-]*$`, and validate IDs as identifiers rather than
paths. Reordering either mapping must leave resource identity and role
unchanged. Report missing and unknown binding IDs separately and in sorted
order; equal counts are no longer a meaningful check.

Account explicitly for recursive configuration merging: a mapping override
merges keys, whereas the current list override replaces the list. An empty
mapping does not erase inherited exclusions. Recommended deletion mechanism:
allow scientific entries to specify `enabled: false`, defaulting to enabled
when omitted. Disabled entries need no physical binding and must never be
staged or consumed; a retained binding for a disabled entry may be ignored so
layered manifests can disable an inherited resource. Reject bindings for IDs
absent from scientific metadata. Test these cases through the actual
Snakemake merge path and document them in the migration instructions.

Build immutable resolved exclusion records containing identifier, kind, role,
and `FileBinding`. Use one record collection in effective references, rather
than preserving `effective_exclusion_metadata` alongside a second binding
collection. Select full-mode PAR and ordinary exclusions from these records;
retain the requirement for exactly one active PAR resource in full mode.

Use stable registry names such as `exclusion_<id>`, with collision checks for
fixed/internal names. Define a stable internal ID for the bundle-derived
excluded-alleles resource and resolve it through the same record type. Bundle
apply must continue using certified bundle scope and exclusions.

Use sorted IDs for deterministic traversal, logs, and staging paths. Confirm
that filtering still takes the union of exclusions. Replace the science-config
comment claiming that source order is scientific: identity must not depend on
ordering. Keep standalone generation's list handling unchanged in this stage.

Migrate science defaults, apply/full examples, fixture manifests, test-created
configurations, and documented examples together. Reject the old list shape
with a migration pointer. Do not silently infer IDs from list position.

Acceptance: independently permuting metadata and bindings produces the same
resolved resources and scientific outputs. Test missing/unknown IDs, disabled
entries, empty active sets, unsafe IDs, role selection, bundle substitution,
and immutable records. Run unit, schema, workflow, and E2E tiers. ID-based
staging paths may change; functional job counts should remain stable.

## Stage 4 — Carry IDs through generation provenance and finish migration docs

Update `generation_context.py`, scope binding dataclasses, the scope adapter,
and `generation_scope.py` so full-mode exclusion identities reach
`resource_identities.json` without being reconstructed by position. Pair ID
and path in one record through adapters; flat Snakemake input lists may still
declare file dependencies but must not become the authoritative association.

For standalone generation's existing list, retain deterministic internal
identifiers without pretending that these are science-config IDs. Document
that narrower exception. A future standalone keyed-input migration can be a
separate interface change.

Review bundle metadata consumers before replacing positional provenance keys.
Preserve reading existing certified bundles; version metadata only if its
contract requires fixed keys. Public lookup and scope filenames should remain
unchanged. Update provenance assertions in `test_full_workflow.py`, registry
tests, and any snapshots deliberately. Expect provenance and staging names to
change; compare scientific tables separately from metadata/checksums.

Add a migration section to `docs/pipeline.md` beside coverage-to-exclusion:

- Show before/after examples using the existing metadata IDs.
- Explain exact key matching, optional indexes, and the distinction between
  the singular `exclusion` resource and the plural exclusion catalog.
- Explain recursive merging and how to disable inherited resources.
- State that PAR roles and the 50% minimum retained fraction are unchanged.
- Note changed prepared filenames and the possible need to regenerate staged
  resources; do not suggest renaming old outputs to bypass provenance checks.
- State which standalone-generation manifests remain valid unchanged.

Check duplicate YAML-key handling before promising uniqueness diagnostics:
JSON Schema only sees parsed mappings, so duplicate source keys may already
have been overwritten by the loader. Reject duplicates at a supported loading
boundary if practical; otherwise document the loader limitation explicitly.

Acceptance: run workflow, E2E, and resilience tiers, including generation
recovery and full-mode bundle certification. Reordering keyed exclusions must
leave normalized provenance associations and scientific outputs unchanged.
Update baselines only for reviewed naming/provenance changes, then run
`git diff --check`.

## Suggested commits

1. Extract schema definitions and characterize composed-schema behavior.
2. Centralize configuration validation and exclusion correspondence checks.
3. Introduce keyed exclusions, resolved records, fixture migration, and the
   user-facing migration note together.
4. Preserve exclusion IDs in generation provenance and complete integration
   coverage/documentation.

Ship the migration note with Stage 3; expand its provenance details in Stage 4.
Keep algorithm changes, new exclusion formats, and unrelated configuration
normalization out of these commits.
