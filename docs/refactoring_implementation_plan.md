# Refactoring implementation plan

Date: 2026-10-12  
Status: Complete; eight iterations implemented and verified

This plan compares proposals L1–L6 and D2–D10 with the current working tree.
Its purpose is to improve readability and consolidate contracts without changing
PSAP calculations, published bundle formats, or resolved output paths.

## Scope and constraints

- Preserve existing uncommitted work, particularly annotation parallelization,
  validation-rule changes, and the accessory annotation script and tests.
- Keep scientific behavior unchanged except for explicitly identified boundary
  fixes and separately reviewed validation changes.
- Preserve generated-bundle filenames, certificate component names, role IDs,
  and compatibility with existing bundles.
- Preserve checkpoint-driven dependency discovery and sparse partition/bucket
  scheduling.
- Do not restore resource hashing.
- Do not retune memory allocations as part of readability changes.
- Do not change the PSAP probability model or scoring method.
- Defer a wholesale validation API migration to a new `ValidationLog` class.

## Findings requiring early attention

### Scope exclusion binding mismatch

At the initial review, `workflow/rules/generate/prepare.smk` declared `input.exclusions` as a list of
paths. `GenerationScopeInputs.exclusions` expected a mapping, and
`workflow/scripts/generate/prepare_generation_scope.py` binds the input directly.
A separate mapping in params does not correct that binding.

A direct adapter-boundary check with a nonempty path list reproduced a
`ValueError`. Empty lists conceal the mismatch because they can be converted to
an empty dictionary.

Scope computation only iterates over the mapping's values. Represent its
exclusions as an ordered tuple of paths and retain scientific IDs in the
workflow plan and provenance mapping.

### Tests omitted from explicit CI tiers

At the initial review, the following files existed but were not assigned to the explicit tiers in
`tests/run_tier.py`:

- `test_fixed_width_bin_contract.py`
- `test_generation_adapter_bindings.py`
- `test_generation_event_assembly.py`
- `test_generation_reference.py`
- `test_generation_statistics.py`
- `test_generation_vectorized_reduction.py`
- `test_generation_window_records.py`
- `test_indel_bundle_contract.py`
- `test_missing_token_contract.py`
- `test_text_and_fai_contracts.py`

These files are now assigned to the unit tier. A new coverage check requires
every Python test file to be assigned exactly once.

## Proposal assessment

| Proposal | Decision | Implementation guidance |
| --- | --- | --- |
| L1: Split plan builders | Proceed | Extract cohesive builders; normalize once; use a named resource-resolution record rather than a positional tuple; consolidate frozen-mapping helpers. |
| L2: Split scope preparation | Proceed with revisions | Remove the shadowed chromosome and positional zip. The unit-interval frame is already reused, so preserve that implementation. |
| L3: Name execution-policy arguments | Proceed without allocation changes | Use existing dataclass field names such as `baseline_mem_mb` and `baseline_runtime_minutes`. Document units; reconcile allocation values separately. |
| L4: One population-state row builder | Already addressed | `population_rows()` already uses a shared columnar `append_rows()` helper. No additional abstraction is needed. |
| L5: Canonical HEM lookup naming | Proceed atomically | Use `lookup_hem_chrX` for internal artifact bindings; retain `hem_chrX` as the role ID and preserve published filenames and certificate keys. |
| L6: Formatters and linters | Proceed in stages | Snakemake lint already runs in CI. Add Python/Snakemake formatting and selected Python/R lint checks, with mechanical formatting isolated from semantic fixes. |
| D2: Raw generation input map | Adapt | Hashing and its adapter no longer exist. Consolidate bindings used by current jobs, with explicit subsets rather than universal dependencies. |
| D3: Optional-input helper | Proceed with distinct conversion directions | Path-to-input-list conversion is not the same operation as `optional_path()` or scalar index lookup. Replace the misleading index-only `as_inputs()` API without conflating these contracts. |
| D4: Forwarding wrappers and aliases | Proceed | Remove static forwarding layers and redundant HEM guards. Rename the scope binding to `effective_mask`; retain checkpoint-dependent callbacks. |
| D5: Registry queries | Proceed with an explicit adapter handoff | Derive resource queries and exclusion specifications from records. IDs alone do not make the rule-local registry available in script processes. |
| D6: Shared bootstrap and binding | Proceed incrementally | Keep a minimal import-path entry point; set Polars thread limits before importing Polars. Put tool discovery in a library module accessible to both libraries and adapters. |
| D7: Shared genomic constants | Split by responsibility | Share contig parsing and named scope constants, preserving distinct validation/generation scopes. Keep model roles beside `GENERATED_ROLES` in `bundle_contract.py`. |
| D8: Policy values and defaults | Prioritize boundary defaults | Remove misleading scientific defaults from adapters. Schemas already validate policy values; retain direct-library safeguards. Type annotations alone do not enforce policy validity. |
| D9: Targets and sample paths | Proceed | Remove the fixed testing-unit wildcard while preserving resolved filenames, logs and benchmarks. Retain genuine wildcard and deferred-rule callbacks. |
| D10: Validation duplicates | Proceed selectively | Share PED/TSV orchestration while preserving format-specific errors. Explicitly test stricter BED behavior; defer the broad `ValidationLog` migration. |

## Implementation sequence

Each phase should be independently reviewable. Commit labels below describe
logical change sets; a phase may need multiple commits. Do not mix unrelated
working-tree changes into them.

### Phase 1 — Coverage and scope adapter boundary

Suggested commit: `fix: align scope exclusion bindings and cover omitted tests`

Tasks:

- [x] Assign all omitted test files to appropriate CI tiers.
- [x] Add a check that each `test_*.py` file has a tier or an explicit documented
  exclusion; register the new coverage test itself.
- [x] Add a regression for a nonempty scope exclusion list.
- [x] Change `GenerationScopeInputs.exclusions` to an ordered tuple of paths.
- [x] Update scope iteration and direct test construction accordingly.
- [x] Remove the unused scope params exclusion mapping.
- [x] Preserve ID-to-path associations in the plan and finalization provenance.

Acceptance criteria:

- Standalone and full-mode scope preparation support multiple exclusions.
- All exclusion files remain declared dependencies.
- Provenance identities remain unchanged.
- Omitted test coverage becomes a CI failure rather than a silent gap.

### Phase 2 — Explicit scientific parameters

Suggested commit: `fix: require scientific settings at generation boundaries`

Tasks:

- [x] Inventory scientific fields already supplied by generation rule params.
- [x] Remove adapter defaults that could silently change scientific behavior,
  particularly indel length limits and `include_indels`.
- [x] Retain legitimate optional-file defaults and avoid unrelated API changes.
- [x] Update dataclass field ordering as needed when defaults are removed.
- [x] Replace tests asserting a default indel length of `1` with explicit-value
  and missing-field tests.
- [x] Require actionable `BindingError` messages for missing rule parameters.
- [x] Keep schema checks and safeguards in directly callable library functions.

Notes:

The plan builder uses indel length defaults of `50`, while several adapter
dataclasses default to `1`. Removing adapter fallbacks prevents a missing rule
parameter from silently changing the candidate universe.

`Literal` annotations may document policy types, but the current binder does
not enforce them. If runtime coercion for literals or enums is introduced,
implement and test it explicitly. Do not remove validation on the assumption
that annotations enforce it.

Acceptance criteria:

- Missing scientific settings fail at the adapter boundary.
- Explicitly configured SNV and indel workflows produce unchanged results.
- Direct library callers still receive clear invalid-policy errors.

### Phase 3 — Plan builders and readable execution policies

Suggested commit: `refactor: separate workflow plan construction stages`

Tasks:

- [x] Let the Snakefile call `build_workflow_plan(config)` directly.
- [x] Keep normalization and its legacy warning inside the plan builder.
- [x] Extract `_model_plan`, `_generation_settings`, `_report_plan`, and
  `_apply_plan`.
- [x] Separate full-mode and standalone generation resource resolution.
- [x] Return a named record with resources, PAR regions, exclusions and explicit
  diagnostic keys; a single prefix cannot express all existing key differences.
- [x] Replace the duplicate frozen-mapping helpers with one generically typed
  helper.
- [x] Convert `ExecutionPolicy` construction to named arguments using existing
  field names.
- [x] Keep all memory/runtime values and retry behavior unchanged.

Acceptance criteria:

- All mode/schema and immutability tests pass.
- Configuration is not mutated and legacy normalization warns only once.
- Configured versus effective reference semantics remain intact.
- Resource requests and retry scaling match their current tests exactly.

### Phase 4 — Generation bindings and canonical names

Suggested commit: `refactor: consolidate resource bindings and lookup names`

Tasks:

- [x] Introduce a dependency-light helper for optional scalar paths to
  Snakemake input lists.
- [x] Keep path-list decoding, scalar index retrieval and path-list creation
  as distinct operations.
- [x] Derive scope, partition and statistics input dictionaries from canonical
  generation resource bindings, selecting only each job's actual dependencies.
- [x] Remove static forwarding wrappers and unnecessary wildcard lambdas.
- [x] Rename the scope output and corresponding typed fields to `effective_mask`.
- [x] Remove the scope output alias dictionary.
- [x] Centralize enabled role selection beside bundle role constants, with the
  corresponding R helper and cross-language tests.
- [x] Rename internal HEM bindings to `lookup_hem_chrX` across rules, contexts,
  dataclasses, Python consumers and R adapters in one coordinated change.
- [x] Remove redundant outer HEM guards while preserving disabled-role empty
  bindings where required by adapters.
- [x] Retain all checkpoint-dependent input functions.

Acceptance criteria:

- Existing bundles remain readable without migration.
- Public filenames, metadata roles and certificate keys are unchanged.
- HEM-enabled and HEM-disabled DAGs resolve correctly.
- No unrelated raw-resource dependencies are added to jobs.
- Checkpoint expansion and sparse partition/bucket edges are unchanged.

### Phase 5 — Scope orchestration

Suggested commit: `refactor: separate generation scope computation and publication`

Tasks:

- [x] Extract resource checks/mask construction, retention and interval
  preparation, manifest construction, and artifact publication.
- [x] Use named intermediate records where several results travel together.
- [x] Make indel-owner expansion preserve `(unit_id, chromosome, start, end)`
  records; update its callers and tests.
- [x] Remove the parallel-list zip and shadowed chromosome variable.
- [x] Continue reusing the existing unit-interval frame.
- [x] Preserve sorting, Parquet row-group sizes, compression and atomic writes.
- [x] Preserve unsupported-contig warning/filter behavior and sparse scheduling.

Acceptance criteria:

- Tests cover upstream anchors, partition boundaries and multiple units sharing
  coordinates without losing their identities.
- Empty or fully excluded scopes retain their defined behavior.
- Manifest order, bucket membership and scientific table contents are unchanged.
- Refactoring does not introduce additional genome-wide materialization.

### Phase 6 — Apply registry, adapters and sample paths

Suggested commits:

- `refactor: expose paired apply resource queries`
- `refactor: share Python adapter bootstrap and tool discovery`
- `refactor: simplify fixed-unit sample targets`

Tasks:

- [x] Add optional-resource and role-based registry queries for declared fields.
- [x] Build exclusion specifications directly from each resource record instead
  of joining separate sequences by position.
- [x] Remove duplicate indel guards where registration already encodes whether
  the resource is enabled.
- [x] Preserve the existing path to validation for missing enabled resources;
  do not turn it into a DAG-construction error.
- [x] Define a serializable resource association for adapters, with IDs mapped
  to declared input fields so paths can be resolved from `snakemake.input`.
- [x] Do not serialize the live registry or rebuild the full plan inside jobs.
- [x] Share bootstrap setup, retaining a minimal import-path entry point in
  adapters and setting `POLARS_MAX_THREADS` before Polars imports.
- [x] Consolidate `require_tool()` in a dependency-light library module.
- [x] Introduce typed apply bindings incrementally, one adapter at a time.
- [x] Define cohort terminal paths once and share per-sample path construction.
- [x] Replace the fixed `{testing_unit}` wildcard with the configured unit in
  output, log and benchmark patterns.
- [x] Remove constant-returning callbacks only where rule-definition ordering
  permits direct references; retain genuine sample-dependent callbacks.

Acceptance criteria:

- External lookup, generated-bundle apply and full-mode DAGs still construct.
- Missing optional resources never produce empty-string file dependencies.
- Custom indexes remain paired with their source files.
- Adapters work through Snakemake's script execution, not only direct imports.
- Resolved output, log and benchmark paths are unchanged.
- Existing annotation parallelization is preserved.

### Phase 7 — Targeted validation consolidation

Suggested commit: `refactor: share validation parsing and metadata orchestration`

Tasks:

- [x] Extract shared PED/TSV orchestration with explicit labels and error
  handling policies.
- [x] Keep meaningful `read_ped()` and `read_samples_tsv()` entry points.
- [x] Preserve row numbering, error ordering, malformed-header behavior,
  phenotype exclusion and sex warnings.
- [x] Centralize VCF contig-header parsing in `genomic_contracts.py`.
- [x] Define distinct named contig scopes: validation currently includes Y,
  while generation does not; canonical naming also supports MT.
- [x] Share ordering definitions without replacing vectorized frame operations
  with per-row Python callbacks or changing manifest order.
- [x] Use the canonical role helper for lookup validation.
- [x] Move testing-unit ID reading onto `iter_bed()` with an explicit diagnostic
  policy and malformed-coordinate regression tests.
- [x] Preserve validation report schemas and existing merge behavior.

Acceptance criteria:

- Existing sample-metadata characterization tests remain unchanged where
  behavior is intended to remain compatible.
- Contig normalization is consistent without broadening scientific scope.
- Stricter BED acceptance/rejection is documented and tested as a deliberate
  validation change, rather than hidden inside a readability refactor.

### Phase 8 — Formatting and CI enforcement

Suggested commits:

- `style: apply configured Python and Snakemake formatting`
- `ci: enforce formatting and selected lint checks`

Tasks:

- [x] Add pinned formatter/linter dependencies and repository configuration.
- [x] Add `ruff format --check`, `ruff check`, and `snakefmt --check`.
- [x] Enable unused-import/variable checks and review Bugbear findings before
  enforcing the selected rule set.
- [x] Configure narrow exceptions for injected Snakemake names and bootstrap
  import ordering; do not suppress findings globally.
- [x] Introduce selected R `lintr` checks with data.table conventions accounted
  for; repair long R expressions without changing semantics.
- [x] Keep mechanical formatting isolated from functional lint fixes.
- [x] Retain existing Snakemake lint tests and their reviewed allowlist.
- [x] Exclude generated outputs and caches from formatting targets.

Acceptance criteria:

- A clean checkout passes formatter and lint checks.
- Formatting produces no scientific output changes.
- Existing tests pass after formatting.
- Allowlist changes are reviewed rather than regenerated blindly.

## Verification strategy

Run targeted tests for each phase, then the full tier sequence in the pinned CI
environment before merging cross-cutting changes:

```bash
conda run --no-capture-output --name portable-psap-ci python tests/run_tier.py unit
conda run --no-capture-output --name portable-psap-ci python tests/run_tier.py schema
conda run --no-capture-output --name portable-psap-ci python tests/run_tier.py workflow
conda run --no-capture-output --name portable-psap-ci python tests/run_tier.py e2e
conda run --no-capture-output --name portable-psap-ci python tests/run_tier.py resilience
```

Required regression dimensions:

| Dimension | Cases |
| --- | --- |
| Workflow mode | Apply, generate, full |
| Lookup source | External legacy tables; generated bundles with existing certificates |
| Model settings | HEM on/off; SNV-only and SNV-plus-indel generation |
| Exclusions | Empty, one and multiple inputs; retained identity/path pairing |
| Resource indexes | Conventional, custom and optional missing indexes |
| Scope boundaries | Partition-spanning indels, upstream anchors, unsupported contigs, empty retained scope |
| Sample metadata | PED/TSV, malformed headers/rows, unknown phenotype, no controls |
| Execution | Fresh DAG, checkpoint-resolved DAG, rerun, partial-output recovery |

Use existing deterministic digest tests for public outputs. For Parquet
artifacts, verify schemas, ordered rows and required row-group properties rather
than imposing new cross-version byte-identity requirements. Any intentional
change to validation acceptance or diagnostics must have its own regression
test and review note.

## Review baseline

During the read-only review, 47 targeted tests passed locally:

- Adapter bindings: 4 tests.
- Execution policies: 5 tests.
- Workflow plan: 33 tests.
- Sample metadata characterization: 5 tests.

The nonempty scope-exclusion failure was reproduced separately at the adapter
boundary. The full workflow and E2E suites were not run for this review. Local
test success is not a substitute for the pinned CI environment verification
required above.


## Implementation record

The eight phases were implemented sequentially, with targeted verification
between structural changes. Formatting was applied in a dedicated tool pass.
The changes remain uncommitted; the suggested commit labels above remain a guide
for grouping the working-tree changes.

| Iteration | Implemented result |
| --- | --- |
| 1 | Scope exclusions bind as ordered paths; all Python tests have CI tiers and a coverage guard. |
| 2 | Scientific generation settings are required at adapter boundaries, including indel length limits. |
| 3 | Plan construction has separate model/settings/report/apply builders and mode-specific resource resolution; execution policies use keyword arguments. |
| 4 | Raw generation inputs are declared from shared bindings; effective-mask and HEM artifact names are consistent; enabled roles share Python/R definitions. |
| 5 | Scope preparation uses named mask/table records and separate computation/publication helpers; indel-owner intervals retain unit identities. |
| 6 | Apply resource associations refer to actual job inputs; Python adapters share bootstrap and typed bindings; executable discovery and fixed-unit sample paths are consolidated. |
| 7 | PED/TSV orchestration and contig parsing are shared; generation/validation scopes remain distinct; BED ID validation reports malformed coordinates. |
| 8 | Ruff, snakefmt and selected R linters are configured and enforced by a new CI job; formatter/linter versions are pinned. |

Additional findings resolved during verification:

- Empty scope frames require an explicit chromosome string type before applying
  the shared ordering map. A regression test covers empty bucket assignment.
- With HEM disabled, population parsing does not supply male AF values. Event
  assembly previously rejected those unused values. It now receives the enabled
  roles explicitly and validates only applicable frequencies; enabled-role
  probability calculations are unchanged. Unit coverage checks both disabled
  HEM acceptance and enabled HEM rejection of missing male frequencies.
- Conda's lintr 3.3 builds require R 4.4 or later. The CI environment pins lintr
  3.1.2, which resolves with the project's R 4.3.3 runtime. The complete updated
  CI environment was installed successfully in a temporary prefix.
- Conditional variant-filter rules now establish default paths before their
  enabled branches, allowing snakefmt to parse the module consistently.

Final verification passed in the updated CI environment (Python 3.11.8,
Snakemake 9.23.1 and R 4.3.3), using the existing declared job environments.
All 217 tests passed across the five tiers:

| Tier | Tests | Result |
| --- | ---: | --- |
| Unit | 167 | Passed |
| Schema | 34 | Passed |
| Workflow | 7 | Passed |
| End-to-end | 5 | Passed |
| Resilience | 4 | Passed |

The external-lookup report digests remained unchanged. Generated SNV/indel
bundles passed certification, full-mode workflows completed, HEM-disabled
generation omitted HEM work products, and recovery recreated identical bundle
digests after removing partition states, statistics shards and lookup buckets.

Ruff checks, Ruff format checks, snakefmt checks, the selected R lint checks and
`git diff --check` also passed. Existing Snakemake lint allowlists and workflow
baselines required no changes.
