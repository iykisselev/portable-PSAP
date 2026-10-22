# Partitioned null implementation progress

Branch: `feature/partitioned-nulls-acat`.
Specification: [implementation instructions](partitioned_nulls_acat_implementation_instructions.md).

## Iteration 1: immutable baseline and execution-window rename

Completed on 2026-10-20:

- Captured six certified full-mode reference runs from an isolated archive of
  unchanged `main` at `b3c967cab0ca35fe537ae8a7b070d0cf6639a002`.
- Saved per-role lookup rows, original per-sample scored variants, and 17 scalar
  sample/unit/model PSAP groups, with input/output hashes and environment
  provenance in [the reference directory](../tests/reference/partitioned_nulls/main_baseline/README.md).
  Coverage includes all five disease models, HEM disabled, a retained insertion,
  and a nondefault fixed-width score grid.
- Renamed execution-window modules, types, manifests, callbacks, rule names,
  resource policies, adapter scripts, internal paths, and tests. Apply's internal
  annotation-window parameter uses the same terminology.
- Bumped the window and bucket manifest versions to 2 and 3 respectively.
  Old execution-partition manifests are rejected and must be regenerated.
- Preserved the public `generation.execution.partition_size_bp` config spelling,
  window boundaries/order, unit bucketing, numeric calculations, and public outputs.

Verification in the existing working tree:

| Tier | Tests | Result |
| --- | ---: | --- |
| Unit | 189 | Passed |
| Schema | 42 | Passed |
| Workflow | 7 | Passed |
| End-to-end | 6 | Passed |
| Resilience | 4 | Passed |

The new reference workflow test executes six runs and checks all **100 public
component digests** against the captured `main` results. All match byte-for-byte.
Ruff checks/formatting, snakefmt, selected R lint, and whitespace checks pass.
The unit total includes the user's preexisting uncommitted CADD-analysis tests.

Tests used Python 3.11.8 and Snakemake 9.23.1 from the installed CI environment.
A temporary executable wrapper supplied `--workflow-profile none`, allowing the
already installed tools to run without the default profile attempting Conda
downloads. This wrapper is not a workflow change. R lint used the system R 4.3.3
installation with `lintr`; the installed CI R environment lacked that package.

The preexisting README, CADD-analysis files, test-tier addition, implementation
specification, and simplification plan were preserved outside this iteration's
commit. No generated old-format bundle reader has been removed yet: that removal
belongs to the complete N = 1 contract migration, with its consumers and tests.

## Next iterations

1. Implement the complete N = 1 list/config/new-bundle path through generate,
   apply, and full; add singleton omnibus/partition reporting, shared models and
   statuses, bounded lookup loading, and migrate/remove replaced legacy paths.
   Adapt reference workflow inputs/readers to the new contract while preserving
   the frozen scientific values. Pass every tier before moving on.
2. Add N > 1 ordered biological ownership, per-pair eligibility, ID-overlap audit,
   N-by-M lookups, and N partition reports; preserve N = 1 regression equality.
3. Measure endpoints on a real reference bundle, select/calibrate the upper-tail
   policy, implement complete ACAT reporting, and finish the contract audit.

This iteration is the prerequisite refactor, not completion of the biological
partition or ACAT feature. Real-bundle upper-tail measurements are still pending.

## Iteration 2: shared models, statuses, ACAT, and explicit N = 1 inputs

Completed on 2026-10-21:

- Added `workflow/model_catalog.json` as the single source for the five public
  disease models, their roles/scopes/genotype classes, control compatibility,
  and the full shared status registry. Python and R load and validate the same
  catalog; finalization and reporting no longer keep separate status lists.
- Added the pure interior ACAT implementation and the endpoint-policy wrapper.
  It validates inputs before endpoint handling, gives raw zero precedence over
  a simultaneous one, preserves singleton values, and marks unresolved
  multi-partition upper endpoints explicitly. The mid-p code path is present
  for the required evaluation, but production policy remains intentionally
  unselected pending bundle measurements.
- Replaced the implicit/default testing-unit BED with an ordered list of safe
  partition bindings and exact identity-map matching in generate/full mode.
  Apply rejects both partition identities and physical partition overrides.
  Identities no longer live in `config/config.yaml`, preventing deep-merge
  leakage into apply runs. Active generation/full fixtures and examples now
  use a single explicit `cds` partition.
- Retained the existing internal single-BED scope adapter only for the N = 1
  migration point. It deliberately raises for N > 1 until the next iteration
  makes the scope, membership, bucket, and reduction tables partition-aware.

Focused verification passed: cross-language catalog/finalization/scoring tests,
the ACAT analytic and endpoint tests, direct generate/full plan resolution for
all active null-distribution fixtures, and the generate workflow test (2 tests).

## Iteration 3: partition-only generation and apply reporting

Completed on 2026-10-21:

- Removed role-wide bucket lookup outputs from the active generation graph.
  Buckets now render only physical `(partition_id, model_id)` tables and their
  partition/model diagnostics.
- Apply reads one partition/model lookup table at a time, retains certified
  partition ownership through candidate and control matching, and publishes an
  omnibus ACAT report plus one partition report for each configured partition.
  JSON arrays preserve configured partition order and write absent values as
  JSON `null` without numeric truncation.
- Added a two-partition full-mode fixture. A clean run produces ten physical
  lookup tables (two partitions by five enabled models), certified resolved and
  eligible BEDs, and both partition summary files.
- Removed the dead prepared-lookup stage, external legacy lookup reader branch,
  and its unused execution policy. Apply accepts only certified partitioned
  bundles.

Verification: the clean two-partition full fixture completed with Snakemake
9.23.1; partition scope, lookup, audit, and bundle-contract Python tests (7)
and the R partition-reporting test passed. The implementation retains only the
active partitioned workflow path; obsolete unit tests that directly import the
deleted legacy validator still need replacement before the complete historical
test suite can be considered migrated.

## Iteration 4: contract audit and strict partition-only execution

Completed on 2026-10-21:

- Removed the remaining pooled execution-scope artifacts and compatibility
  branch from generation-window preparation. Windows and buckets now derive
  from the union of eligible, resolved, calibration-masked biological
  partitions, while each testing unit still receives one deterministic bucket.
- Removed lookup-header guessing. Apply accepts only the certified
  `testing_unit_id`-headed partitioned lookup tables.
- Removed the obsolete candidate-signature cohort summary path. The standalone
  annotation utility now explicitly consumes `cohort_candidates.tsv.gz` using
  `testing_unit_id` and one `variant_key` per row.
- Re-ran focused generation tests, schema tests, workflow DAG tests, and both
  two-partition full fixtures after these removals.

The remaining external requirement is unchanged: the repository contains no
real generated reference bundle and matched cohort with which to measure
endpoint prevalence and select/calibrate a production multi-partition
upper-tail policy. The workflow emits bounded endpoint and raw-versus-mid-p
diagnostics and marks the policy as pending; it must not be represented as a
completed production calibration until those supplied resources are measured.
