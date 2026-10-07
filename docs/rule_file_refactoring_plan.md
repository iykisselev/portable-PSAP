# Rule-file boilerplate and target-interface refactoring plan

## Goal

Make rule definitions shorter and more uniform without changing scientific
inputs, resource-policy values, DAG topology, or generated artifacts. The end
state has one deployment vocabulary, one public default target (`all`) in every
mode, one source of generation target inputs, and compact variant-filter rules.

This plan deliberately separates syntax-only normalization from the larger
Snakemake inheritance change. Each stage should be independently reviewable
and leave the workflow executable.

## Baseline and guardrails

Before the first edit, record the current rule count, deployment-directive
count, job counts, and public targets. The current workflow has three logical
software environments (`validation`, `bioinformatics`, and `r`), but each rule
spells out its Conda, module, and container binding. Generate rules also embed
script and environment paths inline, while apply rules use path constants.

Preserve these invariants throughout:

- `generate` still emits the complete certified null-distribution bundle and
  recoverable partition/bucket outputs.
- `apply` still emits the three cohort deliverables.
- `full` still requires both sets of deliverables in one DAG.
- Every rule retains its existing `threads`, resource-policy callback, Conda
  environment, module set, and container image.
- Optional calibration-region and exclusion resources remain real input-file
  dependencies when enabled.

For every stage, run the unit tier, fresh DAG baseline test, all three
mode-specific lint runs, and `git diff --check`. Re-run the full smoke E2E tier
after stages 3 and 5 because they touch target selection and shell commands.

## Stage 1 — Canonical deployment and path helpers

**Risk:** low. **Expected behavior change:** none.

Add `env(tool)` in `workflow/rules/deployment.smk`, next to `mod(tool)` and
`sif(tool)`. It should return the pinned Conda file under `CONDA_ENV_DIR`, for
example `env("validation") -> .../envs/validation.yaml`. Keep all three
deployment mechanisms present on every rule; this stage only gives them a
consistent name.

Add `GENERATION_SCRIPT_DIR = Path(workflow.basedir) / "scripts" / "generate"`
in `workflow/rules/generate/context.smk`. Convert all generate-rule `conda`
and `script` directives to `env(...)` and `GENERATION_SCRIPT_DIR`, respectively.
Convert apply-rule `conda` directives to `env(...)` as well. This eliminates
the current mixture of inline `Path(...)`, `CONDA_ENV_DIR`, and branch-specific
script path construction.

Acceptance checks:

- No rule contains `str(CONDA_ENV_DIR / ...)` or an inline
  `Path(workflow.basedir) / "scripts"` expression.
- Every rule resolves to exactly the same environment/image/module binding as
  before. Add a focused test for `env`, `mod`, and `sif` if deployment helpers
  are not already directly tested.
- All DAG job counts remain unchanged.

## Stage 2 — Mechanical readability sweep

**Risk:** low. **Expected behavior change:** none.

Apply small, reviewable simplifications while Stage 1 has made deployment lines
uniform:

- Collapse needlessly multi-line `str(...)` log-path expressions in
  `annotation.smk`, `samples.smk`, and `cohort.smk` where the complete path is
  a single readable expression.
- Replace the identity-copy `LOOKUP_TABLES` mapping in `scoring.smk` with direct
  reads from `LOOKUP_SOURCE.tables`; retain the existing HEM-enabled condition.
- Add a short comment in `workflow/shared/bootstrap.R` explaining that all R
  library files are sourced unconditionally to give every adapter a stable
  runtime namespace. Do not change loading behavior in this stage.
- Correct `Coverage resource` in `docs/apply_workflow_architecture.md` to the
  current `calibration regions` / exclusion terminology and verify nearby prose
  describes it as an optional callable-region filter rather than a coverage
  measurement.

Acceptance checks:

- `rg` confirms the copied lookup mapping and obsolete coverage wording are
  gone.
- R syntax tests still pass.
- The diff changes no rule inputs, outputs, shell command arguments, or
  resource-policy names.

## Stage 3 — One default target and one generation output declaration

**Risk:** medium: this intentionally changes named Snakemake targets.

First make the public-interface decision explicit: retire
`null_distributions` and `full_results` in favor of `all`, matching the
existing apply-mode name. Do not retain aliases unless a documented external
automation consumer requires a transition period; aliases hide stale commands
and duplicate target definitions.

Create `generation_outputs(wildcards)` in generation context. It should return
the dynamic recoverable partition artifacts, statistics, lookup shards, and
diagnostics currently declared by `null_distributions`; retain the static
bundle marker explicitly in each target so an initial dry-run can display its
dependency before checkpoint expansion. The helper must keep the HEM lookup
conditional and must use checkpoint-aware helper functions only at DAG
evaluation time.

Then define exactly one `rule all` in each mode branch:

- `generate`: `generation_outputs`.
- `apply`: the three cohort deliverables.
- `full`: the three cohort deliverables plus `generation_outputs`.

This removes the duplicated full-mode generation input block while retaining
named inputs where they make the target output self-describing.

Update the README, `docs/pipeline.md`, and generate architecture documentation
to show `snakemake ... all` for every mode. Update workflow tests to assert
that the retired target names fail with `No rule to produce ...`, and retain a
fresh-DAG assertion that full mode schedules both branches.

Acceptance checks:

- `snakemake --dry-run all` works in apply, generate, and full fixtures.
- `snakemake --dry-run null_distributions` and `full_results` fail after the
  documented interface transition.
- Fresh DAG counts and E2E output digests are unchanged; only the selected
  target rule name changes.

## Stage 4 — Remove duplicated optional-filter shell bodies

**Risk:** medium: command construction changes, although semantics do not.

Refactor `filter_calibration_regions` first. Add one parameter built by a
small, testable Python helper: empty when calibration regions are absent, or a
shell-escaped `--regions-file <path>` argument when present. Use one `bcftools
view` invocation, one log message that states whether filtering is enabled,
and the existing index command.

`filter_excluded_records` cannot be reduced to the identical command shape:
with exclusions it uses `bcftools isec -C -w1`, and without exclusions it uses
`bcftools view`. Keep that scientific operation distinction. Reduce duplicated
logging and indexing around the branch, but do not force the two subcommands
through a synthetic common command. This corrects the suggestion's assumption
that both rules can be halved in the same way.

Use `shlex.quote` (or an equivalently tested Snakemake-safe formatter) when
constructing any optional shell fragment. Add a unit test with a path containing
spaces and a workflow dry-run for both resource-present and resource-absent
fixtures.

Acceptance checks:

- Calibration-enabled and calibration-absent runs produce the same commands
  and output paths as before, apart from harmless log wording.
- Exclusion-enabled and exclusion-absent runs retain their respective
  `isec` and `view` semantics.
- The filter E2E fixture preserves its expected output digest.

## Stage 5 — Evaluate and, if proven safe, adopt base deployment rules

**Risk:** medium. **Prerequisite:** Stage 1 is complete.

Prototype three abstract rules in an isolated branch or a small fixture:
`base_validation`, `base_bioinformatics`, and `base_r`. Each carries only the
deployment triple (`conda: env(...)`, `envmodules: *mod(...)`, and
`container: sif(...)`). Do not put threads or resources in the base rules:
those are rule-specific policy declarations and must remain visible at the
concrete call site.

Use Snakemake's `use rule ... with:` inheritance only after proving that it
preserves all of the following in the installed minimum-supported Snakemake
version:

- Conda, module, and container directives in local, Slurm, and containerized
  execution profiles.
- Rule-local `threads`, `resources`, `benchmark`, `log`, `script`, and
  `wildcard_constraints` declarations.
- Correct rule names and scheduler resource callbacks in dry-run output.

If the prototype is clean, convert one environment family at a time:
validation rules, then bioinformatics rules, then R rules. Keep each family in
a separate commit. If inheritance obscures a directive, stop after Stages 1–4:
the `env` helper already removes path drift with far lower maintenance cost.

Acceptance checks:

- Every concrete rule has the same declared resource-policy name and job
  allocation as before.
- Lint is clean and the workflow baseline changes only if Snakemake reports a
  different synthetic base-rule representation; functional job counts must not
  change.
- Apply, generate, and full E2E fixtures complete successfully.

## Stage 6 — Consolidate, document, and protect the result

**Risk:** low after prior stages.

Remove compatibility comments and temporary target-transition notes only after
downstream automation has moved to `all`. Add a compact contributor note to
the workflow architecture docs:

- choose a logical deployment name through `env`, `mod`, and `sif`;
- choose the appropriate base rule if Stage 5 was adopted;
- add target inputs through the shared generation helper rather than copying a
  checkpoint-driven block; and
- keep optional shell arguments constructed through the tested helper.

Finish by updating the recorded workflow baseline only when a deliberate
target/interface change affects its expected text. Do not update it merely to
mask an unexpected DAG change.

## Recommended commit boundaries

1. Deployment/path helper normalization.
2. Readability and documentation sweep.
3. Unified `all` target plus target/docs/tests migration.
4. Calibration filter command simplification and filter regressions.
5. Base-rule inheritance proof of concept and, if accepted, environment-family
   conversions.
6. Final architecture documentation and baseline reconciliation.
