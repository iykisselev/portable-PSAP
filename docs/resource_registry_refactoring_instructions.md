# Resource Registry Refactoring Instructions

## Decision

Stop copying `ApplyResourceRegistry` values into `WORKFLOW_*` and
`VALIDATION_*` globals.  Rule inputs should name the resource they consume and
the representation required by that rule (`source`, `source_index`,
`staged_data`, or `staged_index`).  This preserves the registry as the single
apply-resource boundary.

This is a behavior-preserving refactor.  It must not change resource staging,
validation inputs, optional-input semantics, wildcard constraints, or the DAG.

The same cleanup is appropriate for the aliases at the top of
`rules/apply/annotation.smk`.  The generation file needs a separate,
follow-on context-object refactor; `GenerationPlan` is not an
`ApplyResourceRegistry`, and several of its current values are derived layouts
rather than registry lookups.

## Important constraints

Do not use a generic `res(name, field="staged_data")` helper.  A default hides
the critical distinction between an external source used by validation or
staging and a staged artifact used downstream.  Prefer these explicit helpers
in `rules/apply/resources.smk`, after `APPLY_RESOURCE_REGISTRY` is created:

```python
def apply_resource(name):
    return APPLY_RESOURCE_REGISTRY.resource(name)


def optional_input(value):
    return [value] if value else []


def optional_resource_field(name, field):
    resource = APPLY_RESOURCE_REGISTRY.resources.get(name)
    return getattr(resource, field) if resource else ""
```

Call sites should use, for example,

```python
cadd_snv=apply_resource("cadd_snv").staged_data,
cadd_snv_index=optional_input(apply_resource("cadd_snv").staged_index),
```

For validation and staging, use `.source` and `.source_index` explicitly.
Keep a condition separate when whether an input exists is governed by model or
bundle state rather than by the value alone.  In particular, do not replace
the following with a bare `optional_input(...)` call without retaining their
guards:

- hemizygous lookup input (`model.hemizygote`);
- optional lookup metadata, validation, and marker inputs, which remain
  value-driven even for a non-bundle lookup source;
- bundle-only eligible units, calibration mask, excluded alleles, and excluded
  units (`LOOKUP_SOURCE.is_bundle`);
- enabled CADD indels (`model.include_indels`);
- plural excluded-allele inputs, which remain `list(LOOKUP_SOURCE.excluded_alleles)`.

`optional_input` is for one optional scalar input only.  It must not be used
to convert an existing iterable into a nested list.

## Stage 0: record the baseline

Before editing, run:

```bash
python tests/run_tier.py unit
python tests/run_tier.py workflow
git diff --check
```

Do not update `tests/workflow_baseline.json` for this refactor.  Initial and
resolved job counts, terminal outputs, and required rules must remain
unchanged.  The lint allowlist is empty for all modes, so new lint findings are
regressions rather than expected output.

## Stage 1: add narrow apply-resource helpers

1. In `workflow/rules/apply/resources.smk`, add `apply_resource(name)` and
   `optional_input(value)` after `APPLY_RESOURCE_REGISTRY` is constructed.
2. Add unit coverage for the pure registry lookups in
   `tests/test_apply_resources.py`.  If helpers remain in a `.smk` file, cover
   their behavior through the workflow tier rather than importing the file as a
   Python module.
3. Do not place `apply_resource` in `deployment.smk` or `paths.smk`: the
   registry does not exist when those files load, and generation mode does not
   own an apply registry.

## Stage 2: reduce `resources.smk` to routing state

Replace field-by-field aliases with direct registry queries in downstream
rules.  Retain only values that are meaningful routing or dynamic-rule state:

- `APPLY_RESOURCE_REGISTRY` itself;
- table and VCF wildcard patterns;
- the eager table-resource policy map and its callback;
- callbacks such as `table_resource_source`, `table_resource_cohort`, and
  `vcf_resource_source` (rewritten to query the registry directly);
- derived paths that have no `ResourceSpec` equivalent, such as the combined
  exclusions VCF path.

Delete aliases only after their callers have moved.  In particular, avoid
keeping both `WORKFLOW_CADD_SNV` and `apply_resource("cadd_snv").staged_data`.

Suggested call-site substitutions:

| Current purpose | Explicit replacement |
| --- | --- |
| Prepared cohort VCF | `apply_resource("cohort").staged_data` |
| Prepared cohort index | `apply_resource("cohort").staged_index` |
| Prepared calibration regions | `apply_resource("calibration_regions").staged_data` |
| Prepared testing units | `apply_resource("testing_units").staged_data` |
| Prepared CADD SNV/indel | `apply_resource("cadd_snv").staged_data` / `apply_resource("cadd_indel").staged_data` |
| Validation source and index | `apply_resource(name).source` / `.source_index` |

For optional registry entries, use `APPLY_RESOURCE_REGISTRY.resources.get(name)`
once at the call site or introduce a clearly named `optional_apply_resource`
or `optional_resource_field(name, field)` helper.  The field must remain an
explicit argument.  Do not make the required-resource helper silently return
`None`.

## Stage 3: migrate apply rules by data-flow boundary

Migrate in these groups, running the relevant tests after each group.

1. `rules/apply/annotation.smk`: remove the CADD and testing-unit aliases at
   the top of the file.  Use staged fields directly in `annotate_variants`.
2. `rules/apply/variants.smk`, `samples.smk`, `scoring.smk`, and `prepare.smk`:
   replace prepared resource globals with staged registry fields.  Preserve
   optional calibration-region and exclusion behavior exactly.
3. `rules/apply/validate.smk`: replace `VALIDATION_*` registry-derived values
   with explicit source fields and `optional_input(...)` only where appropriate.
4. `rules/apply/cohort.smk`: confirm it has no remaining flattened registry
   dependency.

At this point, `validation_helpers.smk` should contain only real shared state:
validation report paths, validation-library paths, validation schema path, and
possibly small callbacks that cannot be expressed directly in a rule.  Delete
the `VALIDATION_*` mirror block and `required_resource`/`optional_resource`
aliases once all callers use `apply_resource`.

## Stage 4: handle exclusions without reintroducing mirrors

Exclusions are plural resources, so they are the main exception to scalar
lookup examples.  Prefer a focused registry-derived helper, for example one
that returns the staged exclusions or source exclusions for a requested
scientific role.  Preserve input order: exclusion order is scientific and is
used when projecting and combining exclusions.

Keep the combined-exclusions derived VCF separate from configured resources;
it is a run-derived registry member with a distinct producer rule.

## Stage 5: refactor generation separately

Do not delete the 22 `GENERATION_*` aliases as part of the apply-registry
change.  First introduce a `GenerationRuleContext` in Python that owns:

- `GenerationPlan`;
- generated scope and execution layouts;
- derived partition, bucket, statistics, log, and bundle paths;
- small helper methods for optional mutation-rate and exclusion inputs.

Then make `rules/generate/context.smk` construct one `GENERATION_CONTEXT` and
delegate its existing wildcard callbacks to it.  Migrate generation rules to
that object gradually.  This keeps derived layout construction in one place
without falsely treating it as a resource registry.

## Tests and acceptance criteria

After each stage, run focused unit tests for the touched Python module.  On
completion, run:

```bash
python tests/run_tier.py unit
python tests/run_tier.py workflow
git diff --check
```

The completed refactor is acceptable only when:

- every apply validation input still reads the original configured source;
- every downstream apply rule still reads the staged artifact;
- optional, bundle-only, and model-enabled inputs retain the same empty-list
  versus populated-list behavior;
- exclusion order and combined-exclusion construction are unchanged;
- all three workflow modes match the recorded DAG baseline and have no lint
  findings;
- no `WORKFLOW_*` or `VALIDATION_*` alias remains solely to mirror a
  `ResourceSpec` field.
