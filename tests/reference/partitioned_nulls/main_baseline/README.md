# Immutable reference for partitioned nulls

Captured on 2026-10-20 from `main` commit
`b3c967cab0ca35fe537ae8a7b070d0cf6639a002`, before changing application code.
The source was an isolated `git archive` extraction. Only the additional test
inputs under `inputs/` were added to that extraction; no runtime code was changed.

Six successful, certified full-mode runs are recorded in `manifest.json`, with
exact configurations, commands, package versions, input hashes, public-output
hashes, and reference-file hashes. The ordinary SNV/HEM, SNV/no-HEM, and indel
fixtures preserve existing coverage. The three `models_*` cases add affected
male/female samples to exercise all five disease models, four models without
HEM, and a retained insertion with positive modeled event mass. The indel cases
use the nondefault 50-point maximum and 0.3 bin width. No external legacy grid
was used.

Each case contains original compressed per-role lookup rows, original scored
variant rows per sample, and `sample_psaps.tsv`: the numerical minimum PSAP per
`(pid, testing_unit_id, model)` derived from those scored rows. Missing groups
are absent, not imputed. These primitive data remain usable after old bundle
readers are deleted; no old runtime implementation is embedded here.

To reproduce a case, extract the pinned revision into an isolated directory,
copy this `inputs/` directory to the same repository-relative path, and run from
that directory with the recorded CI environment on `PATH`:

```bash
snakemake --workflow-profile none --snakefile workflow/Snakefile \
  --configfile /path/to/main_baseline/models_snv_hem/run.yaml \
  --config outdir=/tmp/psap-reference-reproduction --cores 2
```

The `none` profile uses the installed tools; the default workflow profile would
install per-rule Conda environments instead. The absolute temporary paths in the
manifest are historical provenance, not runtime dependencies for these tests.
Outputs were copied byte-for-byte; scalar minima preserve the selected original
numeric text. The reference checksums must not be refreshed from feature-branch
outputs to make a regression pass.

During the mechanical window rename, `test_reference_workflow.py` compares every
recorded public component digest. During the subsequent intentional new-contract
migration, adapt that test's run configurations and output readers to the new
contract and compare the frozen scientific rows/values and model projections.
Keep the original captured files and provenance unchanged. Additional PSAP ties
may add supporting variants later, but must not change the scalar reference PSAPs.
