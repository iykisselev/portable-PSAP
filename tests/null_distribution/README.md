# Null-distribution fixtures

These are independent fixtures for the null-distribution implementation. They
are deliberately small enough to inspect and calculate by hand.

- `probability_oracle.json` contains event probabilities and expected maximum
  distributions. `test_null_distribution_oracle.py` enumerates all Bernoulli
  outcomes rather than using production distribution code.
- `population_af.vcf` selects `AF_CUSTOM` as the overall configured INFO tag,
  and supplies deliberately distinct `AF_FEMALE` and `AF_MALE` tags for X;
  contains a different `AF_DECOY`, a `Number=R` field, a valid AF zero, missing
  AF, a failed-filter record and a scalar INFO field on a multiallelic record.
- `testing_units.bed`, `coverage.bed` and `par_regions.bed` exercise union-base
  coverage with overlapping intervals and a mixed PAR/non-PAR X unit.
- `generation_scores.tsv` covers every possible SNV in the effective
  calibration mask. Removing any required row is a hard preparation failure.
- `fixture_expectations.json` records AF state, interval and pending-model
  expectations for AF extraction and interval arithmetic.

`generated_bundle_v1` captures the fixed-width decimal-safe generated-grid
mapping for regression testing. The HEM workflow produces
`generated_bundle_v2`, which adds certified X-policy metadata and an enabled
`lookup_hem_chrX.tsv.gz` component. Both remain distinct from the existing
`external_legacy_v0` mapping.
