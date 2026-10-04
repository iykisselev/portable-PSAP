# Null-distribution fixtures

These are independent fixtures for the null-distribution implementation. They
are deliberately small enough to inspect and calculate by hand.

- `population_af.vcf` selects `AF_CUSTOM` as the overall configured INFO tag,
  and supplies deliberately distinct `AF_FEMALE` and `AF_MALE` tags for X;
  contains a different `AF_DECOY`, a `Number=R` field, a valid AF zero, missing
  AF, a failed-filter record and a scalar INFO field on a multiallelic record.
- `testing_units.bed`, `exclusion.bed` and `par_regions.bed` exercise retained
  fractions with overlapping intervals and a mixed PAR/non-PAR X unit.
- `generation_scores.tsv` covers every possible SNV in the effective
  calibration mask. Removing any required row is a hard preparation failure.

The workflow produces `generated_bundle_v2`, which uses the fixed-width
decimal-safe generated grid, certified X-policy metadata and an optional
`lookup_hem_chrX.tsv.gz` component. Its bin mapping is deliberately distinct
from the existing `external_legacy_v0` mapping.
