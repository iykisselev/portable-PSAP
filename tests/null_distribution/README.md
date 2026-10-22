# Null-distribution fixtures

These independent fixtures cover the null-distribution implementation. They
are small enough to inspect and calculate by hand.

- `population_af.vcf` selects `AF_CUSTOM` as the overall configured INFO tag,
  and supplies deliberately distinct `AF_FEMALE` and `AF_MALE` tags for X;
  contains a different `AF_DECOY`, a `Number=R` field, a valid AF zero, missing
  AF, a failed-filter record and a scalar INFO field on a multiallelic record.
- `testing_units.bed`, `exclusion.bed` and `par_regions.bed` exercise retained
  fractions with overlapping intervals and a mixed PAR/non-PAR X unit.
- `generation_scores.tsv` covers every possible SNV in the effective
  calibration mask. Removing any required row is a hard preparation failure.

Every fixture produces `partitioned_null_bundle_v1`. Each configured
partition has one certified physical lookup table for every enabled disease
model, alongside resolved and eligible partition BEDs. The bundle records the
ordered overlap policy, candidate audit, and scope-status diagnostics.
