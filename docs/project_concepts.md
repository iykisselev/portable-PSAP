# Portable PSAP: concepts and scientific logic

Portable PSAP prioritizes rare variants in a cohort by comparing their
deleteriousness with a population-derived expectation for the genomic region in
which they occur. It implements the single-case PSAP framework and a
testing-unit (for example, gene) extension. It does not establish that a
variant causes disease: it produces reproducible, model-aware evidence for
prioritization.

## The question the workflow answers

For a candidate variant in a testing unit, the workflow asks: *under a chosen
inheritance model and the population variation represented by the lookup
table, how unusual is an event at least this deleterious in this unit?*

The answer is a PSAP value. Lower values indicate a more unusual event under
the selected null model. A value is meaningful only with the assembly,
testing-unit definition, coverage scope, exclusion policy, CADD release, and
genotype model that were used to construct its lookup table.

The workflow supports three ways to obtain these values:

| Mode | Purpose |
| --- | --- |
| `apply` | Score a cohort using supplied lookup tables or a previously certified bundle. |
| `generate` | Build and certify a lookup bundle from population allele frequencies and score resources. |
| `full` | Generate a bundle and immediately apply that same bundle to a cohort. |

## Scientific objects

### Testing units

A testing unit is the entity over which rare-variant events are evaluated.
The default unit type is a gene, but the workflow treats the BED's fourth
column as an opaque, stable identifier. The same identifiers must appear in
the lookup tables. A variant may belong to more than one unit and is then
evaluated once for every unit membership.

### Scores

CADD-PHRED is used as the variant-level severity score. During application,
the score is clipped to the lookup grid's configured range and mapped to a
fixed score bin. When a sample has several compatible variants in one unit and
model, all variants tied for the highest score are retained. This avoids an
arbitrary choice among equal best variants.

### Inheritance roles

The lookup roles correspond to genotype-event probabilities:

| Role | Compatible genotype event | Null probability at allele frequency `f` |
| --- | --- | --- |
| `het` | Diploid heterozygote | `2f(1-f)` |
| `hom` | Diploid homozygous alternate | `f²` |
| `hem_chrX` | Male non-PAR chromosome-X alternate | `f` |

The hemizygous role is enabled only when `model.hemizygote` is true. It is
restricted to non-pseudoautosomal chromosome X. With the sex-specific
frequency policy, female frequencies serve the X `het`/`hom` roles and male
frequencies serve `hem_chrX`.

### The null distribution

Generation creates a score-bin survival row for each testing unit and role.
For the alleles assigned to a unit, their compatible genotype-event
probabilities are treated as independent. The workflow aggregates the
probability of at least one event in each score-bin tail, conditional on there
being at least one compatible event in the unit. The resulting row maps a
score bin to the null probability of an event at least that severe.

This is the quantity subsequently looked up as `PSAP`; it is not recalculated
from the cohort. Consequently, a cohort run is deterministic for a fixed VCF,
bundle, and configuration.

## Defining the scientific scope

The calibration scope is deliberately narrower than a raw testing-unit BED:

1. The workflow starts with autosomes 1–22 and chromosome X represented by the
   FASTA index.
2. It removes pseudoautosomal intervals and any configured generation
   exclusions.
3. It intersects the remaining sequence with the callable/coverage resource.
4. It measures each unit's covered fraction over its supported, non-excluded
   bases and retains it only when it meets
   `generation_settings.coverage_minimum`.

The generated bundle publishes the surviving unit intervals, the exact
calibration mask, alleles rejected during generation, and units excluded for
coverage or zero event mass. Those artifacts make the lookup's scope
inspectable rather than implicit.

For alleles absent from the population panel, generation uses the explicit
`generation.unobserved_policy`:

- `zero` assigns an event probability of zero; or
- `mutation_rate` substitutes a supplied per-allele mutation rate, subject to
  the configured chromosome-X mutation-rate policy.

The policy is recorded in bundle metadata. It is a scientific choice, not a
performance option.

## Applying a lookup to a cohort

The apply workflow first verifies that cohort metadata, reference sequence,
resource formats, testing-unit IDs, and model settings are compatible. It then
stages resources, filters the VCF to callable regions, removes excluded
records, and attaches CADD and testing-unit annotations.

For every selected sample, the workflow extracts genotypes, standardizes their
representation, and applies sex-aware genotype classification. Compatible
genotypes are joined to testing-unit/CADD annotations and scored against the
role-specific lookup table.

Affected participants produce candidate records. Unaffected participants
produce control-evidence records. The final report compares an affected
candidate only with controls that carry the identical allele in the same
testing unit; dominant heterozygous candidates may use evidence from any
configured control model, while the other models require a model match.

## Why bundles are certified

A generated bundle contains more than lookup tables. It also includes a
metadata record, component checksums, scope artifacts, excluded units and
alleles, a validation report, and a `validation.ok` marker. The marker is
created only after structural and semantic validation succeeds. Apply mode
uses that certificate to reject incomplete or incompatible bundles before
expensive cohort processing begins.

In `full` mode, this guarantee is stronger: the generated lookup tables,
eligible testing units, calibration regions, and excluded alleles become the
effective apply resources. The configured FASTA and CADD score resource remain
the shared physical bindings. This prevents the generation and application
halves from silently using different scopes.

## Important interpretation limits

- PSAP ranks statistical rarity under the chosen population model; it is not a
  clinical classification or a replacement for segregation, phenotype, or
  functional evidence.
- Results assume all resources use the declared assembly and compatible contig
  naming.
- A score is only comparable within a matching lookup contract, binning
  scheme, testing-unit set, and model configuration.
- Coverage and exclusions change the scientific denominator. Replacing either
  requires generating a new compatible bundle rather than reusing an old one.

For implementation-level details, see
[the generation architecture](generate_workflow_architecture.md) and
[the apply architecture](apply_workflow_architecture.md).
