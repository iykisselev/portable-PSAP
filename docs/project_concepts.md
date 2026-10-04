# Portable PSAP: concepts and scientific logic

Portable PSAP helps prioritize rare variants in a cohort. It compares the
severity of an observed variant with what would be expected from population
variation in the same genomic region. The workflow implements the single-case
PSAP framework and extends it to testing units such as genes.

The result is reproducible, model-aware evidence for prioritization—not a
claim that a variant causes disease.

## What the workflow is asking

For each candidate, the workflow asks: *under the chosen inheritance model and
the population variation represented by the lookup table, how unusual is an
event at least this deleterious in this testing unit?*

The answer is a PSAP value. Smaller values mean that the event is less expected
under the selected null model. A value only makes sense alongside the
assembly, testing-unit definition, retained scope, exclusions, CADD release,
and genotype model used to build its lookup table.

There are three ways to use the workflow:

| Mode | Purpose |
| --- | --- |
| `apply` | Score a cohort using supplied lookup tables or a previously certified bundle. |
| `generate` | Build and certify a lookup bundle from population allele frequencies and score resources. |
| `full` | Generate a bundle and immediately apply that same bundle to a cohort. |

## The main scientific pieces

### Testing units

A testing unit is the thing over which rare-variant events are evaluated. The
default unit is a gene, but the workflow treats the BED's fourth column as a
stable identifier and does not assume what it means. Those same identifiers
must appear in the lookup tables. If a variant overlaps several units, it is
evaluated once for each membership.

### Scores

CADD-PHRED supplies the variant-level severity score. During application, the
score is clipped to the lookup grid's range and mapped to a fixed bin. If a
sample has several compatible variants with the same highest score in one unit
and model, all of them are kept. That avoids making an arbitrary choice among
equally strong candidates.

### Inheritance roles

The lookup tables are split into genotype roles. Their null probabilities are:

| Role | Compatible genotype event | Null probability at allele frequency `f` |
| --- | --- | --- |
| `het` | Diploid heterozygote | `2f(1-f)` |
| `hom` | Diploid homozygous alternate | `f²` |
| `hem_chrX` | Male non-PAR chromosome-X alternate | `f` |

The hemizygous role is enabled only when `model.hemizygote` is true, and only
for non-pseudoautosomal chromosome X. With sex-specific frequencies, female
frequencies are used for X `het`/`hom`, while male frequencies are used for
`hem_chrX`.

### The null distribution

Generation builds one score-bin survival row for each testing unit and role.
For the alleles assigned to a unit, compatible genotype-event probabilities are
treated as independent. The workflow then combines them into the probability
of seeing at least one event in each score-bin tail, conditional on seeing at
least one compatible event in the unit. Each row therefore maps a score bin to
the null probability of an event at least that severe.

That conditional tail probability is what application looks up as `PSAP`; it
is not recalculated from the cohort. For a fixed VCF, bundle, and configuration,
the cohort result is therefore deterministic.

## How the calibration scope is defined

The raw testing-unit BED is not used as-is. The calibration scope is built in
four steps:

1. The workflow starts with autosomes 1–22 and chromosome X represented by the
   FASTA index.
2. It removes pseudoautosomal intervals. Each unit's supported, non-PAR bases
   form the denominator for its retained fraction.
3. It subtracts the union of `exclusion` and the other interval exclusions.
   These inputs are regions to remove, whatever their source.
4. It divides the remaining unit bases by the denominator from step 2 and
   requires `generation_settings.minimum_retained_fraction` (default 0.5).
   Exactly 50% passes; a unit with no remaining bases never passes.

The generated bundle records the surviving unit intervals, the exact
calibration mask, rejected alleles, and units excluded for low retention or zero
event mass. These artifacts make the lookup's scope visible and reviewable.

For an allele missing from the population panel, generation follows the
explicit `generation.unobserved_policy`:

- `zero` assigns an event probability of zero; or
- `mutation_rate` substitutes a supplied per-allele mutation rate, subject to
  the configured chromosome-X mutation-rate policy.

The policy is written to the bundle metadata. It is a scientific assumption,
not a performance setting.

## Applying a lookup to a cohort

When applying a lookup, the workflow first checks the cohort metadata, reference
sequence, resource formats, testing-unit IDs, and model settings. It stages the
resources, keeps variants in callable regions, removes excluded records, and
adds CADD and testing-unit annotations.

For each selected sample, it extracts genotypes, standardizes their
representation, and applies sex-aware genotype classification. Compatible
genotypes are joined to their testing-unit and CADD annotations, then scored
against the appropriate role-specific lookup table.

Affected participants produce candidate records; unaffected participants
produce control evidence. The final report compares a candidate only with
controls carrying the identical allele in the same testing unit. Dominant
heterozygous candidates may use evidence from any configured control model,
while the other models require a matching model.

## Why bundles are certified

A generated bundle is more than a set of lookup tables. It also contains
metadata, component checksums, scope artifacts, excluded units and alleles, a
validation report, and a `validation.ok` marker. The marker is written only
after structural and semantic checks pass. Apply mode uses that certificate to
reject incomplete or incompatible bundles before expensive cohort processing
starts.

In `full` mode, the guarantee extends across both halves of the DAG. The
generated lookup tables, eligible units, calibration regions, and excluded
alleles become the effective apply resources. The configured FASTA and CADD
score resource remain the shared physical bindings, so generation and
application cannot silently drift to different scopes.

## Keep these interpretation limits in mind

- PSAP ranks statistical rarity under the chosen population model; it is not a
  clinical classification or a replacement for segregation, phenotype, or
  functional evidence.
- Results assume all resources use the declared assembly and compatible contig
  naming.
- A score is only comparable within a matching lookup contract, binning
  scheme, testing-unit set, and model configuration.
- Coverage and exclusions change the scientific denominator. Replacing either
  requires generating a new compatible bundle rather than reusing an old one.

For an end-to-end operational reference, see the
[pipeline guide](pipeline.md). For implementation-level details, see the
[generation architecture](generate_workflow_architecture.md) and
[apply architecture](apply_workflow_architecture.md).
