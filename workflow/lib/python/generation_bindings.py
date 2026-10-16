"""Typed Snakemake-to-library bindings for generation operations.

The generation library is invoked only by its script adapters, but its inputs
are still a public internal contract: direct tests construct the same objects.
These frozen records make optional paths, ordered shards, and native settings
explicit instead of relying on anonymous attribute bags.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Mapping


PathLike = str | Path


@dataclass(frozen=True)
class GenerationScopeInputs:
    population_vcf: PathLike
    population_vcf_index: PathLike | None
    score_snv: PathLike
    score_snv_index: PathLike | None
    exclusion: PathLike | None
    testing_units: PathLike
    par_regions: PathLike
    exclusions: tuple[PathLike, ...]
    fasta: PathLike
    fasta_fai: PathLike
    mutation_rate: PathLike | None
    mutation_rate_index: PathLike | None
    score_indel: PathLike | None = None
    score_indel_index: PathLike | None = None
    mutation_rate_indel: PathLike | None = None
    mutation_rate_indel_index: PathLike | None = None
    testing_unit_partitions: tuple[PathLike, ...] = ()


@dataclass(frozen=True)
class GenerationScopeOutputs:
    windows: PathLike
    buckets: PathLike
    units: PathLike
    unit_intervals: PathLike
    masked_unit_intervals: PathLike
    effective_mask: PathLike
    resolved_partition_intervals: PathLike | None = None
    partition_statuses: PathLike | None = None


@dataclass(frozen=True)
class GenerationScopeParameters:
    minimum_retained_fraction: float
    window_size: int
    unit_bucket_size: int
    compression: str
    af_tags: Mapping[str, str]
    include_indels: bool
    indel_max_reference_length: int
    partition_ids: tuple[str, ...] = ()


@dataclass(frozen=True)
class GenerationWindowInputs:
    windows: PathLike
    units: PathLike
    unit_intervals: PathLike
    effective_mask: PathLike
    population_vcf: PathLike
    population_vcf_index: PathLike | None
    score_snv: PathLike
    score_snv_index: PathLike | None
    fasta: PathLike
    fasta_fai: PathLike
    score_indel: PathLike | None = None
    score_indel_index: PathLike | None = None


@dataclass(frozen=True)
class GenerationWindowOutputs:
    allele_states: PathLike
    allele_membership: PathLike
    allele_scores: PathLike
    candidate_audit: PathLike


@dataclass(frozen=True)
class GenerationWindowParameters:
    af_tags: Mapping[str, str]
    x_af_policy: str
    hemizygote: bool
    compression: str
    include_indels: bool
    indel_max_reference_length: int
    indel_max_alternate_length: int


@dataclass(frozen=True)
class WindowStatisticsInputs:
    windows: PathLike
    buckets: PathLike
    allele_states: PathLike
    allele_membership: PathLike
    allele_scores: PathLike
    mutation_rate: PathLike | None
    mutation_rate_index: PathLike | None
    mutation_rate_indel: PathLike | None = None
    mutation_rate_indel_index: PathLike | None = None


@dataclass(frozen=True)
class WindowStatisticsOutputs:
    statistics: PathLike


@dataclass(frozen=True)
class WindowStatisticsParameters:
    x_af_policy: str
    unobserved_policy: str
    x_mutation_rate_policy: str
    score_max: float
    bin_width: float
    roles: tuple[str, ...]
    compression: str


@dataclass(frozen=True)
class LookupReductionInputs:
    buckets: PathLike
    window_statistics: tuple[PathLike, ...]


@dataclass(frozen=True)
class LookupReductionOutputs:
    lookup_het: PathLike
    lookup_hom: PathLike
    lookup_hem_chrX: PathLike | None
    diagnostics: PathLike


@dataclass(frozen=True)
class LookupReductionParameters:
    bin_count: int
    roles: tuple[str, ...]


@dataclass(frozen=True)
class BundleFinalizationInputs:
    lookup_het: tuple[PathLike, ...]
    lookup_hom: tuple[PathLike, ...]
    lookup_hem_chrX: tuple[PathLike, ...]
    diagnostics: tuple[PathLike, ...]
    units: PathLike
    masked_unit_intervals: PathLike
    allele_states: tuple[PathLike, ...]
    candidate_audit: tuple[PathLike, ...]
    effective_mask: PathLike


@dataclass(frozen=True)
class BundleFinalizationOutputs:
    lookup_het: PathLike
    lookup_hom: PathLike
    lookup_hem_chrX: PathLike | None
    eligible_bed: PathLike
    calibration_regions: PathLike
    excluded_alleles: PathLike
    candidate_audit: PathLike
    population_af_zero_imputations: PathLike
    excluded_units: PathLike
    metadata: PathLike


@dataclass(frozen=True)
class BundleFinalizationParameters:
    bin_count: int
    minimum_retained_fraction: float
    bin_width: float
    score_max: float
    policy: str
    conversion: str
    assembly: str
    testing_unit_type: str
    af_tags: Mapping[str, str]
    x_af_policy: str
    x_mutation_rate_policy: str
    hemizygote: bool
    variant_classes: tuple[str, ...]
    resource_identities: Mapping[str, Mapping[str, str]]
    indel_candidate_universe: str
    indel_membership: str
    indel_mask_policy: str
    indel_max_reference_length: int
    indel_max_alternate_length: int


@dataclass(frozen=True)
class BundleValidationInputs:
    lookup_het: PathLike
    lookup_hom: PathLike
    lookup_hem_chrX: PathLike | None
    eligible_bed: PathLike
    calibration_regions: PathLike
    excluded_alleles: PathLike
    candidate_audit: PathLike
    population_af_zero_imputations: PathLike
    excluded_units: PathLike
    metadata: PathLike
    schema: PathLike


@dataclass(frozen=True)
class BundleValidationOutputs:
    validation: PathLike
    marker: PathLike
