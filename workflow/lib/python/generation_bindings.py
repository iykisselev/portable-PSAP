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


@dataclass(frozen=True)
class GenerationScopeOutputs:
    partitions: PathLike
    buckets: PathLike
    units: PathLike
    unit_intervals: PathLike
    masked_unit_intervals: PathLike
    mask: PathLike
    resource_identities: PathLike


@dataclass(frozen=True)
class GenerationScopeParameters:
    minimum_retained_fraction: float
    partition_size: int
    unit_bucket_size: int
    compression: str


@dataclass(frozen=True)
class GenerationPartitionInputs:
    partitions: PathLike
    units: PathLike
    unit_intervals: PathLike
    mask: PathLike
    population_vcf: PathLike
    population_vcf_index: PathLike | None
    score_snv: PathLike
    score_snv_index: PathLike | None
    fasta: PathLike
    fasta_fai: PathLike


@dataclass(frozen=True)
class GenerationPartitionOutputs:
    allele_states: PathLike
    allele_membership: PathLike
    allele_scores: PathLike


@dataclass(frozen=True)
class GenerationPartitionParameters:
    af_tags: Mapping[str, str]
    x_af_policy: str
    hemizygote: bool
    compression: str


@dataclass(frozen=True)
class PartitionStatisticsInputs:
    partitions: PathLike
    buckets: PathLike
    allele_states: PathLike
    allele_membership: PathLike
    allele_scores: PathLike
    mutation_rate: PathLike | None
    mutation_rate_index: PathLike | None


@dataclass(frozen=True)
class PartitionStatisticsOutputs:
    statistics: PathLike


@dataclass(frozen=True)
class PartitionStatisticsParameters:
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
    partition_statistics: tuple[PathLike, ...]


@dataclass(frozen=True)
class LookupReductionOutputs:
    lookup_het: PathLike
    lookup_hom: PathLike
    lookup_hem: PathLike | None
    diagnostics: PathLike


@dataclass(frozen=True)
class LookupReductionParameters:
    bin_count: int
    roles: tuple[str, ...]


@dataclass(frozen=True)
class BundleFinalizationInputs:
    lookup_het: tuple[PathLike, ...]
    lookup_hom: tuple[PathLike, ...]
    lookup_hem: tuple[PathLike, ...]
    diagnostics: tuple[PathLike, ...]
    units: PathLike
    masked_unit_intervals: PathLike
    allele_states: tuple[PathLike, ...]
    effective_mask: PathLike
    resource_identities: PathLike


@dataclass(frozen=True)
class BundleFinalizationOutputs:
    lookup_het: PathLike
    lookup_hom: PathLike
    lookup_hem: PathLike | None
    eligible_bed: PathLike
    calibration_regions: PathLike
    excluded_alleles: PathLike
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
    resource_ids: Mapping[str, str]


@dataclass(frozen=True)
class BundleValidationInputs:
    lookup_het: PathLike
    lookup_hom: PathLike
    lookup_hem: PathLike | None
    eligible_bed: PathLike
    calibration_regions: PathLike
    excluded_alleles: PathLike
    excluded_units: PathLike
    metadata: PathLike
    schema: PathLike


@dataclass(frozen=True)
class BundleValidationOutputs:
    validation: PathLike
    marker: PathLike
