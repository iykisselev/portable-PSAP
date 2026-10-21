"""Typed file and parameter contracts for apply script adapters."""

from dataclasses import dataclass
from pathlib import Path
from typing import Mapping

PathLike = str | Path


@dataclass(frozen=True)
class ReportOutput:
    report: PathLike


@dataclass(frozen=True)
class CohortValidationInputs:
    metadata: PathLike
    vcf: PathLike
    vcf_index: PathLike | None
    fasta: PathLike
    fasta_index: PathLike
    fasta_gzi: PathLike | None


@dataclass(frozen=True)
class CohortValidationOutputs:
    report: PathLike
    manifest: PathLike


@dataclass(frozen=True)
class CohortValidationParameters:
    metadata_kind: str
    hemizygote_model: bool


@dataclass(frozen=True)
class ResourceValidationInputs:
    vcf: PathLike


@dataclass(frozen=True)
class ResourceValidationParameters:
    unit: str
    include_indels: bool
    resource_fields: Mapping[str, Mapping[str, str]]


@dataclass(frozen=True)
class ModelValidationInputs:
    lookup_tables: tuple[PathLike, ...]
    lookup_metadata: PathLike
    lookup_validation: PathLike
    lookup_marker: PathLike
    lookup_calibration_mask: PathLike
    lookup_excluded_alleles: tuple[PathLike, ...]
    resolved_partitions: tuple[PathLike, ...]
    eligible_partitions: tuple[PathLike, ...]
    partition_statuses: PathLike
    testing_unit_id_overlap: PathLike


@dataclass(frozen=True)
class ModelValidationParameters:
    score_max: float
    assembly: str
    unit: str
    hemizygote_model: bool
    include_indels: bool
    lookup_contract: str
    lookup_bin_width: float | None
    lookup_component_keys: tuple[str, ...]


@dataclass(frozen=True)
class PreparationInputs:
    source: PathLike
    source_index: PathLike | None
    cohort: PathLike | None = None


@dataclass(frozen=True)
class PreparationOutputs:
    data: PathLike
    index: PathLike


@dataclass(frozen=True)
class PreparationParameters:
    kind: str
    tsv_layout: str


@dataclass(frozen=True)
class ExclusionInputs:
    variants: PathLike


@dataclass(frozen=True)
class ExclusionOutputs:
    vcf: PathLike
    index: PathLike


@dataclass(frozen=True)
class ExclusionParameters:
    sources: tuple[str, ...]
    indexed_tsv_sources: tuple[PathLike, ...]


@dataclass(frozen=True)
class SplitInputs:
    vcf: PathLike


@dataclass(frozen=True)
class SplitOutputs:
    vcf: tuple[PathLike, ...]
    csi: tuple[PathLike, ...]


@dataclass(frozen=True)
class SplitParameters:
    samples: tuple[str, ...]


@dataclass(frozen=True)
class MergeInputs:
    cohort: PathLike
    resources: PathLike
    model: PathLike


@dataclass(frozen=True)
class GateInputs:
    report: PathLike


@dataclass(frozen=True)
class GateOutputs:
    ok: PathLike


@dataclass(frozen=True)
class GateParameters:
    failure_message: str
    include_errors: bool
    success_marker: str


@dataclass(frozen=True)
class AnnotationInputs:
    annotations: PathLike


@dataclass(frozen=True)
class AnnotationParameters:
    include_indels: bool
