"""Immutable resolution of Portable PSAP workflow configuration.

Snakemake parses rule modules in one shared Python namespace.  Resolving the
run mode and physical resource bindings here keeps that parsing layer from
rewriting configuration or maintaining parallel, partially-effective binding
dictionaries.
"""

from __future__ import annotations

from copy import deepcopy
from dataclasses import dataclass, field
import json
import math
from pathlib import Path
import re
from types import MappingProxyType
from typing import Any, Mapping, TypeVar

from generation_resource_policy import check_generation_index_binding
from config_validation import ConfigurationValidationError, validate_configuration
from io_contracts import is_missing_token
from workflow_bindings import as_inputs
from bundle_contract import (
    BUNDLE_ARTIFACT_FILENAMES,
    CURRENT_GENERATED_BUNDLE_CONTRACT as GENERATED_BUNDLE_CONTRACT,
    INDEL_GENERATED_BUNDLE_CONTRACT,
    require_supported_contract,
    fixed_width_bin_count,
)
from partitioned_bundle_contract import (
    PARTITIONED_NULL_BUNDLE_CONTRACT,
    bundle_component_paths,
    require_partition_id,
)


class WorkflowPlanError(ValueError):
    """Raised when a valid-looking configuration cannot form one workflow plan."""


def normalize_run_config(config: Mapping[str, Any]) -> dict[str, Any]:
    """Copy a manifest without applying removed configuration migrations."""
    return deepcopy(dict(config))


def _reject_legacy_exclusion_lists(config: Mapping[str, Any]) -> None:
    """Give users an actionable error for the ID-keyed exclusion migration."""
    reference = config.get("reference")
    if isinstance(reference, Mapping) and isinstance(reference.get("exclusions"), list):
        raise WorkflowPlanError(
            "reference.exclusions must be a mapping keyed by exclusion ID; "
            "see docs/pipeline.md#keyed-reference-exclusions"
        )
    bindings = config.get("bindings")
    if not isinstance(bindings, Mapping):
        return
    bound_reference = bindings.get("reference")
    if isinstance(bound_reference, Mapping) and isinstance(
        bound_reference.get("exclusions"), list
    ):
        raise WorkflowPlanError(
            "bindings.reference.exclusions must be a mapping keyed by exclusion ID; "
            "see docs/pipeline.md#keyed-reference-exclusions"
        )


def _mapping(value: Any, key: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise WorkflowPlanError(f"{key} must be a mapping")
    return value


def _required(value: Any, key: str) -> str:
    if is_missing_token(value):
        raise WorkflowPlanError(f"Required configuration value is missing: {key}")
    path = str(value).strip()
    return str(Path(path).expanduser())


def _optional(value: Any) -> str | None:
    if is_missing_token(value):
        return None
    path = str(value).strip()
    return str(Path(path).expanduser())


def _number(value: Any, key: str, converter: type[float] | type[int]) -> float | int:
    try:
        result = converter(value)
    except (TypeError, ValueError, OverflowError) as exc:
        raise WorkflowPlanError(f"{key} must be a {converter.__name__}") from exc
    if isinstance(result, float) and not math.isfinite(result):
        raise WorkflowPlanError(f"{key} must be finite")
    return result


def _boolean(value: Any, key: str) -> bool:
    """Require a real configuration boolean rather than coercing truthy values."""
    if not isinstance(value, bool):
        raise WorkflowPlanError(f"{key} must be a boolean")
    return value


K = TypeVar("K")
V = TypeVar("V")


def _frozen_mapping(values: Mapping[K, V]) -> Mapping[K, V]:
    return MappingProxyType(dict(values))


def generated_bundle_paths(bundle_dir: Path) -> Mapping[str, str]:
    """Return non-lookup public paths for a partitioned generated bundle."""
    return _frozen_mapping(
        {
            "metadata": str(bundle_dir / "metadata.json"),
            "validation": str(bundle_dir / "validation.json"),
            "marker": str(bundle_dir / "validation.ok"),
            "calibration_regions": str(bundle_dir / "calibration_regions.bed.gz"),
            "excluded_alleles": str(bundle_dir / "excluded_alleles.tsv.gz"),
            "candidate_audit": str(bundle_dir / "candidate_audit.tsv.gz"),
            "population_af_zero_imputations": str(
                bundle_dir / "population_af_zero_imputations.tsv.gz"
            ),
            "unit_partition_status": str(bundle_dir / "unit_partition_status.tsv.gz"),
            "testing_unit_id_overlap": str(bundle_dir / "testing_unit_id_overlap.tsv.gz"),
        }
    )


@dataclass(frozen=True)
class FileBinding:
    """One data file and its optional declared index."""

    data: str | None
    index: str | None = None

    def data_inputs(self) -> list[str]:
        return as_inputs(self.data)

    def index_inputs(self) -> list[str]:
        return as_inputs(self.index)

    @classmethod
    def from_config(
        cls, value: Any, key: str, *, required: bool = True
    ) -> "FileBinding":
        mapping = _mapping(value, key)
        data = (
            _required(mapping.get("data"), f"{key}.data")
            if required
            else _optional(mapping.get("data"))
        )
        return cls(data=data, index=_optional(mapping.get("index")))


@dataclass(frozen=True)
class TestingUnitPartition:
    """One ordered biological partition and its source BED binding."""

    identifier: str
    bed: FileBinding


def _testing_unit_partitions(value: Any, key: str) -> tuple[TestingUnitPartition, ...]:
    if not isinstance(value, list) or not value:
        raise WorkflowPlanError(f"{key} must be a nonempty ordered list")
    partitions: list[TestingUnitPartition] = []
    seen: set[str] = set()
    for index, item in enumerate(value):
        item_key = f"{key}[{index}]"
        mapping = _mapping(item, item_key)
        identifier = _required(mapping.get("id"), f"{item_key}.id")
        if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_-]*", identifier):
            raise WorkflowPlanError(
                f"{item_key}.id must match ^[A-Za-z0-9][A-Za-z0-9_-]*$"
            )
        if identifier in seen:
            raise WorkflowPlanError(f"{key} contains duplicate partition ID {identifier!r}")
        seen.add(identifier)
        partitions.append(
            TestingUnitPartition(
                identifier=identifier,
                bed=FileBinding.from_config(mapping.get("bed"), f"{item_key}.bed"),
            )
        )
    return tuple(partitions)


def _validate_partition_identity_set(
    identities: Mapping[str, str],
    partitions: tuple[TestingUnitPartition, ...],
    *,
    identity_key: str,
    binding_key: str,
) -> None:
    identity_ids, binding_ids = set(identities), {item.identifier for item in partitions}
    missing, extra = sorted(binding_ids - identity_ids), sorted(identity_ids - binding_ids)
    if missing or extra:
        details = []
        if missing:
            details.append("missing identities: " + ", ".join(missing))
        if extra:
            details.append("identities without bindings: " + ", ".join(extra))
        raise WorkflowPlanError(
            f"{identity_key} and {binding_key} must have exactly the same IDs; "
            + "; ".join(details)
        )


@dataclass(frozen=True)
class FastaBinding:
    """One reference FASTA and its explicitly declared companion indexes."""

    data: str
    fai: str
    gzi: str | None

    @classmethod
    def from_config(cls, value: Any, key: str) -> "FastaBinding":
        mapping = _mapping(value, key)
        return cls(
            data=_required(mapping.get("data"), f"{key}.data"),
            fai=_required(mapping.get("fai"), f"{key}.fai"),
            gzi=_optional(mapping.get("gzi")),
        )


@dataclass(frozen=True)
class ExclusionMetadata:
    """Scientific identity and role of one exclusion binding."""

    identifier: str
    kind: str
    role: str
    enabled: bool


@dataclass(frozen=True)
class ExclusionBinding:
    """One active exclusion's stable identity, scientific role, and path."""

    identifier: str
    kind: str
    role: str
    binding: FileBinding


@dataclass(frozen=True)
class ReferenceMetadata:
    """Scientific resource identities that accompany a physical binding set."""

    assembly: str
    fasta_id: str
    exclusion_id: str
    exclusion_kind: str
    testing_unit_partitions: Mapping[str, str]
    cadd_snv_id: str
    cadd_indel_id: str | None
    exclusions: Mapping[str, ExclusionMetadata]


@dataclass(frozen=True)
class ReferenceBindings:
    """Immutable physical bindings used to execute one reference scope."""

    fasta: FastaBinding
    exclusion: FileBinding
    testing_unit_partitions: tuple[TestingUnitPartition, ...]
    exclusions: Mapping[str, FileBinding]
    cadd_snv: FileBinding
    cadd_indel: FileBinding | None
    calibration_regions: FileBinding | None = None

@dataclass(frozen=True)
class LookupSource:
    """Resolved lookup tables and optional generated-bundle certificate files."""

    source_type: str
    contract: str
    tables: Mapping[tuple[str, str], str]
    partition_ids: tuple[str, ...]
    score_max: float
    bin_width: float | None
    metadata: str | None = None
    validation: str | None = None
    marker: str | None = None
    resolved_partitions: Mapping[str, str] = field(default_factory=lambda: MappingProxyType({}))
    eligible_partitions: Mapping[str, str] = field(default_factory=lambda: MappingProxyType({}))
    partition_statuses: str | None = None
    calibration_mask: str | None = None
    excluded_alleles: tuple[str, ...] = ()

    @property
    def is_bundle(self) -> bool:
        return self.source_type == "bundle"


@dataclass(frozen=True)
class ModelPlan:
    unit: str
    hemizygote: bool
    include_indels: bool
    score_max: float


@dataclass(frozen=True)
class GenerationSettings:
    minimum_retained_fraction: float
    score_bin_width: float
    mutation_conversion: str
    variant_classes: tuple[str, ...]
    indel_candidate_universe: str
    indel_membership: str
    indel_mask_policy: str
    indel_max_reference_length: int
    indel_max_alternate_length: int
    x_af_policy: str
    x_mutation_rate_policy: str


@dataclass(frozen=True)
class GenerationPlan:
    """Generation inputs, policy, and paths for generate/full modes."""

    output_dir: Path
    log_dir: Path
    population_vcf: FileBinding
    score_snv: FileBinding
    score_indel: FileBinding | None
    exclusion: FileBinding
    testing_unit_partitions: tuple[TestingUnitPartition, ...]
    fasta: FastaBinding
    par_regions: FileBinding
    exclusions: tuple[ExclusionBinding, ...]
    mutation_rate: FileBinding | None
    mutation_rate_indel: FileBinding | None
    unobserved_policy: str
    af_source: Mapping[str, str]
    window_size_bp: int
    unit_bucket_size: int
    parquet_compression: str

    @property
    def scalable_dir(self) -> Path:
        return self.output_dir / "scalable"

    @property
    def bundle_dir(self) -> Path:
        return self.output_dir / "bundle"


@dataclass(frozen=True)
class InputBindings:
    vcf: FileBinding
    ped: str | None
    samples: str | None


@dataclass(frozen=True)
class ApplyPlan:
    """Apply inputs and the configured/effective reference distinction."""

    output_dir: Path
    log_dir: Path
    inputs: InputBindings
    configured_reference: ReferenceBindings
    effective_reference: ReferenceBindings
    effective_exclusions: tuple[ExclusionBinding, ...]
    lookup_source: LookupSource


@dataclass(frozen=True)
class ReportPlan:
    candidate_id_separator: str
    candidate_id_ordinal_width: int
    deduplicate_controls: bool


@dataclass(frozen=True)
class WorkflowPlan:
    """The complete, immutable interpretation of one merged run configuration."""

    mode: str
    output_root: Path
    model: ModelPlan
    reference: ReferenceMetadata
    generation_settings: GenerationSettings
    report: ReportPlan
    generation: GenerationPlan | None
    apply: ApplyPlan | None


def binding_data(
    binding: FileBinding | FastaBinding | None, *, required: bool = True
) -> str:
    """Return a binding's data path, preserving a clear missing-binding error."""
    if binding is None or binding.data is None:
        if required:
            raise WorkflowPlanError("Required binding has no data path")
        return ""
    return binding.data


def _reference_metadata(config: Mapping[str, Any]) -> ReferenceMetadata:
    reference = _mapping(config.get("reference"), "reference")
    exclusions = _mapping(reference.get("exclusions"), "reference.exclusions")
    exclusion_metadata = {
        identifier: ExclusionMetadata(
            identifier=identifier,
            kind=_required(item.get("kind"), f"reference.exclusions.{identifier}.kind"),
            role=_required(item.get("role"), f"reference.exclusions.{identifier}.role"),
            enabled=_boolean(
                item.get("enabled", True), f"reference.exclusions.{identifier}.enabled"
            ),
        )
        for identifier, value in exclusions.items()
        for item in (_mapping(value, f"reference.exclusions.{identifier}"),)
    }
    cadd = _mapping(reference.get("cadd"), "reference.cadd")
    indel = cadd.get("indel")
    partition_identities = _mapping(
        reference.get("testing_unit_partitions", {}),
        "reference.testing_unit_partitions",
    )
    return ReferenceMetadata(
        assembly=_required(reference.get("assembly"), "reference.assembly"),
        fasta_id=_required(
            _mapping(reference.get("fasta"), "reference.fasta").get("id"),
            "reference.fasta.id",
        ),
        exclusion_id=_required(
            _mapping(reference.get("exclusion"), "reference.exclusion").get("id"),
            "reference.exclusion.id",
        ),
        exclusion_kind=_required(
            _mapping(reference.get("exclusion"), "reference.exclusion").get("kind"),
            "reference.exclusion.kind",
        ),
        testing_unit_partitions=_frozen_mapping(
            {
                identifier: _required(
                    _mapping(value, f"reference.testing_unit_partitions.{identifier}").get(
                        "id"
                    ),
                    f"reference.testing_unit_partitions.{identifier}.id",
                )
                for identifier, value in partition_identities.items()
            }
        ),
        cadd_snv_id=_required(
            _mapping(cadd.get("snv"), "reference.cadd.snv").get("id"),
            "reference.cadd.snv.id",
        ),
        cadd_indel_id=(
            _required(
                _mapping(indel, "reference.cadd.indel").get("id"),
                "reference.cadd.indel.id",
            )
            if indel is not None
            else None
        ),
        exclusions=_frozen_mapping(exclusion_metadata),
    )


def _reference_bindings(
    value: Any, key: str, *, require_lookup: bool
) -> ReferenceBindings:
    reference = _mapping(value, key)
    bundle = (
        require_lookup and reference.get("lookup_source", {}).get("type") == "bundle"
    )
    cadd = _mapping(reference.get("cadd"), f"{key}.cadd")
    indel = cadd.get("indel")
    if require_lookup:
        _mapping(reference.get("lookup_source"), f"{key}.lookup_source")
    return ReferenceBindings(
        fasta=FastaBinding.from_config(reference.get("fasta"), f"{key}.fasta"),
        exclusion=FileBinding.from_config(
            reference.get("exclusion", {"data": None}),
            f"{key}.exclusion",
            required=False,
        ),
        testing_unit_partitions=(
            _testing_unit_partitions(
                reference.get("testing_unit_partitions"),
                f"{key}.testing_unit_partitions",
            )
            if not bundle
            else ()
        ),
        exclusions=_frozen_mapping(
            {
                identifier: FileBinding.from_config(
                    item, f"{key}.exclusions.{identifier}"
                )
                for identifier, value in _mapping(
                    reference.get("exclusions", {}), f"{key}.exclusions"
                ).items()
                for item in (_mapping(value, f"{key}.exclusions.{identifier}"),)
            }
        ),
        cadd_snv=FileBinding.from_config(cadd.get("snv"), f"{key}.cadd.snv"),
        cadd_indel=(
            FileBinding.from_config(indel, f"{key}.cadd.indel", required=False)
            if indel is not None
            else None
        ),
    )


def _resolved_reference_exclusions(
    metadata: ReferenceMetadata, bindings: ReferenceBindings
) -> tuple[ExclusionBinding, ...]:
    """Join active scientific exclusions to paths by stable identifier."""
    configured_ids = set(metadata.exclusions)
    binding_ids = set(bindings.exclusions)
    unknown_ids = sorted(binding_ids - configured_ids)
    if unknown_ids:
        raise WorkflowPlanError(
            "bindings.reference.exclusions has no scientific metadata for: "
            + ", ".join(unknown_ids)
            + "; declare each ID with kind and role in reference.exclusions, "
            "or remove its physical binding"
        )
    active_ids = {
        identifier
        for identifier, exclusion in metadata.exclusions.items()
        if exclusion.enabled
    }
    missing_ids = sorted(active_ids - binding_ids)
    if missing_ids:
        raise WorkflowPlanError(
            "bindings.reference.exclusions is missing active bindings for: "
            + ", ".join(missing_ids)
            + "; supply data/index under bindings.reference.exclusions "
            "or explicitly disable the corresponding scientific exclusion"
        )
    return tuple(
        ExclusionBinding(
            identifier=identifier,
            kind=metadata.exclusions[identifier].kind,
            role=metadata.exclusions[identifier].role,
            binding=bindings.exclusions[identifier],
        )
        for identifier in sorted(active_ids)
    )


def _lookup_source(
    value: Any,
    key: str,
    *,
    configured_reference: ReferenceBindings,
    model: ModelPlan,
    settings: GenerationSettings,
) -> LookupSource:
    source = _mapping(value, key)
    if source.get("type") != "bundle" or "directory" not in source:
        raise WorkflowPlanError(
            f"{key} only accepts {{type: bundle, directory: ...}} with "
            "partitioned_null_bundle_v1"
        )
    lookup = _directory_lookup_source(source, key)
    if model.hemizygote != any(model_id == "X-linked-hem" for _, model_id in lookup.tables):
        raise WorkflowPlanError(f"{key}.directory: bundle HEM model differs from model.hemizygote")
    return lookup


def _directory_lookup_source(source: Mapping[str, Any], key: str) -> LookupSource:
    """Resolve small, unverified declarations; the validation job certifies bytes."""
    directory = Path(_required(source["directory"], f"{key}.directory"))
    metadata_path = directory / BUNDLE_ARTIFACT_FILENAMES["metadata"]
    try:
        metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
        contract = str(metadata["contract"])
        if contract != PARTITIONED_NULL_BUNDLE_CONTRACT:
            raise ValueError("bundle contract must be partitioned_null_bundle_v1")
        grid = metadata["binning"]
        maximum, width = float(grid["score_max"]), float(grid["bin_width"])
        fixed_width_bin_count(maximum, width)
        components = metadata["components"]

        def component(name: str) -> str:
            filename = components[name]["path"]
            path = directory / filename
            if not path.is_file():
                raise ValueError(f"missing component {path}")
            return str(path)

        partitions = metadata["partitions"]
        partition_ids = tuple(entry["id"] for entry in partitions)
        if not partition_ids or len(partition_ids) != len(set(partition_ids)):
            raise ValueError("partitions must have unique IDs")
        for partition_id in partition_ids:
            require_partition_id(partition_id)
        if [entry["order"] for entry in partitions] != list(range(len(partition_ids))):
            raise ValueError("partitions must have contiguous configured order")
        tables = {}
        resolved, eligible = {}, {}
        for entry in partitions:
            partition_id = entry["id"]
            resolved[partition_id] = component(entry["resolved_component"])
            eligible[partition_id] = component(entry["eligible_component"])
            for model_id, component_name in entry["lookups"].items():
                tables[(partition_id, model_id)] = component(component_name)
        certificate = {}
        for name in ("validation", "marker"):
            path = directory / BUNDLE_ARTIFACT_FILENAMES[name]
            if not path.is_file():
                raise ValueError(f"missing certificate file {path}")
            certificate[name] = str(path)
        return LookupSource(
            source_type="bundle",
            contract=contract,
            tables=_frozen_mapping(tables),
            partition_ids=partition_ids,
            score_max=maximum,
            bin_width=width,
            metadata=str(metadata_path),
            **certificate,
            calibration_mask=component("calibration_regions"),
            excluded_alleles=(component("excluded_alleles"),),
            resolved_partitions=_frozen_mapping(resolved),
            eligible_partitions=_frozen_mapping(eligible),
            partition_statuses=component("unit_partition_status"),
        )
    except (OSError, ValueError, TypeError, KeyError, OverflowError) as exc:
        raise WorkflowPlanError(
            f"{key}.directory: cannot resolve bundle {directory}: {exc}. "
            "Point directory to a complete, validated generated bundle."
        ) from exc


@dataclass(frozen=True)
class GenerationResources:
    """Resolved resource scope and its configuration diagnostic keys."""

    resources: ReferenceBindings
    par_regions: FileBinding
    exclusions: tuple[ExclusionBinding, ...]
    snv_key: str
    indel_key: str


def _full_mode_resources(
    configured_reference: ReferenceBindings | None,
    reference_metadata: ReferenceMetadata,
) -> GenerationResources:
    if configured_reference is None:
        raise WorkflowPlanError("full mode requires bindings.reference")
    resolved_exclusions = _resolved_reference_exclusions(
        reference_metadata, configured_reference
    )
    par_exclusions = [
        exclusion
        for exclusion in resolved_exclusions
        if exclusion.role == "par_regions"
    ]
    if len(par_exclusions) != 1:
        raise WorkflowPlanError(
            "full mode requires exactly one bindings.reference exclusion with role par_regions"
        )
    resources = configured_reference
    par_regions = par_exclusions[0].binding
    exclusions = tuple(
        exclusion
        for exclusion in resolved_exclusions
        if exclusion.role == "variants_exclude"
    )
    return GenerationResources(
        resources,
        par_regions,
        exclusions,
        "bindings.reference.cadd.snv",
        "bindings.reference.cadd.indel",
    )


def _standalone_resources(generation: Mapping[str, Any]) -> GenerationResources:
    resources = ReferenceBindings(
        fasta=FastaBinding.from_config(generation.get("fasta"), "generation.fasta"),
        exclusion=FileBinding.from_config(
            generation.get("exclusion"), "generation.exclusion", required=False
        ),
        testing_unit_partitions=_testing_unit_partitions(
            generation.get("testing_unit_partitions"),
            "generation.testing_unit_partitions",
        ),
        exclusions=_frozen_mapping({}),
        cadd_snv=FileBinding.from_config(
            generation.get("score_snv"), "generation.score_snv"
        ),
        cadd_indel=(
            FileBinding.from_config(
                generation.get("score_indel"),
                "generation.score_indel",
                required=False,
            )
            if generation.get("score_indel") is not None
            else None
        ),
    )
    par_regions = FileBinding.from_config(
        generation.get("par_regions"), "generation.par_regions"
    )
    exclusions = tuple(
        ExclusionBinding(
            identifier=f"standalone-{index:03d}",
            kind="bed",
            role="variants_exclude",
            binding=FileBinding.from_config(item, f"generation.exclusions[{index}]"),
        )
        for index, item in enumerate(generation.get("exclusions", []))
    )
    return GenerationResources(
        resources,
        par_regions,
        exclusions,
        "generation.score_snv",
        "generation.score_indel",
    )


def _generation_plan(
    config: Mapping[str, Any],
    *,
    output_root: Path,
    mode: str,
    include_indels: bool,
    configured_reference: ReferenceBindings | None,
    reference_metadata: ReferenceMetadata,
) -> GenerationPlan:
    generation = _mapping(config.get("generation"), "generation")
    run = _mapping(config.get("run"), "run")
    prefix = _required(run.get("generation_prefix"), "run.generation_prefix")
    execution = _mapping(generation.get("execution", {}), "generation.execution")
    resolved = (
        _full_mode_resources(configured_reference, reference_metadata)
        if mode == "full"
        else _standalone_resources(generation)
    )
    resources = resolved.resources
    if not reference_metadata.testing_unit_partitions:
        raise WorkflowPlanError(
            "reference.testing_unit_partitions is required in generate and full modes"
        )
    _validate_partition_identity_set(
        reference_metadata.testing_unit_partitions,
        resources.testing_unit_partitions,
        identity_key="reference.testing_unit_partitions",
        binding_key=(
            "bindings.reference.testing_unit_partitions"
            if mode == "full"
            else "generation.testing_unit_partitions"
        ),
    )
    par_regions = resolved.par_regions
    exclusions = resolved.exclusions
    af_source = _frozen_mapping(
        {
            name: _required(value, f"generation.af_source.{name}")
            for name, value in _mapping(
                generation.get("af_source"), "generation.af_source"
            ).items()
            if value is not None
        }
    )
    if "overall" not in af_source:
        raise WorkflowPlanError("generation.af_source.overall is required")
    population_vcf = FileBinding.from_config(
        generation.get("population_vcf"), "generation.population_vcf"
    )
    mutation_rate = (
        FileBinding.from_config(generation["mutation_rate"], "generation.mutation_rate")
        if generation.get("mutation_rate") is not None
        else None
    )
    mutation_rate_indel = (
        FileBinding.from_config(
            generation["mutation_rate_indel"], "generation.mutation_rate_indel"
        )
        if generation.get("mutation_rate_indel") is not None
        else None
    )
    if include_indels and resources.cadd_indel is None:
        key = resolved.indel_key
        raise WorkflowPlanError(f"Indel generation requires {key}")
    if include_indels and reference_metadata.cadd_indel_id is None:
        raise WorkflowPlanError("Indel generation requires reference.cadd.indel.id")
    allow_unindexed_fixtures = execution.get("allow_unindexed_fixture_inputs", False)
    if not isinstance(allow_unindexed_fixtures, bool):
        raise WorkflowPlanError(
            "generation.execution.allow_unindexed_fixture_inputs must be a boolean"
        )
    window_resources = [
        ("generation.population_vcf", population_vcf),
        (resolved.snv_key, resources.cadd_snv),
    ]
    if include_indels:
        assert resources.cadd_indel is not None
        window_resources.append(
            (
                resolved.indel_key,
                resources.cadd_indel,
            )
        )
    if mutation_rate is not None:
        window_resources.append(("generation.mutation_rate", mutation_rate))
    if mutation_rate_indel is not None:
        window_resources.append(("generation.mutation_rate_indel", mutation_rate_indel))
    for key, resource in window_resources:
        try:
            check_generation_index_binding(
                resource.data,
                resource.index,
                key=key,
                allow_unindexed_fixtures=allow_unindexed_fixtures,
            )
        except ValueError as exc:
            raise WorkflowPlanError(str(exc)) from exc
    return GenerationPlan(
        output_dir=output_root / "generate" / prefix,
        log_dir=output_root / "logs" / "generate" / prefix,
        population_vcf=population_vcf,
        score_snv=resources.cadd_snv,
        score_indel=resources.cadd_indel,
        exclusion=resources.exclusion,
        testing_unit_partitions=resources.testing_unit_partitions,
        fasta=resources.fasta,
        par_regions=par_regions,
        exclusions=exclusions,
        mutation_rate=mutation_rate,
        mutation_rate_indel=mutation_rate_indel,
        unobserved_policy=_required(
            generation.get("unobserved_policy"), "generation.unobserved_policy"
        ),
        af_source=af_source,
        window_size_bp=_number(
            execution.get("partition_size_bp", 5_000_000),
            "generation.execution.partition_size_bp",
            int,
        ),
        unit_bucket_size=_number(
            execution.get("unit_bucket_size", 250),
            "generation.execution.unit_bucket_size",
            int,
        ),
        parquet_compression=str(execution.get("parquet_compression", "zstd")),
    )


def _generated_lookup_source(
    generation: GenerationPlan, model: ModelPlan, settings: GenerationSettings
) -> LookupSource:
    bundle_paths = generated_bundle_paths(generation.bundle_dir)
    partition_ids = tuple(partition.identifier for partition in generation.testing_unit_partitions)
    tables = bundle_component_paths(generation.bundle_dir, partition_ids, model.hemizygote)
    return LookupSource(
        source_type="bundle",
        contract=PARTITIONED_NULL_BUNDLE_CONTRACT,
        tables=tables,
        partition_ids=partition_ids,
        score_max=model.score_max,
        bin_width=settings.score_bin_width,
        metadata=bundle_paths["metadata"],
        validation=bundle_paths["validation"],
        marker=bundle_paths["marker"],
        calibration_mask=bundle_paths["calibration_regions"],
        excluded_alleles=(bundle_paths["excluded_alleles"],),
        resolved_partitions=_frozen_mapping({
            partition_id: str(generation.bundle_dir / f"resolved.{partition_id}.bed.gz")
            for partition_id in partition_ids
        }),
        eligible_partitions=_frozen_mapping({
            partition_id: str(generation.bundle_dir / f"eligible.{partition_id}.bed.gz")
            for partition_id in partition_ids
        }),
        partition_statuses=bundle_paths["unit_partition_status"],
    )


def _effective_bundle_reference(
    configured_reference: ReferenceBindings, lookup_source: LookupSource
) -> tuple[ReferenceBindings, tuple[ExclusionBinding, ...]]:
    bundle_exclusions = tuple(
        ExclusionBinding(
            identifier=(
                "generated-bundle-excluded-alleles"
                if index == 0
                else f"generated-bundle-excluded-alleles-{index:03d}"
            ),
            kind="tsv",
            role="variants_exclude",
            binding=FileBinding(path),
        )
        for index, path in enumerate(lookup_source.excluded_alleles)
    )
    return (
        ReferenceBindings(
            fasta=configured_reference.fasta,
            exclusion=FileBinding(None),
            calibration_regions=FileBinding(lookup_source.calibration_mask),
            testing_unit_partitions=(),
            exclusions=_frozen_mapping(
                {
                    exclusion.identifier: exclusion.binding
                    for exclusion in bundle_exclusions
                }
            ),
            cadd_snv=configured_reference.cadd_snv,
            cadd_indel=configured_reference.cadd_indel,
        ),
        bundle_exclusions,
    )


def _model_plan(config: Mapping[str, Any]) -> ModelPlan:
    model_config = _mapping(config.get("model"), "model")
    return ModelPlan(
        unit=_required(model_config.get("unit"), "model.unit"),
        hemizygote=_boolean(model_config.get("hemizygote"), "model.hemizygote"),
        include_indels=_boolean(
            model_config.get("include_indels"), "model.include_indels"
        ),
        score_max=_number(model_config.get("score_max"), "model.score_max", float),
    )


def _generation_settings(
    config: Mapping[str, Any], mode: str, model: ModelPlan
) -> GenerationSettings:
    generation_settings_config = _mapping(
        config.get("generation_settings"), "generation_settings"
    )
    x_chromosome = _mapping(
        generation_settings_config.get("x_chromosome"),
        "generation_settings.x_chromosome",
    )
    indels = _mapping(
        generation_settings_config.get("indels", {}), "generation_settings.indels"
    )
    settings = GenerationSettings(
        minimum_retained_fraction=_number(
            generation_settings_config.get("minimum_retained_fraction"),
            "generation_settings.minimum_retained_fraction",
            float,
        ),
        score_bin_width=_number(
            generation_settings_config.get("score_bin_width"),
            "generation_settings.score_bin_width",
            float,
        ),
        mutation_conversion=_required(
            generation_settings_config.get("mutation_conversion"),
            "generation_settings.mutation_conversion",
        ),
        variant_classes=tuple(generation_settings_config.get("variant_classes", [])),
        indel_candidate_universe=str(
            indels.get("candidate_universe", "score_resource")
        ),
        indel_membership=str(indels.get("membership", "reference_span_overlap")),
        indel_mask_policy=str(indels.get("mask_policy", "reference_span_contained")),
        indel_max_reference_length=_number(
            indels.get("max_reference_length", 50),
            "generation_settings.indels.max_reference_length",
            int,
        ),
        indel_max_alternate_length=_number(
            indels.get("max_alternate_length", 50),
            "generation_settings.indels.max_alternate_length",
            int,
        ),
        x_af_policy=_required(
            x_chromosome.get("af_policy"), "generation_settings.x_chromosome.af_policy"
        ),
        x_mutation_rate_policy=_required(
            x_chromosome.get("mutation_rate_policy"),
            "generation_settings.x_chromosome.mutation_rate_policy",
        ),
    )
    if mode in {"generate", "full"}:
        expected_classes = ("snv", "indel") if model.include_indels else ("snv",)
        if settings.variant_classes != expected_classes:
            raise WorkflowPlanError(
                "generation_settings.variant_classes must be "
                f"{list(expected_classes)!r} when model.include_indels is "
                f"{model.include_indels!r}"
            )
    if (
        settings.indel_candidate_universe != "score_resource"
        or settings.indel_membership != "reference_span_overlap"
        or settings.indel_mask_policy != "reference_span_contained"
        or settings.indel_max_reference_length < 1
        or settings.indel_max_alternate_length < 1
    ):
        raise WorkflowPlanError("generation_settings.indels has an unsupported policy")
    return settings


def _report_plan(config: Mapping[str, Any]) -> ReportPlan:
    report_config = _mapping(config.get("report"), "report")
    candidate_id = _mapping(report_config.get("candidate_id"), "report.candidate_id")
    return ReportPlan(
        candidate_id_separator=_required(
            candidate_id.get("separator"), "report.candidate_id.separator"
        ),
        candidate_id_ordinal_width=_number(
            candidate_id.get("ordinal_width"), "report.candidate_id.ordinal_width", int
        ),
        deduplicate_controls=bool(report_config.get("deduplicate_controls", True)),
    )


def _apply_plan(
    config: Mapping[str, Any],
    *,
    output_root: Path,
    mode: str,
    configured_reference: ReferenceBindings | None,
    generation: GenerationPlan | None,
    reference_metadata: ReferenceMetadata,
    model: ModelPlan,
    settings: GenerationSettings,
) -> ApplyPlan:
    bindings = _mapping(config.get("bindings"), "bindings")
    if configured_reference is None:
        raise WorkflowPlanError(f"{mode} mode requires bindings.reference")
    run = _mapping(config.get("run"), "run")
    prefix = _required(run.get("apply_prefix"), "run.apply_prefix")
    if mode == "full":
        if generation is None:
            raise WorkflowPlanError("full mode requires generation inputs")
        lookup_source = _generated_lookup_source(generation, model, settings)
        effective_reference, effective_exclusions = _effective_bundle_reference(
            configured_reference, lookup_source
        )
    else:
        source_config = _mapping(bindings.get("reference"), "bindings.reference").get(
            "lookup_source"
        )
        lookup_source = _lookup_source(
            source_config,
            "bindings.reference.lookup_source",
            configured_reference=configured_reference,
            model=model,
            settings=settings,
        )
        if lookup_source.is_bundle:
            effective_reference, effective_exclusions = _effective_bundle_reference(
                configured_reference, lookup_source
            )
        else:
            effective_reference = configured_reference
            effective_exclusions = _resolved_reference_exclusions(
                reference_metadata, configured_reference
            )
    inputs_config = _mapping(config.get("inputs"), "inputs")
    return ApplyPlan(
        output_dir=output_root / "apply" / prefix,
        log_dir=output_root / "logs" / "apply" / prefix,
        inputs=InputBindings(
            vcf=FileBinding.from_config(inputs_config.get("vcf"), "inputs.vcf"),
            ped=_optional(inputs_config.get("ped")),
            samples=_optional(inputs_config.get("samples")),
        ),
        configured_reference=configured_reference,
        effective_reference=effective_reference,
        effective_exclusions=effective_exclusions,
        lookup_source=lookup_source,
    )


def build_workflow_plan(config: Mapping[str, Any]) -> WorkflowPlan:
    """Resolve validated configuration into immutable mode-specific records."""
    normalized = normalize_run_config(config)
    _reject_legacy_exclusion_lists(normalized)
    try:
        validate_configuration(normalized)
    except ConfigurationValidationError as exc:
        raise WorkflowPlanError(str(exc)) from exc
    mode = _required(
        _mapping(normalized.get("workflow"), "workflow").get("mode"), "workflow.mode"
    )
    if mode not in {"apply", "generate", "full"}:
        raise WorkflowPlanError(f"Unsupported workflow.mode: {mode!r}")
    reference_config = _mapping(normalized.get("reference"), "reference")
    if mode == "apply" and "testing_unit_partitions" in reference_config:
        raise WorkflowPlanError(
            "apply mode forbids reference.testing_unit_partitions; use the certified bundle"
        )
    if mode == "apply":
        generation_config = normalized.get("generation")
        bindings_config = normalized.get("bindings")
        bound_reference = (
            bindings_config.get("reference", {})
            if isinstance(bindings_config, Mapping)
            else {}
        )
        if (
            isinstance(generation_config, Mapping)
            and "testing_unit_partitions" in generation_config
        ) or (
            isinstance(bound_reference, Mapping)
            and "testing_unit_partitions" in bound_reference
        ):
            raise WorkflowPlanError(
                "apply mode forbids physical testing-unit partition bindings"
            )
    output_root = Path(_required(normalized.get("outdir"), "outdir"))
    model = _model_plan(normalized)
    settings = _generation_settings(normalized, mode, model)
    reference_metadata = _reference_metadata(normalized)
    report = _report_plan(normalized)

    configured_reference: ReferenceBindings | None = None
    generation: GenerationPlan | None = None
    apply: ApplyPlan | None = None
    if mode in {"apply", "full"}:
        bindings = _mapping(normalized.get("bindings"), "bindings")
        configured_reference = _reference_bindings(
            bindings.get("reference"),
            "bindings.reference",
            require_lookup=mode == "apply",
        )
    if mode in {"generate", "full"}:
        generation = _generation_plan(
            normalized,
            output_root=output_root,
            mode=mode,
            include_indels=model.include_indels,
            configured_reference=configured_reference,
            reference_metadata=reference_metadata,
        )
    if mode in {"apply", "full"}:
        apply = _apply_plan(
            normalized,
            output_root=output_root,
            mode=mode,
            configured_reference=configured_reference,
            generation=generation,
            reference_metadata=reference_metadata,
            model=model,
            settings=settings,
        )
    return WorkflowPlan(
        mode=mode,
        output_root=output_root,
        model=model,
        reference=reference_metadata,
        generation_settings=settings,
        report=report,
        generation=generation,
        apply=apply,
    )
