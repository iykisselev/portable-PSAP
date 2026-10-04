"""Immutable resolution of Portable PSAP workflow configuration.

Snakemake parses rule modules in one shared Python namespace.  Resolving the
run mode and physical resource bindings here keeps that parsing layer from
rewriting configuration or maintaining parallel, partially-effective binding
dictionaries.
"""

from __future__ import annotations

from copy import deepcopy
from dataclasses import dataclass
from pathlib import Path
from types import MappingProxyType
from typing import Any, Mapping
import warnings

from generation_resource_policy import check_generation_index_binding
from config_validation import ConfigurationValidationError, validate_configuration
from bundle_contract import (
    BUNDLE_COMPONENT_FILENAMES,
    CURRENT_GENERATED_BUNDLE_CONTRACT as GENERATED_BUNDLE_CONTRACT,
)


class WorkflowPlanError(ValueError):
    """Raised when a valid-looking configuration cannot form one workflow plan."""


def normalize_run_config(config: Mapping[str, Any]) -> dict[str, Any]:
    """Copy a legacy manifest into the explicit apply-mode configuration shape."""
    normalized = deepcopy(dict(config))
    if "workflow" not in normalized:
        warnings.warn(
            "Run manifests without workflow.mode are deprecated and currently "
            "default to apply mode; add workflow: {mode: apply} before workflow "
            "configuration contract v2, when workflow.mode becomes mandatory.",
            DeprecationWarning,
            stacklevel=2,
        )
        normalized["workflow"] = {"mode": "apply"}
    return normalized


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
    if value is None:
        raise WorkflowPlanError(f"Required configuration value is missing: {key}")
    path = str(value).strip()
    if not path or path.lower() in {"na", "none", "null", ".", "-"}:
        raise WorkflowPlanError(f"Required configuration value is missing: {key}")
    return str(Path(path).expanduser())


def _optional(value: Any) -> str | None:
    if value is None:
        return None
    path = str(value).strip()
    if path.lower() in {"", "na", "none", "null", ".", "-"}:
        return None
    return str(Path(path).expanduser())


def _number(value: Any, key: str, converter: type[float] | type[int]) -> float | int:
    try:
        return converter(value)
    except (TypeError, ValueError) as exc:
        raise WorkflowPlanError(f"{key} must be a {converter.__name__}") from exc


def _boolean(value: Any, key: str) -> bool:
    """Require a real configuration boolean rather than coercing truthy values."""
    if not isinstance(value, bool):
        raise WorkflowPlanError(f"{key} must be a boolean")
    return value


def _frozen_mapping(values: Mapping[str, str]) -> Mapping[str, str]:
    return MappingProxyType(dict(values))


def _frozen_binding_mapping(
    values: Mapping[str, FileBinding]
) -> Mapping[str, FileBinding]:
    return MappingProxyType(dict(values))


def generated_bundle_paths(bundle_dir: Path) -> Mapping[str, str]:
    """Return the public generated-bundle component paths for one output root."""
    return _frozen_mapping(
        {
            component: str(bundle_dir / filename)
            for component, filename in BUNDLE_COMPONENT_FILENAMES.items()
        }
    )


@dataclass(frozen=True)
class FileBinding:
    """One data file and its optional declared index."""

    data: str | None
    index: str | None = None

    @classmethod
    def from_config(
        cls, value: Any, key: str, *, required: bool = True
    ) -> "FileBinding":
        mapping = _mapping(value, key)
        data = _required(mapping.get("data"), f"{key}.data") if required else _optional(mapping.get("data"))
        return cls(data=data, index=_optional(mapping.get("index")))


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
    testing_unit_id: str
    cadd_snv_id: str
    cadd_indel_id: str | None
    exclusions: Mapping[str, ExclusionMetadata]


@dataclass(frozen=True)
class ReferenceBindings:
    """Immutable physical bindings used to execute one reference scope."""

    fasta: FastaBinding
    exclusion: FileBinding
    testing_unit_bed: FileBinding
    exclusions: Mapping[str, FileBinding]
    cadd_snv: FileBinding
    cadd_indel: FileBinding | None
    calibration_regions: FileBinding | None = None


@dataclass(frozen=True)
class LookupSource:
    """Resolved lookup tables and optional generated-bundle certificate files."""

    source_type: str
    contract: str
    tables: Mapping[str, str]
    score_max: float
    bin_width: float | None
    metadata: str | None = None
    validation: str | None = None
    marker: str | None = None
    eligible_units: str | None = None
    calibration_mask: str | None = None
    excluded_alleles: tuple[str, ...] = ()
    excluded_units: str | None = None

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
    x_af_policy: str
    x_mutation_rate_policy: str


@dataclass(frozen=True)
class GenerationPlan:
    """Generation inputs, policy, and paths for generate/full modes."""

    output_dir: Path
    log_dir: Path
    population_vcf: FileBinding
    score_snv: FileBinding
    exclusion: FileBinding
    testing_unit_bed: FileBinding
    fasta: FastaBinding
    par_regions: FileBinding
    exclusions: tuple[ExclusionBinding, ...]
    mutation_rate: FileBinding | None
    unobserved_policy: str
    af_source: Mapping[str, str]
    partition_size_bp: int
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


def binding_data(binding: FileBinding | FastaBinding | None, *, required: bool = True) -> str:
    """Return a binding's data path, preserving a clear missing-binding error."""
    if binding is None or binding.data is None:
        if required:
            raise WorkflowPlanError("Required binding has no data path")
        return ""
    return binding.data


def binding_index(binding: FileBinding | None) -> str:
    """Return only the index declared by the run manifest."""
    return "" if binding is None or binding.index is None else binding.index


def _reference_metadata(config: Mapping[str, Any]) -> ReferenceMetadata:
    reference = _mapping(config.get("reference"), "reference")
    exclusions = _mapping(reference.get("exclusions"), "reference.exclusions")
    exclusion_metadata = {
        identifier: ExclusionMetadata(
            identifier=identifier,
            kind=_required(item.get("kind"), f"reference.exclusions.{identifier}.kind"),
            role=_required(item.get("role"), f"reference.exclusions.{identifier}.role"),
            enabled=_boolean(item.get("enabled", True), f"reference.exclusions.{identifier}.enabled"),
        )
        for identifier, value in exclusions.items()
        for item in (_mapping(value, f"reference.exclusions.{identifier}"),)
    }
    cadd = _mapping(reference.get("cadd"), "reference.cadd")
    indel = cadd.get("indel")
    return ReferenceMetadata(
        assembly=_required(reference.get("assembly"), "reference.assembly"),
        fasta_id=_required(_mapping(reference.get("fasta"), "reference.fasta").get("id"), "reference.fasta.id"),
        exclusion_id=_required(_mapping(reference.get("exclusion"), "reference.exclusion").get("id"), "reference.exclusion.id"),
        exclusion_kind=_required(_mapping(reference.get("exclusion"), "reference.exclusion").get("kind"), "reference.exclusion.kind"),
        testing_unit_id=_required(_mapping(reference.get("testing_unit_bed"), "reference.testing_unit_bed").get("id"), "reference.testing_unit_bed.id"),
        cadd_snv_id=_required(_mapping(cadd.get("snv"), "reference.cadd.snv").get("id"), "reference.cadd.snv.id"),
        cadd_indel_id=(
            _required(_mapping(indel, "reference.cadd.indel").get("id"), "reference.cadd.indel.id")
            if indel is not None
            else None
        ),
        exclusions=_frozen_mapping(exclusion_metadata),
    )


def _reference_bindings(value: Any, key: str, *, require_lookup: bool) -> ReferenceBindings:
    reference = _mapping(value, key)
    cadd = _mapping(reference.get("cadd"), f"{key}.cadd")
    indel = cadd.get("indel")
    if require_lookup:
        _mapping(reference.get("lookup_source"), f"{key}.lookup_source")
    return ReferenceBindings(
        fasta=FastaBinding.from_config(reference.get("fasta"), f"{key}.fasta"),
        exclusion=FileBinding.from_config(
            reference.get("exclusion"), f"{key}.exclusion", required=False
        ),
        testing_unit_bed=FileBinding.from_config(
            reference.get("testing_unit_bed"), f"{key}.testing_unit_bed"
        ),
        exclusions=_frozen_binding_mapping(
            {
                identifier: FileBinding.from_config(item, f"{key}.exclusions.{identifier}")
                for identifier, value in _mapping(
                    reference.get("exclusions"), f"{key}.exclusions"
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
    source_type = _required(source.get("type"), f"{key}.type")
    contract = _required(source.get("contract"), f"{key}.contract")
    tables = _frozen_mapping(
        {
            name: _required(path, f"{key}.tables.{name}")
            for name, path in _mapping(source.get("tables"), f"{key}.tables").items()
        }
    )
    for role in ("het", "hom"):
        if role not in tables:
            raise WorkflowPlanError(f"{key}.tables.{role} is required")
    if model.hemizygote and "hem_chrX" not in tables:
        raise WorkflowPlanError(
            f"{key}.tables.hem_chrX is required when model.hemizygote is true"
        )

    is_bundle = source_type == "bundle"
    metadata = _optional(source.get("metadata"))
    validation = _optional(source.get("validation"))
    marker = _optional(source.get("marker"))
    eligible_units = (
        _optional(source.get("eligible_units"))
        if is_bundle
        else binding_data(configured_reference.testing_unit_bed)
    )
    calibration_mask = (
        _optional(source.get("calibration_mask"))
        if is_bundle
        else None
    )
    excluded_alleles = (
        (_required(source.get("excluded_alleles"), f"{key}.excluded_alleles"),)
        if is_bundle
        else ()
    )
    excluded_units = _optional(source.get("excluded_units")) if is_bundle else None
    if is_bundle and not all(
        (metadata, validation, marker, eligible_units, calibration_mask, excluded_units)
    ):
        raise WorkflowPlanError(
            f"{key} bundle sources require metadata, validation, marker, eligible_units, "
            "calibration_mask, and excluded_units"
        )
    return LookupSource(
        source_type=source_type,
        contract=contract,
        tables=tables,
        score_max=model.score_max,
        bin_width=settings.score_bin_width
        if contract in {"generated_bundle_v1", GENERATED_BUNDLE_CONTRACT}
        else None,
        metadata=metadata,
        validation=validation,
        marker=marker,
        eligible_units=eligible_units,
        calibration_mask=calibration_mask,
        excluded_alleles=excluded_alleles,
        excluded_units=excluded_units,
    )


def _generation_plan(
    config: Mapping[str, Any],
    *,
    output_root: Path,
    mode: str,
    configured_reference: ReferenceBindings | None,
    reference_metadata: ReferenceMetadata,
) -> GenerationPlan:
    generation = _mapping(config.get("generation"), "generation")
    run = _mapping(config.get("run"), "run")
    prefix = _required(run.get("generation_prefix"), "run.generation_prefix")
    execution = _mapping(generation.get("execution", {}), "generation.execution")
    if mode == "full":
        if configured_reference is None:
            raise WorkflowPlanError("full mode requires bindings.reference")
        resolved_exclusions = _resolved_reference_exclusions(
            reference_metadata, configured_reference
        )
        par_exclusions = [
            exclusion for exclusion in resolved_exclusions if exclusion.role == "par_regions"
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
    else:
        resources = ReferenceBindings(
            fasta=FastaBinding.from_config(generation.get("fasta"), "generation.fasta"),
            exclusion=FileBinding.from_config(
                generation.get("exclusion"), "generation.exclusion", required=False
            ),
            testing_unit_bed=FileBinding.from_config(
                generation.get("testing_unit_bed"), "generation.testing_unit_bed"
            ),
            exclusions=_frozen_binding_mapping({}),
            cadd_snv=FileBinding.from_config(generation.get("score_snv"), "generation.score_snv"),
            cadd_indel=None,
        )
        par_regions = FileBinding.from_config(generation.get("par_regions"), "generation.par_regions")
        exclusions = tuple(
            ExclusionBinding(
                identifier=f"standalone-{index:03d}",
                kind="bed",
                role="variants_exclude",
                binding=FileBinding.from_config(
                    item, f"generation.exclusions[{index}]"
                ),
            )
            for index, item in enumerate(generation.get("exclusions", []))
        )
    af_source = _frozen_mapping(
        {
            name: _required(value, f"generation.af_source.{name}")
            for name, value in _mapping(generation.get("af_source"), "generation.af_source").items()
            if value is not None
        }
    )
    if "overall" not in af_source:
        raise WorkflowPlanError("generation.af_source.overall is required")
    population_vcf = FileBinding.from_config(generation.get("population_vcf"), "generation.population_vcf")
    mutation_rate = (
        FileBinding.from_config(generation["mutation_rate"], "generation.mutation_rate")
        if generation.get("mutation_rate") is not None else None
    )
    allow_unindexed_fixtures = execution.get("allow_unindexed_fixture_inputs", False)
    if not isinstance(allow_unindexed_fixtures, bool):
        raise WorkflowPlanError("generation.execution.allow_unindexed_fixture_inputs must be a boolean")
    partition_resources = [
        ("generation.population_vcf", population_vcf),
        ("bindings.reference.cadd.snv" if mode == "full" else "generation.score_snv", resources.cadd_snv),
    ]
    if mutation_rate is not None:
        partition_resources.append(("generation.mutation_rate", mutation_rate))
    for key, resource in partition_resources:
        try:
            check_generation_index_binding(
                resource.data, resource.index, key=key, allow_unindexed_fixtures=allow_unindexed_fixtures,
            )
        except ValueError as exc:
            raise WorkflowPlanError(str(exc)) from exc
    return GenerationPlan(
        output_dir=output_root / "generate" / prefix,
        log_dir=output_root / "logs" / "generate" / prefix,
        population_vcf=population_vcf,
        score_snv=resources.cadd_snv,
        exclusion=resources.exclusion,
        testing_unit_bed=resources.testing_unit_bed,
        fasta=resources.fasta,
        par_regions=par_regions,
        exclusions=exclusions,
        mutation_rate=mutation_rate,
        unobserved_policy=_required(
            generation.get("unobserved_policy"), "generation.unobserved_policy"
        ),
        af_source=af_source,
        partition_size_bp=_number(
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
    tables = {
        "het": bundle_paths["lookup_het"],
        "hom": bundle_paths["lookup_hom"],
    }
    if model.hemizygote:
        tables["hem_chrX"] = bundle_paths["lookup_hem"]
    return LookupSource(
        source_type="bundle",
        contract=GENERATED_BUNDLE_CONTRACT,
        tables=_frozen_mapping(tables),
        score_max=model.score_max,
        bin_width=settings.score_bin_width,
        metadata=bundle_paths["metadata"],
        validation=bundle_paths["validation"],
        marker=bundle_paths["marker"],
        eligible_units=bundle_paths["eligible_bed"],
        calibration_mask=bundle_paths["calibration_regions"],
        excluded_alleles=(bundle_paths["excluded_alleles"],),
        excluded_units=bundle_paths["excluded_units"],
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
            testing_unit_bed=FileBinding(lookup_source.eligible_units),
            exclusions=_frozen_binding_mapping(
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


def build_workflow_plan(config: Mapping[str, Any]) -> WorkflowPlan:
    """Resolve validated configuration into immutable mode-specific records."""
    normalized = normalize_run_config(config)
    _reject_legacy_exclusion_lists(normalized)
    try:
        validate_configuration(normalized)
    except ConfigurationValidationError as exc:
        raise WorkflowPlanError(str(exc)) from exc
    mode = _required(_mapping(normalized.get("workflow"), "workflow").get("mode"), "workflow.mode")
    if mode not in {"apply", "generate", "full"}:
        raise WorkflowPlanError(f"Unsupported workflow.mode: {mode!r}")
    output_root = Path(_required(normalized.get("outdir"), "outdir"))
    model_config = _mapping(normalized.get("model"), "model")
    model = ModelPlan(
        unit=_required(model_config.get("unit"), "model.unit"),
        hemizygote=_boolean(model_config.get("hemizygote"), "model.hemizygote"),
        include_indels=_boolean(
            model_config.get("include_indels"), "model.include_indels"
        ),
        score_max=_number(model_config.get("score_max"), "model.score_max", float),
    )
    generation_settings_config = _mapping(
        normalized.get("generation_settings"), "generation_settings"
    )
    x_chromosome = _mapping(
        generation_settings_config.get("x_chromosome"), "generation_settings.x_chromosome"
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
        x_af_policy=_required(x_chromosome.get("af_policy"), "generation_settings.x_chromosome.af_policy"),
        x_mutation_rate_policy=_required(
            x_chromosome.get("mutation_rate_policy"),
            "generation_settings.x_chromosome.mutation_rate_policy",
        ),
    )
    reference_metadata = _reference_metadata(normalized)
    report_config = _mapping(normalized.get("report"), "report")
    candidate_id = _mapping(report_config.get("candidate_id"), "report.candidate_id")
    report = ReportPlan(
        candidate_id_separator=_required(candidate_id.get("separator"), "report.candidate_id.separator"),
        candidate_id_ordinal_width=_number(
            candidate_id.get("ordinal_width"), "report.candidate_id.ordinal_width", int
        ),
        deduplicate_controls=bool(report_config.get("deduplicate_controls", True)),
    )

    configured_reference: ReferenceBindings | None = None
    generation: GenerationPlan | None = None
    apply: ApplyPlan | None = None
    if mode in {"apply", "full"}:
        bindings = _mapping(normalized.get("bindings"), "bindings")
        configured_reference = _reference_bindings(
            bindings.get("reference"), "bindings.reference", require_lookup=mode == "apply"
        )
    if mode in {"generate", "full"}:
        generation = _generation_plan(
            normalized,
            output_root=output_root,
            mode=mode,
            configured_reference=configured_reference,
            reference_metadata=reference_metadata,
        )
    if mode in {"apply", "full"}:
        if configured_reference is None:
            raise WorkflowPlanError(f"{mode} mode requires bindings.reference")
        run = _mapping(normalized.get("run"), "run")
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
        inputs_config = _mapping(normalized.get("inputs"), "inputs")
        apply = ApplyPlan(
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
