"""Describe the files used to validate, stage, and record one apply run.

One resource catalog is built from the effective ``WorkflowPlan``. That keeps
validation and staging pointed at the same file for each scientific resource.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import re
from types import MappingProxyType
from typing import Mapping

from workflow_plan import (
    FileBinding,
    FastaBinding,
    WorkflowPlan,
    binding_data,
)
from workflow_bindings import as_inputs


@dataclass(frozen=True)
class ResourceSpec:
    """One physical or derived scientific resource with its staging contract."""

    logical_name: str
    source: str
    source_index: str
    kind: str
    scientific_role: str
    requires_chromosome_normalization: bool
    staged_data: str
    staged_index: str
    tsv_layout: str | None
    origin: str
    companion_indexes: tuple[str, ...] = ()
    stage: bool = True


@dataclass(frozen=True)
class ApplyResourceRegistry:
    """Immutable catalog and stable paths for one apply/full workflow run."""

    resources: Mapping[str, ResourceSpec]
    prepared_dir: Path

    def resource(self, logical_name: str) -> ResourceSpec:
        """Return one named resource with a clear error for stale consumers."""
        try:
            return self.resources[logical_name]
        except KeyError as exc:
            raise KeyError(f"Unknown apply resource: {logical_name!r}") from exc

    def inputs(self, name: str, field: str = "source") -> list[str]:
        """Declare one enabled resource field, preserving absent optionals."""
        resource = self.resources.get(name)
        return as_inputs(getattr(resource, field)) if resource else []

    def exclusion_specs(self) -> list[str]:
        return [
            f"{self.resource(name).kind}:{self.resource(name).staged_data}"
            for name in self.exclusion_ids
        ]

    def validation_fields(self) -> dict[str, dict[str, str]]:
        """Serialize resource identities as references to declared input fields."""
        roles = {
            "calibration_regions",
            "variants_exclude",
            "par_regions",
            "cadd_snv",
            "cadd_indel",
        }
        return {
            name: {
                "data": f"resource_{name}",
                "index": f"resource_{name}_index",
                "role": resource.scientific_role,
            }
            for name, resource in self.resources.items()
            if resource.scientific_role in roles
        }

    def validation_inputs(self) -> dict[str, list[str]]:
        result = {}
        for name, fields in self.validation_fields().items():
            result[fields["data"]] = self.inputs(name)
            result[fields["index"]] = self.inputs(name, "source_index")
        return result

    @property
    def prepared_table_dir(self) -> Path:
        return self.prepared_dir / "tables"

    @property
    def table_resources(self) -> Mapping[str, ResourceSpec]:
        return MappingProxyType(
            {
                name: resource_spec
                for name, resource_spec in self.resources.items()
                if resource_spec.stage and resource_spec.kind in {"bed", "tsv"}
            }
        )

    @property
    def table_resource_ids(self) -> tuple[str, ...]:
        return tuple(self.table_resources)

    @property
    def table_resource_pattern(self) -> str:
        return f"(?:{'|'.join(map(re.escape, self.table_resource_ids))})"

    @property
    def exclusion_ids(self) -> tuple[str, ...]:
        return tuple(
            name
            for name, resource_spec in self.table_resources.items()
            if resource_spec.scientific_role in {"variants_exclude", "par_regions"}
        )

    @property
    def workflow_exclusions(self) -> tuple[str, ...]:
        return tuple(self.resource(name).staged_data for name in self.exclusion_ids)

    @property
    def exclusion_sources(self) -> tuple[str, ...]:
        return tuple(self.resource(name).source for name in self.exclusion_ids)

    @property
    def combined_exclusions(self) -> str:
        return self.resource("combined_exclusions").staged_data

    @property
    def combined_exclusions_index(self) -> str:
        return self.resource("combined_exclusions").staged_index


def build_apply_resource_registry(
    plan: WorkflowPlan, prepared_dir: Path
) -> ApplyResourceRegistry:
    """Build the one resource catalog used by apply validation and staging."""
    if plan.apply is None:
        raise ValueError(
            "Apply resource registry requires an apply or full WorkflowPlan"
        )

    apply_plan = plan.apply
    effective_reference = apply_plan.effective_reference
    resource_origin = (
        "generated_bundle" if apply_plan.lookup_source.is_bundle else "configured"
    )
    resource_specs: dict[str, ResourceSpec] = {}
    prepared_table_dir = prepared_dir / "tables"

    def add_table_resource(
        logical_name: str,
        resource_binding: FileBinding | None,
        kind: str,
        scientific_role: str,
        *,
        requires_chromosome_normalization: bool = False,
        tsv_layout: str | None = None,
        origin: str = resource_origin,
    ) -> None:
        if resource_binding is None or not resource_binding.data:
            return
        staged_data = str(prepared_table_dir / f"{logical_name}.{kind}.gz")
        resource_specs[logical_name] = ResourceSpec(
            logical_name=logical_name,
            source=resource_binding.data,
            source_index=resource_binding.index or "",
            kind=kind,
            scientific_role=scientific_role,
            requires_chromosome_normalization=requires_chromosome_normalization,
            staged_data=staged_data,
            staged_index=f"{staged_data}.tbi",
            tsv_layout=tsv_layout,
            origin=origin,
        )

    cohort_vcf_source = binding_data(apply_plan.inputs.vcf)
    staged_cohort_vcf = str(prepared_dir / "cohort.vcf.gz")
    resource_specs["cohort"] = ResourceSpec(
        logical_name="cohort",
        source=cohort_vcf_source,
        source_index=apply_plan.inputs.vcf.index or "",
        kind="vcf",
        scientific_role="cohort_vcf",
        requires_chromosome_normalization=False,
        staged_data=staged_cohort_vcf,
        staged_index=f"{staged_cohort_vcf}.tbi",
        tsv_layout=None,
        origin="input",
    )
    resource_specs["fasta"] = _fasta_spec(effective_reference.fasta)
    add_table_resource(
        "calibration_regions",
        effective_reference.calibration_regions,
        "bed",
        "calibration_regions",
        requires_chromosome_normalization=True,
    )
    add_table_resource(
        "exclusion",
        effective_reference.exclusion,
        plan.reference.exclusion_kind,
        "variants_exclude",
        requires_chromosome_normalization=True,
    )
    # Bundle apply reads certified per-partition BEDs directly.  A pooled
    # testing-unit staging resource would lose the unit/partition relation.
    add_table_resource(
        "cadd_snv",
        effective_reference.cadd_snv,
        "tsv",
        "cadd_snv",
        requires_chromosome_normalization=True,
        tsv_layout="position",
        origin="configured",
    )
    if plan.model.include_indels:
        add_table_resource(
            "cadd_indel",
            effective_reference.cadd_indel,
            "tsv",
            "cadd_indel",
            requires_chromosome_normalization=True,
            tsv_layout="position",
            origin="configured",
        )

    for exclusion in apply_plan.effective_exclusions:
        add_table_resource(
            f"exclusion_{exclusion.identifier}",
            exclusion.binding,
            exclusion.kind,
            exclusion.role,
            origin=resource_origin,
        )

    if any(
        resource_spec.scientific_role in {"variants_exclude", "par_regions"}
        for resource_spec in resource_specs.values()
    ):
        combined_exclusions_vcf = str(prepared_dir / "combined_exclusions.vcf.gz")
        resource_specs["combined_exclusions"] = ResourceSpec(
            logical_name="combined_exclusions",
            source=combined_exclusions_vcf,
            source_index="",
            kind="vcf",
            scientific_role="combined_exclusions",
            requires_chromosome_normalization=False,
            staged_data=combined_exclusions_vcf,
            staged_index=f"{combined_exclusions_vcf}.tbi",
            tsv_layout=None,
            origin="run_derived",
            stage=False,
        )
    return ApplyResourceRegistry(
        resources=MappingProxyType(dict(resource_specs)),
        prepared_dir=prepared_dir,
    )


def _fasta_spec(fasta_binding: FastaBinding) -> ResourceSpec:
    """Keep the configured FASTA and its declared indexes in the catalog."""
    companion_indexes = (fasta_binding.gzi,) if fasta_binding.gzi else ()
    return ResourceSpec(
        logical_name="fasta",
        source=fasta_binding.data,
        source_index=fasta_binding.fai,
        kind="fasta",
        scientific_role="reference_fasta",
        requires_chromosome_normalization=False,
        staged_data=fasta_binding.data,
        staged_index=fasta_binding.fai,
        tsv_layout=None,
        origin="configured",
        companion_indexes=companion_indexes,
        stage=False,
    )
