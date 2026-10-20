"""Derive generation paths and work layout without opening input files."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from types import MappingProxyType
from typing import Mapping

from bundle_contract import fixed_width_bin_count, model_roles
from partitioned_bundle_contract import bundle_component_paths, require_partition_id
from workflow_plan import (
    GenerationPlan,
    GenerationSettings,
    ModelPlan,
    generated_bundle_paths,
)


@dataclass(frozen=True)
class GenerationScopeArtifacts:
    """Prepared-artifact paths, including deferred resource provenance."""

    artifacts: Mapping[str, str]


@dataclass(frozen=True)
class GenerationExecutionLayout:
    """Output directories, score bins, and enabled model roles for one run."""

    calculated_dir: Path
    statistics_dir: Path
    lookup_bucket_dir: Path
    bin_count: int
    roles: tuple[str, ...]
    bundle_paths: Mapping[str, str]
    partition_lookup_paths: Mapping[tuple[str, str], str]
    resolved_partition_paths: Mapping[str, str]
    eligible_partition_paths: Mapping[str, str]


def build_generation_scope_artifacts(plan: GenerationPlan) -> GenerationScopeArtifacts:
    """Choose the scope-checkpoint file paths without reading genomic inputs."""
    scalable_dir = plan.scalable_dir
    prepared_dir = scalable_dir / "prepared"
    return GenerationScopeArtifacts(
        artifacts=MappingProxyType(
            {
                "windows": str(prepared_dir / "windows.json"),
                "buckets": str(prepared_dir / "unit_buckets.json"),
                "units": str(prepared_dir / "provisional_units.parquet"),
                "unit_intervals": str(prepared_dir / "unit_intervals.parquet"),
                "masked_unit_intervals": str(
                    prepared_dir / "masked_unit_intervals.parquet"
                ),
                "effective_mask": str(
                    prepared_dir / "effective_calibration_mask.parquet"
                ),
                "resolved_partition_intervals": str(
                    prepared_dir / "resolved_partition_intervals.parquet"
                ),
                "partition_statuses": str(prepared_dir / "unit_partition_status.parquet"),
                "testing_unit_id_overlap": str(
                    prepared_dir / "testing_unit_id_overlap.parquet"
                ),
            }
        ),
    )


def build_generation_execution_layout(
    plan: GenerationPlan,
    model: ModelPlan,
    settings: GenerationSettings,
) -> GenerationExecutionLayout:
    """Choose reproducible output paths and work buckets from the run plan."""
    calculated_dir = plan.scalable_dir / "calculated"
    statistics_dir = calculated_dir / "window_statistics"
    lookup_bucket_dir = calculated_dir / "lookup_buckets"
    bin_count = fixed_width_bin_count(model.score_max, settings.score_bin_width)
    return GenerationExecutionLayout(
        calculated_dir=calculated_dir,
        statistics_dir=statistics_dir,
        lookup_bucket_dir=lookup_bucket_dir,
        bin_count=bin_count,
        roles=model_roles(model.hemizygote),
        bundle_paths=generated_bundle_paths(plan.bundle_dir),
        partition_lookup_paths=bundle_component_paths(
            plan.bundle_dir,
            tuple(partition.identifier for partition in plan.testing_unit_partitions),
            model.hemizygote,
        ),
        resolved_partition_paths=MappingProxyType(
            {
                partition.identifier: str(
                    plan.bundle_dir / f"resolved.{require_partition_id(partition.identifier)}.bed.gz"
                )
                for partition in plan.testing_unit_partitions
            }
        ),
        eligible_partition_paths=MappingProxyType(
            {
                partition.identifier: str(
                    plan.bundle_dir / f"eligible.{require_partition_id(partition.identifier)}.bed.gz"
                )
                for partition in plan.testing_unit_partitions
            }
        ),
    )
