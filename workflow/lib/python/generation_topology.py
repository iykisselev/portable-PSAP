"""Static generation paths and policy derived without inspecting inputs."""

from __future__ import annotations

from dataclasses import dataclass
from decimal import Decimal, ROUND_CEILING
from pathlib import Path
from types import MappingProxyType
from typing import Mapping

from workflow_plan import (
    GenerationPlan,
    GenerationSettings,
    ModelPlan,
    generated_bundle_paths,
)


@dataclass(frozen=True)
class GenerationScopeArtifacts:
    """Deterministic outputs of the data-dependent scope checkpoint."""

    artifacts: Mapping[str, str]


@dataclass(frozen=True)
class GenerationExecutionLayout:
    """Derived generation paths, buckets, and enabled model roles."""

    calculated_dir: Path
    statistics_dir: Path
    lookup_bucket_dir: Path
    bin_count: int
    roles: tuple[str, ...]
    bundle_paths: Mapping[str, str]


def build_generation_scope_artifacts(plan: GenerationPlan) -> GenerationScopeArtifacts:
    """Return scope-checkpoint paths without reading genomic resource files."""
    scalable_dir = plan.scalable_dir
    prepared_dir = scalable_dir / "prepared"
    return GenerationScopeArtifacts(
        artifacts=MappingProxyType(
            {
                "partitions": str(prepared_dir / "partitions.json"),
                "buckets": str(prepared_dir / "unit_buckets.json"),
                "units": str(prepared_dir / "provisional_units.parquet"),
                "unit_intervals": str(prepared_dir / "unit_intervals.parquet"),
                "masked_unit_intervals": str(
                    prepared_dir / "masked_unit_intervals.parquet"
                ),
                "effective_mask": str(
                    prepared_dir / "effective_calibration_mask.parquet"
                ),
                "resource_identities": str(prepared_dir / "resource_identities.json"),
            }
        ),
    )


def build_generation_execution_layout(
    plan: GenerationPlan,
    model: ModelPlan,
    settings: GenerationSettings,
) -> GenerationExecutionLayout:
    """Derive deterministic paths and work buckets from immutable plans."""
    calculated_dir = plan.scalable_dir / "calculated"
    statistics_dir = calculated_dir / "partition_statistics"
    lookup_bucket_dir = calculated_dir / "lookup_buckets"
    bin_count = int(
        (
            Decimal(str(model.score_max)) / Decimal(str(settings.score_bin_width))
        ).to_integral_value(rounding=ROUND_CEILING)
    )
    return GenerationExecutionLayout(
        calculated_dir=calculated_dir,
        statistics_dir=statistics_dir,
        lookup_bucket_dir=lookup_bucket_dir,
        bin_count=bin_count,
        roles=("het", "hom") + (("hem_chrX",) if model.hemizygote else ()),
        bundle_paths=generated_bundle_paths(plan.bundle_dir),
    )
