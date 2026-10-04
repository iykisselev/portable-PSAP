"""One rule-facing context for generation plans and their derived layouts."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Mapping

from generation_topology import (
    GenerationExecutionLayout,
    GenerationScopeArtifacts,
    build_generation_execution_layout,
    build_generation_scope_artifacts,
)
from workflow_plan import (
    FileBinding,
    GenerationPlan,
    GenerationSettings,
    ModelPlan,
    WorkflowPlan,
)


@dataclass(frozen=True)
class GenerationRuleContext:
    """Expose generation inputs, paths, and derived execution layout together."""

    plan: GenerationPlan
    model: ModelPlan
    settings: GenerationSettings
    scope_layout: GenerationScopeArtifacts
    execution_layout: GenerationExecutionLayout

    @classmethod
    def from_workflow_plan(cls, workflow_plan: WorkflowPlan) -> "GenerationRuleContext":
        """Build the generation-only rule context from the resolved workflow plan."""
        if workflow_plan.generation is None:
            raise ValueError("Generation rules require a generate or full WorkflowPlan")
        generation_plan = workflow_plan.generation
        return cls(
            plan=generation_plan,
            model=workflow_plan.model,
            settings=workflow_plan.generation_settings,
            scope_layout=build_generation_scope_artifacts(generation_plan),
            execution_layout=build_generation_execution_layout(
                generation_plan,
                workflow_plan.model,
                workflow_plan.generation_settings,
            ),
        )

    @property
    def partition_dir(self) -> Path:
        return self.plan.scalable_dir / "partitions"

    @property
    def log_dir(self) -> Path:
        return self.plan.log_dir

    @property
    def scope_paths(self) -> Mapping[str, str]:
        return self.scope_layout.artifacts

    @property
    def statistics_dir(self) -> Path:
        return self.execution_layout.statistics_dir

    @property
    def lookup_bucket_dir(self) -> Path:
        return self.execution_layout.lookup_bucket_dir

    @property
    def bin_count(self) -> int:
        return self.execution_layout.bin_count

    @property
    def roles(self) -> tuple[str, ...]:
        return self.execution_layout.roles

    @property
    def hem_enabled(self) -> bool:
        return self.model.hemizygote

    @property
    def bundle_paths(self) -> Mapping[str, str]:
        return self.execution_layout.bundle_paths

    def index_input(self, binding: FileBinding) -> list[str]:
        """Return a declared index as an optional Snakemake input list."""
        return [binding.index] if binding.index else []

    def mutation_rate_data(self) -> list[str]:
        """Return the optional mutation-rate table without a placeholder path."""
        return [self.plan.mutation_rate.data] if self.plan.mutation_rate else []

    def mutation_rate_index(self) -> list[str]:
        """Return the optional mutation-rate index without a placeholder path."""
        if self.plan.mutation_rate and self.plan.mutation_rate.index:
            return [self.plan.mutation_rate.index]
        return []

    def exclusion_data(self) -> list[str]:
        """Return exclusion dependencies in deterministic identifier order."""
        return [
            exclusion.binding.data
            for exclusion in self.plan.exclusions
            if exclusion.binding.data
        ]

    def exclusion_paths_by_id(self) -> Mapping[str, str]:
        """Return the authoritative exclusion identity-to-path association."""
        return {
            exclusion.identifier: exclusion.binding.data
            for exclusion in self.plan.exclusions
            if exclusion.binding.data
        }

    def partition_artifact(self, partition_id: str, artifact_name: str) -> str:
        """Return one partition artifact path for a known artifact name."""
        suffixes = {
            "allele_states": "allele_states.parquet",
            "allele_membership": "allele_membership.parquet",
            "allele_scores": "allele_scores.parquet",
        }
        return str(self.partition_dir / f"{partition_id}.{suffixes[artifact_name]}")

    def bucket_artifact(self, bucket_id: str, artifact_name: str) -> str:
        """Return one lookup-bucket artifact path for a known artifact name."""
        suffixes = {
            "lookup_het": "het.tsv",
            "lookup_hom": "hom.tsv",
            "lookup_hem": "hem_chrX.tsv",
            "diagnostics": "diagnostics.tsv",
        }
        return str(self.lookup_bucket_dir / f"{bucket_id}.{suffixes[artifact_name]}")

    def statistics_artifact(self, partition_id: str, bucket_id: str) -> str:
        """Return one partition-local statistics artifact path."""
        return str(self.statistics_dir / f"{partition_id}.bucket-{bucket_id}.parquet")
