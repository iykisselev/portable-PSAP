"""One rule-facing context for generation plans and their derived layouts."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from types import MappingProxyType
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
    ReferenceMetadata,
    WorkflowPlan,
)


@dataclass(frozen=True)
class GenerationRuleContext:
    """Expose generation inputs, paths, and derived execution layout together."""

    plan: GenerationPlan
    model: ModelPlan
    settings: GenerationSettings
    reference: ReferenceMetadata
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
            reference=workflow_plan.reference,
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

    def bundle_outputs(self) -> Mapping[str, str | list[str]]:
        """Return finalization outputs, retaining an empty optional hem binding."""
        outputs = {
            name: path
            for name, path in self.bundle_paths.items()
            if name not in {"validation", "marker", "lookup_hem"}
        }
        outputs["lookup_hem"] = (
            self.bundle_paths["lookup_hem"] if self.hem_enabled else []
        )
        return MappingProxyType(outputs)

    @property
    def bundle_resource_ids(self) -> Mapping[str, str]:
        """Return stable scientific identities for generated-bundle provenance."""
        return MappingProxyType(
            {
                "fasta": self.reference.fasta_id,
                "exclusion": self.reference.exclusion_id,
                "testing_units": self.reference.testing_unit_id,
                "score_snv": self.reference.cadd_snv_id,
                "score_indel": self.reference.cadd_indel_id or "score_indel",
                "population_vcf": f"population_vcf:{self.plan.af_source['overall']}",
                "par_regions": "par_regions",
                "mutation_rate": "mutation_rate",
                "mutation_rate_indel": "mutation_rate_indel",
            }
        )

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

    def indel_mutation_rate_data(self) -> list[str]:
        return [self.plan.mutation_rate_indel.data] if self.plan.mutation_rate_indel else []

    def indel_mutation_rate_index(self) -> list[str]:
        if self.plan.mutation_rate_indel and self.plan.mutation_rate_indel.index:
            return [self.plan.mutation_rate_indel.index]
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
            "candidate_audit": "candidate_audit.parquet",
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
