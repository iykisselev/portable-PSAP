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
    def window_dir(self) -> Path:
        return self.plan.scalable_dir / "windows"

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
    def bundle_paths(self) -> Mapping[str, str]:
        return self.execution_layout.bundle_paths

    def bundle_outputs(self) -> Mapping[str, str | list[str]]:
        """Return the non-lookup outputs owned by bundle finalization."""
        return MappingProxyType(
            {
                (
                    "partition_statuses"
                    if name == "unit_partition_status"
                    else "testing_unit_id_overlap"
                    if name == "testing_unit_id_overlap"
                    else name
                ): path
                for name, path in self.bundle_paths.items()
                if name not in {"validation", "marker"}
            }
        )

    @property
    def bundle_resource_identities(self) -> Mapping[str, Mapping[str, str]]:
        """Return declared scientific identities and source paths without hashing inputs."""
        identities = {
            "fasta": self.reference.fasta_id,
            "exclusion": self.reference.exclusion_id,
            "score_snv": self.reference.cadd_snv_id,
            "score_indel": self.reference.cadd_indel_id or "score_indel",
            "population_vcf": f"population_vcf:{self.plan.af_source['overall']}",
            "par_regions": "par_regions",
            "mutation_rate": "mutation_rate",
            "mutation_rate_indel": "mutation_rate_indel",
        }
        paths = {name: binding.data for name, binding in self.raw_bindings().items()}
        for partition in self.plan.testing_unit_partitions:
            name = f"testing_unit_partition.{partition.identifier}"
            identities[name] = self.reference.testing_unit_partitions.get(
                partition.identifier, partition.identifier
            )
            paths[name] = partition.bed.data
        paths.update(self.exclusion_paths_by_id())
        return MappingProxyType(
            {
                name: MappingProxyType(
                    {"identity": identities.get(name, name), "path": str(path)}
                )
                for name, path in paths.items()
                if path
            }
        )

    def raw_bindings(self) -> Mapping[str, FileBinding]:
        """Keep each raw data path associated with its declared index."""
        return {
            "population_vcf": self.plan.population_vcf,
            "score_snv": self.plan.score_snv,
            "score_indel": self.plan.score_indel or FileBinding(None),
            "exclusion": self.plan.exclusion,
            "testing_units": self.plan.testing_unit_partitions[0].bed,
            "par_regions": self.plan.par_regions,
            "fasta": FileBinding(self.plan.fasta.data),
            "mutation_rate": self.plan.mutation_rate or FileBinding(None),
            "mutation_rate_indel": self.plan.mutation_rate_indel or FileBinding(None),
        }

    def raw_inputs(self, operation: str) -> dict[str, list[str] | str]:
        """Declare only the raw dependencies consumed by one generation step."""
        names = {
            "scope": (
                "population_vcf",
                "score_snv",
                "score_indel",
                "exclusion",
                "testing_units",
                "par_regions",
                "fasta",
                "mutation_rate",
                "mutation_rate_indel",
            ),
            "window": ("population_vcf", "score_snv", "score_indel", "fasta"),
            "statistics": ("mutation_rate", "mutation_rate_indel"),
        }[operation]
        bindings = self.raw_bindings()
        result: dict[str, list[str] | str] = {}
        indexed = {
            "population_vcf",
            "score_snv",
            "score_indel",
            "mutation_rate",
            "mutation_rate_indel",
        }
        for name in names:
            binding = bindings[name]
            if name == "score_indel" and not self.model.include_indels:
                binding = FileBinding(None)
            result[name] = binding.data_inputs()
            if name in indexed:
                result[f"{name}_index"] = binding.index_inputs()
        if operation == "scope":
            result["testing_unit_partitions"] = [
                partition.bed.data
                for partition in self.plan.testing_unit_partitions
                if partition.bed.data
            ]
        if operation in {"scope", "window"}:
            result["fasta_fai"] = self.plan.fasta.fai
        if operation == "scope":
            result["exclusions"] = self.exclusion_data()
        return result

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

    def window_artifact(self, window_id: str, artifact_name: str) -> str:
        """Return one window artifact path for a known artifact name."""
        suffixes = {
            "allele_states": "allele_states.parquet",
            "allele_membership": "allele_membership.parquet",
            "allele_scores": "allele_scores.parquet",
            "candidate_audit": "candidate_audit.parquet",
        }
        return str(self.window_dir / f"{window_id}.{suffixes[artifact_name]}")

    def bucket_artifact(self, bucket_id: str, artifact_name: str) -> str:
        """Return one lookup-bucket artifact path for a known artifact name."""
        suffixes = {
            "lookup_het": "het.tsv",
            "lookup_hom": "hom.tsv",
            "lookup_hem_chrX": "hem_chrX.tsv",
            "diagnostics": "diagnostics.tsv",
            "partition_model_diagnostics": "partition-model-diagnostics.tsv",
        }
        return str(self.lookup_bucket_dir / f"{bucket_id}.{suffixes[artifact_name]}")

    def bucket_partition_model_artifact(
        self, bucket_id: str, partition_id: str, model_id: str
    ) -> str:
        return str(self.lookup_bucket_dir / f"{bucket_id}.{partition_id}.{model_id}.tsv")

    def statistics_artifact(self, window_id: str, bucket_id: str) -> str:
        """Return one window-local statistics artifact path."""
        return str(self.statistics_dir / f"{window_id}.bucket-{bucket_id}.parquet")
