"""Scheduler-independent, named execution policies for workflow rules.

Rules retain explicit scheduler and software declarations. They select only a
logical policy name here, so baseline resources and size-based retry scaling
are tested in one place rather than duplicated as numeric multipliers.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from pathlib import Path
from types import MappingProxyType
from typing import Callable, Iterable, Mapping


@dataclass(frozen=True)
class ExecutionPolicy:
    """Scheduler mem_mb requests and runtime minutes; input scaling uses GiB.

    Existing allocation values are retained independently of byte-unit naming.
    """

    baseline_mem_mb: int
    baseline_runtime_minutes: int
    mem_mb_per_gib: int = 0
    runtime_minutes_per_gib: int = 0
    scale_with_attempt: bool = False


@dataclass(frozen=True)
class RuleResourceCallbacks:
    """Snakemake callbacks that resolve both resources from one policy."""

    mem_mb: Callable[[object, Iterable[object], int], int]
    runtime: Callable[[object, Iterable[object], int], int]


RULE_POLICY_NAMES = frozenset(
    {
        "prepare_generation_scope",
        "prepare_generation_window",
        "calculate_window_null_distribution",
        "reduce_null_distribution_bucket",
        "finalize_null_distribution_bundle",
        "validate_null_distribution_bundle",
        "validate_cohort",
        "validate_resources",
        "validate_model",
        "merge_validation",
        "require_valid_inputs",
        "split_genotypes",
        "extract_gt",
        "prepare_sample_genotypes",
        "annotate_variants",
        "annotate_partition_memberships",
        "summarize_annotations",
        "require_valid_annotations",
        "filter_calibration_regions",
        "filter_excluded_records",
        "build_sample_candidates",
        "build_control_reference",
        "build_cohort_report",
        "prepare_vcf_resource",
        "combine_exclusions",
        "apply_psap",
    }
)


TABLE_RESOURCE_POLICY_NAMES: Mapping[str, str] = MappingProxyType(
    {
        "calibration_regions": "stage_calibration_regions",
        "testing_units": "stage_testing_units",
        "cadd_snv": "stage_cadd_snv",
        "cadd_indel": "stage_cadd_indel",
        "variants_exclude": "stage_variants_exclude",
        "par_regions": "stage_par_regions",
    }
)


EXECUTION_POLICIES: Mapping[str, ExecutionPolicy] = MappingProxyType(
    {
        "prepare_generation_scope": ExecutionPolicy(
            baseline_mem_mb=4096, baseline_runtime_minutes=120
        ),
        "prepare_generation_window": ExecutionPolicy(
            baseline_mem_mb=16384, baseline_runtime_minutes=360
        ),
        "calculate_window_null_distribution": ExecutionPolicy(
            baseline_mem_mb=16384, baseline_runtime_minutes=240
        ),
        "reduce_null_distribution_bucket": ExecutionPolicy(
            baseline_mem_mb=4096, baseline_runtime_minutes=120
        ),
        "finalize_null_distribution_bundle": ExecutionPolicy(
            baseline_mem_mb=8192, baseline_runtime_minutes=120
        ),
        "validate_null_distribution_bundle": ExecutionPolicy(
            baseline_mem_mb=1024, baseline_runtime_minutes=120
        ),
        "validate_cohort": ExecutionPolicy(
            baseline_mem_mb=4096, baseline_runtime_minutes=30
        ),
        "validate_resources": ExecutionPolicy(
            baseline_mem_mb=4096, baseline_runtime_minutes=30
        ),
        "validate_model": ExecutionPolicy(
            baseline_mem_mb=512, baseline_runtime_minutes=10
        ),
        "merge_validation": ExecutionPolicy(
            baseline_mem_mb=512, baseline_runtime_minutes=5
        ),
        "require_valid_inputs": ExecutionPolicy(
            baseline_mem_mb=512, baseline_runtime_minutes=5
        ),
        "split_genotypes": ExecutionPolicy(
            baseline_mem_mb=4096, baseline_runtime_minutes=120
        ),
        "extract_gt": ExecutionPolicy(
            baseline_mem_mb=4096, baseline_runtime_minutes=120
        ),
        "prepare_sample_genotypes": ExecutionPolicy(
            baseline_mem_mb=16000, baseline_runtime_minutes=180
        ),
        "annotate_variants": ExecutionPolicy(
            baseline_mem_mb=32000,
            baseline_runtime_minutes=720,
            mem_mb_per_gib=512,
            runtime_minutes_per_gib=30,
            scale_with_attempt=True,
        ),
        "annotate_partition_memberships": ExecutionPolicy(
            baseline_mem_mb=4096, baseline_runtime_minutes=60
        ),
        "summarize_annotations": ExecutionPolicy(
            baseline_mem_mb=1024, baseline_runtime_minutes=30
        ),
        "require_valid_annotations": ExecutionPolicy(
            baseline_mem_mb=512, baseline_runtime_minutes=5
        ),
        "filter_calibration_regions": ExecutionPolicy(
            baseline_mem_mb=4096, baseline_runtime_minutes=180
        ),
        "filter_excluded_records": ExecutionPolicy(
            baseline_mem_mb=4096, baseline_runtime_minutes=180
        ),
        "build_sample_candidates": ExecutionPolicy(
            baseline_mem_mb=4096, baseline_runtime_minutes=60
        ),
        "build_control_reference": ExecutionPolicy(
            baseline_mem_mb=16000,
            baseline_runtime_minutes=360,
            mem_mb_per_gib=256,
            runtime_minutes_per_gib=20,
            scale_with_attempt=True,
        ),
        "build_cohort_report": ExecutionPolicy(
            baseline_mem_mb=16000,
            baseline_runtime_minutes=360,
            mem_mb_per_gib=256,
            runtime_minutes_per_gib=20,
            scale_with_attempt=True,
        ),
        "prepare_vcf_resource": ExecutionPolicy(
            baseline_mem_mb=4096,
            baseline_runtime_minutes=180,
            scale_with_attempt=True,
        ),
        "combine_exclusions": ExecutionPolicy(
            baseline_mem_mb=12000, baseline_runtime_minutes=180
        ),
        "apply_psap": ExecutionPolicy(
            baseline_mem_mb=16000, baseline_runtime_minutes=360
        ),
        "stage_calibration_regions": ExecutionPolicy(
            baseline_mem_mb=2048,
            baseline_runtime_minutes=60,
            scale_with_attempt=True,
        ),
        "stage_testing_units": ExecutionPolicy(
            baseline_mem_mb=2048,
            baseline_runtime_minutes=60,
            scale_with_attempt=True,
        ),
        "stage_cadd_snv": ExecutionPolicy(
            baseline_mem_mb=2048,
            baseline_runtime_minutes=360,
            mem_mb_per_gib=256,
            runtime_minutes_per_gib=20,
            scale_with_attempt=True,
        ),
        "stage_cadd_indel": ExecutionPolicy(
            baseline_mem_mb=2048,
            baseline_runtime_minutes=360,
            mem_mb_per_gib=256,
            runtime_minutes_per_gib=20,
            scale_with_attempt=True,
        ),
        "stage_variants_exclude": ExecutionPolicy(
            baseline_mem_mb=4096,
            baseline_runtime_minutes=180,
            scale_with_attempt=True,
        ),
        "stage_par_regions": ExecutionPolicy(
            baseline_mem_mb=4096,
            baseline_runtime_minutes=180,
            scale_with_attempt=True,
        ),
    }
)


def input_size_mb(input_files: Iterable[object]) -> float:
    """Return declared local-input bytes in MiB, ignoring unavailable paths."""
    size_mb = 0.0
    for item in input_files:
        try:
            size_mb += Path(str(item)).stat().st_size / (1024 * 1024)
        except OSError:
            pass
    return size_mb


def execution_policy(name: str) -> ExecutionPolicy:
    """Resolve a named policy and reject an unregistered resource request."""
    try:
        return EXECUTION_POLICIES[name]
    except KeyError:
        raise ValueError(f"Unknown execution policy: {name!r}") from None


def table_resource_policy(scientific_role: str) -> str:
    """Return the explicit staging policy for one supported table-resource role."""
    try:
        policy_name = TABLE_RESOURCE_POLICY_NAMES[scientific_role]
    except KeyError:
        raise ValueError(
            f"No staging policy is registered for table-resource role: "
            f"{scientific_role!r}"
        ) from None
    execution_policy(policy_name)
    return policy_name


def scaled_value(
    baseline: int,
    per_gib: int,
    input_files: Iterable[object],
    attempt: int,
    *,
    scale_with_attempt: bool,
) -> int:
    """Apply whole-GiB input scaling and optional Snakemake retry scaling."""
    return int(
        (baseline + math.ceil(input_size_mb(input_files) / 1024) * per_gib)
        * (attempt if scale_with_attempt else 1)
    )


def resource_mem_mb(name: str, input_files: Iterable[object], attempt: int) -> int:
    """Return memory requested for one logical execution policy."""
    policy = execution_policy(name)
    return scaled_value(
        policy.baseline_mem_mb,
        policy.mem_mb_per_gib,
        input_files,
        attempt,
        scale_with_attempt=policy.scale_with_attempt,
    )


def resource_runtime_minutes(
    name: str, input_files: Iterable[object], attempt: int
) -> int:
    """Return wall time requested for one logical execution policy."""
    policy = execution_policy(name)
    return scaled_value(
        policy.baseline_runtime_minutes,
        policy.runtime_minutes_per_gib,
        input_files,
        attempt,
        scale_with_attempt=policy.scale_with_attempt,
    )


def policy_resources(
    policy_name: str | Callable[[object], str],
) -> RuleResourceCallbacks:
    """Build Snakemake resource callbacks from one static or dynamic policy name."""
    if isinstance(policy_name, str):
        execution_policy(policy_name)

        def resolve_policy(_wildcards: object) -> str:
            return policy_name

    else:
        resolve_policy = policy_name

    def memory(wildcards: object, input: Iterable[object], attempt: int) -> int:
        return resource_mem_mb(resolve_policy(wildcards), input, attempt)

    def runtime(wildcards: object, input: Iterable[object], attempt: int) -> int:
        return resource_runtime_minutes(resolve_policy(wildcards), input, attempt)

    return RuleResourceCallbacks(mem_mb=memory, runtime=runtime)
