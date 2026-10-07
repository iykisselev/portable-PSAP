"""Partition-manifest adapters for bounded null-distribution generation."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from generation_manifests import read_partition_manifest
from genomic_contracts import read_fai


@dataclass(frozen=True)
class GenomePartition:
    """One half-open genomic region with a stable output identifier."""

    identifier: str
    chromosome: str
    start: int
    end: int
    ordinal: int


def read_fasta_index(path: str | Path) -> dict[str, tuple[str, int]]:
    """Compatibility view of the shared FAI records for partitioning."""
    return {chromosome: (record.name, record.length) for chromosome, record in read_fai(path).items()}


def build_partitions(
    fasta_index: str | Path,
    testing_unit_intervals: list[tuple[str, int, int]],
    partition_size: int,
    *,
    candidate_owner_intervals: list[tuple[str, int, int]] | None = None,
) -> list[GenomePartition]:
    """Return fixed-size windows that intersect eligible candidate intervals.

    ``candidate_owner_intervals`` is an optional indel anchor halo.  It is
    separate from the ordinary callable-unit intervals so downstream code can
    preserve exact span membership while still scheduling a preceding owner
    partition for a boundary-crossing deletion.
    """
    if partition_size <= 0:
        raise ValueError("generation.execution.partition_size_bp must be positive")
    fasta_contigs = read_fasta_index(fasta_index)
    candidate_intervals = list(testing_unit_intervals)
    if candidate_owner_intervals is not None:
        candidate_intervals.extend(candidate_owner_intervals)
    chromosomes = {chromosome for chromosome, _, _ in candidate_intervals}
    missing = sorted(chromosomes - set(fasta_contigs))
    if missing:
        raise ValueError("Testing-unit contigs are absent from the FASTA index: " + ", ".join(missing))

    selected_windows: dict[str, set[int]] = {}
    for chromosome, interval_start, interval_end in candidate_intervals:
        _, contig_length = fasta_contigs[chromosome]
        if interval_start < 0 or interval_end <= interval_start:
            raise ValueError(f"Testing-unit interval on {chromosome} is invalid")
        if interval_end > contig_length:
            raise ValueError(
                f"Testing-unit interval on {chromosome} extends beyond the FASTA contig"
            )
        first_window = interval_start // partition_size
        last_window = (interval_end - 1) // partition_size
        selected_windows.setdefault(chromosome, set()).update(
            range(first_window, last_window + 1)
        )

    partitions: list[GenomePartition] = []
    ordinal = 0
    for chromosome in sorted(selected_windows):
        _, length = fasta_contigs[chromosome]
        for window in sorted(selected_windows[chromosome]):
            start = window * partition_size
            end = min(start + partition_size, length)
            partitions.append(
                GenomePartition(
                    identifier=f"{chromosome}_{start:012d}_{end:012d}",
                    chromosome=chromosome,
                    start=start,
                    end=end,
                    ordinal=ordinal,
                )
            )
            ordinal += 1
    return partitions


def indel_candidate_owner_intervals(
    testing_unit_intervals: list[tuple[str, int, int]],
    max_reference_length: int,
) -> list[tuple[str, int, int]]:
    """Expand unit intervals leftward for indel anchor-owner discovery.

    An indel with a REF of at most ``max_reference_length`` can overlap a unit
    while its anchor is up to ``max_reference_length - 1`` bases upstream.
    The returned intervals are not membership intervals: exact span overlap is
    still checked by :func:`generation_indels.indel_candidate_memberships`.
    """
    if max_reference_length < 1:
        raise ValueError("Indel maximum reference length must be positive")
    return [
        (chromosome, max(0, start - max_reference_length + 1), end)
        for chromosome, start, end in testing_unit_intervals
    ]


def read_genome_partitions(path: str | Path) -> list[GenomePartition]:
    """Load validated scope partitions into the partition-local record type."""
    return [
        GenomePartition(
            identifier=entry.identifier,
            chromosome=entry.chromosome,
            start=entry.start,
            end=entry.end,
            ordinal=entry.ordinal,
        )
        for entry in read_partition_manifest(path)
    ]


def partition_from_manifest(path: str | Path, identifier: str) -> GenomePartition:
    """Resolve one requested wildcard without reconstructing coordinates from text."""
    for partition in read_genome_partitions(path):
        if partition.identifier == identifier:
            return partition
    raise ValueError(f"Unknown generation partition {identifier!r}")
