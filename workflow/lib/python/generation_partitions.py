"""Partition-manifest adapters for bounded null-distribution generation."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from generation_manifests import read_partition_manifest as _read_partition_manifest
from genomic_contracts import canonical_chromosome


@dataclass(frozen=True)
class GenomePartition:
    """One half-open genomic region with a stable output identifier."""

    identifier: str
    chromosome: str
    start: int
    end: int
    ordinal: int


def read_fasta_index(path: str | Path) -> dict[str, tuple[str, int]]:
    """Return canonical contigs, source FASTA names, and lengths from a FAI."""
    contigs: dict[str, tuple[str, int]] = {}
    with Path(path).open(encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, start=1):
            fields = raw.rstrip("\n").split("\t")
            if len(fields) < 2:
                raise ValueError(f"FAI line {line_number} needs a name and length")
            chromosome = canonical_chromosome(fields[0])
            if chromosome in contigs:
                raise ValueError(f"FAI has duplicate canonical contig {chromosome!r}")
            length = int(fields[1])
            if length <= 0:
                raise ValueError(f"FAI contig {fields[0]!r} has a non-positive length")
            contigs[chromosome] = (fields[0], length)
    return contigs


def build_partitions(
    fasta_index: str | Path,
    chromosomes: set[str],
    partition_size: int,
) -> list[GenomePartition]:
    """Divide testing-unit contigs into stable, non-overlapping windows."""
    if partition_size <= 0:
        raise ValueError("generation.execution.partition_size_bp must be positive")
    fasta_contigs = read_fasta_index(fasta_index)
    missing = sorted(chromosomes - set(fasta_contigs))
    if missing:
        raise ValueError("Testing-unit contigs are absent from the FASTA index: " + ", ".join(missing))

    partitions: list[GenomePartition] = []
    ordinal = 0
    for chromosome in sorted(chromosomes):
        _, length = fasta_contigs[chromosome]
        for start in range(0, length, partition_size):
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


def read_partition_manifest(path: str | Path) -> list[GenomePartition]:
    """Load validated scope partitions into the partition-local record type."""
    return [
        GenomePartition(
            identifier=entry.identifier,
            chromosome=entry.chromosome,
            start=entry.start,
            end=entry.end,
            ordinal=entry.ordinal,
        )
        for entry in _read_partition_manifest(path)
    ]


def partition_from_manifest(path: str | Path, identifier: str) -> GenomePartition:
    """Resolve one requested wildcard without reconstructing coordinates from text."""
    for partition in read_partition_manifest(path):
        if partition.identifier == identifier:
            return partition
    raise ValueError(f"Unknown generation partition {identifier!r}")
