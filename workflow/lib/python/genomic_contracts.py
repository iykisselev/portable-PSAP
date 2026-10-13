"""Canonical genomic identifiers shared by generation and validation code.

The generated-bundle contract uses chromosome names without a ``chr`` prefix,
``MT`` for mitochondrial sequence, one-based VCF positions, and uppercase
alleles.  Keep that conversion in this deliberately small module so resource
readers do not develop subtly different keys.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


DNA_ALLELES = frozenset({"A", "C", "G", "T"})
AUTOSOMES = tuple(str(value) for value in range(1, 23))
GENERATION_CONTIGS = AUTOSOMES + ("X",)
VALIDATION_CONTIGS = GENERATION_CONTIGS + ("Y",)
CONTIG_ORDER = {
    name: ordinal for ordinal, name in enumerate((*VALIDATION_CONTIGS, "MT"))
}


def contig_sort_key(chromosome: str) -> tuple[int, str]:
    canonical = canonical_chromosome(chromosome)
    return CONTIG_ORDER.get(canonical, len(CONTIG_ORDER)), canonical


def parse_contig_id(header_line: str) -> str | None:
    """Read the declared VCF contig ID without changing its source spelling."""
    prefix = "##contig=<ID="
    if not header_line.startswith(prefix):
        return None
    identifier = header_line[len(prefix) :].split(",", 1)[0].split(">", 1)[0]
    return identifier.strip() or None


class AlleleContractError(ValueError):
    """Raised when an allele cannot satisfy the shared VCF identity contract."""


@dataclass(frozen=True)
class FaiRecord:
    """One canonicalized FASTA-index record with byte-addressing layout."""

    name: str
    length: int
    offset: int
    line_bases: int
    line_bytes: int


def read_fai(path: str | Path) -> dict[str, FaiRecord]:
    """Read a strict FAI into canonical contig keys.

    A five-column FAI is required even when a caller needs only contig lengths:
    accepting an incomplete index in one stage and rejecting it in another
    produces failures much later than the resource boundary.
    """
    records: dict[str, FaiRecord] = {}
    with Path(path).open(encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, start=1):
            fields = raw.rstrip("\n").split("\t")
            if len(fields) < 5:
                raise ValueError(f"FAI line {line_number} needs five columns")
            name = fields[0].strip()
            if not name:
                raise ValueError(f"FAI line {line_number} has an empty contig name")
            try:
                length, offset, line_bases, line_bytes = map(int, fields[1:5])
            except ValueError as exc:
                raise ValueError(
                    f"FAI line {line_number} has non-integer layout values"
                ) from exc
            if length <= 0 or offset < 0 or line_bases <= 0 or line_bytes < line_bases:
                raise ValueError(f"FAI line {line_number} has an invalid contig layout")
            chromosome = canonical_chromosome(name)
            if chromosome in records:
                raise ValueError(f"FAI has duplicate canonical contig {chromosome!r}")
            records[chromosome] = FaiRecord(
                name, length, offset, line_bases, line_bytes
            )
    return records


def canonical_chromosome(chromosome: str) -> str:
    """Return the bundle's canonical chromosome spelling.

    This normalizes only representation, not scientific scope: callers remain
    responsible for deciding which contigs their contract permits.
    """
    value = str(chromosome).strip()
    if value.lower().startswith("chr"):
        value = value[3:]
    value = value.upper()
    return "MT" if value == "M" else value


def allele_key(
    chromosome: str, position: int, reference: str, alternate: str
) -> tuple[str, int, str, str]:
    """Return the canonical, VCF-style identity for one alternate allele."""
    return (
        canonical_chromosome(chromosome),
        int(position),
        str(reference).upper(),
        str(alternate).upper(),
    )


def variant_class(reference: str, alternate: str) -> str:
    """Classify one normalized biallelic VCF allele for null generation.

    This is deliberately a sequence-level classifier: reference-aware
    left-alignment validation belongs to the resource-validation stage.
    """
    reference = str(reference).upper()
    alternate = str(alternate).upper()
    if not reference or not alternate:
        return "unsupported"
    if set(reference) - DNA_ALLELES or set(alternate) - DNA_ALLELES:
        return "unsupported"
    if len(reference) == len(alternate) == 1 and reference != alternate:
        return "snv"
    if len(reference) != len(alternate):
        return "indel"
    return "unsupported"


def reference_span(position: int, reference: str) -> tuple[int, int]:
    """Return the zero-based, half-open reference span of one VCF allele."""
    try:
        position = int(position)
    except (TypeError, ValueError) as exc:
        raise AlleleContractError("VCF position must be an integer") from exc
    reference = str(reference).upper()
    if position < 1:
        raise AlleleContractError("VCF position must be positive")
    if not reference or set(reference) - DNA_ALLELES:
        raise AlleleContractError("REF must be a nonempty A/C/G/T sequence")
    start = position - 1
    return start, start + len(reference)


def is_minimal_allele(reference: str, alternate: str) -> bool:
    """Return whether a supported allele has no trim-able common bases.

    VCF minimization may only trim a common prefix or suffix while both
    alleles retain at least one base.  This deliberately does not establish
    left alignment; use :func:`is_left_aligned_allele` for that independent
    normalization property.
    """
    reference = str(reference).upper()
    alternate = str(alternate).upper()
    if variant_class(reference, alternate) == "unsupported":
        return False
    return not (
        len(reference) > 1
        and len(alternate) > 1
        and (reference[0] == alternate[0] or reference[-1] == alternate[-1])
    )


def is_left_aligned_allele(position: int, reference: str, alternate: str) -> bool:
    """Return whether a supported, minimal allele cannot shift one base left.

    For an indel, a left shift is possible exactly when the terminal REF and
    ALT bases agree.  The shifted representation replaces the preceding
    reference base plus the truncated REF with that base plus the truncated
    ALT.  Position one therefore has no possible left shift.  The caller is
    responsible for validating REF against the configured FASTA.
    """
    try:
        position = int(position)
    except (TypeError, ValueError) as exc:
        raise AlleleContractError("VCF position must be an integer") from exc
    if position < 1:
        raise AlleleContractError("VCF position must be positive")
    reference = str(reference).upper()
    alternate = str(alternate).upper()
    if variant_class(reference, alternate) == "unsupported":
        return False
    return position == 1 or reference[-1] != alternate[-1]


def validated_allele_key(
    chromosome: str, position: int, reference: str, alternate: str
) -> tuple[str, int, str, str]:
    """Return a canonical key after enforcing supported sequence syntax."""
    try:
        key = allele_key(chromosome, position, reference, alternate)
    except (TypeError, ValueError) as exc:
        raise AlleleContractError("VCF position must be an integer") from exc
    if not key[0]:
        raise AlleleContractError("Chromosome must not be empty")
    reference_span(key[1], key[2])
    if variant_class(key[2], key[3]) == "unsupported":
        raise AlleleContractError(
            "Allele must be an SNV or length-changing A/C/G/T indel"
        )
    return key


def allele_key_text(key: tuple[str, int, str, str]) -> str:
    """Render a canonical alternate-allele key for diagnostics and joins."""
    return f"{key[0]}:{key[1]}:{key[2]}:{key[3]}"
