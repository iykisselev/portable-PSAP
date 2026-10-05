"""Canonical genomic identifiers shared by generation and validation code.

The generated-bundle contract uses chromosome names without a ``chr`` prefix,
``MT`` for mitochondrial sequence, one-based VCF positions, and uppercase
alleles.  Keep that conversion in this deliberately small module so resource
readers do not develop subtly different keys.
"""

from __future__ import annotations


DNA_ALLELES = frozenset({"A", "C", "G", "T"})
VARIANT_CLASSES = frozenset({"snv", "indel", "unsupported"})


class AlleleContractError(ValueError):
    """Raised when an allele cannot satisfy the shared VCF identity contract."""


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
