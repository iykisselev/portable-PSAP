"""Canonical genomic identifiers shared by generation and validation code.

The generated-bundle contract uses chromosome names without a ``chr`` prefix,
``MT`` for mitochondrial sequence, one-based VCF positions, and uppercase
alleles.  Keep that conversion in this deliberately small module so resource
readers do not develop subtly different keys.
"""

from __future__ import annotations


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


def allele_key_text(key: tuple[str, int, str, str]) -> str:
    """Render a canonical alternate-allele key for diagnostics and joins."""
    return f"{key[0]}:{key[1]}:{key[2]}:{key[3]}"
