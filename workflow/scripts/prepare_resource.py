"""Stage one VCF or table resource through injected Snakemake I/O."""

from __future__ import annotations

from pathlib import Path

from apply.prepare_genomic_resource import prepare_table, prepare_vcf, vcf_chromosome_map


def main() -> None:
    source = Path(snakemake.input.source)
    output = Path(snakemake.output.data)
    kind = str(snakemake.params.kind)
    if kind == "vcf":
        prepare_vcf(source, output)
        return
    chromosome_map = vcf_chromosome_map(Path(snakemake.input.cohort[0])) if snakemake.input.cohort else None
    prepare_table(source, output, kind, str(snakemake.params.tsv_layout), chromosome_map)


if __name__ == "__main__":
    main()
