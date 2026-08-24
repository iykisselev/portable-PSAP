"""Stage one VCF or table resource through injected Snakemake I/O."""

from pathlib import Path

from _library import add_library_path

add_library_path(snakemake.scriptdir)

from resources import prepare_table, prepare_vcf, vcf_chromosome_map


def main() -> None:
    source = Path(snakemake.input.source)
    output = Path(snakemake.output.data)
    kind = str(snakemake.params.kind)
    if kind == "vcf":
        prepare_vcf(source, output)
        return
    chromosome_map = vcf_chromosome_map(Path(snakemake.input.cohort[0])) if snakemake.input.cohort else None
    source_index = (
        Path(snakemake.input.source_index[0])
        if snakemake.input.source_index
        else None
    )
    prepare_table(
        source,
        output,
        kind,
        str(snakemake.params.tsv_layout),
        chromosome_map,
        source_index,
    )


if __name__ == "__main__":
    main()
