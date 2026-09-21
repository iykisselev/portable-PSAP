"""Stage one VCF or table resource through injected Snakemake I/O."""

from pathlib import Path
import sys

SHARED_SCRIPT_DIR = Path(snakemake.scriptdir).parents[1] / "shared"
sys.path.insert(0, str(SHARED_SCRIPT_DIR))
from _library import add_library_path

add_library_path(snakemake.scriptdir)

from resources import prepare_table, prepare_vcf, vcf_chromosome_map


def main() -> None:
    """Stage one configured VCF or tabular resource for workflow consumption.

    Reads the selected Snakemake resource and writes its normalized staged copy.
    VCF resources retain their own format; tabular resources may use the cohort
    VCF to harmonize chromosome labels and an optional source index.

    Raises:
        OSError: If a configured input or output cannot be accessed.
        ValueError: If the resource layout is unsupported or malformed.
    """
    source = Path(snakemake.input.source)
    output = Path(snakemake.output.data)
    kind = str(snakemake.params.kind)
    if kind == "vcf":
        source_index = (
            Path(snakemake.input.source_index[0])
            if snakemake.input.source_index
            else None
        )
        prepare_vcf(source, output, source_index)
        return
    # Tabular chromosome labels must agree with the cohort before joins, whereas
    # a staged VCF preserves its native contig metadata.
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
