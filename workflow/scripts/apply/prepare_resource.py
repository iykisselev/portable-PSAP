"""Stage one VCF or table resource through injected Snakemake I/O."""

from pathlib import Path
import sys

sys.path.insert(0, str(Path(snakemake.scriptdir).parents[1] / "shared"))
from _bootstrap import bootstrap

bootstrap(snakemake)
from apply_bindings import PreparationInputs, PreparationOutputs, PreparationParameters
from _library import bind


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
    inputs = bind(PreparationInputs, snakemake.input)
    outputs = bind(PreparationOutputs, snakemake.output)
    params = bind(PreparationParameters, snakemake.params)
    source = Path(inputs.source)
    output = Path(outputs.data)
    source_index = Path(inputs.source_index) if inputs.source_index else None
    if params.kind == "vcf":
        prepare_vcf(source, output, source_index)
        return
    chromosome_map = vcf_chromosome_map(Path(inputs.cohort)) if inputs.cohort else None
    prepare_table(
        source, output, params.kind, params.tsv_layout, chromosome_map, source_index
    )


if __name__ == "__main__":
    main()
