"""Split non-reference genotypes into one compact VCF per planned sample."""

import subprocess
import tempfile
from pathlib import Path


import sys

sys.path.insert(0, str(Path(snakemake.scriptdir).parents[1] / "shared"))
from _bootstrap import bootstrap

bootstrap(snakemake)
from apply_bindings import SplitInputs, SplitOutputs, SplitParameters
from _library import bind


from execution_tools import require_tool


def write_sample_map(path: Path, samples: list[str]) -> None:
    """Write the ``+split`` map, retaining names and choosing output bases.

    The second field keeps each source sample name.  Supplying the third field
    avoids the plugin's default filename sanitization and makes its outputs
    match the declared Snakemake paths.
    """
    path.write_text(
        "".join(f"{sample}\t-\t{sample}\n" for sample in samples),
        encoding="utf-8",
    )


def main() -> None:
    """Run one genotype-aware cohort pass and verify every declared output."""
    inputs = bind(SplitInputs, snakemake.input)
    outputs = bind(SplitOutputs, snakemake.output)
    params = bind(SplitParameters, snakemake.params)
    samples = list(params.samples)
    vcf_outputs = [Path(path) for path in outputs.vcf]
    indices = [Path(path) for path in outputs.csi]
    if len(samples) != len(vcf_outputs) or len(samples) != len(indices):
        raise RuntimeError("Split sample map and declared VCF outputs differ in length")

    if not vcf_outputs:
        return
    output_dir = vcf_outputs[0].parent
    output_dir.mkdir(parents=True, exist_ok=True)
    bcftools = require_tool("bcftools")

    with tempfile.TemporaryDirectory(prefix="psap-split-genotypes-") as temporary:
        sample_map = Path(temporary) / "samples.tsv"
        write_sample_map(sample_map, samples)
        command = [
            bcftools,
            "+split",
            str(inputs.vcf),
            "-S",
            str(sample_map),
            "-k",
            "INFO,FMT/GT",
            "-i",
            'GT="alt"',
            "-Oz",
            "--write-index",
            "-o",
            str(output_dir),
        ]
        with Path(snakemake.log[0]).open("wt", encoding="utf-8") as log:
            subprocess.run(command, check=True, stdout=log, stderr=subprocess.STDOUT)

    missing = [str(path) for path in (*vcf_outputs, *indices) if not path.is_file()]
    if missing:
        raise RuntimeError(
            "bcftools +split did not create declared output(s): " + ", ".join(missing)
        )


if __name__ == "__main__":
    main()
