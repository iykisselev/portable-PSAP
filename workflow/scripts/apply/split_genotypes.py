"""Split non-reference genotypes into one compact VCF per planned sample."""

import shutil
import subprocess
import tempfile
from pathlib import Path


def require_tool(name: str) -> str:
    """Return a required executable with an actionable failure."""
    executable = shutil.which(name)
    if executable is None:
        raise RuntimeError(f"Required executable is unavailable: {name}")
    return executable


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
    samples = list(snakemake.params.samples)
    outputs = [Path(path) for path in snakemake.output.vcf]
    indices = [Path(path) for path in snakemake.output.csi]
    if len(samples) != len(outputs) or len(samples) != len(indices):
        raise RuntimeError("Split sample map and declared VCF outputs differ in length")

    if not outputs:
        return
    output_dir = outputs[0].parent
    output_dir.mkdir(parents=True, exist_ok=True)
    bcftools = require_tool("bcftools")

    with tempfile.TemporaryDirectory(prefix="psap-split-genotypes-") as temporary:
        sample_map = Path(temporary) / "samples.tsv"
        write_sample_map(sample_map, samples)
        command = [
            bcftools,
            "+split",
            str(snakemake.input.vcf),
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

    missing = [str(path) for path in (*outputs, *indices) if not path.is_file()]
    if missing:
        raise RuntimeError(
            "bcftools +split did not create declared output(s): " + ", ".join(missing)
        )


if __name__ == "__main__":
    main()
