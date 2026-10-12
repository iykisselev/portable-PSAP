"""Project configured exclusions from a sites-only filtered BCF stream."""

import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(snakemake.scriptdir).parents[1] / "shared"))
from _bootstrap import bootstrap

bootstrap(snakemake)
from apply_bindings import ExclusionInputs, ExclusionOutputs, ExclusionParameters
from _library import bind


from exclusions import collect_exclusions, project_vcf_records


from execution_tools import require_tool


def main() -> None:
    """Project exclusions from the calibration-filtered BCF coordinate space.

    A sites-only ``bcftools view`` stream avoids reading cohort genotype
    columns. The matched records are written directly as an indexed BGZF VCF.
    """
    inputs = bind(ExclusionInputs, snakemake.input)
    outputs = bind(ExclusionOutputs, snakemake.output)
    params = bind(ExclusionParameters, snakemake.params)
    intervals, variants, indexed_tsv_sources = collect_exclusions(
        params.sources,
        indexed_tsv_paths=params.indexed_tsv_sources,
    )
    bcftools = require_tool("bcftools")
    bgzip = require_tool("bgzip")
    source_path = Path(inputs.variants)
    output_path = Path(outputs.vcf)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    reader = subprocess.Popen(
        [bcftools, "view", "-G", "-Ov", str(source_path)],
        stdout=subprocess.PIPE,
        text=True,
    )
    with output_path.open("wb") as output_handle:
        writer = subprocess.Popen(
            [bgzip, "--threads", str(snakemake.threads), "--stdout"],
            stdin=subprocess.PIPE,
            stdout=output_handle,
            text=True,
        )
        try:
            assert reader.stdout is not None
            assert writer.stdin is not None
            project_vcf_records(
                reader.stdout,
                writer.stdin,
                intervals,
                variants,
                str(source_path),
                indexed_tsv_paths=indexed_tsv_sources,
            )
        finally:
            if reader.stdout is not None:
                reader.stdout.close()
            if writer.stdin is not None:
                writer.stdin.close()
        if reader.wait() != 0:
            raise subprocess.CalledProcessError(reader.returncode, reader.args)
        if writer.wait() != 0:
            raise subprocess.CalledProcessError(writer.returncode, writer.args)
    subprocess.run(
        [bcftools, "index", "--force", "--tbi", str(output_path)],
        check=True,
    )


if __name__ == "__main__":
    main()
