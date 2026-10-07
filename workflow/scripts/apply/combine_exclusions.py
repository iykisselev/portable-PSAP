"""Project configured exclusions from a sites-only filtered BCF stream."""

import shutil
import subprocess
import sys
from pathlib import Path

SHARED_SCRIPT_DIR = Path(snakemake.scriptdir).parents[1] / "shared"
sys.path.insert(0, str(SHARED_SCRIPT_DIR))

from _library import add_library_path
add_library_path(snakemake.scriptdir)

from exclusions import collect_exclusions, project_vcf_records


def require_tool(name: str) -> str:
    executable = shutil.which(name)
    if executable is None:
        raise RuntimeError(f"Required executable is unavailable: {name}")
    return executable


def main() -> None:
    """Project exclusions from the calibration-filtered BCF coordinate space.

    A sites-only ``bcftools view`` stream avoids reading cohort genotype
    columns. The matched records are written directly as an indexed BGZF VCF.
    """
    intervals, variants = collect_exclusions(snakemake.params.sources)
    bcftools = require_tool("bcftools")
    bgzip = require_tool("bgzip")
    source_path = Path(snakemake.input.variants)
    output_path = Path(snakemake.output.vcf)
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
                reader.stdout, writer.stdin, intervals, variants, str(source_path)
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
