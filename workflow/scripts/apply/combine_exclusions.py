"""Project configured exclusions through injected Snakemake I/O."""

from pathlib import Path
import sys

SHARED_SCRIPT_DIR = Path(snakemake.scriptdir).parents[1] / "shared"
sys.path.insert(0, str(SHARED_SCRIPT_DIR))

from _library import add_library_path
add_library_path(snakemake.scriptdir)

from exclusions import collect_exclusions, write_projected_vcf


def main() -> None:
    """Project configured exclusion sources onto the cohort VCF coordinate space.

    Reads source specifications and the cohort VCF from Snakemake, then writes
    the merged projected VCF used by downstream exclusion filtering. Invalid
    source records propagate the validation errors raised by ``exclusions``.
    """
    intervals, variants = collect_exclusions(snakemake.params.sources)
    write_projected_vcf(Path(snakemake.input.cohort), Path(snakemake.output.vcf), intervals, variants)


if __name__ == "__main__":
    main()
