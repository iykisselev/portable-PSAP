#!/usr/bin/env python3
"""Create a BGZF-compressed, indexed workflow copy of a genomic resource."""

from __future__ import annotations

import argparse
import gzip
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path


BGZF_EOF_BLOCK = bytes.fromhex(
    "1f8b08040000000000ff0600424302001b0003000000000000000000"
)

def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--input", required=True)
    result.add_argument("--output", required=True)
    result.add_argument("--kind", choices=("vcf", "bed", "tsv"), required=True)
    result.add_argument("--tsv-layout", choices=("auto", "position"), default="auto")
    result.add_argument(
        "--chromosome-map-vcf",
        help="Rewrite table chromosome names to the conventions in this VCF.",
    )
    return result


def is_gzip(path: Path) -> bool:
    with path.open("rb") as handle:
        return handle.read(2) == b"\x1f\x8b"


def is_bgzip(path: Path) -> bool:
    try:
        with path.open("rb") as handle:
            header = handle.read(18)
            if len(header) < 18 or header[:2] != b"\x1f\x8b" or not header[3] & 4:
                return False
            extra_length = int.from_bytes(header[10:12], "little")
            extra = header[12:] + handle.read(max(0, extra_length - 6))
    except OSError:
        return False

    offset = 0
    while offset + 4 <= len(extra):
        length = int.from_bytes(extra[offset + 2 : offset + 4], "little")
        if extra[offset : offset + 2] == b"BC" and length == 2:
            return True
        offset += 4 + length
    return False


def source_has_index(path: Path) -> bool:
    return Path(f"{path}.tbi").is_file() or Path(f"{path}.csi").is_file()


def has_bgzf_eof(path: Path) -> bool:
    try:
        with path.open("rb") as handle:
            handle.seek(-len(BGZF_EOF_BLOCK), 2)
            return handle.read() == BGZF_EOF_BLOCK
    except OSError:
        return False


def require_tool(name: str) -> str:
    executable = shutil.which(name)
    if executable is None:
        raise RuntimeError(f"Required executable is unavailable: {name}")
    return executable


def run(command: list[str], **kwargs) -> None:
    rendered = " ".join(command)
    print(f"RUN: {rendered}", file=sys.stderr)
    subprocess.run(command, check=True, **kwargs)


def prepare_vcf(source: Path, output: Path) -> None:
    bcftools = require_tool("bcftools")
    run(
        [
            bcftools,
            "sort",
            "--output-type",
            "z",
            "--output-file",
            str(output),
            str(source),
        ]
    )
    run([bcftools, "index", "--force", "--tbi", str(output)])


def chromosome_key(name: str) -> str:
    normalized = name.strip()
    if normalized.lower().startswith("chr"):
        normalized = normalized[3:]
    normalized = normalized.upper()
    return "MT" if normalized == "M" else normalized


def vcf_chromosome_map(path: Path) -> dict[str, str]:
    bcftools = require_tool("bcftools")
    result = subprocess.run(
        [bcftools, "view", "--header-only", str(path)],
        text=True,
        stdout=subprocess.PIPE,
        check=True,
    )
    names: list[str] = []
    prefix = "##contig=<ID="
    for line in result.stdout.splitlines():
        if line.startswith(prefix):
            name = line[len(prefix) :].split(",", 1)[0].split(">", 1)[0]
            if name:
                names.append(name)

    # Indexed record statistics cover VCFs whose headers omit contig metadata.
    if not names:
        stats = subprocess.run(
            [bcftools, "index", "--stats", str(path)],
            text=True,
            stdout=subprocess.PIPE,
            check=True,
        )
        names.extend(
            line.split("\t", 1)[0]
            for line in stats.stdout.splitlines()
            if line.strip()
        )

    chromosome_map: dict[str, str] = {}
    for name in names:
        key = chromosome_key(name)
        previous = chromosome_map.get(key)
        if previous is not None and previous != name:
            raise ValueError(
                f"VCF uses conflicting chromosome aliases {previous!r} and {name!r}"
            )
        chromosome_map[key] = name
    return chromosome_map


def copy_as_bgzip(
    source: Path,
    output: Path,
    chromosome_map: dict[str, str] | None = None,
) -> None:
    bgzip = require_tool("bgzip")
    opener = gzip.open if is_gzip(source) else open
    sort_tool = require_tool("sort")

    with tempfile.TemporaryDirectory(prefix="prepare_genomic_resource-") as temp_dir:
        temp_path = Path(temp_dir)
        header_path = temp_path / "header"
        rows_path = temp_path / "rows"
        sorted_rows_path = temp_path / "rows.sorted"
        rewritten = 0

        with (
            opener(source, "rt", encoding="utf-8", errors="replace") as source_handle,
            header_path.open("w", encoding="utf-8") as header_handle,
            rows_path.open("w", encoding="utf-8") as rows_handle,
        ):
            for line in source_handle:
                if (
                    not line.strip()
                    or line.startswith(("#", "track ", "browser "))
                ):
                    header_handle.write(line)
                    continue

                if chromosome_map:
                    chromosome, separator, remainder = line.partition("\t")
                    replacement = chromosome_map.get(chromosome_key(chromosome))
                    if separator and replacement and replacement != chromosome:
                        line = replacement + separator + remainder
                        rewritten += 1
                rows_handle.write(line)

        run(
            [
                sort_tool,
                "--stable",
                "--field-separator=\t",
                "--key=1,1",
                "--key=2,2n",
                "--output",
                str(sorted_rows_path),
                str(rows_path),
            ]
        )

        with output.open("wb") as output_handle:
            process = subprocess.Popen(
                [bgzip, "--stdout"],
                stdin=subprocess.PIPE,
                stdout=output_handle,
            )
            assert process.stdin is not None
            with header_path.open("rb") as header_handle:
                shutil.copyfileobj(header_handle, process.stdin)
            with sorted_rows_path.open("rb") as rows_handle:
                shutil.copyfileobj(rows_handle, process.stdin)
            process.stdin.close()
            return_code = process.wait()
    if return_code != 0:
        raise RuntimeError(f"bgzip exited with status {return_code}")
    if chromosome_map is not None:
        print(f"Rewrote chromosome names in {rewritten} data row(s)", file=sys.stderr)


def indexed_chromosome_names(path: Path, tabix: str) -> set[str]:
    result = subprocess.run(
        [tabix, "--list-chroms", str(path)],
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    if result.returncode != 0:
        return set()
    return {line.strip() for line in result.stdout.splitlines() if line.strip()}


def chromosome_rewrites_required(
    names: set[str],
    chromosome_map: dict[str, str],
) -> bool:
    return any(
        replacement is not None and replacement != name
        for name in names
        if (replacement := chromosome_map.get(chromosome_key(name))) is not None
    )


def link_compatible_table(source: Path, output: Path) -> None:
    if output.exists() or output.is_symlink():
        output.unlink()
    output.symlink_to(source.resolve())
    print(f"Reused chromosome-compatible indexed resource: {source}", file=sys.stderr)


def copy_indexed_as_bgzip(
    source: Path,
    output: Path,
    chromosome_map: dict[str, str] | None,
) -> None:
    """Rewrite an already sorted/indexed table without materializing or sorting it."""
    bgzip = require_tool("bgzip")
    opener = gzip.open if is_gzip(source) else open
    rewritten = 0

    with (
        opener(source, "rt", encoding="utf-8", errors="replace") as source_handle,
        output.open("wb") as output_handle,
    ):
        process = subprocess.Popen(
            [bgzip, "--stdout"],
            stdin=subprocess.PIPE,
            stdout=output_handle,
        )
        assert process.stdin is not None
        for line in source_handle:
            if chromosome_map and line.strip() and not line.startswith(
                ("#", "track ", "browser ")
            ):
                chromosome, separator, remainder = line.partition("\t")
                replacement = chromosome_map.get(chromosome_key(chromosome))
                if separator and replacement and replacement != chromosome:
                    line = replacement + separator + remainder
                    rewritten += 1
            process.stdin.write(line.encode("utf-8"))
        process.stdin.close()
        return_code = process.wait()

    if return_code != 0:
        raise RuntimeError(f"bgzip exited with status {return_code}")
    if chromosome_map is not None:
        print(f"Rewrote chromosome names in {rewritten} data row(s)", file=sys.stderr)


def first_data_columns(path: Path) -> list[str]:
    opener = gzip.open if is_gzip(path) else open
    with opener(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            return stripped.split("\t")
    return []


def prepare_table(
    source: Path,
    output: Path,
    kind: str,
    tsv_layout: str,
    chromosome_map: dict[str, str] | None = None,
) -> None:
    tabix = require_tool("tabix")
    if kind == "bed":
        index_args = ["--preset", "bed"]
    elif tsv_layout == "position":
        index_args = ["--sequence", "1", "--begin", "2", "--end", "2"]
    else:
        columns = first_data_columns(source)
        interval_layout = False
        if len(columns) >= 3:
            try:
                int(columns[1])
                int(columns[2])
                interval_layout = True
            except ValueError:
                pass
        index_args = [
            "--sequence",
            "1",
            "--begin",
            "2",
            "--end",
            "3" if interval_layout else "2",
        ]

    source_ready = (
        is_bgzip(source) and has_bgzf_eof(source) and source_has_index(source)
    )
    if source_ready:
        source_chromosomes = indexed_chromosome_names(source, tabix)
        compatible = (
            chromosome_map is None
            or (
                bool(source_chromosomes)
                and not chromosome_rewrites_required(
                    source_chromosomes,
                    chromosome_map,
                )
            )
        )
        if compatible:
            link_compatible_table(source, output)
        else:
            copy_indexed_as_bgzip(source, output, chromosome_map)
    else:
        copy_as_bgzip(source, output, chromosome_map)

    run([tabix, "--force", *index_args, str(output)])


def main() -> int:
    args = parser().parse_args()
    source = Path(args.input)
    output = Path(args.output)
    if not source.is_file():
        raise FileNotFoundError(f"Input resource does not exist: {source}")

    reasons = []
    if not is_bgzip(source):
        reasons.append("not BGZF-compressed")
    elif not has_bgzf_eof(source):
        reasons.append("BGZF EOF marker is absent")
    if not source_has_index(source):
        reasons.append("not indexed")
    if args.chromosome_map_vcf:
        reasons.append("chromosome names will be harmonized to the cohort VCF")
    print(f"WARNING: preparing {source}: " + "; ".join(reasons), file=sys.stderr)

    output.parent.mkdir(parents=True, exist_ok=True)
    if args.kind == "vcf":
        prepare_vcf(source, output)
    else:
        chromosome_map = (
            vcf_chromosome_map(Path(args.chromosome_map_vcf))
            if args.chromosome_map_vcf
            else None
        )
        prepare_table(
            source,
            output,
            args.kind,
            args.tsv_layout,
            chromosome_map,
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
