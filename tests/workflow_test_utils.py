"""Shared subprocess helpers for workflow-level regression tests."""

from __future__ import annotations

import os
from pathlib import Path
import re
import shutil
import subprocess
from typing import Sequence


ROOT = Path(__file__).resolve().parents[1]
SNAKEMAKE = shutil.which("snakemake")


def workflow_command(
    manifest: str | Path,
    outdir: Path,
    *,
    dry_run: bool = False,
    lint: bool = False,
    targets: Sequence[str] = (),
) -> list[str]:
    """Build one reproducible local Snakemake invocation."""
    if SNAKEMAKE is None:
        raise RuntimeError("Snakemake is not available")
    manifest_path = Path(manifest)
    if not manifest_path.is_absolute():
        manifest_path = ROOT / manifest_path
    command = [
        SNAKEMAKE,
        "--snakefile",
        str(ROOT / "workflow" / "Snakefile"),
        "--configfile",
        str(manifest_path),
        "--config",
        f"outdir={outdir}",
        "--cores",
        "2",
    ]
    if dry_run:
        command.append("--dry-run")
    if lint:
        command.append("--lint")
    command.extend(targets)
    return command


def run_workflow(
    manifest: str | Path,
    outdir: Path,
    *,
    dry_run: bool = False,
    lint: bool = False,
    targets: Sequence[str] = (),
    timeout: int = 180,
) -> subprocess.CompletedProcess[str]:
    """Run the workflow with an isolated Snakemake source cache."""
    cache_dir = outdir.parent / "cache"
    cache_dir.mkdir(parents=True, exist_ok=True)
    environment = os.environ.copy()
    environment["XDG_CACHE_HOME"] = str(cache_dir)
    return subprocess.run(
        workflow_command(
            manifest,
            outdir,
            dry_run=dry_run,
            lint=lint,
            targets=targets,
        ),
        cwd=ROOT,
        env=environment,
        capture_output=True,
        text=True,
        timeout=timeout,
    )


def command_output(result: subprocess.CompletedProcess[str]) -> str:
    """Return both streams for assertions and useful failure messages."""
    return result.stdout + "\n" + result.stderr


def initial_job_count(result: subprocess.CompletedProcess[str]) -> int:
    """Read the final dry-run job-statistics total emitted by Snakemake."""
    totals = re.findall(r"^total\s+(\d+)\s*$", command_output(result), re.MULTILINE)
    if not totals:
        raise ValueError("Snakemake did not emit a Job stats total")
    return int(totals[-1])


def planned_job_count(result: subprocess.CompletedProcess[str]) -> int:
    """Read planned jobs, treating an up-to-date dry-run as zero work."""
    if "Nothing to be done" in command_output(result):
        return 0
    return initial_job_count(result)


def completed_step_count(result: subprocess.CompletedProcess[str]) -> int:
    """Read the completed-step total from a successful execution transcript."""
    steps = re.findall(r"(\d+) of \1 steps \(100%\) done", command_output(result))
    if not steps:
        raise ValueError("Snakemake did not emit a completed-step total")
    return int(steps[-1])
