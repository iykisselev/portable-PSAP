"""Restart coverage for independently recoverable generated topology artifacts."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import tempfile
import unittest

from workflow_test_utils import SNAKEMAKE, command_output, run_workflow


GENERATE_MANIFEST = "tests/null_distribution/generate.run.yaml"
NO_HEM_MANIFEST = "tests/null_distribution/generate_no_hem.run.yaml"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def bundle_digests(bundle: Path) -> dict[str, str]:
    return {
        path.name: sha256(path) for path in sorted(bundle.iterdir()) if path.is_file()
    }


@unittest.skipUnless(SNAKEMAKE, "Snakemake is required for generation recovery tests")
class GenerationRecoveryTests(unittest.TestCase):
    def run_workflow_assert(
        self, manifest: str, outdir: Path, *, dry_run: bool = False
    ):
        result = run_workflow(manifest, outdir, dry_run=dry_run)
        output = command_output(result)
        self.assertEqual(0, result.returncode, output[-6000:])
        return output

    def test_missing_window_state_statistics_shard_and_lookup_bucket_recover(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory(
            prefix="portable-psap-generation-recovery-"
        ) as temporary:
            outdir = Path(temporary) / "output"
            self.run_workflow_assert(GENERATE_MANIFEST, outdir)

            generation = outdir / "generate" / "fixture"
            bundle = generation / "bundle"
            expected_bundle = bundle_digests(bundle)
            prepared = generation / "scalable" / "prepared"
            windows = json.loads(
                (prepared / "windows.json").read_text(encoding="utf-8")
            )
            buckets = json.loads(
                (prepared / "unit_buckets.json").read_text(encoding="utf-8")
            )
            window_ids = [row["window_id"] for row in windows["windows"]]
            first_bucket = buckets["buckets"][0]
            bucket_id = first_bucket["bucket_id"]
            statistics_window_id = first_bucket["window_ids"][0]
            self.assertGreaterEqual(len(window_ids), 2)

            state = (
                generation
                / "scalable"
                / "windows"
                / f"{window_ids[0]}.allele_states.parquet"
            )
            state.unlink()
            state_plan = self.run_workflow_assert(
                GENERATE_MANIFEST, outdir, dry_run=True
            )
            self.assertIn("rule prepare_generation_window:", state_plan)
            self.assertEqual(1, state_plan.count("rule prepare_generation_window:"))
            self.assertIn(str(state), state_plan)
            self.assertNotIn(f"wildcards: window={window_ids[1]}", state_plan)
            self.run_workflow_assert(GENERATE_MANIFEST, outdir)
            self.assertEqual(expected_bundle, bundle_digests(bundle))

            statistics = (
                generation
                / "scalable"
                / "calculated"
                / "window_statistics"
                / f"{statistics_window_id}.bucket-{bucket_id}.parquet"
            )
            statistics.unlink()
            statistics_plan = self.run_workflow_assert(
                GENERATE_MANIFEST, outdir, dry_run=True
            )
            self.assertIn("rule calculate_window_null_distribution:", statistics_plan)
            self.assertNotIn("rule prepare_generation_window:", statistics_plan)
            self.run_workflow_assert(GENERATE_MANIFEST, outdir)
            self.assertEqual(expected_bundle, bundle_digests(bundle))

            lookup = (
                generation
                / "scalable"
                / "calculated"
                / "lookup_buckets"
                / f"{bucket_id}.cds.DOM-het.tsv"
            )
            lookup.unlink()
            lookup_plan = self.run_workflow_assert(
                GENERATE_MANIFEST, outdir, dry_run=True
            )
            self.assertIn("rule reduce_null_distribution_bucket:", lookup_plan)
            self.assertNotIn("rule calculate_window_null_distribution:", lookup_plan)
            self.assertNotIn("rule prepare_generation_window:", lookup_plan)
            self.run_workflow_assert(GENERATE_MANIFEST, outdir)
            self.assertEqual(expected_bundle, bundle_digests(bundle))

    def test_hem_disabled_bundle_omits_hem_work_products(self) -> None:
        with tempfile.TemporaryDirectory(
            prefix="portable-psap-generation-no-hem-"
        ) as temporary:
            outdir = Path(temporary) / "output"
            self.run_workflow_assert(NO_HEM_MANIFEST, outdir)

            bundle = outdir / "generate" / "fixture-no-hem" / "bundle"
            metadata = json.loads(
                (bundle / "metadata.json").read_text(encoding="utf-8")
            )
            lookup_buckets = (
                outdir
                / "generate"
                / "fixture-no-hem"
                / "scalable"
                / "calculated"
                / "lookup_buckets"
            )
            self.assertNotIn("hem_chrX", metadata["model_assumptions"]["roles"])
            self.assertNotIn("lookup.cds.X-linked-hem", metadata["components"])
            self.assertFalse((bundle / "lookup.cds.X-linked-hem.tsv.gz").exists())
            self.assertFalse(any(lookup_buckets.glob("*.X-linked-hem.tsv")))


if __name__ == "__main__":
    unittest.main()
