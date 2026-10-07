"""Focused probability aggregation contracts independent of Snakemake jobs."""

from __future__ import annotations

import math
import sys
import unittest
from pathlib import Path

import polars as pl


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "workflow" / "lib" / "python"))
sys.path.insert(0, str(ROOT / "tests"))

from generation_probability import role_probability_expression  # noqa: E402
from oracles import distribution  # noqa: E402


class GenerationProbabilityTests(unittest.TestCase):
    def test_additive_statistics_match_the_probability_oracle_case(self) -> None:
        status, event_mass, event_count, survival = distribution(
            {
                3: (0, math.log1p(-0.18), 1),
                5: (0, math.log1p(-0.32), 1),
            },
            14,
        )
        self.assertEqual(status, "ok")
        self.assertEqual(event_count, 2)
        self.assertAlmostEqual(event_mass, 0.4424)
        assert survival is not None
        self.assertAlmostEqual(survival[4], 0.32 / 0.4424)
        self.assertTrue(all(left >= right for left, right in zip(survival, survival[1:])))

    def test_event_models_preserve_probability_bounds_and_hem_identity(self) -> None:
        values = pl.DataFrame({"effective_af": [0.0, 0.2, 0.5, 1.0]}).select(
            role_probability_expression("het").alias("het"),
            role_probability_expression("hom").alias("hom"),
            role_probability_expression("hem_chrX").alias("hem"),
            "effective_af",
        )
        for row in values.iter_rows(named=True):
            self.assertTrue(0 <= row["het"] <= 1)
            self.assertTrue(0 <= row["hom"] <= 1)
            self.assertEqual(row["hem"], row["effective_af"])


if __name__ == "__main__":
    unittest.main()
