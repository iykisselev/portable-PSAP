"""Contracts for the production generation event-frequency adapter."""

from __future__ import annotations

from pathlib import Path
import sys
import unittest

import polars as pl


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "workflow" / "lib" / "python"))

from generation_event_assembly import (  # noqa: E402
    EventAssemblyError,
    assign_effective_frequencies,
    role_event_frames,
)


class GenerationEventAssemblyTests(unittest.TestCase):
    def events(self) -> pl.DataFrame:
        return pl.DataFrame(
            {
                "testing_unit_id": ["unit", "unit", "unit", "rejected"],
                "CHROM": ["1", "X", "1", "1"],
                "POS": [10, 20, 30, 40],
                "REF": ["A", "C", "T", "G"],
                "ALT": ["AT", "CT", "A", "GA"],
                "variant_class": ["indel", "indel", "snv", "indel"],
                "score": [10.0, 20.0, 30.0, 40.0],
                "allele_state": ["observed", "observed", "absent", "rejected"],
                "af_overall": [0.2, 0.3, None, None],
                "af_female": [None, 0.25, None, None],
                "af_male": [None, 0.4, None, None],
                "mu": [None, None, 0.01, None],
                "mu_indel": [None, None, None, None],
                "unit_bucket": [0, 0, 0, 0],
            }
        )

    def test_uses_sex_specific_observed_and_class_specific_absence_frequencies(
        self,
    ) -> None:
        events = assign_effective_frequencies(
            self.events(),
            unobserved_policy="mutation_rate",
            x_af_policy="sex_specific",
            x_mutation_rate_policy="sex_averaged_proxy",
            roles=("het", "hom", "hem_chrX"),
        )
        by_position = {row["POS"]: row for row in events.iter_rows(named=True)}
        self.assertEqual(3, events.height)
        self.assertEqual(
            (
                by_position[20]["effective_af_standard"],
                by_position[20]["effective_af_hem"],
            ),
            (0.25, 0.4),
        )
        self.assertEqual(
            (
                by_position[30]["effective_af_standard"],
                by_position[30]["assignment_source"],
            ),
            (0.01, "absent_mutation_rate"),
        )
        roles = role_event_frames(events, ("het", "hom", "hem_chrX"))
        self.assertEqual(3, roles["het"].height)
        self.assertEqual(1, roles["hem_chrX"].height)
        self.assertEqual([0.4], roles["hem_chrX"]["effective_af"].to_list())

    def test_rejects_absent_x_mutation_rate_fallback_without_proxy_policy(self) -> None:
        events = (
            self.events()
            .with_row_index("row")
            .with_columns(
                pl.when(pl.col("row") == 2)
                .then(pl.lit("X"))
                .otherwise(pl.col("CHROM"))
                .alias("CHROM")
            )
            .drop("row")
        )
        with self.assertRaisesRegex(EventAssemblyError, "absent chrX alleles"):
            assign_effective_frequencies(
                events,
                unobserved_policy="mutation_rate",
                x_af_policy="overall_proxy",
                x_mutation_rate_policy="reject",
                roles=("het", "hom", "hem_chrX"),
            )

    def test_missing_male_af_is_required_only_when_hem_is_enabled(self):
        events = self.events().with_columns(
            pl.lit(None, dtype=pl.Float64).alias("af_male")
        )
        settings = dict(
            unobserved_policy="zero",
            x_af_policy="sex_specific",
            x_mutation_rate_policy="reject",
        )
        assembled = assign_effective_frequencies(
            events, roles=("het", "hom"), **settings
        )
        self.assertEqual([0.2, 0.25, 0.0], assembled["effective_af_standard"].to_list())
        with self.assertRaisesRegex(EventAssemblyError, "missing or invalid"):
            assign_effective_frequencies(
                events, roles=("het", "hom", "hem_chrX"), **settings
            )


if __name__ == "__main__":
    unittest.main()
