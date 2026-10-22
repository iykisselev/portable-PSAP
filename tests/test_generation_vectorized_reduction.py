"""Property coverage for the vectorized lookup reducer against its reference."""

from __future__ import annotations

from collections import defaultdict
import math
from pathlib import Path
import random
import sys
from tempfile import TemporaryDirectory
import unittest

import polars as pl


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "workflow" / "lib" / "python"))
sys.path.insert(0, str(ROOT / "tests"))

from generation_probability import (  # noqa: E402
    aggregate_statistics_vectorized,
    vectorized_lookup_rows,
)
from oracles import distribution, merge_statistics  # noqa: E402


STATISTICS_COLUMNS = [
    "testing_unit_id",
    "partition_id",
    "role",
    "bin",
    "certain_count",
    "log_none_sum",
    "event_count",
    "snv_event_count",
    "indel_event_count",
]


def reference_rows(
    paths: list[Path], units: list[str], roles: tuple[str, ...], bin_count: int
):
    """Render the pre-vectorization reduction for property comparison."""
    grouped = merge_statistics(paths)
    lookups: dict[str, dict[str, list[float]]] = defaultdict(dict)
    diagnostics = []
    for unit in units:
        for role in roles:
            values = grouped[(unit, role)]
            status, event_mass, event_count, survival = distribution(
                {bin_number: entry[:3] for bin_number, entry in values.items()},
                bin_count,
            )
            diagnostics.append(
                (
                    unit,
                    role,
                    event_mass,
                    status,
                    event_count,
                    sum(entry[3] for entry in values.values()),
                    sum(entry[4] for entry in values.values()),
                )
            )
            if survival is not None:
                lookups[role][unit] = survival
    return lookups, diagnostics


class VectorizedReductionTests(unittest.TestCase):
    def test_seeded_random_window_statistics_match_the_reference_within_one_ulp(
        self,
    ) -> None:
        units = [f"UNIT_{index:02d}" for index in range(6)]
        roles = ("het", "hom", "hem_chrX")
        bin_count = 17
        with TemporaryDirectory() as temporary:
            directory = Path(temporary)
            for seed in range(10):
                randomizer = random.Random(seed)
                paths = []
                for shard in range(5):
                    rows = []
                    for unit in units:
                        for role in roles:
                            for bin_number in range(1, bin_count + 1):
                                if randomizer.random() >= 0.35:
                                    continue
                                event_count = randomizer.randint(1, 12)
                                snv_event_count = randomizer.randint(0, event_count)
                                certain_count = int(randomizer.random() < 0.08)
                                log_none_sum = (
                                    0.0
                                    if certain_count
                                    else math.fsum(
                                        math.log1p(-randomizer.uniform(1e-8, 0.5))
                                        for _ in range(event_count)
                                    )
                                )
                                rows.append(
                                    (
                                        unit,
                                        "cds",
                                        role,
                                        bin_number,
                                        certain_count,
                                        log_none_sum,
                                        event_count,
                                        snv_event_count,
                                        event_count - snv_event_count,
                                    )
                                )
                    path = directory / f"seed-{seed}-shard-{shard}.parquet"
                    pl.DataFrame(
                        rows, schema=STATISTICS_COLUMNS, orient="row"
                    ).write_parquet(path)
                    paths.append(path)

                expected_lookups, expected_diagnostics = reference_rows(
                    paths, units, roles, bin_count
                )
                actual_lookups, actual_diagnostics = vectorized_lookup_rows(
                    aggregate_statistics_vectorized(paths).drop("partition_id"), units, roles, bin_count
                )
                for expected, actual in zip(
                    expected_diagnostics, actual_diagnostics, strict=True
                ):
                    self.assertEqual(expected[:2], actual[:2])
                    self.assertEqual(expected[3:], actual[3:])
                    self.assertLessEqual(
                        abs(expected[2] - actual[2]),
                        max(math.ulp(expected[2]), math.ulp(actual[2])),
                    )
                for role in roles:
                    actual_by_unit = {
                        str(row[0]): row[1:] for row in actual_lookups[role]
                    }
                    self.assertEqual(set(expected_lookups[role]), set(actual_by_unit))
                    for unit, expected_survival in expected_lookups[role].items():
                        for expected, actual in zip(
                            expected_survival, actual_by_unit[unit], strict=True
                        ):
                            self.assertLessEqual(
                                abs(expected - actual),
                                max(math.ulp(expected), math.ulp(actual)),
                            )


if __name__ == "__main__":
    unittest.main()
