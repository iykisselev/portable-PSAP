"""Shared indexed-window allele-record reader contracts."""

from __future__ import annotations

from pathlib import Path
import sys
from tempfile import TemporaryDirectory
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "workflow" / "lib" / "python"))

from generation_population import window_records  # noqa: E402


class WindowRecordTests(unittest.TestCase):
    def test_skips_non_data_and_keeps_only_canonical_owner_rows(self) -> None:
        with TemporaryDirectory() as temporary:
            resource = Path(temporary) / "records.tsv"
            resource.write_text(
                "# comment\n\n"
                "chr1\t1\tA\tC\tleft-boundary\n"
                "chr1\t2\tA\tG\towned\n"
                "2\t3\tA\tT\tother-chromosome\n"
                "1\t5\tA\tC\tright-boundary\n"
                "1\t6\tA\tG\tpast-window\n",
                encoding="utf-8",
            )
            records = list(window_records(resource, None, "1", 1, 5, min_columns=5))

        self.assertEqual(
            [
                (4, ["chr1", "2", "A", "G", "owned"], ("1", 2, "A", "G")),
                (6, ["1", "5", "A", "C", "right-boundary"], ("1", 5, "A", "C")),
            ],
            records,
        )

    def test_rejects_short_noncomment_data_rows_with_a_source_location(self) -> None:
        with TemporaryDirectory() as temporary:
            resource = Path(temporary) / "records.tsv"
            resource.write_text("# comment\n1\t2\tA\n", encoding="utf-8")
            with self.assertRaisesRegex(
                ValueError, r"records\.tsv: line 2 needs 5 columns"
            ):
                list(window_records(resource, None, "1", 0, 10, min_columns=5))


if __name__ == "__main__":
    unittest.main()
