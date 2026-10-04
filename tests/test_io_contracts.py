"""Focused checks for atomic generation artifact publication helpers."""

from __future__ import annotations

import gzip
import sys
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "workflow" / "lib" / "python"))

from io_contracts import (  # noqa: E402
    component_inventory,
    write_gzip_tsv_atomic,
    write_text_atomic,
)


class IOContractTests(unittest.TestCase):
    def test_text_writer_returns_destination_and_controls_parent_creation(self) -> None:
        with TemporaryDirectory() as temporary:
            root = Path(temporary)
            destination = root / "created" / "artifact.txt"
            self.assertEqual(
                destination,
                write_text_atomic(destination, "complete\n"),
            )
            self.assertEqual("complete\n", destination.read_text(encoding="utf-8"))

            missing_parent = root / "missing" / "artifact.txt"
            with self.assertRaises(FileNotFoundError):
                write_text_atomic(missing_parent, "not written", create_parents=False)

    def test_gzip_component_is_deterministic_and_has_a_complete_inventory(self) -> None:
        with TemporaryDirectory() as temporary:
            path = Path(temporary) / "component.tsv.gz"
            count, first = write_gzip_tsv_atomic(path, ["id", "value"], [("a", 1), ("b", 2)])
            first_bytes = path.read_bytes()
            _, second = write_gzip_tsv_atomic(path, ["id", "value"], [("a", 1), ("b", 2)])
            self.assertEqual(count, 2)
            self.assertEqual(first_bytes, path.read_bytes())
            self.assertEqual(first, second)
            self.assertEqual(first, component_inventory(path))
            with gzip.open(path, "rt", encoding="utf-8") as handle:
                self.assertEqual(handle.read(), "id\tvalue\na\t1\nb\t2\n")


if __name__ == "__main__":
    unittest.main()
