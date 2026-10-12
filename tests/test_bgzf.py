"""Focused contracts for shared BGZF framing checks."""

from __future__ import annotations

import gzip
from pathlib import Path
import sys
from tempfile import TemporaryDirectory
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "workflow" / "lib" / "python"))

from bgzf import BGZF_EOF_BLOCK, has_bgzf_eof, is_bgzf  # noqa: E402


def gzip_header_with_extra(extra: bytes) -> bytes:
    return b"\x1f\x8b\x08\x04" + b"\x00" * 6 + len(extra).to_bytes(2, "little") + extra


class BgzfTests(unittest.TestCase):
    def test_detects_a_complete_bgzf_extra_field_and_eof_block(self) -> None:
        with TemporaryDirectory() as temporary:
            path = Path(temporary) / "complete.bgz"
            path.write_bytes(
                gzip_header_with_extra(b"BC\x02\x00\x00\x00") + BGZF_EOF_BLOCK
            )
            self.assertTrue(is_bgzf(path))
            self.assertTrue(has_bgzf_eof(path))

    def test_rejects_plain_gzip_and_truncated_extra_fields(self) -> None:
        with TemporaryDirectory() as temporary:
            directory = Path(temporary)
            plain = directory / "plain.gz"
            plain.write_bytes(gzip.compress(b"not bgzf"))
            truncated = directory / "truncated.bgz"
            truncated.write_bytes(
                b"\x1f\x8b\x08\x04" + b"\x00" * 6 + b"\x06\x00BC\x02\x00"
            )
            self.assertFalse(is_bgzf(plain))
            self.assertFalse(is_bgzf(truncated))
            self.assertFalse(has_bgzf_eof(plain))
            self.assertFalse(has_bgzf_eof(truncated))
