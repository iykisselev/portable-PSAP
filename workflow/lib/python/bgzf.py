"""Small, shared checks for BGZF framing used by validation and staging."""

from __future__ import annotations

from pathlib import Path


BGZF_EOF_BLOCK = bytes.fromhex(
    "1f8b08040000000000ff0600424302001b0003000000000000000000"
)


def is_bgzf(path: str | Path) -> bool:
    """Return whether a file has a complete gzip ``BC`` extra subfield."""
    try:
        with Path(path).open("rb") as handle:
            header = handle.read(12)
            if len(header) != 12 or header[:2] != b"\x1f\x8b" or not header[3] & 0x04:
                return False
            extra_length = int.from_bytes(header[10:12], byteorder="little")
            extra = handle.read(extra_length)
    except OSError:
        return False

    if len(extra) != extra_length:
        return False
    offset = 0
    while offset + 4 <= len(extra):
        subfield_length = int.from_bytes(extra[offset + 2 : offset + 4], "little")
        subfield_end = offset + 4 + subfield_length
        if subfield_end > len(extra):
            return False
        if extra[offset : offset + 2] == b"BC" and subfield_length == 2:
            return True
        offset = subfield_end
    return False


def has_bgzf_eof(path: str | Path) -> bool:
    """Return whether a file ends with BGZF's mandatory empty-block marker."""
    try:
        with Path(path).open("rb") as handle:
            handle.seek(-len(BGZF_EOF_BLOCK), 2)
            return handle.read() == BGZF_EOF_BLOCK
    except OSError:
        return False
