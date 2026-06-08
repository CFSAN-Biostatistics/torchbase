"""Shared utilities for torch converters."""

import logging as _logging
from pathlib import Path

_log = _logging.getLogger("torchbase.conversions.fetch")


def fetch_file(url: str, dest: Path, *, progress: bool = True) -> Path:
    """Download url to dest, streaming with an optional progress indicator.

    Returns dest.
    """
    import requests

    dest = Path(dest)
    dest.parent.mkdir(parents=True, exist_ok=True)

    if progress:
        _log.info("Downloading %s ...", dest.name)

    resp = requests.get(url, timeout=120, stream=True)
    resp.raise_for_status()
    with open(dest, "wb") as f:
        for chunk in resp.iter_content(chunk_size=1 << 16):
            f.write(chunk)
    return dest
