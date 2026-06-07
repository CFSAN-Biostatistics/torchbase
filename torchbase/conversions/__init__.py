"""Shared utilities for torch converters."""

from pathlib import Path


def fetch_file(url: str, dest: Path, *, progress: bool = True) -> Path:
    """Download url to dest, streaming with an optional progress indicator.

    Returns dest.
    """
    import requests

    dest = Path(dest)
    dest.parent.mkdir(parents=True, exist_ok=True)

    if progress:
        print(f"  Downloading {dest.name} ...", flush=True)

    resp = requests.get(url, timeout=120, stream=True)
    resp.raise_for_status()
    with open(dest, "wb") as f:
        for chunk in resp.iter_content(chunk_size=1 << 16):
            f.write(chunk)
    return dest
