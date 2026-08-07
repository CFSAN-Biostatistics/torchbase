"""Shared utilities for torch converters."""

import logging as _logging
from pathlib import Path
from typing import Union

_log = _logging.getLogger("torchbase.conversions.fetch")


def fetch_file(
    url: str,
    dest: Path,
    *,
    progress: bool = True,
    verify: Union[bool, str] = True,
    retries: int = 5,
) -> Path:
    """Download url to dest, streaming with progress, retry, and resume support.

    Downloads to a .part temporary file; renames to dest on completion.
    If a .part file already exists, resumes from the current offset via an
    HTTP Range request.

    Returns dest.
    """
    import requests
    from requests.adapters import HTTPAdapter
    from urllib3.util.retry import Retry

    dest = Path(dest)
    dest.parent.mkdir(parents=True, exist_ok=True)
    part = dest.with_suffix(dest.suffix + ".part")

    session = requests.Session()
    retry_strategy = Retry(
        total=retries,
        backoff_factor=1.5,
        status_forcelist=[429, 500, 502, 503, 504],
        allowed_methods=["GET"],
    )
    adapter = HTTPAdapter(max_retries=retry_strategy)
    session.mount("https://", adapter)
    session.mount("http://", adapter)

    resume_from = part.stat().st_size if part.exists() else 0
    headers = {}
    if resume_from:
        headers["Range"] = f"bytes={resume_from}-"
        _log.debug("Resuming %s from byte %d", dest.name, resume_from)

    resp = session.get(url, timeout=(10, 300), stream=True, verify=verify, headers=headers)

    # 416 = Range Not Satisfiable: offset is past EOF, file is already complete
    if resp.status_code == 416:
        _log.debug("%s already complete (server returned 416)", dest.name)
        part.rename(dest)
        return dest

    resp.raise_for_status()

    total = int(resp.headers.get("content-length", 0)) + resume_from
    mode = "ab" if resume_from else "wb"

    with open(part, mode) as f:
        if progress:
            try:
                from tqdm import tqdm
                bar = tqdm(
                    total=total or None,
                    initial=resume_from,
                    unit="B",
                    unit_scale=True,
                    unit_divisor=1024,
                    desc=dest.name,
                    leave=False,
                )
                for chunk in resp.iter_content(chunk_size=1 << 16):
                    f.write(chunk)
                    bar.update(len(chunk))
                bar.close()
            except ImportError:
                _log.info("Downloading %s ...", dest.name)
                for chunk in resp.iter_content(chunk_size=1 << 16):
                    f.write(chunk)
        else:
            for chunk in resp.iter_content(chunk_size=1 << 16):
                f.write(chunk)

    part.rename(dest)
    return dest
