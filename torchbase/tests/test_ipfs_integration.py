"""IPFS round-trip integration tests using the Kubo HTTP API.

Tests are skipped when a Kubo daemon is not reachable at localhost:5001.
When IPFS is available, tests exercise real download, pin, and list behavior.
"""

from pathlib import Path

import pytest
import requests


def _ipfs_available():
    try:
        r = requests.post("http://localhost:5001/api/v0/id", timeout=2)
        return r.status_code == 200
    except Exception:
        return False


ipfs_available = pytest.mark.skipif(
    not _ipfs_available(),
    reason="Kubo IPFS daemon not available at localhost:5001",
)

# A small well-known CID (empty UnixFS directory, stable across Kubo versions)
KNOWN_CID = "QmUNLLsPACCz1vLxQVkXqqLX5R1X345qqfHbsf67hvA3Nn"


class TestIPFSDownload:
    @ipfs_available
    def test_download_torch_returns_path(self):
        """download_torch returns a local Path."""
        from torchbase.torchfs import download_torch

        result = download_torch(KNOWN_CID)
        assert isinstance(result, Path)

    @ipfs_available
    def test_download_torch_path_exists(self):
        """download_torch result path exists on disk after download."""
        from torchbase.torchfs import download_torch

        result = download_torch(KNOWN_CID)
        assert result.exists()

    @ipfs_available
    def test_download_torch_is_idempotent(self):
        """Calling download_torch twice for same CID returns same path."""
        from torchbase.torchfs import download_torch

        p1 = download_torch(KNOWN_CID)
        p2 = download_torch(KNOWN_CID)
        assert p1 == p2

    @ipfs_available
    def test_exists_returns_true_after_download(self):
        """exists() returns True for a CID that has been downloaded."""
        from torchbase.torchfs import download_torch, exists

        download_torch(KNOWN_CID)
        assert exists({}, {"cid": KNOWN_CID}) is True

    @ipfs_available
    def test_exists_returns_false_for_uncached_cid(self):
        """exists() returns False for a CID not in local cache."""
        from torchbase.torchfs import exists

        fake_cid = "QmFakeNotRealCidThatWillNeverExistInCache12345678"
        assert exists({}, {"cid": fake_cid}) is False


class TestIPFSUnavailable:
    def test_download_raises_on_connection_error(self, monkeypatch):
        """download_torch raises when Kubo is unreachable."""
        import torchbase.torchfs as torchfs

        def bad_post(*args, **kwargs):
            raise requests.ConnectionError("refused")

        monkeypatch.setattr(torchfs.requests, "post", bad_post)

        import uuid
        unique_cid = f"QmTest{uuid.uuid4().hex}"
        with pytest.raises(Exception):
            torchfs.download_torch(unique_cid)

    def test_exists_returns_false_for_missing_cid(self):
        """exists() returns False when cache path does not exist."""
        from torchbase.torchfs import exists

        result = exists({}, {"cid": "QmDefinitelyNotCachedXXXX"})
        assert result is False

    def test_exists_returns_false_for_empty_cid(self):
        """exists() returns False when torch_entry has no cid field."""
        from torchbase.torchfs import exists

        assert exists({}, {}) is False
