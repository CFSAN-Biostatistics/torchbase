"""Tests for torchbase.chain — IPLD commit chain primitives."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Dict, List
from unittest.mock import MagicMock, patch

import pytest
import toml

from torchbase.chain import (
    _entries_toml_string,
    block_signature_message,
    block_to_toml,
    get_chain_head,
    make_genesis_block,
    make_update_block,
    reconstruct_manifest,
    verify_block,
    walk_chain,
)
from torchbase.signing import (
    FileKeySigner,
    _b64u,
    _b64u_decode,
    generate_software_keypair,
)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture()
def key_pair(tmp_path):
    priv, pub = generate_software_keypair("testns", tmp_path)
    signer = FileKeySigner(priv)
    pub_b64 = pub.read_text().strip()
    return signer, pub_b64


# ---------------------------------------------------------------------------
# Block construction tests
# ---------------------------------------------------------------------------

class TestGenesisBlock:
    def test_fields_present(self, key_pair):
        signer, _ = key_pair
        block = make_genesis_block("testns", signer)
        assert block["type"] == "genesis"
        assert block["namespace"] == "testns"
        assert "public_key" in block
        assert "timestamp" in block
        assert "signature" in block

    def test_signature_valid(self, key_pair):
        signer, pub_b64 = key_pair
        block = make_genesis_block("testns", signer)
        assert verify_block(block, pub_b64)

    def test_tampered_namespace_fails(self, key_pair):
        signer, pub_b64 = key_pair
        block = make_genesis_block("testns", signer)
        block["namespace"] = "evil"
        assert not verify_block(block, pub_b64)

    def test_tampered_public_key_fails(self, key_pair):
        signer, pub_b64 = key_pair
        block = make_genesis_block("testns", signer)
        block["public_key"] = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA="
        assert not verify_block(block, pub_b64)


class TestUpdateBlock:
    def _entries(self):
        return {
            "testns/seqsero2": {
                "2.0.0": "QmTorchCIDabc123",
                "latest": "QmTorchCIDabc123",
                "signatures": {"2.0.0": "dummysig"},
            }
        }

    def test_fields_present(self, key_pair):
        signer, _ = key_pair
        block = make_update_block("testns", "QmPrev123", self._entries(), signer)
        assert block["type"] == "update"
        assert block["namespace"] == "testns"
        assert block["previous"] == "QmPrev123"
        assert "timestamp" in block
        assert "signature" in block
        # Entries should be in block
        assert "testns/seqsero2" in block

    def test_signature_valid(self, key_pair):
        signer, pub_b64 = key_pair
        block = make_update_block("testns", "QmPrev123", self._entries(), signer)
        assert verify_block(block, pub_b64)

    def test_tampered_previous_fails(self, key_pair):
        signer, pub_b64 = key_pair
        block = make_update_block("testns", "QmPrev123", self._entries(), signer)
        block["previous"] = "QmEvil999"
        assert not verify_block(block, pub_b64)

    def test_tampered_entry_fails(self, key_pair):
        signer, pub_b64 = key_pair
        block = make_update_block("testns", "QmPrev123", self._entries(), signer)
        block["testns/seqsero2"]["2.0.0"] = "QmEvilCID"
        assert not verify_block(block, pub_b64)


class TestBlockSignatureMessage:
    def test_genesis_message_format(self, key_pair):
        signer, _ = key_pair
        block = make_genesis_block("testns", signer)
        msg = block_signature_message(block).decode()
        expected_prefix = f"genesis:testns:{block['public_key']}:{block['timestamp']}"
        assert msg == expected_prefix

    def test_update_message_includes_entries_hash(self, key_pair):
        signer, _ = key_pair
        entries = {"testns/tool": {"1.0.0": "QmABC"}}
        block = make_update_block("testns", "QmPrev", entries, signer)
        msg = block_signature_message(block).decode()
        assert msg.startswith("update:testns:QmPrev:")
        # Last token is a 64-char SHA-256 hex digest (rsplit to avoid timestamp colons)
        sha256_part = msg.rsplit(":", 1)[-1]
        assert len(sha256_part) == 64
        int(sha256_part, 16)  # must be valid hex

    def test_unknown_type_raises(self):
        with pytest.raises(ValueError):
            block_signature_message({"type": "unknown", "namespace": "x"})


# ---------------------------------------------------------------------------
# Serialisation tests
# ---------------------------------------------------------------------------

class TestBlockToToml:
    def test_genesis_round_trip(self, key_pair):
        signer, _ = key_pair
        block = make_genesis_block("testns", signer)
        toml_str = block_to_toml(block)
        parsed = toml.loads(toml_str)
        assert parsed["type"] == "genesis"
        assert parsed["namespace"] == "testns"
        assert parsed["public_key"] == block["public_key"]
        assert parsed["signature"] == block["signature"]

    def test_update_round_trip_preserves_entries(self, key_pair):
        signer, _ = key_pair
        entries = {
            "testns/seqsero2": {
                "2.0.0": "QmTorch",
                "latest": "QmTorch",
            }
        }
        block = make_update_block("testns", "QmPrev", entries, signer)
        toml_str = block_to_toml(block)
        parsed = toml.loads(toml_str)
        assert parsed["type"] == "update"
        assert parsed["testns/seqsero2"]["2.0.0"] == "QmTorch"

    def test_header_fields_come_before_entries(self, key_pair):
        signer, _ = key_pair
        entries = {"testns/tool": {"1.0.0": "QmXYZ"}}
        block = make_update_block("testns", "QmPrev", entries, signer)
        toml_str = block_to_toml(block)
        lines = [l for l in toml_str.splitlines() if l.strip()]
        # First non-empty line should be 'type'
        assert lines[0].startswith("type")


# ---------------------------------------------------------------------------
# Chain walk tests (Kubo mocked)
# ---------------------------------------------------------------------------

class TestWalkChain:
    def _make_chain(self, key_pair, n_updates: int = 2):
        """Build an in-memory chain, return (blocks, cids) genesis-first."""
        signer, pub_b64 = key_pair
        genesis = make_genesis_block("testns", signer)
        genesis_toml = block_to_toml(genesis)
        genesis_cid = "QmGenesis"

        blocks = [genesis]
        cids = [genesis_cid]
        toml_strings = {genesis_cid: genesis_toml}

        prev_cid = genesis_cid
        for i in range(n_updates):
            entries = {f"testns/tool{i}": {f"{i}.0.0": f"QmTorch{i}"}}
            upd = make_update_block("testns", prev_cid, entries, signer)
            upd_toml = block_to_toml(upd)
            upd_cid = f"QmUpdate{i}"
            toml_strings[upd_cid] = upd_toml
            blocks.append(upd)
            cids.append(upd_cid)
            prev_cid = upd_cid

        return blocks, cids, toml_strings

    def _mock_fetch(self, toml_strings):
        def _fetch(cid, node="127.0.0.1", port=5001):
            if cid not in toml_strings:
                raise KeyError(f"CID not found: {cid}")
            return toml.loads(toml_strings[cid])
        return _fetch

    def test_walk_returns_genesis_first(self, key_pair):
        blocks, cids, toml_strings = self._make_chain(key_pair, n_updates=2)
        head_cid = cids[-1]
        with patch("torchbase.chain.fetch_block", self._mock_fetch(toml_strings)):
            chain = walk_chain(head_cid)
        assert chain[0]["type"] == "genesis"
        assert len(chain) == 3

    def test_walk_verifies_all_sigs(self, key_pair):
        blocks, cids, toml_strings = self._make_chain(key_pair, n_updates=1)
        head_cid = cids[-1]
        with patch("torchbase.chain.fetch_block", self._mock_fetch(toml_strings)):
            chain = walk_chain(head_cid)
        assert len(chain) == 2

    def test_walk_detects_tampered_signature(self, key_pair):
        """Tampering a block's signature field is caught during verification."""
        blocks, cids, toml_strings = self._make_chain(key_pair, n_updates=1)
        # Corrupt the update block's signature (keep previous intact so walk completes)
        upd_cid = cids[-1]
        parsed = toml.loads(toml_strings[upd_cid])
        parsed["signature"] = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA="
        toml_strings[upd_cid] = toml.dumps(parsed)
        with patch("torchbase.chain.fetch_block", self._mock_fetch(toml_strings)):
            with pytest.raises(ValueError, match="invalid signature"):
                walk_chain(upd_cid)

    def test_walk_detects_broken_previous_link(self, key_pair):
        """A previous pointer to a missing CID raises ValueError."""
        blocks, cids, toml_strings = self._make_chain(key_pair, n_updates=1)
        # Tamper the previous pointer so it points to a non-existent CID
        upd_cid = cids[-1]
        parsed = toml.loads(toml_strings[upd_cid])
        parsed["previous"] = "QmDoesNotExist"
        toml_strings[upd_cid] = toml.dumps(parsed)
        with patch("torchbase.chain.fetch_block", self._mock_fetch(toml_strings)):
            with pytest.raises(ValueError, match="Could not fetch block"):
                walk_chain(upd_cid)

    def test_walk_detects_missing_previous(self, key_pair):
        signer, _ = key_pair
        # Manually build an update block with no 'previous' field
        upd = {"type": "update", "namespace": "testns", "timestamp": "x", "signature": "sig"}
        toml_strings = {"QmHead": toml.dumps(upd)}
        with patch("torchbase.chain.fetch_block", self._mock_fetch(toml_strings)):
            with pytest.raises(ValueError, match="no 'previous' link"):
                walk_chain("QmHead")

    def test_walk_genesis_only(self, key_pair):
        signer, pub_b64 = key_pair
        genesis = make_genesis_block("testns", signer)
        toml_strings = {"QmGen": block_to_toml(genesis)}
        with patch("torchbase.chain.fetch_block", self._mock_fetch(toml_strings)):
            chain = walk_chain("QmGen")
        assert len(chain) == 1
        assert chain[0]["type"] == "genesis"


# ---------------------------------------------------------------------------
# Manifest reconstruction tests
# ---------------------------------------------------------------------------

class TestReconstructManifest:
    def test_empty_chain(self, key_pair):
        signer, _ = key_pair
        genesis = make_genesis_block("testns", signer)
        manifest = reconstruct_manifest([genesis])
        assert manifest == {}

    def test_single_update(self, key_pair):
        signer, _ = key_pair
        genesis = make_genesis_block("testns", signer)
        entries = {"testns/seqsero2": {"1.0.0": "QmABC", "latest": "QmABC"}}
        upd = make_update_block("testns", "QmGen", entries, signer)
        manifest = reconstruct_manifest([genesis, upd])
        assert manifest["testns/seqsero2"]["1.0.0"] == "QmABC"

    def test_newer_entry_wins(self, key_pair):
        signer, _ = key_pair
        genesis = make_genesis_block("testns", signer)
        e1 = {"testns/seqsero2": {"1.0.0": "QmOld", "latest": "QmOld"}}
        e2 = {"testns/seqsero2": {"2.0.0": "QmNew", "latest": "QmNew"}}
        upd1 = make_update_block("testns", "QmGen", e1, signer)
        upd2 = make_update_block("testns", "QmUpd1", e2, signer)
        manifest = reconstruct_manifest([genesis, upd1, upd2])
        assert manifest["testns/seqsero2"]["1.0.0"] == "QmOld"
        assert manifest["testns/seqsero2"]["2.0.0"] == "QmNew"
        assert manifest["testns/seqsero2"]["latest"] == "QmNew"

    def test_multiple_torches(self, key_pair):
        signer, _ = key_pair
        genesis = make_genesis_block("testns", signer)
        e1 = {"testns/toolA": {"1.0.0": "QmA"}}
        e2 = {"testns/toolB": {"1.0.0": "QmB"}}
        upd1 = make_update_block("testns", "QmGen", e1, signer)
        upd2 = make_update_block("testns", "QmU1", e2, signer)
        manifest = reconstruct_manifest([genesis, upd1, upd2])
        assert "testns/toolA" in manifest
        assert "testns/toolB" in manifest


# ---------------------------------------------------------------------------
# upload_directory multipart format tests
# ---------------------------------------------------------------------------

class TestUploadDirectory:
    def test_multipart_has_one_part_per_file(self, tmp_path):
        from torchbase.chain import upload_directory

        torch_dir = tmp_path / "1.0.0.torch"
        torch_dir.mkdir()
        (torch_dir / "metadata.toml").write_text("[meta]\n")
        resources = torch_dir / "_resources"
        resources.mkdir()
        (resources / "locus.fasta").write_text(">seq\nATCG\n")

        captured = {}

        def mock_post(url, params=None, files=None, timeout=None):
            captured["files"] = files
            resp = MagicMock()
            resp.raise_for_status = MagicMock()
            resp.json.return_value = {"Hash": "QmMockRoot"}
            return resp

        with patch("requests.post", mock_post):
            cid = upload_directory(torch_dir)

        assert cid == "QmMockRoot"
        # Should have exactly 2 file parts
        assert len(captured["files"]) == 2

    def test_relative_path_from_parent(self, tmp_path):
        from torchbase.chain import upload_directory

        torch_dir = tmp_path / "1.0.0.torch"
        torch_dir.mkdir()
        (torch_dir / "profiles.tsv").write_text("ST\t1\n")

        captured_names = []

        def mock_post(url, params=None, files=None, timeout=None):
            for _, (name, _, _) in files:
                captured_names.append(name)
            resp = MagicMock()
            resp.raise_for_status = MagicMock()
            resp.json.return_value = {"Hash": "QmR"}
            return resp

        with patch("requests.post", mock_post):
            upload_directory(torch_dir)

        # All relative paths should start with the torch dir name, not parent
        for name in captured_names:
            assert name.startswith("1.0.0.torch/"), f"Bad rel path: {name}"


# ---------------------------------------------------------------------------
# get_chain_head tests
# ---------------------------------------------------------------------------

class TestGetChainHead:
    def test_resolves_cid(self):
        def mock_post(url, params=None, timeout=None):
            resp = MagicMock()
            resp.raise_for_status = MagicMock()
            resp.json.return_value = {"Path": "/ipfs/QmHead123"}
            return resp

        with patch("requests.post", mock_post):
            cid = get_chain_head("/ipns/k51qtest")
        assert cid == "QmHead123"

    def test_returns_none_on_error(self):
        def mock_post(url, params=None, timeout=None):
            raise Exception("no route to host")

        with patch("requests.post", mock_post):
            result = get_chain_head("/ipns/k51qtest")
        assert result is None


# ---------------------------------------------------------------------------
# RegistryManager chain resolution tests
# ---------------------------------------------------------------------------

class TestRegistryManagerChainResolution:
    def test_resolves_from_chain_when_namespace_configured(self, key_pair, tmp_path):
        from torchbase.registry import RegistryManager
        from torchbase.config import RegistryConfig

        signer, pub_b64 = key_pair
        genesis = make_genesis_block("testns", signer)
        entries = {
            "testns/seqsero2": {"2.0.0": "QmTorchCID", "latest": "QmTorchCID"}
        }
        upd = make_update_block("testns", "QmGenCID", entries, signer)

        chain = [genesis, upd]
        cids = ["QmGenCID", "QmHeadCID"]
        toml_map = dict(zip(cids, [block_to_toml(b) for b in chain]))

        def mock_fetch_block(cid, node="127.0.0.1", port=5001):
            return toml.loads(toml_map[cid])

        def mock_get_head(ipns, node="127.0.0.1", port=5001):
            return "QmHeadCID"

        config = RegistryConfig()
        config.namespaces = {"testns": "/ipns/k51qtest"}

        manager = RegistryManager(config)

        with patch("torchbase.chain.fetch_block", mock_fetch_block):
            with patch("torchbase.chain.get_chain_head", mock_get_head):
                cid = manager.resolve("testns/seqsero2")

        assert cid == "QmTorchCID"

    def test_falls_through_to_http_registry_when_chain_empty(self):
        from torchbase.registry import RegistryManager
        from torchbase.config import RegistryConfig

        config = RegistryConfig()
        config.namespaces = {"testns": "/ipns/k51qtest"}
        config.default_registry = "https://example.com/registry.toml"

        manifest_data = toml.dumps({
            "testns/seqsero2": {"1.0.0": "QmHTTPCID", "latest": "QmHTTPCID"}
        })

        def mock_get_head(ipns, node="127.0.0.1", port=5001):
            return None  # chain not yet published

        def mock_get(url, timeout=None):
            resp = MagicMock()
            resp.raise_for_status = MagicMock()
            resp.text = manifest_data
            return resp

        manager = RegistryManager(config)
        with patch("torchbase.chain.get_chain_head", mock_get_head):
            with patch("requests.get", mock_get):
                cid = manager.resolve("testns/seqsero2")

        assert cid == "QmHTTPCID"
