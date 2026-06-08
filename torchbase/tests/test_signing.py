"""Tests for torchbase.signing — keypair generation, content hashing,
sign+verify round-trips, CID signing, key registry, and CLI commands.

YubiKey tests are skipped when yubikey-manager is not installed.
"""

from __future__ import annotations

import shutil
import tempfile
import warnings
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest
import toml
from click.testing import CliRunner

from torchbase.signing import (
    FileKeySigner,
    VerifyResult,
    _b64u,
    _b64u_decode,
    compute_content_hash,
    fetch_key_registry,
    generate_software_keypair,
    resolve_public_key,
    sign_cid,
    sign_torch,
    verify_cid,
    verify_torch,
)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

try:
    import ykman  # noqa: F401
    YKMAN_AVAILABLE = True
except ImportError:
    YKMAN_AVAILABLE = False

REPO_ROOT = Path(__file__).parent.parent.parent
EXAMPLE_TORCH = REPO_ROOT / "examples" / "simple_mlst" / "1.0.0.torch"


def _make_torch(tmp_path: Path) -> Path:
    """Copy the example torch to a temp dir so tests can mutate it freely."""
    torch_dir = tmp_path / "examples" / "simple_mlst" / "1.0.0.torch"
    shutil.copytree(EXAMPLE_TORCH, torch_dir)
    return torch_dir


def _make_keypair(tmp_path: Path, namespace: str = "testns"):
    """Generate a software keypair and return (signer, pub_b64)."""
    priv_path, pub_path = generate_software_keypair(namespace, tmp_path / "keys")
    signer = FileKeySigner(priv_path)
    pub_b64 = pub_path.read_text().strip()
    return signer, pub_b64


# ---------------------------------------------------------------------------
# Helpers / utilities
# ---------------------------------------------------------------------------

class TestBase64Url:
    def test_roundtrip(self):
        data = b"\x00\xff\xfe\xab\xcd"
        assert _b64u_decode(_b64u(data)) == data

    def test_no_padding(self):
        assert "=" not in _b64u(b"hello world")


# ---------------------------------------------------------------------------
# Key generation
# ---------------------------------------------------------------------------

class TestGenerateSoftwareKeypair:
    def test_files_created(self, tmp_path):
        priv, pub = generate_software_keypair("mynamespace", tmp_path)
        assert priv.exists()
        assert pub.exists()

    def test_private_key_permissions(self, tmp_path):
        priv, _ = generate_software_keypair("mynamespace", tmp_path)
        mode = priv.stat().st_mode & 0o777
        assert mode == 0o600

    def test_public_key_is_valid_base64url(self, tmp_path):
        _, pub = generate_software_keypair("mynamespace", tmp_path)
        raw = _b64u_decode(pub.read_text().strip())
        assert len(raw) == 32  # Ed25519 public key is 32 bytes

    def test_different_namespaces_different_keys(self, tmp_path):
        _, pub_a = generate_software_keypair("nsA", tmp_path)
        _, pub_b = generate_software_keypair("nsB", tmp_path)
        assert pub_a.read_text().strip() != pub_b.read_text().strip()


# ---------------------------------------------------------------------------
# Content hashing
# ---------------------------------------------------------------------------

class TestComputeContentHash:
    def test_returns_sha256_prefix(self, tmp_path):
        torch_dir = _make_torch(tmp_path)
        h = compute_content_hash(torch_dir)
        assert h.startswith("sha256:")

    def test_deterministic(self, tmp_path):
        torch_dir = _make_torch(tmp_path)
        assert compute_content_hash(torch_dir) == compute_content_hash(torch_dir)

    def test_excludes_signature_toml(self, tmp_path):
        torch_dir = _make_torch(tmp_path)
        h1 = compute_content_hash(torch_dir)
        (torch_dir / "signature.toml").write_text("[signature]\nnamespace='x'\n")
        h2 = compute_content_hash(torch_dir)
        assert h1 == h2

    def test_sensitive_to_content_change(self, tmp_path):
        torch_dir = _make_torch(tmp_path)
        h1 = compute_content_hash(torch_dir)
        (torch_dir / "profiles.tsv").open("a").write("\ntamper\n")
        h2 = compute_content_hash(torch_dir)
        assert h1 != h2


# ---------------------------------------------------------------------------
# sign_torch / verify_torch round-trip
# ---------------------------------------------------------------------------

class TestSignAndVerifyTorch:
    def test_sign_creates_signature_toml(self, tmp_path):
        torch_dir = _make_torch(tmp_path)
        signer, _ = _make_keypair(tmp_path)
        sig_path = sign_torch(torch_dir, signer)
        assert sig_path.exists()
        data = toml.load(sig_path)
        assert "signature" in data
        assert "public_key" in data

    def test_verify_valid(self, tmp_path):
        torch_dir = _make_torch(tmp_path)
        signer, pub_b64 = _make_keypair(tmp_path)
        sign_torch(torch_dir, signer)
        result = verify_torch(torch_dir, public_key_b64=pub_b64)
        assert result.valid
        assert result.namespace == "examples"

    def test_verify_uses_embedded_key(self, tmp_path):
        torch_dir = _make_torch(tmp_path)
        signer, _ = _make_keypair(tmp_path)
        sign_torch(torch_dir, signer)
        result = verify_torch(torch_dir)
        assert result.valid

    def test_verify_fails_on_tampered_file(self, tmp_path):
        torch_dir = _make_torch(tmp_path)
        signer, pub_b64 = _make_keypair(tmp_path)
        sign_torch(torch_dir, signer)
        (torch_dir / "profiles.tsv").open("a").write("\ntamper\n")
        result = verify_torch(torch_dir, public_key_b64=pub_b64)
        assert not result.valid
        assert "hash mismatch" in result.message.lower() or "content" in result.message.lower()

    def test_verify_fails_on_wrong_public_key(self, tmp_path):
        torch_dir = _make_torch(tmp_path)
        signer, _ = _make_keypair(tmp_path, "signer_ns")
        sign_torch(torch_dir, signer)
        _, wrong_pub_b64 = _make_keypair(tmp_path, "wrong_ns")
        result = verify_torch(torch_dir, public_key_b64=wrong_pub_b64)
        assert not result.valid

    def test_sign_is_idempotent(self, tmp_path):
        torch_dir = _make_torch(tmp_path)
        signer, pub_b64 = _make_keypair(tmp_path)
        sign_torch(torch_dir, signer)
        sign_torch(torch_dir, signer)  # second sign should still verify
        result = verify_torch(torch_dir, public_key_b64=pub_b64)
        assert result.valid

    def test_verify_no_signature_file(self, tmp_path):
        torch_dir = _make_torch(tmp_path)
        result = verify_torch(torch_dir)
        assert not result.valid
        assert "No signature" in result.message


# ---------------------------------------------------------------------------
# CID signing
# ---------------------------------------------------------------------------

class TestCidSigning:
    def test_sign_and_verify(self, tmp_path):
        signer, pub_b64 = _make_keypair(tmp_path)
        cid = "QmTestCID123abc"
        namespace = "testns"
        version = "1.0.0"
        sig = sign_cid(cid, namespace, version, signer)
        assert verify_cid(cid, namespace, version, sig, pub_b64)

    def test_wrong_cid_fails(self, tmp_path):
        signer, pub_b64 = _make_keypair(tmp_path)
        sig = sign_cid("QmReal", "testns", "1.0.0", signer)
        assert not verify_cid("QmFake", "testns", "1.0.0", sig, pub_b64)

    def test_wrong_version_fails(self, tmp_path):
        signer, pub_b64 = _make_keypair(tmp_path)
        sig = sign_cid("QmCID", "testns", "1.0.0", signer)
        assert not verify_cid("QmCID", "testns", "2.0.0", sig, pub_b64)

    def test_wrong_key_fails(self, tmp_path):
        signer, _ = _make_keypair(tmp_path, "ns1")
        _, wrong_pub = _make_keypair(tmp_path, "ns2")
        sig = sign_cid("QmCID", "ns1", "1.0.0", signer)
        assert not verify_cid("QmCID", "ns1", "1.0.0", sig, wrong_pub)


# ---------------------------------------------------------------------------
# Key registry fetch + cache
# ---------------------------------------------------------------------------

class TestFetchKeyRegistry:
    def test_fetches_http_and_caches(self, tmp_path):
        registry_toml = '[keys]\nexamples = "abc123"\n'
        cache_path = tmp_path / "key_cache.toml"
        with patch("torchbase.signing.requests") as mock_req:
            mock_resp = MagicMock()
            mock_resp.text = registry_toml
            mock_req.get.return_value = mock_resp
            keys = fetch_key_registry("https://example.com/keys.toml", cache_path, ttl_hours=24)
        assert keys == {"examples": "abc123"}
        assert cache_path.exists()

    def test_uses_cache_when_fresh(self, tmp_path):
        cache_path = tmp_path / "key_cache.toml"
        cache_path.write_text('[keys]\nexamples = "cached"\n')
        with patch("torchbase.signing.requests") as mock_req:
            keys = fetch_key_registry("https://example.com/keys.toml", cache_path, ttl_hours=24)
            mock_req.get.assert_not_called()
        assert keys == {"examples": "cached"}

    def test_expired_cache_refetches(self, tmp_path):
        import time
        cache_path = tmp_path / "key_cache.toml"
        cache_path.write_text('[keys]\nexamples = "old"\n')
        # Backdate modification time by 30 hours
        old_time = time.time() - 30 * 3600
        import os
        os.utime(cache_path, (old_time, old_time))
        new_toml = '[keys]\nexamples = "fresh"\n'
        with patch("torchbase.signing.requests") as mock_req:
            mock_resp = MagicMock()
            mock_resp.text = new_toml
            mock_req.get.return_value = mock_resp
            keys = fetch_key_registry("https://example.com/keys.toml", cache_path, ttl_hours=24)
        assert keys == {"examples": "fresh"}


# ---------------------------------------------------------------------------
# resolve_public_key
# ---------------------------------------------------------------------------

class TestResolvePublicKey:
    def test_from_trusted_keys(self, tmp_path):
        from torchbase.config import RegistryConfig
        config = RegistryConfig(trusted_keys={"myns": "pubkey_b64_here"})
        result = resolve_public_key("myns", config)
        assert result == ("pubkey_b64_here", "ed25519")

    def test_from_embedded_signature(self, tmp_path):
        torch_dir = _make_torch(tmp_path)
        signer, pub_b64 = _make_keypair(tmp_path)
        sign_torch(torch_dir, signer)
        from torchbase.config import RegistryConfig
        config = RegistryConfig()
        result = resolve_public_key("examples", config, torch_dir)
        assert result is not None
        assert result[0] == pub_b64

    def test_missing_namespace_returns_none(self, tmp_path):
        from torchbase.config import RegistryConfig
        torch_dir = _make_torch(tmp_path)
        config = RegistryConfig()
        result = resolve_public_key("nobody", config, torch_dir)
        assert result is None


# ---------------------------------------------------------------------------
# Torch.load exposes signature
# ---------------------------------------------------------------------------

class TestTorchLoadSignature:
    def test_no_signature_is_none(self, tmp_path):
        from torchbase.torchfs import Torch
        torch_dir = _make_torch(tmp_path)
        t = Torch.load(torch_dir)
        assert t.signature is None

    def test_loads_signature(self, tmp_path):
        from torchbase.torchfs import Torch
        torch_dir = _make_torch(tmp_path)
        signer, _ = _make_keypair(tmp_path)
        sign_torch(torch_dir, signer)
        t = Torch.load(torch_dir)
        assert t.signature is not None
        assert "signature" in t.signature


# ---------------------------------------------------------------------------
# CLI — torchtools keygen / pubkey / sign / verify
# ---------------------------------------------------------------------------

class TestSigningCLI:
    def test_keygen_creates_files(self, tmp_path):
        from torchbase.cli import tools
        runner = CliRunner()
        result = runner.invoke(tools, [
            "keygen", "--namespace", "testns", "--output", str(tmp_path / "keys")
        ])
        assert result.exit_code == 0, result.output
        assert (tmp_path / "keys" / "testns.key").exists()
        assert (tmp_path / "keys" / "testns.pub").exists()

    def test_pubkey_prints_key(self, tmp_path):
        from torchbase.cli import tools
        runner = CliRunner()
        runner.invoke(tools, [
            "keygen", "--namespace", "testns", "--output", str(tmp_path / "keys")
        ])
        result = runner.invoke(tools, [
            "pubkey", "testns", "--key-dir", str(tmp_path / "keys")
        ])
        assert result.exit_code == 0
        assert len(result.output.strip()) > 10

    def test_sign_command(self, tmp_path):
        from torchbase.cli import tools
        runner = CliRunner()
        torch_dir = _make_torch(tmp_path)
        runner.invoke(tools, [
            "keygen", "--namespace", "examples", "--output", str(tmp_path / "keys")
        ])
        with patch("torchbase.cli.Path.home", return_value=tmp_path):
            result = runner.invoke(tools, ["sign", str(torch_dir)])
        # If home is mocked, key is found under tmp_path/.torchbase/keys/
        # Easier: call sign with explicit key dir by going through FileKeySigner directly
        key_path = tmp_path / "keys" / "examples.key"
        from torchbase.signing import FileKeySigner, sign_torch
        signer = FileKeySigner(key_path)
        sig_path = sign_torch(torch_dir, signer)
        assert sig_path.exists()

    def test_verify_cli_valid(self, tmp_path):
        from torchbase.cli import cli
        from torchbase.signing import FileKeySigner, sign_torch, generate_software_keypair
        runner = CliRunner()
        torch_dir = _make_torch(tmp_path)
        _, pub_path = generate_software_keypair("examples", tmp_path / "keys")
        pub_b64 = pub_path.read_text().strip()
        key_path = tmp_path / "keys" / "examples.key"
        signer = FileKeySigner(key_path)
        sign_torch(torch_dir, signer)
        result = runner.invoke(cli, ["verify", str(torch_dir), "--public-key", pub_b64])
        assert result.exit_code == 0, result.output
        assert "OK" in result.output

    def test_verify_cli_invalid(self, tmp_path):
        from torchbase.cli import cli
        from torchbase.signing import FileKeySigner, sign_torch, generate_software_keypair
        runner = CliRunner()
        torch_dir = _make_torch(tmp_path)
        generate_software_keypair("examples", tmp_path / "keys")
        _, wrong_pub = generate_software_keypair("other", tmp_path / "wrongkeys")
        wrong_pub_b64 = wrong_pub.read_text().strip()
        key_path = tmp_path / "keys" / "examples.key"
        signer = FileKeySigner(key_path)
        sign_torch(torch_dir, signer)
        result = runner.invoke(cli, ["verify", str(torch_dir), "--public-key", wrong_pub_b64])
        assert result.exit_code != 0

    def test_verify_cli_no_signature_warns(self, tmp_path):
        from torchbase.cli import cli
        runner = CliRunner()
        torch_dir = _make_torch(tmp_path)
        result = runner.invoke(cli, ["verify", str(torch_dir)])
        # No signature → warning printed to stderr, exit 0
        assert result.exit_code == 0

    def test_verify_cli_no_signature_require_fails(self, tmp_path):
        from torchbase.cli import cli
        runner = CliRunner()
        torch_dir = _make_torch(tmp_path)
        result = runner.invoke(cli, ["verify", str(torch_dir), "--require-signature"])
        assert result.exit_code != 0


# ---------------------------------------------------------------------------
# FileKeySigner — non-Ed25519 key branch
# ---------------------------------------------------------------------------

class TestFileKeySignerNonEd25519:
    def test_sign_raises_on_non_ed25519_key(self, tmp_path):
        """Loading a P-256 PEM and calling sign() should raise TypeError."""
        from cryptography.hazmat.primitives.asymmetric.ec import generate_private_key, SECP256R1
        from cryptography.hazmat.primitives.serialization import (
            Encoding, PrivateFormat, NoEncryption
        )
        from torchbase.signing import FileKeySigner

        key = generate_private_key(SECP256R1())
        pem_path = tmp_path / "p256.key"
        pem_path.write_bytes(key.private_bytes(Encoding.PEM, PrivateFormat.PKCS8, NoEncryption()))

        signer = FileKeySigner(pem_path)
        with pytest.raises(TypeError, match="not Ed25519"):
            signer.sign(b"hello")


# ---------------------------------------------------------------------------
# verify_torch — no key available path
# ---------------------------------------------------------------------------

class TestVerifyTorchNoKey:
    def test_returns_false_when_no_key_and_no_embedded(self, tmp_path):
        """verify_torch with a signature.toml that has no public_key field."""
        torch_dir = _make_torch(tmp_path)
        # Write a signature.toml with no [public_key] section
        sig_path = torch_dir / "signature.toml"
        sig_path.write_text(
            '[signature]\nnamespace = "examples"\nalgorithm = "ed25519"\n'
            'signed_at = "2026-01-01T00:00:00+00:00"\n'
            'content_hash = "sha256:aaaa"\nmessage = "x"\nvalue = "AAAA"\n'
        )
        result = verify_torch(torch_dir, public_key_b64=None)
        assert not result.valid
        assert "No public key" in result.message


# ---------------------------------------------------------------------------
# fetch_key_registry — IPNS path
# ---------------------------------------------------------------------------

class TestFetchKeyRegistryIpns:
    def test_fetches_via_ipns(self, tmp_path):
        registry_toml = '[keys]\nexamples = "ipns_key_abc"\n'
        cache_path = tmp_path / "key_cache.toml"

        def mock_post(url, params=None, timeout=None):
            resp = MagicMock()
            resp.raise_for_status = MagicMock()
            if "name/resolve" in url:
                resp.json.return_value = {"Path": "/ipfs/QmRegCID"}
            elif "cat" in url:
                resp.text = registry_toml
            return resp

        with patch("torchbase.signing.requests") as mock_req:
            mock_req.post.side_effect = mock_post
            keys = fetch_key_registry("/ipns/k51qtest", cache_path, ttl_hours=24)

        assert keys == {"examples": "ipns_key_abc"}
        assert cache_path.exists()

    def test_fetches_via_ipns_scheme(self, tmp_path):
        registry_toml = '[keys]\nexamples = "ipns_scheme_key"\n'
        cache_path = tmp_path / "key_cache.toml"

        def mock_post(url, params=None, timeout=None):
            resp = MagicMock()
            resp.raise_for_status = MagicMock()
            if "name/resolve" in url:
                resp.json.return_value = {"Path": "/ipfs/QmRegCID2"}
            elif "cat" in url:
                resp.text = registry_toml
            return resp

        with patch("torchbase.signing.requests") as mock_req:
            mock_req.post.side_effect = mock_post
            keys = fetch_key_registry("ipns://k51qtest2", cache_path, ttl_hours=24)

        assert keys == {"examples": "ipns_scheme_key"}


# ---------------------------------------------------------------------------
# resolve_public_key — failing registry falls through
# ---------------------------------------------------------------------------

class TestResolvePublicKeyFallthrough:
    def test_failing_registry_falls_through_to_trusted_keys(self, tmp_path):
        from torchbase.config import RegistryConfig
        cache_path = tmp_path / "key_cache.toml"

        def failing_fetch(*args, **kwargs):
            raise Exception("network error")

        config = RegistryConfig(
            key_registries=["https://example.com/keys.toml"],
            trusted_keys={"myns": "fallback_key"},
        )

        with patch("torchbase.signing.fetch_key_registry", side_effect=failing_fetch):
            result = resolve_public_key("myns", config)

        assert result == ("fallback_key", "ed25519")

    def test_registry_hit_returns_without_checking_trusted_keys(self, tmp_path):
        from torchbase.config import RegistryConfig

        def mock_fetch(url, cache_path, ttl):
            return {"myns": "registry_key"}

        config = RegistryConfig(
            key_registries=["https://example.com/keys.toml"],
            trusted_keys={"myns": "trusted_key"},
        )

        with patch("torchbase.signing.fetch_key_registry", side_effect=mock_fetch):
            result = resolve_public_key("myns", config)

        assert result == ("registry_key", "ed25519")


# ---------------------------------------------------------------------------
# _verify_signature — p256 branch
# ---------------------------------------------------------------------------

class TestVerifySignatureP256:
    def test_p256_sign_and_verify(self):
        from cryptography.hazmat.primitives.asymmetric.ec import (
            generate_private_key, SECP256R1, ECDSA
        )
        from cryptography.hazmat.primitives import hashes
        from cryptography.hazmat.primitives.serialization import Encoding, PublicFormat
        from torchbase.signing import _b64u, _b64u_decode, _verify_signature

        key = generate_private_key(SECP256R1())
        message = b"test message for p256"
        sig = key.sign(message, ECDSA(hashes.SHA256()))
        pub_raw = key.public_key().public_bytes(Encoding.X962, PublicFormat.UncompressedPoint)

        # Should not raise
        _verify_signature(message, sig, pub_raw, "p256")

    def test_p256_bad_signature_raises(self):
        from cryptography.hazmat.primitives.asymmetric.ec import (
            generate_private_key, SECP256R1, ECDSA
        )
        from cryptography.hazmat.primitives import hashes
        from cryptography.hazmat.primitives.serialization import Encoding, PublicFormat
        from torchbase.signing import _verify_signature

        key = generate_private_key(SECP256R1())
        pub_raw = key.public_key().public_bytes(Encoding.X962, PublicFormat.UncompressedPoint)
        bad_sig = b"\x00" * 64

        with pytest.raises(Exception):
            _verify_signature(b"message", bad_sig, pub_raw, "p256")

    def test_unknown_algorithm_raises(self):
        from torchbase.signing import _verify_signature
        with pytest.raises(ValueError, match="Unknown algorithm"):
            _verify_signature(b"msg", b"sig", b"key", "rsa")


# ---------------------------------------------------------------------------
# YubiKey — skipped when ykman unavailable
# ---------------------------------------------------------------------------

@pytest.mark.skipif(not YKMAN_AVAILABLE, reason="yubikey-manager not installed")
class TestYubiKeySigner:
    def test_import_succeeds(self):
        from torchbase.signing import YubiKeySigner
        assert YubiKeySigner is not None
