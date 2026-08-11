"""Torch cryptographic signing and verification.

Two signing backends:
  FileKeySigner  — Ed25519 private key stored as a PEM file
  YubiKeySigner  — signing performed on-device via YubiKey PIV (key never exported)

The content hash covers every file in the torch directory except signature.toml
itself, so signing is idempotent and the signature file can be added without
changing what was signed.
"""

from __future__ import annotations

import base64
import hashlib
import os
import warnings
from collections import namedtuple
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Optional, TYPE_CHECKING

import requests
import toml

if TYPE_CHECKING:
    pass


# ---------------------------------------------------------------------------
# Result type
# ---------------------------------------------------------------------------

VerifyResult = namedtuple("VerifyResult", ["valid", "namespace", "algorithm", "message"])


# ---------------------------------------------------------------------------
# Signer protocol (structural, no ABC overhead)
# ---------------------------------------------------------------------------

class FileKeySigner:
    """Ed25519 signer backed by a PEM private key file."""

    algorithm = "ed25519"

    def __init__(self, key_path: Path):
        from cryptography.hazmat.primitives.serialization import load_pem_private_key
        key_path = Path(key_path)
        raw = key_path.read_bytes()
        self._private_key = load_pem_private_key(raw, password=None)

    def sign(self, message: bytes) -> bytes:
        from cryptography.hazmat.primitives.asymmetric.ed25519 import Ed25519PrivateKey
        if not isinstance(self._private_key, Ed25519PrivateKey):
            raise TypeError("Key is not Ed25519")
        return self._private_key.sign(message)

    def public_key_bytes(self) -> bytes:
        from cryptography.hazmat.primitives.serialization import (
            Encoding, PublicFormat
        )
        return self._private_key.public_key().public_bytes(
            Encoding.Raw, PublicFormat.Raw
        )


class YubiKeySigner:
    """On-device signer via YubiKey PIV. Requires yubikey-manager."""

    def __init__(self, slot: str = "9c", pin: Optional[str] = None):
        try:
            from ykman.device import connect_to_device
            from yubikit.piv import PivSession, SLOT
        except ImportError:
            raise ImportError(
                "yubikey-manager is required for YubiKey support. "
                "Install with: pip install 'torchbase[yubikey]'"
            )
        self._slot_str = slot
        self._pin = pin
        self._slot_const = _piv_slot(slot)
        self._connect = connect_to_device
        self._PivSession = PivSession
        self.algorithm = self._detect_algorithm()

    def _session(self):
        conn, _, _ = self._connect(None)
        return conn, self._PivSession(conn)

    def _detect_algorithm(self) -> str:
        conn, piv = self._session()
        try:
            meta = piv.get_slot_metadata(self._slot_const)
            kt = str(meta.key_type)
            return "ed25519" if "ED25519" in kt.upper() else "p256"
        finally:
            conn.close()

    def sign(self, message: bytes) -> bytes:
        from yubikit.piv import KEY_TYPE
        conn, piv = self._session()
        try:
            pin = self._pin
            if pin is None:
                import click
                pin = click.prompt("YubiKey PIV PIN", hide_input=True)
            piv.verify_pin(pin)
            if self.algorithm == "ed25519":
                from cryptography.hazmat.primitives.asymmetric.ed25519 import Ed25519PublicKey
                return piv.sign(self._slot_const, KEY_TYPE.ED25519, message)
            else:
                from cryptography.hazmat.primitives.asymmetric.ec import ECDSA
                from cryptography.hazmat.primitives import hashes
                return piv.sign(
                    self._slot_const, KEY_TYPE.ECCP256, message,
                    ECDSA(hashes.SHA256())
                )
        finally:
            conn.close()

    def public_key_bytes(self) -> bytes:
        from cryptography.hazmat.primitives.serialization import Encoding, PublicFormat
        conn, piv = self._session()
        try:
            pub = piv.get_slot_metadata(self._slot_const).public_key
            if self.algorithm == "ed25519":
                return pub.public_bytes(Encoding.Raw, PublicFormat.Raw)
            else:
                return pub.public_bytes(Encoding.X962, PublicFormat.UncompressedPoint)
        finally:
            conn.close()


def _piv_slot(slot_str: str):
    from yubikit.piv import SLOT
    return {
        "9a": SLOT.AUTHENTICATION,
        "9c": SLOT.SIGNATURE,
        "9d": SLOT.KEY_MANAGEMENT,
        "9e": SLOT.CARD_AUTH,
    }[slot_str]


# ---------------------------------------------------------------------------
# Key generation
# ---------------------------------------------------------------------------

def generate_software_keypair(namespace: str, output_dir: Path) -> tuple:
    """Generate an Ed25519 keypair for a namespace.

    Returns (private_key_path, public_key_path).
    Private key written with mode 0600.
    """
    from cryptography.hazmat.primitives.asymmetric.ed25519 import Ed25519PrivateKey
    from cryptography.hazmat.primitives.serialization import (
        Encoding, PrivateFormat, PublicFormat, NoEncryption
    )
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    private_key = Ed25519PrivateKey.generate()

    priv_path = output_dir / f"{namespace}.key"
    priv_path.write_bytes(
        private_key.private_bytes(Encoding.PEM, PrivateFormat.PKCS8, NoEncryption())
    )
    os.chmod(priv_path, 0o600)

    pub_bytes = private_key.public_key().public_bytes(Encoding.Raw, PublicFormat.Raw)
    pub_path = output_dir / f"{namespace}.pub"
    pub_path.write_text(_b64u(pub_bytes) + "\n", encoding="utf-8")

    return priv_path, pub_path


def setup_yubikey_slot(
    slot: str,
    namespace: str,
    output_dir: Path,
    pin: Optional[str] = None,
) -> Path:
    """Generate a key on a YubiKey PIV slot and export the public key.

    Returns path to the written <namespace>.pub file.
    """
    try:
        from ykman.device import connect_to_device
        from yubikit.piv import PivSession, KEY_TYPE, MANAGEMENT_KEY_TYPE
    except ImportError:
        raise ImportError(
            "yubikey-manager is required. Install with: pip install 'torchbase[yubikey]'"
        )
    from cryptography.hazmat.primitives.serialization import Encoding, PublicFormat

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    slot_const = _piv_slot(slot)
    conn, _, _ = connect_to_device(None)
    piv = PivSession(conn)
    try:
        if pin is None:
            import click
            pin = click.prompt("YubiKey PIV PIN", hide_input=True)
        piv.verify_pin(pin)
        pub = piv.generate_key(slot_const, KEY_TYPE.ECCP256)
        pub_bytes = pub.public_bytes(Encoding.X962, PublicFormat.UncompressedPoint)
    finally:
        conn.close()

    pub_path = output_dir / f"{namespace}.pub"
    pub_path.write_text(_b64u(pub_bytes) + "\n", encoding="utf-8")
    return pub_path


# ---------------------------------------------------------------------------
# Content hashing
# ---------------------------------------------------------------------------

def compute_content_hash(torch_path: Path) -> str:
    """Canonical SHA-256 of all torch files except signature.toml.

    Returns "sha256:<hex>".
    """
    torch_path = Path(torch_path)
    lines = []
    for f in sorted(torch_path.rglob("*")):
        if not f.is_file():
            continue
        rel = f.relative_to(torch_path)
        if rel.parts[0] == "signature.toml":
            continue
        file_hash = hashlib.sha256(f.read_bytes()).hexdigest()
        lines.append(f"{file_hash}  {rel.as_posix()}\n")
    digest = hashlib.sha256("".join(sorted(lines)).encode()).hexdigest()
    return f"sha256:{digest}"


# ---------------------------------------------------------------------------
# Torch signing and verification
# ---------------------------------------------------------------------------

def sign_torch(torch_path: Path, signer) -> Path:
    """Sign a torch directory, writing signature.toml.

    Returns path to signature.toml.
    """
    torch_path = Path(torch_path)
    meta_path = torch_path / "metadata.toml"
    meta = toml.load(meta_path)
    namespace = meta["namespace"]
    name = meta["name"]
    version = str(meta["version"])

    content_hash = compute_content_hash(torch_path)
    message_str = f"{namespace}:{name}:{version}:{content_hash}"
    sig_bytes = signer.sign(message_str.encode())
    pub_bytes = signer.public_key_bytes()

    sig_data = {
        "signature": {
            "namespace": namespace,
            "algorithm": signer.algorithm,
            "signed_at": datetime.now(timezone.utc).isoformat(),
            "content_hash": content_hash,
            "message": message_str,
            "value": _b64u(sig_bytes),
        },
        "public_key": {
            "key": _b64u(pub_bytes),
        },
    }
    sig_path = torch_path / "signature.toml"
    with open(sig_path, "w", encoding="utf-8") as f:
        toml.dump(sig_data, f)
    return sig_path


def verify_torch(
    torch_path: Path,
    public_key_b64: Optional[str] = None,
    algorithm: Optional[str] = None,
) -> VerifyResult:
    """Verify the content signature of a torch.

    If public_key_b64 is None, uses the embedded key in signature.toml.
    Returns VerifyResult(valid, namespace, algorithm, message).
    """
    torch_path = Path(torch_path)
    sig_path = torch_path / "signature.toml"

    if not sig_path.exists():
        return VerifyResult(False, "", "", "No signature.toml found")

    sig_data = toml.load(sig_path)
    sig_section = sig_data.get("signature", {})
    pub_section = sig_data.get("public_key", {})

    namespace = sig_section.get("namespace", "")
    algo = algorithm or sig_section.get("algorithm", "ed25519")
    stored_hash = sig_section.get("content_hash", "")
    message_str = sig_section.get("message", "")
    sig_b64 = sig_section.get("value", "")
    embedded_key_b64 = pub_section.get("key", "")

    key_b64 = public_key_b64 or embedded_key_b64
    if not key_b64:
        return VerifyResult(False, namespace, algo, "No public key available for verification")

    # Recompute content hash
    actual_hash = compute_content_hash(torch_path)
    if actual_hash != stored_hash:
        return VerifyResult(
            False, namespace, algo,
            f"Content hash mismatch: stored={stored_hash} actual={actual_hash}"
        )

    # Verify signature
    try:
        _verify_signature(message_str.encode(), _b64u_decode(sig_b64), _b64u_decode(key_b64), algo)
    except Exception as e:
        return VerifyResult(False, namespace, algo, f"Signature invalid: {e}")

    return VerifyResult(True, namespace, algo, f"Valid signature from '{namespace}'")


# ---------------------------------------------------------------------------
# CID signing and verification
# ---------------------------------------------------------------------------

def _cid_message(cid: str, namespace: str, version: str) -> bytes:
    return f"{namespace}:{version}:{cid}".encode()


def sign_cid(cid: str, namespace: str, version: str, signer) -> str:
    """Sign a CID; return base64url-encoded signature."""
    return _b64u(signer.sign(_cid_message(cid, namespace, version)))


def verify_cid(
    cid: str,
    namespace: str,
    version: str,
    sig_b64: str,
    public_key_b64: str,
    algorithm: str = "ed25519",
) -> bool:
    """Verify a CID signature. Returns True if valid."""
    try:
        _verify_signature(
            _cid_message(cid, namespace, version),
            _b64u_decode(sig_b64),
            _b64u_decode(public_key_b64),
            algorithm,
        )
        return True
    except Exception:
        return False


# ---------------------------------------------------------------------------
# Key registry
# ---------------------------------------------------------------------------

def fetch_key_registry(
    url: str,
    cache_path: Path,
    ttl_hours: int = 24,
) -> Dict[str, str]:
    """Fetch a key registry TOML from url, caching at cache_path.

    Returns dict mapping namespace → base64url public key.
    """
    from torchbase.torchfs import _kubo_url, node, port

    cache_path = Path(cache_path)

    # Serve from cache if fresh
    if cache_path.exists():
        mtime = cache_path.stat().st_mtime
        age_hours = (datetime.now().timestamp() - mtime) / 3600
        if age_hours < ttl_hours:
            return toml.load(cache_path).get("keys", {})

    # Resolve IPNS if needed
    if url.startswith("/ipns/") or url.startswith("ipns://"):
        ipns_name = url.replace("/ipns/", "").replace("ipns://", "")
        resolve = requests.post(
            f"{_kubo_url(node, port)}/name/resolve",
            params={"arg": ipns_name}, timeout=30,
        )
        resolve.raise_for_status()
        cid = resolve.json()["Path"].lstrip("/ipfs/")
        cat = requests.post(
            f"{_kubo_url(node, port)}/cat",
            params={"arg": cid}, timeout=30,
        )
        cat.raise_for_status()
        registry = toml.loads(cat.text)
    else:
        response = requests.get(url, timeout=30)
        response.raise_for_status()
        registry = toml.loads(response.text)

    cache_path.parent.mkdir(parents=True, exist_ok=True)
    with open(cache_path, "w", encoding="utf-8") as f:
        toml.dump(registry, f)

    return registry.get("keys", {})


def resolve_public_key(
    namespace: str,
    config,
    torch_path: Optional[Path] = None,
) -> Optional[tuple]:
    """Find the public key and algorithm for a namespace.

    Search order: key registries in config → trusted_keys in config → embedded in signature.toml.
    Returns (public_key_b64, algorithm) or None if not found.
    """
    cache_path = Path.home() / ".torchbase" / "key_cache.toml"
    ttl = getattr(config, "key_cache_ttl_hours", 24)

    for registry_url in getattr(config, "key_registries", []):
        try:
            keys = fetch_key_registry(registry_url, cache_path, ttl)
            if namespace in keys:
                return keys[namespace], "ed25519"
        except Exception:
            pass

    trusted = getattr(config, "trusted_keys", {})
    if namespace in trusted:
        return trusted[namespace], "ed25519"

    if torch_path is not None:
        sig_path = Path(torch_path) / "signature.toml"
        if sig_path.exists():
            sig_data = toml.load(sig_path)
            key = sig_data.get("public_key", {}).get("key")
            algo = sig_data.get("signature", {}).get("algorithm", "ed25519")
            if key:
                return key, algo

    return None


# ---------------------------------------------------------------------------
# Internal crypto helpers
# ---------------------------------------------------------------------------

def _b64u(data: bytes) -> str:
    return base64.urlsafe_b64encode(data).rstrip(b"=").decode()


def _b64u_decode(s: str) -> bytes:
    pad = 4 - len(s) % 4
    return base64.urlsafe_b64decode(s + "=" * (pad % 4))


def _verify_signature(message: bytes, signature: bytes, public_key_raw: bytes, algorithm: str) -> None:
    """Verify signature; raises on failure."""
    if algorithm == "ed25519":
        from cryptography.hazmat.primitives.asymmetric.ed25519 import Ed25519PublicKey
        pub = Ed25519PublicKey.from_public_bytes(public_key_raw)
        pub.verify(signature, message)
    elif algorithm == "p256":
        from cryptography.hazmat.primitives.asymmetric.ec import (
            ECDSA, EllipticCurvePublicKey, SECP256R1
        )
        from cryptography.hazmat.primitives import hashes
        from cryptography.hazmat.primitives.asymmetric.utils import decode_dss_signature
        from cryptography.hazmat.backends import default_backend
        from cryptography.hazmat.primitives.asymmetric.ec import (
            EllipticCurvePublicNumbers, SECP256R1
        )
        # public_key_raw is uncompressed point (65 bytes: 04 || x || y)
        pub = EllipticCurvePublicKey.from_encoded_point(SECP256R1(), public_key_raw)
        pub.verify(signature, message, ECDSA(hashes.SHA256()))
    else:
        raise ValueError(f"Unknown algorithm: {algorithm}")
