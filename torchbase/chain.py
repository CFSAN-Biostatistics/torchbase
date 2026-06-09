"""IPLD commit chain for namespace registry.

Each block is a TOML document stored as an IPFS object.  Blocks form a
singly-linked list (each references the previous by CID); IPNS tracks the
head.  The chain is tamper-evident: altering any block changes its CID,
which breaks every forward link.

Block formats
-------------
Genesis block (namespace claim)::

    type       = "genesis"
    namespace  = "us-fda-hfp"
    public_key = "<base64url-ed25519-pubkey>"
    timestamp  = "2026-06-07T12:00:00+00:00"
    signature  = "<base64url>"
    # signs: "genesis:{namespace}:{public_key}:{timestamp}"

Update block (torch publication)::

    type      = "update"
    namespace = "us-fda-hfp"
    previous  = "<CID of previous block>"
    timestamp = "2026-06-07T12:00:00+00:00"
    signature = "<base64url>"
    # signs: "update:{namespace}:{previous}:{timestamp}:{sha256(entries_toml)}"

    ["us-fda-hfp/seqsero2"]
    "2.0.0" = "<torch CID>"
    latest   = "<torch CID>"

    ["us-fda-hfp/seqsero2".signatures]
    "2.0.0" = "<base64url CID sig>"

The entries section is the same format as the flat manifest TOML that
RegistryManager already parses.  Replaying blocks oldest-first reconstructs a
full manifest with no changes to the resolution layer.
"""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Optional, TYPE_CHECKING

import toml

if TYPE_CHECKING:
    pass


# ---------------------------------------------------------------------------
# Canonical Kubo base URL
# ---------------------------------------------------------------------------

def _kubo_url(node: str, port: int) -> str:
    return f"http://{node}:{port}/api/v0"


# ---------------------------------------------------------------------------
# Block construction
# ---------------------------------------------------------------------------

def block_signature_message(block: Dict) -> bytes:
    """Return the canonical bytes that the block signature covers."""
    if block["type"] == "genesis":
        msg = (
            f"genesis:{block['namespace']}:{block['public_key']}:{block['timestamp']}"
        )
    elif block["type"] == "update":
        entries_toml = _entries_toml_string(block)
        entries_hash = hashlib.sha256(entries_toml.encode()).hexdigest()
        msg = (
            f"update:{block['namespace']}:{block['previous']}"
            f":{block['timestamp']}:{entries_hash}"
        )
    else:
        raise ValueError(f"Unknown block type: {block['type']!r}")
    return msg.encode()


def make_genesis_block(namespace: str, signer) -> Dict:
    """Build and sign a genesis block for *namespace*.

    Args:
        namespace: Namespace string (e.g. "us-fda-hfp").
        signer: A signer object with .sign(bytes)->bytes and
            .public_key_bytes()->bytes attributes.

    Returns:
        Signed genesis block dict ready for serialisation.
    """
    from torchbase.signing import _b64u

    pub_b64 = _b64u(signer.public_key_bytes())
    timestamp = datetime.now(timezone.utc).isoformat()

    block = {
        "type": "genesis",
        "namespace": namespace,
        "public_key": pub_b64,
        "timestamp": timestamp,
    }
    sig_bytes = signer.sign(block_signature_message(block))
    block["signature"] = _b64u(sig_bytes)
    return block


def make_update_block(
    namespace: str,
    previous_cid: str,
    entries: Dict,
    signer,
) -> Dict:
    """Build and sign an update block for *namespace*.

    Args:
        namespace: Namespace string.
        previous_cid: CID of the immediately preceding block.
        entries: Manifest-format dict mapping "namespace/name" to version
            dicts (same structure RegistryManager consumes).
        signer: Signer object.

    Returns:
        Signed update block dict.
    """
    from torchbase.signing import _b64u

    timestamp = datetime.now(timezone.utc).isoformat()

    block: Dict = {
        "type": "update",
        "namespace": namespace,
        "previous": previous_cid,
        "timestamp": timestamp,
        **entries,
    }
    sig_bytes = signer.sign(block_signature_message(block))
    block["signature"] = _b64u(sig_bytes)
    return block


def _entries_toml_string(block: Dict) -> str:
    """Extract the entries portion of an update block as a TOML string.

    Header fields are excluded; only the manifest-format entries remain.
    """
    header_keys = {"type", "namespace", "previous", "timestamp", "signature"}
    entries = {k: v for k, v in block.items() if k not in header_keys}
    return toml.dumps(entries)


# ---------------------------------------------------------------------------
# Block serialisation
# ---------------------------------------------------------------------------

def block_to_toml(block: Dict) -> str:
    """Serialise a block dict to a TOML string.

    Header fields come first (type, namespace, previous/public_key, timestamp,
    signature), then the entries section so the file is human-readable.
    """
    header_keys = ["type", "namespace", "public_key", "previous", "timestamp", "signature"]
    lines = []
    for key in header_keys:
        if key in block:
            lines.append(f'{key} = {_toml_scalar(block[key])}')

    # Entries (update block only): everything not in header_keys
    entry_keys = {k: v for k, v in block.items() if k not in header_keys}
    if entry_keys:
        lines.append("")
        lines.append(toml.dumps(entry_keys).rstrip())

    return "\n".join(lines) + "\n"


def _toml_scalar(value) -> str:
    """Render a Python scalar as a TOML inline value."""
    return json.dumps(value)  # JSON scalars are valid TOML scalars


# ---------------------------------------------------------------------------
# IPFS / Kubo primitives
# ---------------------------------------------------------------------------

def upload_block(
    block_toml: str,
    node: str = "127.0.0.1",
    port: int = 5001,
) -> str:
    """Upload a TOML string to Kubo and return its CID.

    Args:
        block_toml: TOML-encoded block text.
        node: Kubo API host.
        port: Kubo API port.

    Returns:
        CID string of the uploaded object.
    """
    import requests as _req

    resp = _req.post(
        f"{_kubo_url(node, port)}/add",
        params={"quieter": "true"},
        files={"file": ("block.toml", block_toml.encode(), "application/octet-stream")},
        timeout=60,
    )
    resp.raise_for_status()
    return resp.json()["Hash"]


def fetch_block(
    cid: str,
    node: str = "127.0.0.1",
    port: int = 5001,
) -> Dict:
    """Fetch and parse a block from Kubo.

    Args:
        cid: CID of the block.
        node: Kubo API host.
        port: Kubo API port.

    Returns:
        Parsed block dict.
    """
    import requests as _req

    resp = _req.post(
        f"{_kubo_url(node, port)}/cat",
        params={"arg": cid},
        timeout=30,
    )
    resp.raise_for_status()
    return toml.loads(resp.text)


def pin_cid(
    cid: str,
    node: str = "127.0.0.1",
    port: int = 5001,
) -> None:
    """Pin a CID in Kubo so it won't be garbage-collected."""
    import requests as _req

    resp = _req.post(
        f"{_kubo_url(node, port)}/pin/add",
        params={"arg": cid},
        timeout=60,
    )
    resp.raise_for_status()


def upload_directory(
    torch_path: Path,
    node: str = "127.0.0.1",
    port: int = 5001,
) -> str:
    """Upload a torch directory to Kubo and return the root CID.

    Kubo /add with recursive=true requires one multipart part per file,
    with the relative path from the *parent* of the torch directory in
    the Content-Disposition filename field.

    Args:
        torch_path: Path to the torch directory (e.g. ".../1.0.0.torch").
        node: Kubo API host.
        port: Kubo API port.

    Returns:
        Root CID of the uploaded directory.
    """
    import requests as _req

    torch_path = Path(torch_path)
    file_parts = []
    for f in sorted(torch_path.rglob("*")):
        if not f.is_file():
            continue
        rel = str(f.relative_to(torch_path.parent))
        file_parts.append(
            ("file", (rel, open(f, "rb"), "application/octet-stream"))
        )

    resp = _req.post(
        f"{_kubo_url(node, port)}/add",
        params={"recursive": "true", "quieter": "true"},
        files=file_parts,
        timeout=300,
    )
    resp.raise_for_status()
    return resp.json()["Hash"]


def get_chain_head(
    ipns_address: str,
    node: str = "127.0.0.1",
    port: int = 5001,
) -> Optional[str]:
    """Resolve an IPNS address to its current CID (chain head).

    Returns None if the name has not been published yet.
    """
    import requests as _req

    try:
        resp = _req.post(
            f"{_kubo_url(node, port)}/name/resolve",
            params={"arg": ipns_address},
            timeout=30,
        )
        resp.raise_for_status()
        path = resp.json()["Path"]
        return path.lstrip("/ipfs/").split("/")[0]
    except Exception:
        return None


def publish_chain_head(
    cid: str,
    key_name: str,
    node: str = "127.0.0.1",
    port: int = 5001,
) -> None:
    """Publish *cid* as the IPNS head for *key_name* in Kubo's keystore."""
    import requests as _req

    resp = _req.post(
        f"{_kubo_url(node, port)}/name/publish",
        params={"arg": cid, "key": key_name},
        timeout=120,
    )
    resp.raise_for_status()


def import_key_to_kubo(
    namespace: str,
    pem_path: Path,
    node: str = "127.0.0.1",
    port: int = 5001,
) -> str:
    """Import an Ed25519 PEM private key into Kubo's keystore.

    Returns the IPNS address (e.g. "/ipns/k51q...") that the key controls.
    """
    import requests as _req

    pem_path = Path(pem_path)
    resp = _req.post(
        f"{_kubo_url(node, port)}/key/import",
        params={"arg": namespace},
        files={"file": open(pem_path, "rb")},
        timeout=30,
    )
    resp.raise_for_status()
    return "/ipns/" + resp.json()["Id"]


# ---------------------------------------------------------------------------
# Chain walking and verification
# ---------------------------------------------------------------------------

def verify_block(block: Dict, public_key_b64: str) -> bool:
    """Verify a block's signature against *public_key_b64*.

    Returns True if valid, False otherwise.
    """
    from torchbase.signing import _b64u_decode, _verify_signature

    sig_b64 = block.get("signature", "")
    try:
        _verify_signature(
            block_signature_message(block),
            _b64u_decode(sig_b64),
            _b64u_decode(public_key_b64),
            "ed25519",
        )
        return True
    except Exception:
        return False


def walk_chain(
    head_cid: str,
    node: str = "127.0.0.1",
    port: int = 5001,
) -> List[Dict]:
    """Fetch and verify the entire chain from genesis to head.

    Walks backwards from *head_cid* to the genesis block by following
    'previous' links, then reverses to return blocks oldest-first.

    Each non-genesis block's signature is verified against the public key
    declared in the genesis block.

    Args:
        head_cid: CID of the most recent block.
        node: Kubo API host.
        port: Kubo API port.

    Returns:
        List of block dicts, oldest (genesis) first.

    Raises:
        ValueError: If a signature is invalid or a block is malformed.
    """
    blocks_with_cids: List[tuple] = []
    current_cid = head_cid

    # Walk backwards collecting blocks
    while True:
        try:
            block = fetch_block(current_cid, node, port)
        except Exception as e:
            raise ValueError(f"Could not fetch block {current_cid}: {e}") from e
        blocks_with_cids.append((current_cid, block))
        if block.get("type") == "genesis":
            break
        previous = block.get("previous")
        if not previous:
            raise ValueError(
                f"Block {current_cid} has no 'previous' link but is not a genesis block"
            )
        current_cid = previous

    # Reverse to get genesis-first order
    blocks_with_cids.reverse()
    blocks = [b for _, b in blocks_with_cids]

    # Extract public key from genesis block
    genesis = blocks[0]
    if genesis.get("type") != "genesis":
        raise ValueError("First block is not a genesis block")
    public_key_b64 = genesis.get("public_key", "")
    if not public_key_b64:
        raise ValueError("Genesis block missing public_key")

    # Verify every block's signature
    for cid, block in blocks_with_cids:
        if not verify_block(block, public_key_b64):
            raise ValueError(f"Block {cid} has an invalid signature")

    return blocks


# ---------------------------------------------------------------------------
# Manifest reconstruction
# ---------------------------------------------------------------------------

def reconstruct_manifest(chain: List[Dict]) -> Dict:
    """Replay update blocks to reconstruct a manifest dict.

    Blocks are replayed oldest-first; later entries for the same key win.

    Args:
        chain: List of block dicts as returned by walk_chain() (genesis first).

    Returns:
        Manifest dict in RegistryManager format.
    """
    header_keys = {"type", "namespace", "public_key", "previous", "timestamp", "signature"}
    manifest: Dict = {}
    for block in chain:
        if block.get("type") != "update":
            continue
        for key, value in block.items():
            if key in header_keys:
                continue
            if isinstance(value, dict):
                if key not in manifest:
                    manifest[key] = {}
                manifest[key].update(value)
    return manifest


# ---------------------------------------------------------------------------
# Log submission
# ---------------------------------------------------------------------------

def submit_to_log(
    genesis_cid: str,
    namespace: str,
    ipns_address: str,
    log_url: str,
) -> None:
    """Submit a genesis CID to a CT-style log operator.

    The log operator records the submission; it does not approve it.
    This is a fire-and-forget operation — errors are raised but the
    caller should treat them as non-fatal for the publish workflow.

    Args:
        genesis_cid: CID of the genesis block.
        namespace: Namespace string.
        ipns_address: IPNS address the namespace controls.
        log_url: HTTP(S) endpoint of the log operator.
    """
    import requests as _req

    payload = {
        "genesis_cid": genesis_cid,
        "namespace": namespace,
        "ipns_address": ipns_address,
    }
    resp = _req.post(log_url, json=payload, timeout=30)
    resp.raise_for_status()
