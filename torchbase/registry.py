"""Registry torch resolution and CID fetching."""

import warnings
from pathlib import Path
from typing import Optional, Dict
import toml
import tempfile
import shutil

try:
    import requests
except ImportError:
    requests = None

from torchbase.config import RegistryConfig


class RegistryManager:
    """Resolves torch references to IPFS CIDs via IPNS registries.

    Supports:
    - Multiple registry fallback (default, then additional)
    - Version constraints (latest, explicit, pinned)
    - Pin precedence over version arguments
    - Mocked IPFS operations (CID → local path mapping)

    Attributes:
        config: RegistryConfig with registry URLs and pins
    """

    def __init__(self, config: Optional[RegistryConfig] = None):
        """Initialize RegistryManager.

        Args:
            config: RegistryConfig with registries and pins. If None,
                   creates empty config.
        """
        if config is None:
            config = RegistryConfig()
        self.config = config

    def resolve(
        self,
        torch_name: str,
        version: Optional[str] = None,
        require_signature: bool = False,
    ) -> str:
        """Resolve torch reference to CID.

        Resolves "namespace/name" to a CID by querying registries in order:
        1. Check if version is pinned (pins take precedence)
        2. Try default registry
        3. Try additional registries in order
        4. Verify CID signature (warn if missing, raise if invalid)
        5. Return CID for resolved version or raise error

        Args:
            torch_name: Torch reference as "namespace/name"
            version: Version constraint:
                   - None: resolve to latest
                   - "X.Y.Z": resolve to explicit version
                   Overridden by pins if present
            require_signature: If True, missing signatures are treated as errors

        Returns:
            CID as string

        Raises:
            ValueError: If torch not found, version not found, or signature invalid
        """
        if "/" not in torch_name:
            raise ValueError(f"Torch name must include namespace: {torch_name}")

        # Check if version is pinned (pins take precedence)
        if torch_name in self.config.pins:
            version = self.config.pins[torch_name]

        # Check chain-based namespace registries first
        namespace = torch_name.split("/")[0]
        chain_ipns = getattr(self.config, "namespaces", {}).get(namespace)
        if chain_ipns:
            try:
                manifest = self._fetch_chain_manifest(namespace, chain_ipns)
                if torch_name in manifest:
                    torch_versions = manifest[torch_name]
                    if version is None:
                        if "latest" in torch_versions:
                            cid = torch_versions["latest"]
                            resolved_version = next(
                                (v for v, c in torch_versions.items()
                                 if v not in ("latest", "signatures", "workflow")
                                 and not isinstance(c, dict) and c == cid),
                                None,
                            )
                        elif torch_versions:
                            cid = next(
                                c for k, c in torch_versions.items()
                                if k not in ("signatures", "workflow")
                                and not isinstance(c, dict)
                            )
                            resolved_version = next(
                                k for k, c in torch_versions.items()
                                if k not in ("signatures", "workflow")
                                and not isinstance(c, dict)
                            )
                        else:
                            resolved_version = None
                            cid = None
                    else:
                        if version not in torch_versions:
                            raise ValueError(
                                f"Version {version} not found for {torch_name}"
                            )
                        cid = torch_versions[version]
                        resolved_version = version

                    if cid:
                        self._verify_cid_signature(
                            torch_name, resolved_version, cid, torch_versions, require_signature
                        )
                        return cid
            except ValueError:
                raise
            except Exception:
                pass

        # Try registries in order
        registries_to_try = []
        if self.config.default_registry:
            registries_to_try.append(self.config.default_registry)
        registries_to_try.extend(self.config.additional_registries)

        for registry_url in registries_to_try:
            manifest = self._fetch_manifest(registry_url)

            if torch_name not in manifest:
                continue

            torch_versions = manifest[torch_name]

            # Determine which version and CID to use
            if version is None:
                if "latest" in torch_versions:
                    cid = torch_versions["latest"]
                    resolved_version = next(
                        (v for v, c in torch_versions.items()
                         if v not in ("latest", "signatures", "workflow")
                         and not isinstance(c, dict) and c == cid),
                        None,
                    )
                elif torch_versions:
                    cid = next(
                        c for k, c in torch_versions.items()
                        if k not in ("signatures", "workflow")
                        and not isinstance(c, dict)
                    )
                    resolved_version = next(
                        k for k, c in torch_versions.items()
                        if k not in ("signatures", "workflow")
                        and not isinstance(c, dict)
                    )
                else:
                    continue
            else:
                if version not in torch_versions:
                    raise ValueError(
                        f"Version {version} not found for {torch_name}"
                    )
                cid = torch_versions[version]
                resolved_version = version

            self._verify_cid_signature(
                torch_name, resolved_version, cid, torch_versions, require_signature
            )
            return cid

        # Not found in any registry
        raise ValueError(f"Torch {torch_name} not found in any registry")

    def _verify_cid_signature(
        self,
        torch_name: str,
        version: Optional[str],
        cid: str,
        torch_versions: Dict,
        require_signature: bool,
    ) -> None:
        """Check manifest CID signature and warn or raise as appropriate."""
        from torchbase.signing import verify_cid, resolve_public_key

        namespace = torch_name.split("/")[0]
        signatures = torch_versions.get("signatures", {})

        if not isinstance(signatures, dict) or version not in signatures:
            msg = f"No CID signature found for {torch_name} {version}"
            if require_signature:
                raise ValueError(msg)
            warnings.warn(msg, stacklevel=3)
            return

        sig_b64 = signatures[version]
        result = resolve_public_key(namespace, self.config)
        if result is None:
            msg = (
                f"No public key found for namespace '{namespace}'; "
                f"cannot verify {torch_name} {version}"
            )
            if require_signature:
                raise ValueError(msg)
            warnings.warn(msg, stacklevel=3)
            return

        pub_key_b64, algorithm = result
        valid = verify_cid(cid, namespace, version or "", sig_b64, pub_key_b64, algorithm)
        if not valid:
            raise ValueError(
                f"CID signature verification failed for {torch_name} {version}"
            )

    def fetch_torch(
        self,
        torch_name: str,
        version: Optional[str] = None,
        pin: bool = False,
        require_signature: bool = False,
    ) -> Path:
        """Fetch torch and return local path.

        Main entry point for resolving torch to a local path.
        Pins take precedence over version arguments.

        Args:
            torch_name: Torch reference as "namespace/name"
            version: Version constraint (None for latest)
            pin: Unused in mock implementation (for API compatibility)
            require_signature: If True, missing signatures are treated as errors

        Returns:
            Path to torch (mocked local path for now)

        Raises:
            ValueError: If torch not found or resolution fails
        """
        cid = self.resolve(torch_name, version=version, require_signature=require_signature)
        return self._cid_to_local_path(cid)

    def _fetch_chain_manifest(self, namespace: str, ipns_address: str) -> Dict:
        """Fetch and reconstruct a manifest from an IPLD commit chain.

        Resolves *ipns_address* to the chain head CID, walks the chain,
        verifies all block signatures, and returns a manifest dict in the
        same format that _fetch_manifest() produces.

        Args:
            namespace: Namespace string (used for error messages).
            ipns_address: IPNS address the namespace controls.

        Returns:
            Manifest dict (may be empty if no update blocks exist yet).
        """
        from torchbase.chain import get_chain_head, walk_chain, reconstruct_manifest
        from torchbase.torchfs import node, port

        head_cid = get_chain_head(ipns_address, node, port)
        if head_cid is None:
            return {}

        chain = walk_chain(head_cid, node, port)
        return reconstruct_manifest(chain)

    def _fetch_manifest(self, registry_url: str) -> Dict:
        """Fetch manifest from registry.

        Supports two URL forms:
        - /ipns/<name> or ipns://<name>  → resolved via Kubo API
        - http(s)://...                  → fetched directly

        Args:
            registry_url: Registry URL or IPNS path

        Returns:
            Dictionary with torch references and CIDs

        Raises:
            ImportError: If requests not available
            Exception: If fetch or parsing fails
        """
        if requests is None:
            raise ImportError("requests library required for registry fetching")

        from torchbase.torchfs import _kubo_url, node, port

        # Detect IPNS reference and resolve via Kubo
        ipns_name = None
        if registry_url.startswith("/ipns/"):
            ipns_name = registry_url[len("/ipns/"):]
        elif registry_url.startswith("ipns://"):
            ipns_name = registry_url[len("ipns://"):]

        if ipns_name:
            # Step 1: resolve IPNS name → CID
            resolve_resp = requests.post(
                f"{_kubo_url(node, port)}/name/resolve",
                params={"arg": ipns_name},
                timeout=30,
            )
            resolve_resp.raise_for_status()
            cid = resolve_resp.json()["Path"].lstrip("/ipfs/")

            # Step 2: fetch manifest TOML from CID
            cat_resp = requests.post(
                f"{_kubo_url(node, port)}/cat",
                params={"arg": cid},
                timeout=30,
            )
            cat_resp.raise_for_status()
            return toml.loads(cat_resp.text)

        # Plain HTTP registry
        response = requests.get(registry_url, timeout=30)
        response.raise_for_status()
        return toml.loads(response.text)

    def _cid_to_local_path(self, cid: str) -> Path:
        """Map CID to local path.

        For now, this is a mocked IPFS implementation that returns
        a predictable local path. In production, this would fetch
        from IPFS.

        Args:
            cid: IPFS CID

        Returns:
            Path object for the torch

        Note:
            Currently returns a mock path. Real IPFS integration would
            be added here.
        """
        cache_path = Path.home() / ".torchbase" / "cache" / cid
        if not cache_path.exists():
            from torchbase.torchfs import download_torch
            try:
                download_torch(cid)
            except Exception:
                pass  # Return path even if download fails; caller checks exists()
        return cache_path

    def pin_torch(
        self,
        torch_name: str,
        version: Optional[str] = None,
        config_path: Optional[Path] = None
    ) -> None:
        """Pin a torch version to config file.

        Pins take effect on subsequent resolve() calls. First call fetches
        latest (or specified version), subsequent calls are no-op.

        If torch has workflow dependency, recursively pins that too.

        Args:
            torch_name: Torch reference as "namespace/name"
            version: Version to pin (None for latest). If already pinned,
                    this is ignored (idempotent).
            config_path: Path to config file to update. Defaults to
                        ~/.torchbase/config.toml

        Raises:
            ValueError: If torch not found or version not found
            Exception: If config update fails (atomic - no partial writes)
        """
        if config_path is None:
            config_path = Path.home() / ".torchbase" / "config.toml"

        config_path = Path(config_path)

        # Load existing config
        if config_path.exists():
            existing_config = toml.load(config_path)
        else:
            existing_config = {}

        # Ensure pins section exists
        if "pins" not in existing_config:
            existing_config["pins"] = {}

        # Check if already pinned (idempotent)
        if torch_name in existing_config["pins"]:
            # Already pinned - no-op
            return

        try:
            # Resolve to get version if not specified
            if version is None:
                # Fetch manifest to find latest version
                registries_to_try = []
                if self.config.default_registry:
                    registries_to_try.append(self.config.default_registry)
                registries_to_try.extend(self.config.additional_registries)

                manifest = None
                for registry_url in registries_to_try:
                    manifest = self._fetch_manifest(registry_url)
                    if torch_name in manifest:
                        break

                if manifest is None or torch_name not in manifest:
                    raise ValueError(f"Torch {torch_name} not found in any registry")

                torch_versions = manifest[torch_name]

                # Find the actual version for "latest"
                if "latest" in torch_versions:
                    latest_cid = torch_versions["latest"]
                    # Find which version corresponds to this CID
                    for ver, cid in torch_versions.items():
                        if ver != "latest" and ver != "workflow" and cid == latest_cid:
                            version = ver
                            break

                if version is None:
                    raise ValueError(f"Could not determine version for {torch_name}")

                # Check for workflow dependency
                workflow_torch = torch_versions.get("workflow")
            else:
                # Verify explicit version exists
                registries_to_try = []
                if self.config.default_registry:
                    registries_to_try.append(self.config.default_registry)
                registries_to_try.extend(self.config.additional_registries)

                manifest = None
                for registry_url in registries_to_try:
                    manifest = self._fetch_manifest(registry_url)
                    if torch_name in manifest:
                        break

                if manifest is None or torch_name not in manifest:
                    raise ValueError(f"Torch {torch_name} not found in any registry")

                torch_versions = manifest[torch_name]

                if version not in torch_versions:
                    raise ValueError(f"Version {version} not found for {torch_name}")

                # Check for workflow dependency
                workflow_torch = torch_versions.get("workflow")

            # Update config with pin
            existing_config["pins"][torch_name] = version

            # Recursively pin workflow dependency if present
            if workflow_torch and workflow_torch not in existing_config["pins"]:
                # Recursively pin the workflow torch (latest version)
                self._pin_workflow_torch(workflow_torch, existing_config)

            # Atomic write: write to temp file then rename
            config_path.parent.mkdir(parents=True, exist_ok=True)
            fd, temp_path = tempfile.mkstemp(
                dir=config_path.parent,
                suffix=".toml.tmp"
            )
            try:
                with open(fd, 'w', encoding="utf-8") as f:
                    toml.dump(existing_config, f)
                # Atomic rename
                shutil.move(temp_path, config_path)
            except Exception:
                # Clean up temp file on error
                try:
                    Path(temp_path).unlink(missing_ok=True)
                except Exception:
                    pass
                raise

        except Exception:
            # Re-raise to ensure no partial writes
            raise

    def _pin_workflow_torch(
        self,
        workflow_torch: str,
        config_data: Dict
    ) -> None:
        """Pin workflow torch dependency.

        Helper for pin_torch to recursively pin workflow dependencies.

        Args:
            workflow_torch: Workflow torch reference as "namespace/name"
            config_data: Config dictionary to update in place

        Raises:
            ValueError: If workflow torch not found
        """
        # Fetch manifest for workflow torch
        registries_to_try = []
        if self.config.default_registry:
            registries_to_try.append(self.config.default_registry)
        registries_to_try.extend(self.config.additional_registries)

        manifest = None
        for registry_url in registries_to_try:
            manifest = self._fetch_manifest(registry_url)
            if workflow_torch in manifest:
                break

        if manifest is None or workflow_torch not in manifest:
            raise ValueError(
                f"Workflow torch {workflow_torch} not found in any registry"
            )

        torch_versions = manifest[workflow_torch]

        # Find the actual version for "latest"
        version = None
        if "latest" in torch_versions:
            latest_cid = torch_versions["latest"]
            # Find which version corresponds to this CID
            for ver, cid in torch_versions.items():
                if ver != "latest" and ver != "workflow" and cid == latest_cid:
                    version = ver
                    break

        if version is None:
            raise ValueError(
                f"Could not determine version for {workflow_torch}"
            )

        # Add to pins
        config_data["pins"][workflow_torch] = version
