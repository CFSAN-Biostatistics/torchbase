"""Registry torch resolution and CID fetching."""

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

    def resolve(self, torch_name: str, version: Optional[str] = None) -> str:
        """Resolve torch reference to CID.

        Resolves "namespace/name" to a CID by querying registries in order:
        1. Check if version is pinned (pins take precedence)
        2. Try default registry
        3. Try additional registries in order
        4. Return CID for resolved version or raise error

        Args:
            torch_name: Torch reference as "namespace/name"
            version: Version constraint:
                   - None: resolve to latest
                   - "X.Y.Z": resolve to explicit version
                   Overridden by pins if present

        Returns:
            CID as string

        Raises:
            ValueError: If torch not found or version not found
        """
        if "/" not in torch_name:
            raise ValueError(f"Torch name must include namespace: {torch_name}")

        # Check if version is pinned (pins take precedence)
        if torch_name in self.config.pins:
            version = self.config.pins[torch_name]

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

            # Determine which version to use
            if version is None:
                # Use latest if available
                if "latest" in torch_versions:
                    return torch_versions["latest"]
                # Otherwise use first available version
                if torch_versions:
                    return list(torch_versions.values())[0]

            else:
                # Use explicit version
                if version not in torch_versions:
                    raise ValueError(
                        f"Version {version} not found for {torch_name}"
                    )
                return torch_versions[version]

        # Not found in any registry
        raise ValueError(f"Torch {torch_name} not found in any registry")

    def fetch_torch(
        self,
        torch_name: str,
        version: Optional[str] = None,
        pin: bool = False
    ) -> Path:
        """Fetch torch and return local path.

        Main entry point for resolving torch to a local path.
        Pins take precedence over version arguments.

        Args:
            torch_name: Torch reference as "namespace/name"
            version: Version constraint (None for latest)
            pin: Unused in mock implementation (for API compatibility)

        Returns:
            Path to torch (mocked local path for now)

        Raises:
            ValueError: If torch not found or resolution fails
        """
        cid = self.resolve(torch_name, version=version)
        return self._cid_to_local_path(cid)

    def _fetch_manifest(self, registry_url: str) -> Dict:
        """Fetch manifest from registry.

        Args:
            registry_url: URL of the registry

        Returns:
            Dictionary with torch references and CIDs

        Raises:
            Exception: If fetch or parsing fails
        """
        if requests is None:
            raise ImportError("requests library required for registry fetching")

        response = requests.get(registry_url)
        response.raise_for_status()

        # Parse TOML manifest
        manifest_data = toml.loads(response.text)
        return manifest_data

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
        # Mock: return a predictable path based on CID
        ipfs_cache_dir = Path("/tmp/ipfs")
        return ipfs_cache_dir / cid

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
                with open(fd, 'w') as f:
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
