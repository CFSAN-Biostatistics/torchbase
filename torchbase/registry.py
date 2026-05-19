"""Registry torch resolution and CID fetching."""

from pathlib import Path
from typing import Optional, Dict
import toml

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
