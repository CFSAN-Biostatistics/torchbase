"""Tests for registry torch resolution."""

import pytest
from pathlib import Path
from unittest.mock import patch, MagicMock
from torchbase.registry import RegistryManager
from torchbase.config import RegistryConfig


class TestRegistryManagerBasics:
    """Test basic RegistryManager initialization and setup."""

    def test_init_with_config(self):
        """Test RegistryManager initialization with config."""
        config = RegistryConfig(
            default_registry="https://registry.example.com",
            additional_registries=["https://registry2.example.com"]
        )
        manager = RegistryManager(config)
        assert manager.config == config

    def test_init_with_none_creates_empty_config(self):
        """Test RegistryManager initialization with None config."""
        manager = RegistryManager(None)
        assert manager.config is not None
        assert isinstance(manager.config, RegistryConfig)


class TestRegistryResolution:
    """Test torch reference resolution to CID."""

    def test_resolve_with_default_registry(self):
        """Test resolving torch reference with default registry."""
        config = RegistryConfig(
            default_registry="https://registry.example.com"
        )
        manager = RegistryManager(config)

        # Mock the manifest fetch
        mock_manifest = {
            "namespace/torchname": {
                "latest": "QmHash1",
                "1.0.0": "QmHash2"
            }
        }

        with patch.object(manager, '_fetch_manifest', return_value=mock_manifest):
            result = manager.resolve("namespace/torchname", version=None)
            assert result == "QmHash1"

    def test_resolve_with_explicit_version(self):
        """Test resolving torch reference with explicit version."""
        config = RegistryConfig(
            default_registry="https://registry.example.com"
        )
        manager = RegistryManager(config)

        mock_manifest = {
            "namespace/torchname": {
                "latest": "QmHash1",
                "1.0.0": "QmHash2",
                "2.0.0": "QmHash3"
            }
        }

        with patch.object(manager, '_fetch_manifest', return_value=mock_manifest):
            result = manager.resolve("namespace/torchname", version="1.0.0")
            assert result == "QmHash2"

    def test_resolve_with_pinned_version(self):
        """Test that pins take precedence over version argument."""
        config = RegistryConfig(
            default_registry="https://registry.example.com",
            pins={"namespace/torchname": "1.5.0"}
        )
        manager = RegistryManager(config)

        mock_manifest = {
            "namespace/torchname": {
                "latest": "QmHash1",
                "1.0.0": "QmHash2",
                "1.5.0": "QmHash5",
                "2.0.0": "QmHash3"
            }
        }

        with patch.object(manager, '_fetch_manifest', return_value=mock_manifest):
            # Even if we pass version=2.0.0, pin should win
            result = manager.resolve("namespace/torchname", version="2.0.0")
            assert result == "QmHash5"

    def test_resolve_latest_when_no_version(self):
        """Test resolving to latest when no version specified."""
        config = RegistryConfig(
            default_registry="https://registry.example.com"
        )
        manager = RegistryManager(config)

        mock_manifest = {
            "namespace/torchname": {
                "latest": "QmHashLatest",
                "1.0.0": "QmHash1"
            }
        }

        with patch.object(manager, '_fetch_manifest', return_value=mock_manifest):
            result = manager.resolve("namespace/torchname", version=None)
            assert result == "QmHashLatest"

    def test_resolve_with_multi_registry_fallback(self):
        """Test multi-registry fallback: try default, then additional."""
        config = RegistryConfig(
            default_registry="https://registry1.example.com",
            additional_registries=[
                "https://registry2.example.com",
                "https://registry3.example.com"
            ]
        )
        manager = RegistryManager(config)

        # Mock: first registry doesn't have torch, second does
        def mock_fetch(registry_url):
            if registry_url == "https://registry1.example.com":
                return {}  # Not found in first
            elif registry_url == "https://registry2.example.com":
                return {
                    "namespace/torchname": {
                        "latest": "QmHashFromSecond"
                    }
                }
            return {}

        with patch.object(manager, '_fetch_manifest', side_effect=mock_fetch):
            result = manager.resolve("namespace/torchname", version=None)
            assert result == "QmHashFromSecond"

    def test_resolve_raises_not_found(self):
        """Test that resolution raises when torch not found in any registry."""
        config = RegistryConfig(
            default_registry="https://registry.example.com"
        )
        manager = RegistryManager(config)

        with patch.object(manager, '_fetch_manifest', return_value={}):
            with pytest.raises(ValueError, match="not found"):
                manager.resolve("namespace/unknown", version=None)

    def test_resolve_raises_version_not_found(self):
        """Test that resolution raises when specific version not found."""
        config = RegistryConfig(
            default_registry="https://registry.example.com"
        )
        manager = RegistryManager(config)

        mock_manifest = {
            "namespace/torchname": {
                "latest": "QmHashLatest",
                "1.0.0": "QmHash1"
            }
        }

        with patch.object(manager, '_fetch_manifest', return_value=mock_manifest):
            with pytest.raises(ValueError, match="Version 9.0.0 not found"):
                manager.resolve("namespace/torchname", version="9.0.0")


class TestFetchTorch:
    """Test the main fetch_torch method with local path mocking."""

    def test_fetch_torch_returns_local_path(self):
        """Test that fetch_torch returns a local path for mocked IPFS."""
        config = RegistryConfig(
            default_registry="https://registry.example.com"
        )
        manager = RegistryManager(config)

        mock_manifest = {
            "namespace/torchname": {
                "latest": "QmHash1"
            }
        }

        expected_path = Path("/tmp/ipfs/QmHash1")

        with patch.object(manager, '_fetch_manifest', return_value=mock_manifest):
            with patch.object(manager, '_cid_to_local_path', return_value=expected_path):
                result = manager.fetch_torch("namespace/torchname", version=None)
                assert result == expected_path

    def test_fetch_torch_with_pin(self):
        """Test fetch_torch with pinned version."""
        config = RegistryConfig(
            default_registry="https://registry.example.com",
            pins={"namespace/torchname": "1.5.0"}
        )
        manager = RegistryManager(config)

        mock_manifest = {
            "namespace/torchname": {
                "latest": "QmHash1",
                "1.5.0": "QmHash5"
            }
        }

        expected_path = Path("/tmp/ipfs/QmHash5")

        with patch.object(manager, '_fetch_manifest', return_value=mock_manifest):
            with patch.object(manager, '_cid_to_local_path', return_value=expected_path):
                result = manager.fetch_torch("namespace/torchname", version="2.0.0")
                assert result == expected_path

    def test_fetch_torch_integration(self):
        """Test full fetch_torch with mocked manifest and IPFS."""
        config = RegistryConfig(
            default_registry="https://registry.example.com"
        )
        manager = RegistryManager(config)

        mock_manifest = {
            "test/myapp": {
                "latest": "QmAbcd1234",
                "1.0.0": "QmAbcd5678"
            }
        }

        def mock_fetch(registry_url):
            return mock_manifest

        expected_path = Path("/mnt/cache/ipfs/QmAbcd1234")

        with patch.object(manager, '_fetch_manifest', side_effect=mock_fetch):
            with patch.object(manager, '_cid_to_local_path', return_value=expected_path):
                result = manager.fetch_torch("test/myapp")
                assert result == expected_path


class TestManifestFetching:
    """Test manifest fetching and parsing."""

    def test_fetch_manifest_network_error(self):
        """Test handling of network errors when fetching manifest."""
        config = RegistryConfig(
            default_registry="https://registry.example.com"
        )
        manager = RegistryManager(config)

        with patch('torchbase.registry.requests.get') as mock_get:
            mock_get.side_effect = Exception("Network error")
            with pytest.raises(Exception, match="Network error"):
                manager._fetch_manifest("https://registry.example.com")

    def test_fetch_manifest_malformed_toml(self):
        """Test handling of malformed TOML in manifest."""
        config = RegistryConfig(
            default_registry="https://registry.example.com"
        )
        manager = RegistryManager(config)

        with patch('torchbase.registry.requests.get') as mock_get:
            mock_response = MagicMock()
            mock_response.text = "invalid toml {{{"
            mock_get.return_value = mock_response

            with pytest.raises(Exception):
                manager._fetch_manifest("https://registry.example.com")


class TestCIDToLocalPath:
    """Test CID to local path mapping (mocked IPFS)."""

    def test_cid_to_local_path_mapping(self):
        """Test that CID is mapped to expected local path."""
        config = RegistryConfig()
        manager = RegistryManager(config)

        cid = "QmHash1234"
        result = manager._cid_to_local_path(cid)

        # Should return a path containing the CID
        assert cid in str(result)
        assert isinstance(result, Path)

    def test_cid_to_local_path_consistency(self):
        """Test that same CID maps to same path."""
        config = RegistryConfig()
        manager = RegistryManager(config)

        cid = "QmHash5678"
        path1 = manager._cid_to_local_path(cid)
        path2 = manager._cid_to_local_path(cid)

        assert path1 == path2


class TestEdgeCases:
    """Test edge cases and error handling."""

    def test_torch_name_without_namespace(self):
        """Test that torch name without namespace raises error."""
        config = RegistryConfig(
            default_registry="https://registry.example.com"
        )
        manager = RegistryManager(config)

        with pytest.raises(ValueError, match="namespace"):
            manager.resolve("justname", version=None)

    def test_empty_registry_config(self):
        """Test behavior with empty registry config."""
        config = RegistryConfig()
        manager = RegistryManager(config)

        assert manager.config.default_registry is None
        assert manager.config.additional_registries == []

    def test_no_additional_registries_still_works(self):
        """Test that manager works with only default registry."""
        config = RegistryConfig(
            default_registry="https://registry.example.com"
        )
        manager = RegistryManager(config)

        mock_manifest = {
            "namespace/torchname": {
                "latest": "QmHash1"
            }
        }

        with patch.object(manager, '_fetch_manifest', return_value=mock_manifest):
            result = manager.resolve("namespace/torchname")
            assert result == "QmHash1"
