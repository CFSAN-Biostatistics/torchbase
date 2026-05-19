#!/usr/bin/env python

"""Tests for registry configuration system."""

import pytest
import tempfile
from pathlib import Path
from torchbase.config import RegistryConfig


class TestRegistryConfigTracer:
    """Tracer bullet: Can load user config from home directory"""

    def test_load_user_config_basic(self):
        """Load basic user config from ~/.torchbase/config.toml"""
        with tempfile.TemporaryDirectory() as tmpdir:
            # Create a user config in home/.torchbase/config.toml
            home_dir = Path(tmpdir)
            config_path = home_dir / ".torchbase" / "config.toml"
            config_path.parent.mkdir(parents=True, exist_ok=True)
            config_path.write_text("""
[registries]
default = "/ipns/registry.torchbase.org"
additional = ["/ipns/myorg.torches.org"]

[pins]
"pubmlst/mlst-database" = "2024-03-15"
""")

            # Load with custom home directory
            config = RegistryConfig.load(user_config_dir=home_dir)

            # Verify basic properties
            assert config.default_registry == "/ipns/registry.torchbase.org"
            assert config.additional_registries == ["/ipns/myorg.torches.org"]
            assert config.pins["pubmlst/mlst-database"] == "2024-03-15"


class TestProjectConfig:
    """Load project-specific config from .torchbase.toml"""

    def test_load_project_config(self):
        """Load project config from current directory"""
        with tempfile.TemporaryDirectory() as tmpdir:
            project_dir = Path(tmpdir)
            config_path = project_dir / ".torchbase.toml"
            config_path.write_text("""
[registries]
default = "/ipns/project.local"

[pins]
"project/custom-torch" = "1.0.0"
""")

            config = RegistryConfig.load(project_config_dir=project_dir)

            assert config.default_registry == "/ipns/project.local"
            assert config.pins["project/custom-torch"] == "1.0.0"


class TestHierarchicalOverride:
    """Directory config overrides user config"""

    def test_project_overrides_user_registries(self):
        """Project config overrides user config for registries"""
        with tempfile.TemporaryDirectory() as tmpdir1:
            with tempfile.TemporaryDirectory() as tmpdir2:
                # Create user config
                home_dir = Path(tmpdir1)
                user_config = home_dir / ".torchbase" / "config.toml"
                user_config.parent.mkdir(parents=True, exist_ok=True)
                user_config.write_text("""
[registries]
default = "/ipns/user.registry"
additional = ["/ipns/user.extra"]

[pins]
"shared/torch" = "1.0.0"
""")

                # Create project config that overrides
                project_dir = Path(tmpdir2)
                project_config = project_dir / ".torchbase.toml"
                project_config.write_text("""
[registries]
default = "/ipns/project.registry"
""")

                config = RegistryConfig.load(
                    user_config_dir=home_dir,
                    project_config_dir=project_dir
                )

                # Project should override default, but not affect other settings
                assert config.default_registry == "/ipns/project.registry"
                assert config.additional_registries == ["/ipns/user.extra"]
                # User pins should still be there
                assert config.pins["shared/torch"] == "1.0.0"

    def test_project_overrides_user_pins(self):
        """Project config overrides user pins"""
        with tempfile.TemporaryDirectory() as tmpdir1:
            with tempfile.TemporaryDirectory() as tmpdir2:
                # Create user config
                home_dir = Path(tmpdir1)
                user_config = home_dir / ".torchbase" / "config.toml"
                user_config.parent.mkdir(parents=True, exist_ok=True)
                user_config.write_text("""
[pins]
"shared/torch" = "1.0.0"
"other/torch" = "2.0.0"
""")

                # Create project config with override
                project_dir = Path(tmpdir2)
                project_config = project_dir / ".torchbase.toml"
                project_config.write_text("""
[pins]
"shared/torch" = "1.5.0"
""")

                config = RegistryConfig.load(
                    user_config_dir=home_dir,
                    project_config_dir=project_dir
                )

                # Project version should override user version
                assert config.pins["shared/torch"] == "1.5.0"
                # Other user pins should remain
                assert config.pins["other/torch"] == "2.0.0"


class TestDefaultBehavior:
    """Graceful defaults when no config files exist"""

    def test_no_config_files(self):
        """Return sensible defaults when no config files exist"""
        with tempfile.TemporaryDirectory() as tmpdir1:
            with tempfile.TemporaryDirectory() as tmpdir2:
                config = RegistryConfig.load(
                    user_config_dir=Path(tmpdir1),
                    project_config_dir=Path(tmpdir2)
                )

                assert config.default_registry is None
                assert config.additional_registries == []
                assert config.pins == {}

    def test_only_additional_registries(self):
        """Config with only additional registries (no default)"""
        with tempfile.TemporaryDirectory() as tmpdir:
            home_dir = Path(tmpdir)
            config_path = home_dir / ".torchbase" / "config.toml"
            config_path.parent.mkdir(parents=True, exist_ok=True)
            config_path.write_text("""
[registries]
additional = ["/ipns/registry1", "/ipns/registry2"]
""")

            config = RegistryConfig.load(user_config_dir=home_dir)

            assert config.default_registry is None
            assert config.additional_registries == ["/ipns/registry1", "/ipns/registry2"]

    def test_empty_registries_section(self):
        """Config with empty registries section"""
        with tempfile.TemporaryDirectory() as tmpdir:
            home_dir = Path(tmpdir)
            config_path = home_dir / ".torchbase" / "config.toml"
            config_path.parent.mkdir(parents=True, exist_ok=True)
            config_path.write_text("""
[registries]
""")

            config = RegistryConfig.load(user_config_dir=home_dir)

            assert config.default_registry is None
            assert config.additional_registries == []


class TestErrorHandling:
    """Invalid TOML raises clear error messages"""

    def test_invalid_user_toml(self):
        """Invalid TOML in user config raises ValueError"""
        with tempfile.TemporaryDirectory() as tmpdir:
            home_dir = Path(tmpdir)
            config_path = home_dir / ".torchbase" / "config.toml"
            config_path.parent.mkdir(parents=True, exist_ok=True)
            config_path.write_text("invalid toml [[[")

            with pytest.raises(ValueError, match="Invalid TOML"):
                RegistryConfig.load(user_config_dir=home_dir)

    def test_invalid_project_toml(self):
        """Invalid TOML in project config raises ValueError"""
        with tempfile.TemporaryDirectory() as tmpdir:
            project_dir = Path(tmpdir)
            config_path = project_dir / ".torchbase.toml"
            config_path.write_text("invalid toml ]]]")

            with pytest.raises(ValueError, match="Invalid TOML"):
                RegistryConfig.load(project_config_dir=project_dir)
