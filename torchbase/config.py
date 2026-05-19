"""Registry configuration parsing and management."""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional
import toml


@dataclass
class RegistryConfig:
    """Container for registry and version pin configuration."""

    default_registry: Optional[str] = None
    additional_registries: List[str] = field(default_factory=list)
    pins: Dict[str, str] = field(default_factory=dict)

    @staticmethod
    def load(user_config_dir: Optional[Path] = None,
             project_config_dir: Optional[Path] = None) -> "RegistryConfig":
        """
        Load registry configuration with hierarchical override.

        Configuration hierarchy (highest to lowest priority):
        1. .torchbase.toml in project_config_dir (current working directory)
        2. ~/.torchbase/config.toml in user_config_dir (user home)
        3. Sensible defaults

        Args:
            user_config_dir: Directory containing .torchbase/config.toml
                           (defaults to user home)
            project_config_dir: Directory containing .torchbase.toml
                              (defaults to current working directory)

        Returns:
            RegistryConfig with merged configuration

        Raises:
            ValueError: If TOML files are invalid
        """
        if user_config_dir is None:
            user_config_dir = Path.home()
        if project_config_dir is None:
            project_config_dir = Path.cwd()

        # Start with defaults
        config = RegistryConfig()

        # Load user config
        user_config_path = Path(user_config_dir) / ".torchbase" / "config.toml"
        if user_config_path.exists():
            try:
                user_data = toml.load(user_config_path)
                config = RegistryConfig._merge_config(config, user_data)
            except toml.TomlDecodeError as e:
                raise ValueError(f"Invalid TOML in {user_config_path}: {e}")

        # Load project config (overrides user config)
        project_config_path = Path(project_config_dir) / ".torchbase.toml"
        if project_config_path.exists():
            try:
                project_data = toml.load(project_config_path)
                config = RegistryConfig._merge_config(config, project_data)
            except toml.TomlDecodeError as e:
                raise ValueError(f"Invalid TOML in {project_config_path}: {e}")

        return config

    @staticmethod
    def _merge_config(config: "RegistryConfig",
                      data: Dict) -> "RegistryConfig":
        """
        Merge loaded TOML data into existing config.

        Later data overrides earlier data (project overrides user).

        Args:
            config: Existing RegistryConfig
            data: Dictionary from toml.load()

        Returns:
            Updated RegistryConfig
        """
        # Merge registries section
        if "registries" in data:
            registries = data["registries"]
            if "default" in registries:
                config.default_registry = registries["default"]
            if "additional" in registries:
                config.additional_registries = registries["additional"]

        # Merge pins section
        if "pins" in data:
            config.pins.update(data["pins"])

        return config
