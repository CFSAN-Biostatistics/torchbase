"""Tests for CLI pin mechanism."""

from pathlib import Path
from click.testing import CliRunner
from unittest.mock import patch, MagicMock
from torchbase.cli import cli


class TestCLIPinMechanism:
    """Test CLI --pin flag for torch pull command."""

    def test_pull_with_pin_flag(self, tmp_path):
        """Test torch pull --pin creates config with pinned version."""
        runner = CliRunner()

        with runner.isolated_filesystem(temp_dir=tmp_path):
            # Mock RegistryManager and config loading at the point of import
            with patch('torchbase.config.RegistryConfig.load') as mock_config_load:
                with patch('torchbase.registry.RegistryManager') as mock_manager_class:
                    mock_config = MagicMock()
                    mock_config_load.return_value = mock_config

                    mock_manager = MagicMock()
                    mock_manager_class.return_value = mock_manager

                    result = runner.invoke(cli, ['pull', '--pin', 'namespace/torch'])

                    # Verify pin_torch was called
                    mock_manager.pin_torch.assert_called_once()
                    call_args = mock_manager.pin_torch.call_args
                    assert call_args[0][0] == 'namespace/torch'
                    assert call_args[1]['version'] is None

                    assert result.exit_code == 0

    def test_pull_with_pin_and_explicit_version(self, tmp_path):
        """Test torch pull --pin --version pins specific version."""
        runner = CliRunner()

        with runner.isolated_filesystem(temp_dir=tmp_path):
            # Mock RegistryManager and config loading
            with patch('torchbase.config.RegistryConfig.load') as mock_config_load:
                with patch('torchbase.registry.RegistryManager') as mock_manager_class:
                    mock_config = MagicMock()
                    mock_config_load.return_value = mock_config

                    mock_manager = MagicMock()
                    mock_manager_class.return_value = mock_manager

                    result = runner.invoke(
                        cli,
                        ['pull', '--pin', '--version', '1.5.0', 'namespace/torch']
                    )

                    # Verify pin_torch was called with explicit version
                    mock_manager.pin_torch.assert_called_once()
                    call_args = mock_manager.pin_torch.call_args
                    assert call_args[0][0] == 'namespace/torch'
                    assert call_args[1]['version'] == '1.5.0'

                    assert result.exit_code == 0

    def test_pull_without_pin_flag(self, tmp_path):
        """Test torch pull without --pin performs normal pull."""
        runner = CliRunner()

        with runner.isolated_filesystem(temp_dir=tmp_path):
            # Mock RegistryManager and config loading
            with patch('torchbase.config.RegistryConfig.load') as mock_config_load:
                with patch('torchbase.registry.RegistryManager') as mock_manager_class:
                    mock_config = MagicMock()
                    mock_config_load.return_value = mock_config

                    mock_manager = MagicMock()
                    mock_manager.fetch_torch.return_value = Path("/tmp/ipfs/QmHash")
                    mock_manager_class.return_value = mock_manager

                    result = runner.invoke(cli, ['pull', 'namespace/torch'])

                    # Verify fetch_torch was called, not pin_torch
                    mock_manager.fetch_torch.assert_called_once()
                    mock_manager.pin_torch.assert_not_called()

                    assert result.exit_code == 0

    def test_pull_pin_error_handling(self, tmp_path):
        """Test that CLI handles pin errors gracefully."""
        runner = CliRunner()

        with runner.isolated_filesystem(temp_dir=tmp_path):
            # Mock RegistryManager to raise error
            with patch('torchbase.config.RegistryConfig.load') as mock_config_load:
                with patch('torchbase.registry.RegistryManager') as mock_manager_class:
                    mock_config = MagicMock()
                    mock_config_load.return_value = mock_config

                    mock_manager = MagicMock()
                    mock_manager.pin_torch.side_effect = ValueError("Torch not found")
                    mock_manager_class.return_value = mock_manager

                    result = runner.invoke(cli, ['pull', '--pin', 'namespace/unknown'])

                    # Verify error is handled
                    assert result.exit_code != 0
                    assert "Torch not found" in result.output

    def test_pull_pin_uses_directory_config(self, tmp_path):
        """Test that --pin writes to directory config when in a project."""
        runner = CliRunner()

        with runner.isolated_filesystem(temp_dir=tmp_path):
            # Create a project directory structure
            project_dir = Path.cwd()
            project_dir.mkdir(exist_ok=True)

            # Mock RegistryManager
            with patch('torchbase.config.RegistryConfig.load') as mock_config_load:
                with patch('torchbase.registry.RegistryManager') as mock_manager_class:
                    mock_config = MagicMock()
                    mock_config_load.return_value = mock_config

                    mock_manager = MagicMock()
                    mock_manager_class.return_value = mock_manager

                    result = runner.invoke(cli, ['pull', '--pin', 'namespace/torch'])

                    # Verify config_path argument points to project config
                    call_args = mock_manager.pin_torch.call_args
                    config_path = call_args[1]['config_path']
                    assert config_path.name == ".torchbase.toml"
                    assert config_path.parent == project_dir

                    assert result.exit_code == 0
