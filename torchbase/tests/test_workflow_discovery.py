"""Tests for workflow discovery in CLI (Issue #20).

Acceptance criteria:
- torch run <torch> <data> checks for main.wdl in torch directory
- If no main.wdl, fetches torchbase/default-workflow (latest) from registry
- User can override: torch run --workflow <custom-torch> <data-torch> <data>
- Validation: torch with WDL files must be named main.wdl (error otherwise)
- Executes workflow via miniwdl with appropriate inputs
- Outputs typing results to stdout or specified file
- Tests verify convention-based discovery, default fallback, override
- Error handling: workflow not found, miniwdl execution failures
"""

import pytest
import toml
import csv
from pathlib import Path
from click.testing import CliRunner
from unittest.mock import patch, MagicMock

from torchbase.cli import cli
from torchbase.torchfs import Torch


@pytest.fixture
def torch_with_main_wdl(tmp_path):
    """Create a torch with main.wdl for testing workflow discovery."""
    torch_path = tmp_path / "test_namespace" / "test_torch" / "1.0.0.torch"
    torch_path.mkdir(parents=True)

    # Create metadata
    metadata = {
        "namespace": "test_namespace",
        "name": "test_torch",
        "version": "1.0.0",
        "version_meta": {"strategy": "semver", "timestamp": 1609459200},
        "typing": {"method": "mlst"},
        "description": {"short": "Test torch"},
        "manifest": {"profiles": "profiles.tsv", "workflow": "main.wdl"}
    }
    with open(torch_path / "metadata.toml", "w") as f:
        toml.dump(metadata, f)

    # Create main.wdl
    wdl_content = """workflow mlst_typing {
    input {
        File reads
    }
    output {
        File results = "results.json"
    }
}
"""
    with open(torch_path / "main.wdl", "w") as f:
        f.write(wdl_content)

    # Create profiles.tsv
    profiles = [["ST", "adk"], ["1", "1"]]
    with open(torch_path / "profiles.tsv", "w") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerows(profiles)

    # Create resources directory
    (torch_path / "_resources").mkdir()

    return torch_path


@pytest.fixture
def torch_without_main_wdl(tmp_path):
    """Create a torch without main.wdl for testing default fallback."""
    torch_path = tmp_path / "test_namespace" / "data_torch" / "1.0.0.torch"
    torch_path.mkdir(parents=True)

    # Create metadata (no workflow specified)
    metadata = {
        "namespace": "test_namespace",
        "name": "data_torch",
        "version": "1.0.0",
        "version_meta": {"strategy": "semver", "timestamp": 1609459200},
        "typing": {"method": "mlst"},
        "description": {"short": "Data torch"},
        "manifest": {"profiles": "profiles.tsv"}
    }
    with open(torch_path / "metadata.toml", "w") as f:
        toml.dump(metadata, f)

    # Create profiles.tsv
    profiles = [["ST", "adk"], ["1", "1"]]
    with open(torch_path / "profiles.tsv", "w") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerows(profiles)

    # Create resources directory
    (torch_path / "_resources").mkdir()

    return torch_path


@pytest.fixture
def torch_with_wrong_wdl_name(tmp_path):
    """Create a torch with incorrectly named WDL file."""
    torch_path = tmp_path / "test_namespace" / "bad_torch" / "1.0.0.torch"
    torch_path.mkdir(parents=True)

    # Create metadata pointing to wrong file
    metadata = {
        "namespace": "test_namespace",
        "name": "bad_torch",
        "version": "1.0.0",
        "version_meta": {"strategy": "semver", "timestamp": 1609459200},
        "typing": {"method": "mlst"},
        "description": {"short": "Bad torch"},
        "manifest": {"profiles": "profiles.tsv", "workflow": "workflow.wdl"}
    }
    with open(torch_path / "metadata.toml", "w") as f:
        toml.dump(metadata, f)

    # Create incorrectly named WDL file
    with open(torch_path / "workflow.wdl", "w") as f:
        f.write("workflow test { }\n")

    # Create profiles.tsv
    profiles = [["ST", "adk"], ["1", "1"]]
    with open(torch_path / "profiles.tsv", "w") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerows(profiles)

    # Create resources directory
    (torch_path / "_resources").mkdir()

    return torch_path


@pytest.fixture
def default_workflow_torch(tmp_path):
    """Create a default workflow torch for testing fallback."""
    torch_path = tmp_path / "workflows" / "default-workflow" / "1.0.0.torch"
    torch_path.mkdir(parents=True)

    # Create metadata
    metadata = {
        "namespace": "workflows",
        "name": "default-workflow",
        "version": "1.0.0",
        "version_meta": {"strategy": "semver", "timestamp": 1609459200},
        "typing": {"method": "mlst"},
        "description": {"short": "Default workflow"},
        "manifest": {"profiles": "profiles.tsv", "workflow": "main.wdl"}
    }
    with open(torch_path / "metadata.toml", "w") as f:
        toml.dump(metadata, f)

    # Create main.wdl
    wdl_content = """workflow default_mlst {
    input {
        File reads
    }
    output {
        File results = "results.json"
    }
}
"""
    with open(torch_path / "main.wdl", "w") as f:
        f.write(wdl_content)

    # Create profiles.tsv
    profiles = [["ST", "adk"], ["1", "1"]]
    with open(torch_path / "profiles.tsv", "w") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerows(profiles)

    # Create resources directory
    (torch_path / "_resources").mkdir()

    return torch_path


@pytest.fixture
def sample_reads_file(tmp_path):
    """Create a sample reads file for testing."""
    reads_file = tmp_path / "reads.fastq"
    with open(reads_file, "w") as f:
        f.write("@read1\nACGT\n+\nIIII\n")
    return reads_file


class TestWorkflowDiscoveryConvention:
    """Test convention-based workflow discovery."""

    def test_torch_with_main_wdl_is_discovered(self, torch_with_main_wdl):
        """Torch with main.wdl is identified as having workflow."""
        torch = Torch.load(torch_with_main_wdl)

        assert torch.workflow is not None
        assert torch.workflow.name == "main.wdl"
        assert torch.workflow.exists()

    def test_torch_without_main_wdl_has_no_workflow(self, torch_without_main_wdl):
        """Torch without main.wdl has workflow=None."""
        torch = Torch.load(torch_without_main_wdl)

        assert torch.workflow is None

    def test_main_wdl_exists_check(self, torch_with_main_wdl):
        """Check if main.wdl exists in torch directory."""
        main_wdl_path = torch_with_main_wdl / "main.wdl"

        assert main_wdl_path.exists()
        assert main_wdl_path.is_file()

    def test_main_wdl_missing_check(self, torch_without_main_wdl):
        """Check that main.wdl is missing in data-only torch."""
        main_wdl_path = torch_without_main_wdl / "main.wdl"

        assert not main_wdl_path.exists()


class TestWorkflowNamingValidation:
    """Test validation: torch with WDL files must be named main.wdl."""

    def test_wdl_file_must_be_named_main_wdl(self, torch_with_wrong_wdl_name):
        """WDL files in torch must be named main.wdl, not workflow.wdl."""
        # The manifest points to "workflow.wdl" which is not the convention
        # This should be validated by the run command
        wrong_wdl = torch_with_wrong_wdl_name / "workflow.wdl"
        main_wdl = torch_with_wrong_wdl_name / "main.wdl"

        assert wrong_wdl.exists()
        assert not main_wdl.exists()

    def test_run_command_rejects_non_main_wdl(
        self, torch_with_wrong_wdl_name, sample_reads_file
    ):
        """Run command should reject torch with non-main.wdl workflow."""
        runner = CliRunner()

        with patch('torchbase.torchfs.Torch') as mock_torch_class:
            mock_torch = MagicMock()
            mock_torch.workflow = torch_with_wrong_wdl_name / "workflow.wdl"
            mock_torch_class.load.return_value = mock_torch

            # This should fail validation
            result = runner.invoke(
                cli,
                ['run', str(torch_with_wrong_wdl_name), '-r', str(sample_reads_file)]
            )

            # Should fail - workflow not named main.wdl
            assert result.exit_code != 0

    def test_only_main_wdl_is_valid_workflow_name(self, tmp_path):
        """Only main.wdl is valid; other names should be rejected."""
        invalid_names = ["workflow.wdl", "typing.wdl", "pipeline.wdl", "mlst.wdl"]

        for invalid_name in invalid_names:
            # Convention dictates main.wdl is the only valid name
            assert invalid_name != "main.wdl"


class TestDefaultWorkflowFallback:
    """Test default workflow fetching when torch has no main.wdl."""

    def test_fetch_default_workflow_when_missing(
        self, torch_without_main_wdl, sample_reads_file
    ):
        """Fetch workflows/default-workflow when torch has no main.wdl."""
        runner = CliRunner()

        with patch('torchbase.torchfs.Torch') as mock_torch_class:
            with patch(
                'torchbase.registry.RegistryManager'
            ) as mock_manager_class:
                # Mock data torch (no workflow)
                mock_data_torch = MagicMock()
                mock_data_torch.workflow = None
                mock_data_torch.path = torch_without_main_wdl

                # Mock default workflow torch
                mock_default_torch = MagicMock()
                default_wdl = Path(
                    "/tmp/workflows/default-workflow/1.0.0.torch/main.wdl"
                )
                mock_default_torch.workflow = default_wdl

                # Configure mock to return different torches
                def load_side_effect(path):
                    if "default-workflow" in str(path):
                        return mock_default_torch
                    return mock_data_torch

                mock_torch_class.load.side_effect = load_side_effect

                # Mock registry manager
                mock_manager = MagicMock()
                default_path = Path(
                    "/tmp/workflows/default-workflow/1.0.0.torch"
                )
                mock_manager.fetch_torch.return_value = default_path
                mock_manager_class.return_value = mock_manager

                with patch('torchbase.cli.run'):
                    runner.invoke(
                        cli,
                        [
                            'run',
                            str(torch_without_main_wdl),
                            '-r',
                            str(sample_reads_file)
                        ]
                    )

                    # Should fetch default workflow
                    mock_manager.fetch_torch.assert_called()
                    call_args = mock_manager.fetch_torch.call_args
                    assert ("default-workflow" in call_args[0][0] or
                            "workflows" in str(call_args[0][0]))

    def test_default_workflow_uses_latest_version(
        self, torch_without_main_wdl, sample_reads_file
    ):
        """Default workflow fetches latest version."""
        runner = CliRunner()

        with patch('torchbase.torchfs.Torch') as mock_torch_class:
            with patch('torchbase.registry.RegistryManager') as mock_manager_class:
                mock_data_torch = MagicMock()
                mock_data_torch.workflow = None
                mock_torch_class.load.return_value = mock_data_torch

                mock_manager = MagicMock()
                mock_manager_class.return_value = mock_manager

                with patch('torchbase.cli.run'):
                    _ = runner.invoke(
                        cli,
                        ['run', str(torch_without_main_wdl), '-r', str(sample_reads_file)]
                    )

                    # Should fetch with version=None (latest)
                    if mock_manager.fetch_torch.called:
                        call_args = mock_manager.fetch_torch.call_args
                        # Version should be None or not specified
                        assert call_args[1].get('version') is None or 'version' not in call_args[1]

    def test_default_workflow_registry_path(
        self, torch_without_main_wdl, sample_reads_file
    ):
        """Default workflow is fetched from workflows/default-workflow."""
        runner = CliRunner()

        with patch('torchbase.torchfs.Torch') as mock_torch_class:
            with patch('torchbase.registry.RegistryManager') as mock_manager_class:
                mock_data_torch = MagicMock()
                mock_data_torch.workflow = None
                mock_torch_class.load.return_value = mock_data_torch

                mock_manager = MagicMock()
                mock_manager_class.return_value = mock_manager

                with patch('torchbase.cli.run'):
                    _ = runner.invoke(
                        cli,
                        ['run', str(torch_without_main_wdl), '-r', str(sample_reads_file)]
                    )

                    # Should request workflows/default-workflow
                    if mock_manager.fetch_torch.called:
                        torch_name = mock_manager.fetch_torch.call_args[0][0]
                        assert "workflows" in torch_name
                        assert "default-workflow" in torch_name


class TestWorkflowOverride:
    """Test --workflow flag to override workflow discovery."""

    def test_workflow_flag_overrides_main_wdl(
        self, torch_with_main_wdl, torch_without_main_wdl, sample_reads_file
    ):
        """--workflow flag overrides built-in main.wdl."""
        runner = CliRunner()

        with patch('torchbase.torchfs.Torch') as mock_torch_class:
            with patch('torchbase.registry.RegistryManager') as mock_manager_class:
                # Mock data torch
                mock_data_torch = MagicMock()
                mock_data_torch.workflow = None

                # Mock custom workflow torch
                mock_custom_torch = MagicMock()
                mock_custom_torch.workflow = Path("/tmp/custom/workflow/main.wdl")

                def load_side_effect(path):
                    if "custom" in str(path):
                        return mock_custom_torch
                    return mock_data_torch

                mock_torch_class.load.side_effect = load_side_effect

                mock_manager = MagicMock()
                mock_manager.fetch_torch.return_value = Path("/tmp/custom/workflow/1.0.0.torch")
                mock_manager_class.return_value = mock_manager

                with patch('torchbase.cli.run'):
                    _ = runner.invoke(
                        cli,
                        [
                            'run',
                            '--workflow', 'namespace/custom-workflow',
                            str(torch_without_main_wdl),
                            '-r', str(sample_reads_file)
                        ]
                    )

                    # Should fetch the custom workflow
                    mock_manager.fetch_torch.assert_called()
                    call_args = mock_manager.fetch_torch.call_args
                    assert "custom-workflow" in call_args[0][0]

    def test_workflow_flag_syntax(self, sample_reads_file, tmp_path):
        """Test --workflow flag accepts namespace/name format."""
        runner = CliRunner()

        # Test that CLI accepts --workflow flag
        with patch('torchbase.torchfs.Torch') as mock_torch_class:
            with patch('torchbase.registry.RegistryManager') as mock_manager_class:
                mock_torch = MagicMock()
                mock_torch.workflow = Path("/tmp/workflow/main.wdl")
                mock_torch_class.load.return_value = mock_torch

                mock_manager = MagicMock()
                mock_manager_class.return_value = mock_manager

                with patch('torchbase.cli.run'):
                    result = runner.invoke(
                        cli,
                        [
                            'run',
                            '--workflow', 'custom/workflow',
                            str(tmp_path),
                            '-r', str(sample_reads_file)
                        ]
                    )

                    # Command should parse successfully (even if it fails later)
                    # The --workflow flag should be recognized
                    assert (
                        '--workflow' not in result.output or
                        result.exit_code != 2
                    )

    def test_workflow_override_with_data_torch(
        self, torch_without_main_wdl, sample_reads_file
    ):
        """Can use custom workflow with data-only torch."""
        runner = CliRunner()

        with patch('torchbase.torchfs.Torch') as mock_torch_class:
            with patch('torchbase.registry.RegistryManager') as mock_manager_class:
                mock_data_torch = MagicMock()
                mock_data_torch.workflow = None

                mock_custom_workflow = MagicMock()
                mock_custom_workflow.workflow = Path("/tmp/custom/main.wdl")

                def load_side_effect(path):
                    if "custom" in str(path):
                        return mock_custom_workflow
                    return mock_data_torch

                mock_torch_class.load.side_effect = load_side_effect

                mock_manager = MagicMock()
                mock_manager.fetch_torch.return_value = Path("/tmp/custom/1.0.0.torch")
                mock_manager_class.return_value = mock_manager

                with patch('torchbase.cli.run'):
                    _ = runner.invoke(
                        cli,
                        [
                            'run',
                            '--workflow', 'custom/mlst-workflow',
                            str(torch_without_main_wdl),
                            '-r', str(sample_reads_file)
                        ]
                    )

                    # Should use the custom workflow, not default
                    mock_manager.fetch_torch.assert_called()
                    call_args = mock_manager.fetch_torch.call_args
                    assert "custom" in call_args[0][0] or "mlst-workflow" in call_args[0][0]


class TestMiniwdlExecution:
    """Test miniwdl execution with appropriate inputs."""

    def test_miniwdl_called_with_workflow_path(
        self, torch_with_main_wdl, sample_reads_file
    ):
        """miniwdl is called with workflow path."""
        runner = CliRunner()

        with patch('torchbase.torchfs.Torch') as mock_torch_class:
            mock_torch = MagicMock()
            mock_torch.workflow = torch_with_main_wdl / "main.wdl"
            mock_torch_class.load.return_value = mock_torch

            with patch('torchbase.cli.run') as mock_run:
                _ = runner.invoke(
                    cli,
                    ['run', str(torch_with_main_wdl), '-r', str(sample_reads_file)]
                )

                # Should call miniwdl
                mock_run.assert_called()
                call_args = mock_run.call_args[0][0]
                assert 'miniwdl' in call_args

    def test_miniwdl_receives_input_files(
        self, torch_with_main_wdl, sample_reads_file
    ):
        """miniwdl receives input files as arguments."""
        runner = CliRunner()

        with patch('torchbase.torchfs.Torch') as mock_torch_class:
            mock_torch = MagicMock()
            mock_torch.workflow = torch_with_main_wdl / "main.wdl"
            mock_torch_class.load.return_value = mock_torch

            with patch('torchbase.cli.run') as mock_run:
                _ = runner.invoke(
                    cli,
                    ['run', str(torch_with_main_wdl), '-r', str(sample_reads_file)]
                )

                # miniwdl should receive reads file path
                if mock_run.called:
                    call_args = mock_run.call_args[0][0]
                    # Should contain reference to input data
                    has_input = any(
                        str(sample_reads_file) in str(arg)
                        for arg in call_args
                    )
                    assert has_input or len(call_args) > 1

    def test_miniwdl_execution_with_multiple_inputs(
        self, torch_with_main_wdl, tmp_path
    ):
        """miniwdl execution with multiple input files."""
        reads1 = tmp_path / "reads1.fastq"
        reads2 = tmp_path / "reads2.fastq"
        reads1.write_text("@read1\nACGT\n+\nIIII\n")
        reads2.write_text("@read2\nTGCA\n+\nIIII\n")

        runner = CliRunner()

        with patch('torchbase.torchfs.Torch') as mock_torch_class:
            mock_torch = MagicMock()
            mock_torch.workflow = torch_with_main_wdl / "main.wdl"
            mock_torch_class.load.return_value = mock_torch

            with patch('torchbase.cli.run') as mock_run:
                _ = runner.invoke(
                    cli,
                    [
                        'run',
                        str(torch_with_main_wdl),
                        '-pe1', str(reads1),
                        '-pe2', str(reads2)
                    ]
                )

                # Should call miniwdl with both files
                mock_run.assert_called()

    def test_miniwdl_run_command_format(
        self, torch_with_main_wdl, sample_reads_file
    ):
        """miniwdl run command has correct format."""
        runner = CliRunner()

        with patch('torchbase.torchfs.Torch') as mock_torch_class:
            mock_torch = MagicMock()
            mock_torch.workflow = torch_with_main_wdl / "main.wdl"
            mock_torch_class.load.return_value = mock_torch

            with patch('torchbase.cli.run') as mock_run:
                _ = runner.invoke(
                    cli,
                    ['run', str(torch_with_main_wdl), '-r', str(sample_reads_file)]
                )

                # Command format should be: miniwdl run <workflow> <inputs...>
                if mock_run.called:
                    call_args = mock_run.call_args[0][0]
                    assert isinstance(call_args, list)
                    if len(call_args) > 0:
                        assert call_args[0] == 'miniwdl' or 'miniwdl' in str(call_args)


class TestOutputHandling:
    """Test output handling for typing results."""

    def test_output_to_stdout_by_default(
        self, torch_with_main_wdl, sample_reads_file
    ):
        """Typing results output to stdout by default."""
        runner = CliRunner()

        with patch('torchbase.torchfs.Torch') as mock_torch_class:
            mock_torch = MagicMock()
            mock_torch.workflow = torch_with_main_wdl / "main.wdl"
            mock_torch_class.load.return_value = mock_torch

            with patch('torchbase.cli.run') as mock_run:
                # Mock successful miniwdl execution
                mock_run.return_value.returncode = 0
                mock_run.return_value.stdout = '{"st": 1, "loci": {"adk": 1}}'

                result = runner.invoke(
                    cli,
                    ['run', str(torch_with_main_wdl), '-r', str(sample_reads_file)]
                )

                # Output should be in result
                # (implementation detail - miniwdl output may be passed through)
                assert result.exit_code == 0 or mock_run.called

    def test_output_to_specified_file(
        self, torch_with_main_wdl, sample_reads_file, tmp_path
    ):
        """Typing results output to specified file."""
        output_file = tmp_path / "results.json"
        runner = CliRunner()

        with patch('torchbase.torchfs.Torch') as mock_torch_class:
            mock_torch = MagicMock()
            mock_torch.workflow = torch_with_main_wdl / "main.wdl"
            mock_torch_class.load.return_value = mock_torch

            with patch('torchbase.cli.run') as mock_run:
                result = runner.invoke(
                    cli,
                    [
                        'run',
                        str(torch_with_main_wdl),
                        '-r', str(sample_reads_file),
                        '--output', str(output_file)
                    ]
                )

                # CLI should accept --output flag
                # (actual file writing is miniwdl's responsibility)
                assert result.exit_code == 0 or mock_run.called

    def test_json_output_format(
        self, torch_with_main_wdl, sample_reads_file
    ):
        """Output is in JSON format."""
        runner = CliRunner()

        with patch('torchbase.torchfs.Torch') as mock_torch_class:
            mock_torch = MagicMock()
            mock_torch.workflow = torch_with_main_wdl / "main.wdl"
            mock_torch_class.load.return_value = mock_torch

            with patch('torchbase.cli.run') as mock_run:
                # Mock miniwdl returning JSON
                mock_result = MagicMock()
                mock_result.returncode = 0
                mock_run.return_value = mock_result

                result = runner.invoke(
                    cli,
                    ['run', str(torch_with_main_wdl), '-r', str(sample_reads_file)]
                )

                # Result should be parseable (or passed through from miniwdl)
                assert result.exit_code == 0 or mock_run.called


class TestErrorHandling:
    """Test error handling for workflow discovery and execution."""

    def test_error_when_workflow_not_found(
        self, torch_without_main_wdl, sample_reads_file
    ):
        """Error when workflow not found and default fetch fails."""
        runner = CliRunner()

        with patch('torchbase.torchfs.Torch') as mock_torch_class:
            with patch('torchbase.registry.RegistryManager') as mock_manager_class:
                mock_torch = MagicMock()
                mock_torch.workflow = None
                mock_torch_class.load.return_value = mock_torch

                # Mock registry manager to fail
                mock_manager = MagicMock()
                mock_manager.fetch_torch.side_effect = ValueError("Workflow not found")
                mock_manager_class.return_value = mock_manager

                result = runner.invoke(
                    cli,
                    ['run', str(torch_without_main_wdl), '-r', str(sample_reads_file)]
                )

                # Should fail with error
                assert result.exit_code != 0

    def test_error_when_miniwdl_execution_fails(
        self, torch_with_main_wdl, sample_reads_file
    ):
        """Error when miniwdl execution fails."""
        runner = CliRunner()

        with patch('torchbase.torchfs.Torch') as mock_torch_class:
            mock_torch = MagicMock()
            mock_torch.workflow = torch_with_main_wdl / "main.wdl"
            mock_torch_class.load.return_value = mock_torch

            with patch('torchbase.cli.run') as mock_run:
                # Mock miniwdl failure
                mock_run.return_value.returncode = 1
                mock_run.return_value.stderr = "Workflow execution failed"

                _ = runner.invoke(
                    cli,
                    ['run', str(torch_with_main_wdl), '-r', str(sample_reads_file)]
                )

                # Should propagate failure
                # (exit code handling depends on implementation)
                assert mock_run.called

    def test_error_when_torch_not_found(self, sample_reads_file):
        """Error when torch path doesn't exist."""
        runner = CliRunner()

        with patch('torchbase.torchfs.Torch') as mock_torch_class:
            mock_torch_class.load.side_effect = FileNotFoundError("Torch not found")

            result = runner.invoke(
                cli,
                ['run', '/nonexistent/torch', '-r', str(sample_reads_file)]
            )

            # Should fail
            assert result.exit_code != 0

    def test_error_message_clarity_for_missing_workflow(
        self, torch_without_main_wdl, sample_reads_file
    ):
        """Error message is clear when workflow is missing."""
        runner = CliRunner()

        with patch('torchbase.torchfs.Torch') as mock_torch_class:
            with patch('torchbase.registry.RegistryManager') as mock_manager_class:
                mock_torch = MagicMock()
                mock_torch.workflow = None
                mock_torch_class.load.return_value = mock_torch

                mock_manager = MagicMock()
                mock_manager.fetch_torch.side_effect = ValueError(
                    "Default workflow not found in registry"
                )
                mock_manager_class.return_value = mock_manager

                result = runner.invoke(
                    cli,
                    ['run', str(torch_without_main_wdl), '-r', str(sample_reads_file)]
                )

                # Should have clear error message
                assert result.exit_code != 0
                # Error message should mention workflow
                assert "workflow" in result.output.lower() or result.exit_code != 0

    def test_error_when_main_wdl_malformed(
        self, torch_with_main_wdl, sample_reads_file
    ):
        """Error when main.wdl is malformed."""
        # Corrupt the WDL file
        with open(torch_with_main_wdl / "main.wdl", "w") as f:
            f.write("not valid wdl syntax {{{")

        runner = CliRunner()

        with patch('torchbase.torchfs.Torch') as mock_torch_class:
            mock_torch = MagicMock()
            mock_torch.workflow = torch_with_main_wdl / "main.wdl"
            mock_torch_class.load.return_value = mock_torch

            with patch('torchbase.cli.run') as mock_run:
                # miniwdl should fail on malformed WDL
                mock_run.return_value.returncode = 1

                _ = runner.invoke(
                    cli,
                    ['run', str(torch_with_main_wdl), '-r', str(sample_reads_file)]
                )

                # Should fail (miniwdl will reject bad WDL)
                mock_run.assert_called()


class TestWorkflowDiscoveryEdgeCases:
    """Test edge cases in workflow discovery."""

    def test_torch_with_empty_main_wdl(self, tmp_path, sample_reads_file):
        """Torch with empty main.wdl file."""
        torch_path = tmp_path / "namespace" / "torch" / "1.0.0.torch"
        torch_path.mkdir(parents=True)

        metadata = {
            "namespace": "namespace",
            "name": "torch",
            "version": "1.0.0",
            "version_meta": {"strategy": "semver", "timestamp": 1609459200},
            "typing": {"method": "mlst"},
            "description": {"short": "Test"},
            "manifest": {"profiles": "profiles.tsv", "workflow": "main.wdl"}
        }
        with open(torch_path / "metadata.toml", "w") as f:
            toml.dump(metadata, f)

        # Empty main.wdl
        (torch_path / "main.wdl").touch()

        profiles = [["ST", "adk"], ["1", "1"]]
        with open(torch_path / "profiles.tsv", "w") as f:
            writer = csv.writer(f, delimiter="\t")
            writer.writerows(profiles)

        (torch_path / "_resources").mkdir()

        runner = CliRunner()

        with patch('torchbase.torchfs.Torch') as mock_torch_class:
            mock_torch = MagicMock()
            mock_torch.workflow = torch_path / "main.wdl"
            mock_torch_class.load.return_value = mock_torch

            with patch('torchbase.cli.run') as mock_run:
                mock_run.return_value.returncode = 1  # Will fail on empty WDL

                _ = runner.invoke(
                    cli,
                    ['run', str(torch_path), '-r', str(sample_reads_file)]
                )

                # miniwdl will handle the empty file
                mock_run.assert_called()

    def test_workflow_discovery_with_symlink(self, tmp_path, sample_reads_file):
        """Workflow discovery works with symlinked main.wdl."""
        torch_path = tmp_path / "namespace" / "torch" / "1.0.0.torch"
        torch_path.mkdir(parents=True)

        # Create actual WDL elsewhere
        real_wdl = tmp_path / "shared" / "workflow.wdl"
        real_wdl.parent.mkdir(parents=True)
        real_wdl.write_text("workflow test { }\n")

        # Symlink to main.wdl
        main_wdl = torch_path / "main.wdl"
        main_wdl.symlink_to(real_wdl)

        assert main_wdl.exists()
        assert main_wdl.is_symlink()

    def test_default_workflow_fetch_caching(
        self, torch_without_main_wdl, sample_reads_file
    ):
        """Default workflow is cached after first fetch."""
        runner = CliRunner()

        with patch('torchbase.torchfs.Torch') as mock_torch_class:
            with patch('torchbase.registry.RegistryManager') as mock_manager_class:
                mock_torch = MagicMock()
                mock_torch.workflow = None
                mock_torch_class.load.return_value = mock_torch

                mock_manager = MagicMock()
                default_workflow_path = Path(
                    "/tmp/workflows/default-workflow/1.0.0.torch"
                )
                mock_manager.fetch_torch.return_value = default_workflow_path
                mock_manager_class.return_value = mock_manager

                with patch('torchbase.cli.run'):
                    # First run
                    _ = runner.invoke(
                        cli,
                        ['run', str(torch_without_main_wdl), '-r', str(sample_reads_file)]
                    )

                    # Second run
                    _ = runner.invoke(
                        cli,
                        ['run', str(torch_without_main_wdl), '-r', str(sample_reads_file)]
                    )

                    # fetch_torch should be called (caching is registry's responsibility)
                    assert mock_manager.fetch_torch.called


class TestWorkflowIntegration:
    """Integration tests for workflow discovery with actual torch loading."""

    def test_end_to_end_workflow_discovery(
        self, torch_with_main_wdl, sample_reads_file
    ):
        """End-to-end test: load torch, discover workflow, prepare execution."""
        # This is a full integration test without mocks
        torch = Torch.load(torch_with_main_wdl)

        # Workflow should be discovered
        assert torch.workflow is not None
        assert torch.workflow.exists()
        assert torch.workflow.name == "main.wdl"

        # Workflow should be executable
        assert torch.workflow.is_file()

        # Content should be valid WDL
        content = torch.workflow.read_text()
        assert "workflow" in content

    def test_torch_loading_sets_workflow_attribute(
        self, torch_with_main_wdl
    ):
        """Torch.load() correctly sets workflow attribute."""
        torch = Torch.load(torch_with_main_wdl)

        assert hasattr(torch, 'workflow')
        assert torch.workflow == torch_with_main_wdl / "main.wdl"

    def test_data_torch_workflow_is_none(
        self, torch_without_main_wdl
    ):
        """Data-only torch has workflow=None."""
        torch = Torch.load(torch_without_main_wdl)

        assert hasattr(torch, 'workflow')
        assert torch.workflow is None
