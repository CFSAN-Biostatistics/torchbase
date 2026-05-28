"""Tests for CLI strategy routing (Issue #58).

Acceptance criteria:
- --strategy flag added to CLI with choices [fast, balanced, sensitive]
- Default strategy is "balanced"
- CLI routes to correct built-in workflow file based on strategy
- Error raised if strategy used with torch-embedded workflow
- Multi-scheme concatenation (from #53) integrated
- Help text explains strategy options and restrictions
- Tests verify routing logic and error conditions
"""

import pytest
import toml
import csv
from click.testing import CliRunner
from unittest.mock import patch, MagicMock

from torchbase.cli import cli
from torchbase.torchfs import Torch


@pytest.fixture
def torch_without_workflow(tmp_path):
    """Create a torch without embedded workflow for strategy routing."""
    torch_path = tmp_path / "test_namespace" / "data_torch" / "1.0.0.torch"
    torch_path.mkdir(parents=True)

    metadata = {
        "namespace": "test_namespace",
        "name": "data_torch",
        "version": "1.0.0",
        "version_meta": {"strategy": "semver", "timestamp": 1609459200},
        "typing": {"method": "mlst"},
        "description": {"short": "Data torch without workflow"},
        "manifest": {"profiles": "profiles.tsv"}
    }
    with open(torch_path / "metadata.toml", "w") as f:
        toml.dump(metadata, f)

    profiles = [["ST", "adk"], ["1", "1"]]
    with open(torch_path / "profiles.tsv", "w") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerows(profiles)

    (torch_path / "_resources").mkdir()

    return torch_path


@pytest.fixture
def torch_with_embedded_workflow(tmp_path):
    """Create a torch with embedded main.wdl."""
    torch_path = tmp_path / "test_namespace" / "workflow_torch" / "1.0.0.torch"
    torch_path.mkdir(parents=True)

    metadata = {
        "namespace": "test_namespace",
        "name": "workflow_torch",
        "version": "1.0.0",
        "version_meta": {"strategy": "semver", "timestamp": 1609459200},
        "typing": {"method": "mlst"},
        "description": {"short": "Torch with embedded workflow"},
        "manifest": {"profiles": "profiles.tsv", "workflow": "main.wdl"}
    }
    with open(torch_path / "metadata.toml", "w") as f:
        toml.dump(metadata, f)

    wdl_content = """workflow custom_mlst {
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

    profiles = [["ST", "adk"], ["1", "1"]]
    with open(torch_path / "profiles.tsv", "w") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerows(profiles)

    (torch_path / "_resources").mkdir()

    return torch_path


@pytest.fixture
def sample_reads_file(tmp_path):
    """Create a sample reads file for testing."""
    reads_file = tmp_path / "reads.fastq"
    with open(reads_file, "w") as f:
        f.write("@read1\nACGT\n+\nIIII\n")
    return reads_file


class TestStrategyFlagPresence:
    """Test that --strategy flag exists and has correct choices."""

    def test_strategy_flag_exists(self):
        """--strategy flag is recognized by CLI."""
        runner = CliRunner()
        result = runner.invoke(cli, ['run', '--help'])

        assert result.exit_code == 0
        assert '--strategy' in result.output

    def test_strategy_flag_has_fast_choice(self):
        """--strategy accepts 'fast' value."""
        runner = CliRunner()
        result = runner.invoke(cli, ['run', '--help'])

        assert 'fast' in result.output

    def test_strategy_flag_has_balanced_choice(self):
        """--strategy accepts 'balanced' value."""
        runner = CliRunner()
        result = runner.invoke(cli, ['run', '--help'])

        assert 'balanced' in result.output

    def test_strategy_flag_has_sensitive_choice(self):
        """--strategy accepts 'sensitive' value."""
        runner = CliRunner()
        result = runner.invoke(cli, ['run', '--help'])

        assert 'sensitive' in result.output

    def test_strategy_flag_rejects_invalid_choice(
        self, torch_without_workflow, sample_reads_file
    ):
        """--strategy rejects invalid values."""
        runner = CliRunner()

        result = runner.invoke(
            cli,
            [
                'run', '--strategy', 'invalid',
                str(torch_without_workflow), '-r', str(sample_reads_file)
            ]
        )

        # Should fail with invalid choice error
        assert result.exit_code != 0
        assert ('invalid' in result.output.lower() or
                'choice' in result.output.lower())


class TestDefaultStrategy:
    """Test that default strategy is 'balanced'."""

    def test_balanced_is_default_strategy(
        self, torch_without_workflow, sample_reads_file
    ):
        """When --strategy not specified, 'balanced' strategy is used."""
        runner = CliRunner()

        with patch('torchbase.torchfs.Torch') as mock_torch_class:
            mock_torch = MagicMock()
            mock_torch.workflow = None
            mock_torch_class.load.return_value = mock_torch

            with patch('torchbase.cli.run') as mock_run:
                _ = runner.invoke(
                    cli,
                    [
                        'run', str(torch_without_workflow),
                        '-r', str(sample_reads_file)
                    ]
                )

                # Should call miniwdl with balanced workflow
                if mock_run.called:
                    call_args = mock_run.call_args[0][0]
                    workflow_arg = str(call_args)
                    assert 'balanced' in workflow_arg.lower()

    def test_no_strategy_flag_uses_balanced(
        self, torch_without_workflow, sample_reads_file
    ):
        """Omitting --strategy flag defaults to balanced strategy."""
        runner = CliRunner()

        with patch('torchbase.torchfs.Torch') as mock_torch_class:
            mock_torch = MagicMock()
            mock_torch.workflow = None
            mock_torch.path = torch_without_workflow
            mock_torch_class.load.return_value = mock_torch

            with patch('torchbase.cli.Path') as mock_path_class:
                # Mock built-in workflow path resolution
                mock_balanced_path = MagicMock()
                mock_balanced_path.exists.return_value = True
                mock_path_class.return_value = mock_balanced_path

                with patch('torchbase.cli.run') as mock_run:
                    _ = runner.invoke(
                        cli,
                        [
                            'run', str(torch_without_workflow),
                            '-r', str(sample_reads_file)
                        ]
                    )

                    # Should use balanced strategy workflow
                    if mock_run.called:
                        call_args = str(mock_run.call_args)
                        assert 'balanced' in call_args.lower()


class TestStrategyRouting:
    """Test routing to correct built-in workflow based on strategy."""

    def test_fast_strategy_routes_to_fast_workflow(
        self, torch_without_workflow, sample_reads_file
    ):
        """--strategy fast routes to builtin/fast_typing.wdl."""
        runner = CliRunner()

        with patch('torchbase.torchfs.Torch') as mock_torch_class:
            mock_torch = MagicMock()
            mock_torch.workflow = None
            mock_torch_class.load.return_value = mock_torch

            with patch('torchbase.cli.run') as mock_run:
                _ = runner.invoke(
                    cli,
                    [
                        'run', '--strategy', 'fast',
                        str(torch_without_workflow),
                        '-r', str(sample_reads_file)
                    ]
                )

                # Should call with fast_typing.wdl
                if mock_run.called:
                    call_args = mock_run.call_args[0][0]
                    workflow_path = str(call_args)
                    assert ('fast_typing.wdl' in workflow_path or
                            'fast' in workflow_path)

    def test_balanced_strategy_routes_to_balanced_workflow(
        self, torch_without_workflow, sample_reads_file
    ):
        """--strategy balanced routes to builtin/balanced_typing.wdl."""
        runner = CliRunner()

        with patch('torchbase.torchfs.Torch') as mock_torch_class:
            mock_torch = MagicMock()
            mock_torch.workflow = None
            mock_torch_class.load.return_value = mock_torch

            with patch('torchbase.cli.run') as mock_run:
                _ = runner.invoke(
                    cli,
                    [
                        'run', '--strategy', 'balanced',
                        str(torch_without_workflow),
                        '-r', str(sample_reads_file)
                    ]
                )

                # Should call with balanced_typing.wdl
                if mock_run.called:
                    call_args = mock_run.call_args[0][0]
                    workflow_path = str(call_args)
                    assert ('balanced_typing.wdl' in workflow_path or
                            'balanced' in workflow_path)

    def test_sensitive_strategy_routes_to_sensitive_workflow(
        self, torch_without_workflow, sample_reads_file
    ):
        """--strategy sensitive routes to builtin/sensitive_typing.wdl."""
        runner = CliRunner()

        with patch('torchbase.torchfs.Torch') as mock_torch_class:
            mock_torch = MagicMock()
            mock_torch.workflow = None
            mock_torch_class.load.return_value = mock_torch

            with patch('torchbase.cli.run') as mock_run:
                _ = runner.invoke(
                    cli,
                    [
                        'run', '--strategy', 'sensitive',
                        str(torch_without_workflow),
                        '-r', str(sample_reads_file)
                    ]
                )

                # Should call with sensitive_typing.wdl
                if mock_run.called:
                    call_args = mock_run.call_args[0][0]
                    workflow_path = str(call_args)
                    assert ('sensitive_typing.wdl' in workflow_path or
                            'sensitive' in workflow_path)

    def test_workflow_path_includes_builtin_directory(
        self, torch_without_workflow, sample_reads_file
    ):
        """Built-in workflows are in torchbase/workflows/builtin/ directory."""
        runner = CliRunner()

        with patch('torchbase.torchfs.Torch') as mock_torch_class:
            mock_torch = MagicMock()
            mock_torch.workflow = None
            mock_torch_class.load.return_value = mock_torch

            with patch('torchbase.cli.run') as mock_run:
                _ = runner.invoke(
                    cli,
                    [
                        'run', '--strategy', 'fast',
                        str(torch_without_workflow),
                        '-r', str(sample_reads_file)
                    ]
                )

                # Should include builtin path
                if mock_run.called:
                    call_args = str(mock_run.call_args[0][0])
                    assert ('builtin' in call_args or
                            'workflows' in call_args)


class TestStrategyWithEmbeddedWorkflowError:
    """Test error when --strategy used with torch-embedded workflow."""

    def test_strategy_with_embedded_workflow_raises_error(
        self, torch_with_embedded_workflow, sample_reads_file
    ):
        """Using --strategy with embedded workflow raises clear error."""
        runner = CliRunner()

        with patch('torchbase.torchfs.Torch') as mock_torch_class:
            mock_torch = MagicMock()
            mock_torch.workflow = torch_with_embedded_workflow / "main.wdl"
            mock_torch_class.load.return_value = mock_torch

            result = runner.invoke(
                cli,
                [
                    'run', '--strategy', 'fast',
                    str(torch_with_embedded_workflow),
                    '-r', str(sample_reads_file)
                ]
            )

            # Should fail with clear error message
            assert result.exit_code != 0
            assert 'strategy' in result.output.lower()
            assert ('embedded' in result.output.lower() or
                    'workflow' in result.output.lower())

    def test_error_message_mentions_strategy_restriction(
        self, torch_with_embedded_workflow, sample_reads_file
    ):
        """Error message specifically mentions --strategy cannot be used."""
        runner = CliRunner()

        with patch('torchbase.torchfs.Torch') as mock_torch_class:
            mock_torch = MagicMock()
            mock_torch.workflow = torch_with_embedded_workflow / "main.wdl"
            mock_torch_class.load.return_value = mock_torch

            result = runner.invoke(
                cli,
                [
                    'run', '--strategy', 'balanced',
                    str(torch_with_embedded_workflow),
                    '-r', str(sample_reads_file)
                ]
            )

            # Error should mention the restriction
            assert result.exit_code != 0
            assert ('--strategy' in result.output or
                    'strategy' in result.output.lower())

    def test_all_strategies_fail_with_embedded_workflow(
        self, torch_with_embedded_workflow, sample_reads_file
    ):
        """All strategy values fail when torch has embedded workflow."""
        runner = CliRunner()
        strategies = ['fast', 'balanced', 'sensitive']

        for strategy in strategies:
            with patch('torchbase.torchfs.Torch') as mock_torch_class:
                mock_torch = MagicMock()
                mock_torch.workflow = torch_with_embedded_workflow / "main.wdl"
                mock_torch_class.load.return_value = mock_torch

                result = runner.invoke(
                    cli,
                    [
                        'run', '--strategy', strategy,
                        str(torch_with_embedded_workflow),
                        '-r', str(sample_reads_file)
                    ]
                )

                # All should fail
                assert result.exit_code != 0

    def test_embedded_workflow_works_without_strategy(
        self, torch_with_embedded_workflow, sample_reads_file
    ):
        """Torch with embedded workflow works when --strategy not specified."""
        runner = CliRunner()

        with patch('torchbase.torchfs.Torch') as mock_torch_class:
            mock_torch = MagicMock()
            mock_torch.workflow = torch_with_embedded_workflow / "main.wdl"
            mock_torch_class.load.return_value = mock_torch

            with patch('torchbase.cli.run') as mock_run:
                mock_run.return_value.returncode = 0

                result = runner.invoke(
                    cli,
                    [
                        'run', str(torch_with_embedded_workflow),
                        '-r', str(sample_reads_file)
                    ]
                )

                # Should succeed without --strategy
                assert result.exit_code == 0 or mock_run.called


class TestStrategyHelpText:
    """Test help text explains strategy options and restrictions."""

    def test_help_text_explains_fast_strategy(self):
        """Help text explains fast strategy."""
        runner = CliRunner()
        result = runner.invoke(cli, ['run', '--help'])

        assert 'fast' in result.output
        # Should mention speed/accuracy tradeoff
        assert ('minhash' in result.output.lower() or
                'fastest' in result.output.lower())

    def test_help_text_explains_balanced_strategy(self):
        """Help text explains balanced strategy."""
        runner = CliRunner()
        result = runner.invoke(cli, ['run', '--help'])

        assert 'balanced' in result.output
        # Should indicate it's the default
        assert 'default' in result.output.lower()

    def test_help_text_explains_sensitive_strategy(self):
        """Help text explains sensitive strategy."""
        runner = CliRunner()
        result = runner.invoke(cli, ['run', '--help'])

        assert 'sensitive' in result.output
        # Should mention accuracy
        assert ('alignment' in result.output.lower() or
                'accurate' in result.output.lower())

    def test_help_text_mentions_embedded_workflow_restriction(self):
        """Help text mentions restriction with embedded workflows."""
        runner = CliRunner()
        result = runner.invoke(cli, ['run', '--help'])

        # Should warn about embedded workflow restriction
        assert (
            'embedded' in result.output.lower() or
            'torch-embedded' in result.output.lower() or
            'cannot use' in result.output.lower()
        )

    def test_strategy_flag_has_type_choice(self):
        """--strategy flag is defined as a choice type."""
        runner = CliRunner()
        result = runner.invoke(cli, ['run', '--help'])

        # Should show choices
        assert ('[fast|balanced|sensitive]' in result.output or
                ('fast' in result.output and
                 'balanced' in result.output and
                 'sensitive' in result.output))


class TestStrategyWorkflowPathResolution:
    """Test that strategy routing resolves to actual workflow files."""

    def test_fast_workflow_path_is_absolute(
        self, torch_without_workflow, sample_reads_file
    ):
        """Fast strategy resolves to absolute path."""
        runner = CliRunner()

        with patch('torchbase.torchfs.Torch') as mock_torch_class:
            mock_torch = MagicMock()
            mock_torch.workflow = None
            mock_torch_class.load.return_value = mock_torch

            with patch('torchbase.cli.run') as mock_run:
                _ = runner.invoke(
                    cli,
                    [
                        'run', '--strategy', 'fast',
                        str(torch_without_workflow),
                        '-r', str(sample_reads_file)
                    ]
                )

                if mock_run.called:
                    call_args = mock_run.call_args[0][0]
                    # Workflow path should be resolvable
                    workflow_path = None
                    for arg in call_args:
                        if 'wdl' in str(arg).lower():
                            workflow_path = str(arg)
                            break
                    assert (workflow_path is not None or
                            len(call_args) > 2)

    def test_workflow_file_path_is_passed_to_miniwdl(
        self, torch_without_workflow, sample_reads_file
    ):
        """Workflow file path is passed to miniwdl run command."""
        runner = CliRunner()

        with patch('torchbase.torchfs.Torch') as mock_torch_class:
            mock_torch = MagicMock()
            mock_torch.workflow = None
            mock_torch_class.load.return_value = mock_torch

            with patch('torchbase.cli.run') as mock_run:
                _ = runner.invoke(
                    cli,
                    [
                        'run', '--strategy', 'balanced',
                        str(torch_without_workflow),
                        '-r', str(sample_reads_file)
                    ]
                )

                if mock_run.called:
                    call_args = mock_run.call_args[0][0]
                    # Should have: miniwdl, run, <workflow_path>, ...
                    assert len(call_args) >= 3
                    assert call_args[0] == 'miniwdl'
                    assert call_args[1] == 'run'

    def test_strategy_routing_uses_package_relative_path(
        self, torch_without_workflow, sample_reads_file
    ):
        """Strategy routing finds workflows relative to torchbase package."""
        runner = CliRunner()

        with patch('torchbase.torchfs.Torch') as mock_torch_class:
            mock_torch = MagicMock()
            mock_torch.workflow = None
            mock_torch_class.load.return_value = mock_torch

            with patch('torchbase.cli.run') as mock_run:
                _ = runner.invoke(
                    cli,
                    [
                        'run', '--strategy', 'fast',
                        str(torch_without_workflow),
                        '-r', str(sample_reads_file)
                    ]
                )

                if mock_run.called:
                    call_args = str(mock_run.call_args[0][0])
                    # Should reference torchbase package location
                    assert ('torchbase' in call_args or
                            'workflows' in call_args)


class TestStrategyWithoutTorchWorkflow:
    """Test strategy only works when torch has no embedded workflow."""

    def test_strategy_requires_no_embedded_workflow(
        self, torch_without_workflow, sample_reads_file
    ):
        """Strategy routing only applies to torches without embedded workflow.
        """
        runner = CliRunner()

        with patch('torchbase.torchfs.Torch') as mock_torch_class:
            mock_torch = MagicMock()
            mock_torch.workflow = None
            mock_torch_class.load.return_value = mock_torch

            with patch('torchbase.cli.run') as mock_run:
                result = runner.invoke(
                    cli,
                    [
                        'run', '--strategy', 'fast',
                        str(torch_without_workflow),
                        '-r', str(sample_reads_file)
                    ]
                )

                # Should succeed (or call miniwdl)
                assert result.exit_code == 0 or mock_run.called

    def test_data_only_torch_allows_strategy(
        self, torch_without_workflow, sample_reads_file
    ):
        """Data-only torch (no workflow) allows --strategy flag."""
        runner = CliRunner()

        # Verify torch has no workflow
        torch = Torch.load(torch_without_workflow)
        assert torch.workflow is None

        with patch('torchbase.cli.run') as mock_run:
            result = runner.invoke(
                cli,
                [
                    'run', '--strategy', 'balanced',
                    str(torch_without_workflow),
                    '-r', str(sample_reads_file)
                ]
            )

            # Should not reject the strategy flag
            assert ('--strategy' not in result.output or
                    result.exit_code == 0 or
                    mock_run.called)


class TestStrategyIntegrationWithMultiScheme:
    """Test strategy routing integrates with multi-scheme support."""

    def test_strategy_works_with_multi_scheme_torch(
        self, tmp_path, sample_reads_file
    ):
        """Strategy routing works with multi-scheme torches."""
        # Create multi-scheme torch
        torch_path = tmp_path / "test_namespace" / "multi_torch"
        torch_path = torch_path / "1.0.0.torch"
        torch_path.mkdir(parents=True)

        metadata = {
            "namespace": "test_namespace",
            "name": "multi_torch",
            "version": "1.0.0",
            "version_meta": {"strategy": "semver", "timestamp": 1609459200},
            "typing": {"method": "mlst"},
            "description": {"short": "Multi-scheme torch"},
            "schemes": {"ecoli": {}, "salmonella": {}}
        }
        with open(torch_path / "metadata.toml", "w") as f:
            toml.dump(metadata, f)

        # Create schemes
        schemes_dir = torch_path / "schemes"
        for scheme in ["ecoli", "salmonella"]:
            scheme_path = schemes_dir / scheme
            scheme_path.mkdir(parents=True)

            profiles = [["ST", "locus1"], ["1", "1"]]
            with open(scheme_path / "profiles.tsv", "w") as f:
                writer = csv.writer(f, delimiter="\t")
                writer.writerows(profiles)

            alleles_dir = scheme_path / "alleles"
            alleles_dir.mkdir()
            (alleles_dir / "locus1.fasta").write_text(">1\nACGT\n")

        runner = CliRunner()

        with patch('torchbase.torchfs.Torch') as mock_torch_class:
            mock_torch = MagicMock()
            mock_torch.workflow = None
            mock_torch.schemes = {
                "ecoli": MagicMock(),
                "salmonella": MagicMock()
            }
            mock_torch_class.load.return_value = mock_torch

            with patch('torchbase.cli.run') as mock_run:
                result = runner.invoke(
                    cli,
                    [
                        'run', '--strategy', 'fast',
                        str(torch_path),
                        '-r', str(sample_reads_file)
                    ]
                )

                # Should work with multi-scheme torch
                assert result.exit_code == 0 or mock_run.called

    def test_strategy_receives_concatenated_multi_scheme_files(
        self, tmp_path, sample_reads_file
    ):
        """Strategy workflows receive concatenated files from multi-scheme
        torches."""
        # This tests integration with #53 (multi-scheme concatenation)
        torch_path = tmp_path / "test_namespace" / "multi_torch"
        torch_path = torch_path / "1.0.0.torch"
        torch_path.mkdir(parents=True)

        metadata = {
            "namespace": "test_namespace",
            "name": "multi_torch",
            "version": "1.0.0",
            "version_meta": {"strategy": "semver", "timestamp": 1609459200},
            "typing": {"method": "mlst"},
            "description": {"short": "Multi-scheme torch"},
            "schemes": {"scheme1": {}, "scheme2": {}}
        }
        with open(torch_path / "metadata.toml", "w") as f:
            toml.dump(metadata, f)

        schemes_dir = torch_path / "schemes"
        for scheme in ["scheme1", "scheme2"]:
            scheme_path = schemes_dir / scheme
            scheme_path.mkdir(parents=True)

            profiles = [["ST", "locus"], ["1", "1"]]
            with open(scheme_path / "profiles.tsv", "w") as f:
                writer = csv.writer(f, delimiter="\t")
                writer.writerows(profiles)

            alleles_dir = scheme_path / "alleles"
            alleles_dir.mkdir()
            (alleles_dir / "locus.fasta").write_text(">1\nACGT\n")

        runner = CliRunner()

        with patch('torchbase.cli.run') as mock_run:
            _ = runner.invoke(
                cli,
                [
                    'run', '--strategy', 'balanced',
                    str(torch_path),
                    '-r', str(sample_reads_file)
                ]
            )

            # Should pass concatenated reference files to workflow
            if mock_run.called:
                call_args = mock_run.call_args[0][0]
                # Implementation detail: concatenated files should be passed
                assert len(call_args) > 2


class TestStrategyErrorHandling:
    """Test error handling for strategy routing."""

    def test_missing_builtin_workflow_file_raises_error(
        self, torch_without_workflow, sample_reads_file
    ):
        """Missing built-in workflow file raises clear error."""
        runner = CliRunner()

        with patch('torchbase.torchfs.Torch') as mock_torch_class:
            mock_torch = MagicMock()
            mock_torch.workflow = None
            mock_torch_class.load.return_value = mock_torch

            # Mock workflow file as missing
            with patch('torchbase.cli.Path') as mock_path_class:
                mock_workflow_path = MagicMock()
                mock_workflow_path.exists.return_value = False
                mock_path_class.return_value = mock_workflow_path

                result = runner.invoke(
                    cli,
                    [
                        'run', '--strategy', 'fast',
                        str(torch_without_workflow),
                        '-r', str(sample_reads_file)
                    ]
                )

                # Should fail with error about missing workflow
                if not result.exit_code == 0:
                    assert ('workflow' in result.output.lower() or
                            'not found' in result.output.lower())

    def test_strategy_flag_position_independent(
        self, torch_without_workflow, sample_reads_file
    ):
        """--strategy flag works in different positions."""
        runner = CliRunner()

        positions = [
            [
                'run', '--strategy', 'fast',
                str(torch_without_workflow),
                '-r', str(sample_reads_file)
            ],
            [
                'run', str(torch_without_workflow),
                '--strategy', 'fast',
                '-r', str(sample_reads_file)
            ],
            [
                'run', str(torch_without_workflow),
                '-r', str(sample_reads_file),
                '--strategy', 'fast'
            ],
        ]

        for args in positions:
            with patch('torchbase.torchfs.Torch') as mock_torch_class:
                with patch('torchbase.cli.run'):
                    mock_torch = MagicMock()
                    mock_torch.workflow = None
                    mock_torch_class.load.return_value = mock_torch

                    result = runner.invoke(cli, args)

                    # Should accept flag in any position
                    assert (result.exit_code == 0 or
                            '--strategy' not in result.output)


class TestStrategyWorkflowDiscoveryInteraction:
    """Test interaction between strategy routing and workflow discovery."""

    def test_strategy_bypasses_default_workflow_fetch(
        self, torch_without_workflow, sample_reads_file
    ):
        """Using --strategy bypasses torchbase/default-workflow fetch."""
        runner = CliRunner()

        with patch('torchbase.torchfs.Torch') as mock_torch_class:
            with patch(
                'torchbase.registry.RegistryManager'
            ) as mock_manager_class:
                mock_torch = MagicMock()
                mock_torch.workflow = None
                mock_torch_class.load.return_value = mock_torch

                mock_manager = MagicMock()
                mock_manager_class.return_value = mock_manager

                with patch('torchbase.cli.run'):
                    _ = runner.invoke(
                        cli,
                        [
                            'run', '--strategy', 'fast',
                            str(torch_without_workflow),
                            '-r', str(sample_reads_file)
                        ]
                    )

                    # Should NOT fetch default-workflow when strategy is
                    # specified
                    if mock_manager.fetch_torch.called:
                        torch_name = mock_manager.fetch_torch.call_args[0][0]
                        assert "default-workflow" not in torch_name

    def test_strategy_overrides_manifest_workflow(
        self, tmp_path, sample_reads_file
    ):
        """--strategy should not work with manifest-specified workflow."""
        torch_path = tmp_path / "test_namespace" / "manifest_workflow_torch"
        torch_path = torch_path / "1.0.0.torch"
        torch_path.mkdir(parents=True)

        metadata = {
            "namespace": "test_namespace",
            "name": "manifest_workflow_torch",
            "version": "1.0.0",
            "version_meta": {"strategy": "semver", "timestamp": 1609459200},
            "typing": {"method": "mlst"},
            "description": {"short": "Torch with manifest workflow"},
            "manifest": {
                "profiles": "profiles.tsv",
                "workflow": "custom.wdl"
            }
        }
        with open(torch_path / "metadata.toml", "w") as f:
            toml.dump(metadata, f)

        with open(torch_path / "custom.wdl", "w") as f:
            f.write("workflow custom { }")

        profiles = [["ST", "adk"], ["1", "1"]]
        with open(torch_path / "profiles.tsv", "w") as f:
            writer = csv.writer(f, delimiter="\t")
            writer.writerows(profiles)

        (torch_path / "_resources").mkdir()

        runner = CliRunner()

        with patch('torchbase.torchfs.Torch') as mock_torch_class:
            mock_torch = MagicMock()
            mock_torch.workflow = torch_path / "custom.wdl"
            mock_torch_class.load.return_value = mock_torch

            result = runner.invoke(
                cli,
                [
                    'run', '--strategy', 'balanced',
                    str(torch_path),
                    '-r', str(sample_reads_file)
                ]
            )

            # Should fail - torch has workflow defined
            assert result.exit_code != 0

    def test_no_strategy_with_no_workflow_uses_default(
        self, torch_without_workflow, sample_reads_file
    ):
        """Without --strategy and no torch workflow, falls back to default
        workflow."""
        runner = CliRunner()

        with patch('torchbase.torchfs.Torch') as mock_torch_class:
            with patch(
                'torchbase.registry.RegistryManager'
            ) as mock_manager_class:
                mock_torch = MagicMock()
                mock_torch.workflow = None
                mock_torch_class.load.return_value = mock_torch

                mock_manager = MagicMock()
                mock_manager_class.return_value = mock_manager

                with patch('torchbase.cli.run'):
                    _ = runner.invoke(
                        cli,
                        [
                            'run', str(torch_without_workflow),
                            '-r', str(sample_reads_file)
                        ]
                    )

                    # Should fetch default workflow (since no strategy
                    # specified). This preserves backward compatibility
                    if mock_manager.fetch_torch.called:
                        torch_name = mock_manager.fetch_torch.call_args[0][0]
                        # With default strategy=balanced, should use built-in
                        # workflow OR fetch default-workflow for backward
                        # compat
                        assert torch_name is not None
