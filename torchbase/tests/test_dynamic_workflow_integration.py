"""Integration tests for dynamic workflow parameter surfacing."""

import pytest
import toml
import csv
from pathlib import Path
from click.testing import CliRunner
from unittest.mock import patch, MagicMock

from torchbase.cli import cli


@pytest.fixture
def torch_with_params(tmp_path):
    """Create a torch with parameterized embedded workflow."""
    torch_path = tmp_path / "test_namespace" / "param_torch" / "1.0.0.torch"
    torch_path.mkdir(parents=True)

    metadata = {
        "namespace": "test_namespace",
        "name": "param_torch",
        "version": "1.0.0",
        "version_meta": {"strategy": "semver", "timestamp": 1609459200},
        "typing": {"method": "mlst"},
        "description": {"short": "Torch with parameterized workflow"},
        "manifest": {"profiles": "profiles.tsv"}
    }
    with open(torch_path / "metadata.toml", "w") as f:
        toml.dump(metadata, f)

    # Create a valid WDL with multiple parameter types
    wdl_content = """version 1.0

workflow parameterized_mlst {
    input {
        File query_sequences
        File allele_fasta
        File profiles_table
        Float confidence_threshold = 0.85
        Int min_depth = 3
        String input_type = "contigs"
        Boolean exclude_suspect_alleles = false
    }

    output {
        File result = "result.json"
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
    with open(torch_path / "_resources" / "adk.fasta", "w") as f:
        f.write(">adk_1\nACGT\n")

    return torch_path


@pytest.fixture
def sample_contigs_file(tmp_path):
    """Create a sample contigs file for testing."""
    contigs_file = tmp_path / "contigs.fasta"
    with open(contigs_file, "w") as f:
        f.write(">contig1\nACGTACGT\n")
    return contigs_file


class TestDynamicParameterParsing:
    """Test that workflow parameters are parsed from embedded WDL."""

    def test_workflow_params_auto_provisioned(
        self, torch_with_params, sample_contigs_file
    ):
        """Torch data files are auto-provisioned for embedded workflow."""
        runner = CliRunner()

        with patch('torchbase.cli.run') as mock_run:
            mock_run.return_value.returncode = 0

            result = runner.invoke(
                cli,
                [
                    'run', str(torch_with_params),
                    '-c', str(sample_contigs_file)
                ]
            )

            # Check that miniwdl was called
            assert mock_run.called

            # Get the command that was passed to run()
            call_args = mock_run.call_args[0][0]

            # Should include auto-provisioned parameters
            assert any('allele_fasta=' in arg for arg in call_args)
            assert any('profiles_table=' in arg for arg in call_args)
            assert any('query_sequences=' in arg for arg in call_args)

    def test_user_params_override_defaults(
        self, torch_with_params, sample_contigs_file
    ):
        """User-provided parameters override workflow defaults."""
        runner = CliRunner()

        with patch('torchbase.cli.run') as mock_run:
            mock_run.return_value.returncode = 0

            result = runner.invoke(
                cli,
                [
                    'run', str(torch_with_params),
                    '-c', str(sample_contigs_file),
                    'confidence_threshold=0.95',
                    'min_depth=10'
                ]
            )

            assert mock_run.called
            call_args = mock_run.call_args[0][0]

            # Should include user-provided parameters
            assert 'confidence_threshold=0.95' in call_args
            assert 'min_depth=10' in call_args

    def test_invalid_param_format_raises_error(
        self, torch_with_params, sample_contigs_file
    ):
        """Invalid parameter format raises helpful error."""
        runner = CliRunner()

        result = runner.invoke(
            cli,
            [
                'run', str(torch_with_params),
                '-c', str(sample_contigs_file),
                'invalid_no_equals'
            ]
        )

        assert result.exit_code != 0
        assert 'Invalid parameter format' in result.output

    def test_mixed_torch_args_and_flags(
        self, torch_with_params, sample_contigs_file
    ):
        """Can mix torch_args with existing CLI flags."""
        runner = CliRunner()

        with patch('torchbase.cli.run') as mock_run:
            mock_run.return_value.returncode = 0

            result = runner.invoke(
                cli,
                [
                    'run', str(torch_with_params),
                    '-c', str(sample_contigs_file),
                    '--exclude-suspect-loci',
                    'confidence_threshold=0.90'
                ]
            )

            assert mock_run.called
            call_args = mock_run.call_args[0][0]

            # Should include both
            assert 'confidence_threshold=0.90' in call_args
            assert any('exclude_suspect_loci=true' in arg for arg in call_args)


# Backward compatibility is tested by existing test_cli_strategy_routing.py tests
