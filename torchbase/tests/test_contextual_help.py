"""Tests for context-aware help in torchbase run command.

RED TESTS - NOT YET IMPLEMENTED
================================
These tests define the expected behavior for contextual help, where
`torchbase run my_torch --help` shows torch-specific parameters
instead of generic help text.

Current status: Tests fail (RED phase of TDD)
See GitHub issue for implementation tracking.

Expected behavior:
- When user runs `torchbase run my_torch --help`, inspect the torch's
  embedded workflow and show available parameters in CLI format
- Show parameter types, defaults, and whether required/optional
- Show which parameters are auto-provisioned by the torch
- Include a tip to run `torchbase workflow inspect` for full diagram
- For data-only torches (no embedded workflow), show standard help
"""

import pytest
import toml
import csv
from pathlib import Path
from click.testing import CliRunner

from torchbase.cli import cli

# Skip all tests in this file until feature is implemented
pytestmark = pytest.mark.skip(reason="RED tests for future feature - not yet implemented")


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

    # Create a valid WDL with various parameter types
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
        File? quality_json
    }

    output {
        File result = "result.json"
    }
}
"""
    with open(torch_path / "main.wdl", "w") as f:
        f.write(wdl_content)

    profiles = [["ST", "adk"], ["1", "1"]]
    with open(torch_path / "profiles.tsv", "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerows(profiles)

    (torch_path / "_resources").mkdir()
    with open(torch_path / "_resources" / "adk.fasta", "w") as f:
        f.write(">adk_1\nACGT\n")

    return torch_path


@pytest.fixture
def torch_without_workflow(tmp_path):
    """Create a torch without embedded workflow."""
    torch_path = tmp_path / "test_namespace" / "data_torch" / "1.0.0.torch"
    torch_path.mkdir(parents=True)

    metadata = {
        "namespace": "test_namespace",
        "name": "data_torch",
        "version": "1.0.0",
        "version_meta": {"strategy": "semver", "timestamp": 1609459200},
        "typing": {"method": "mlst"},
        "description": {"short": "Data-only torch"},
        "manifest": {"profiles": "profiles.tsv"}
    }
    with open(torch_path / "metadata.toml", "w") as f:
        toml.dump(metadata, f)

    profiles = [["ST", "adk"], ["1", "1"]]
    with open(torch_path / "profiles.tsv", "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerows(profiles)

    (torch_path / "_resources").mkdir()
    with open(torch_path / "_resources" / "adk.fasta", "w") as f:
        f.write(">adk_1\nACGT\n")

    return torch_path


class TestContextualHelp:
    """Test that --help shows torch-specific parameters."""

    def test_help_with_embedded_workflow_shows_parameters(self, torch_with_params):
        """torchbase run my_torch --help shows torch-specific parameters."""
        runner = CliRunner()
        result = runner.invoke(cli, ['run', str(torch_with_params), '--help'])

        assert result.exit_code == 0

        # Should show torch name/path
        assert 'param_torch' in result.output or str(torch_with_params) in result.output

        # Should show workflow parameters
        assert 'confidence_threshold' in result.output
        assert 'min_depth' in result.output
        assert 'input_type' in result.output

        # Should show defaults
        assert '0.85' in result.output  # default for confidence_threshold
        assert '3' in result.output     # default for min_depth
        assert 'contigs' in result.output  # default for input_type

    def test_help_shows_parameter_types(self, torch_with_params):
        """Help output shows parameter types (Float, Int, String, Boolean)."""
        runner = CliRunner()
        result = runner.invoke(cli, ['run', str(torch_with_params), '--help'])

        assert result.exit_code == 0

        # Should indicate parameter types
        # Could be shown as "FLOAT", "INT", "STRING", etc.
        output_lower = result.output.lower()
        assert any(word in output_lower for word in ['float', 'int', 'string', 'boolean'])

    def test_help_distinguishes_required_vs_optional(self, torch_with_params):
        """Help output distinguishes required from optional parameters."""
        runner = CliRunner()
        result = runner.invoke(cli, ['run', str(torch_with_params), '--help'])

        assert result.exit_code == 0

        # Should have sections or markers for required vs optional
        output_lower = result.output.lower()
        assert 'required' in output_lower or 'optional' in output_lower

    def test_help_shows_standard_cli_options(self, torch_with_params):
        """Help output still shows standard CLI options like -c/--contigs."""
        runner = CliRunner()
        result = runner.invoke(cli, ['run', str(torch_with_params), '--help'])

        assert result.exit_code == 0

        # Should show standard input options
        assert '--contigs' in result.output or '-c' in result.output
        assert '--reads' in result.output or '-r' in result.output

    def test_help_mentions_workflow_inspect(self, torch_with_params):
        """Help output mentions how to get full workflow diagram."""
        runner = CliRunner()
        result = runner.invoke(cli, ['run', str(torch_with_params), '--help'])

        assert result.exit_code == 0

        # Should mention workflow inspect command for detailed view
        assert 'workflow inspect' in result.output

    def test_help_without_torch_shows_generic_help(self):
        """torchbase run --help (no torch) shows generic help."""
        runner = CliRunner()
        result = runner.invoke(cli, ['run', '--help'])

        assert result.exit_code == 0

        # Should show generic usage without torch-specific parameters
        assert 'Usage:' in result.output
        assert 'TORCH' in result.output  # Shows TORCH as an argument

        # Should NOT show workflow-specific parameters
        assert 'confidence_threshold' not in result.output

    def test_help_with_data_only_torch_shows_strategy_info(self, torch_without_workflow):
        """Help for data-only torch shows --strategy options."""
        runner = CliRunner()
        result = runner.invoke(cli, ['run', str(torch_without_workflow), '--help'])

        assert result.exit_code == 0

        # Should mention strategy flag since torch has no embedded workflow
        assert '--strategy' in result.output
        assert any(s in result.output for s in ['fast', 'balanced', 'sensitive'])

    def test_help_format_is_cli_style_not_diagram(self, torch_with_params):
        """Help output is traditional CLI style, not ASCII box diagram."""
        runner = CliRunner()
        result = runner.invoke(cli, ['run', str(torch_with_params), '--help'])

        assert result.exit_code == 0

        # Should NOT be the ASCII box diagram from workflow inspect
        assert '┌' not in result.output  # Box drawing characters
        assert '│' not in result.output
        assert '└' not in result.output

        # Should use traditional CLI help format
        assert 'Usage:' in result.output

    def test_help_shows_parameter_format(self, torch_with_params):
        """Help shows that parameters use key=value format."""
        runner = CliRunner()
        result = runner.invoke(cli, ['run', str(torch_with_params), '--help'])

        assert result.exit_code == 0

        # Should indicate key=value format for parameters
        assert '=' in result.output

    def test_help_explains_auto_provisioned_parameters(self, torch_with_params):
        """Help explains which parameters are auto-provisioned by torch."""
        runner = CliRunner()
        result = runner.invoke(cli, ['run', str(torch_with_params), '--help'])

        assert result.exit_code == 0

        # Should indicate that allele_fasta and profiles_table are provided automatically
        # This could be shown as notes or different formatting
        output_lower = result.output.lower()

        # Either explicitly mentions auto-provision or groups parameters differently
        assert ('auto' in output_lower or
                'provided' in output_lower or
                'allele_fasta' in output_lower)


class TestHelpEdgeCases:
    """Test edge cases for contextual help."""

    def test_help_with_invalid_torch_path(self, tmp_path):
        """Help with non-existent torch path shows error."""
        runner = CliRunner()
        nonexistent = tmp_path / "nonexistent" / "torch"
        result = runner.invoke(cli, ['run', str(nonexistent), '--help'])

        # Could either show error or fall back to generic help
        # Implementation decides which behavior is better
        assert result.exit_code != 0 or 'Usage:' in result.output

    def test_help_position_independent(self, torch_with_params):
        """--help works regardless of position in arguments."""
        runner = CliRunner()

        # help before torch
        result1 = runner.invoke(cli, ['run', '--help', str(torch_with_params)])

        # help after torch
        result2 = runner.invoke(cli, ['run', str(torch_with_params), '--help'])

        # Both should show help (implementation may choose which style)
        assert result1.exit_code == 0 or result2.exit_code == 0
        assert 'Usage:' in result1.output or 'Usage:' in result2.output

    def test_help_with_workflow_parse_error(self, tmp_path):
        """Help gracefully handles workflows that fail to parse."""
        torch_path = tmp_path / "test" / "broken" / "1.0.0.torch"
        torch_path.mkdir(parents=True)

        metadata = {
            "namespace": "test",
            "name": "broken",
            "version": "1.0.0",
            "manifest": {"profiles": "profiles.tsv"}
        }
        with open(torch_path / "metadata.toml", "w") as f:
            toml.dump(metadata, f)

        # Invalid WDL (no version declaration)
        with open(torch_path / "main.wdl", "w") as f:
            f.write("workflow broken { input { File x } }")

        profiles = [["ST", "adk"], ["1", "1"]]
        with open(torch_path / "profiles.tsv", "w", newline="") as f:
            writer = csv.writer(f, delimiter="\t")
            writer.writerows(profiles)

        (torch_path / "_resources").mkdir()

        runner = CliRunner()
        result = runner.invoke(cli, ['run', str(torch_path), '--help'])

        # Should show help even if workflow parsing fails
        assert result.exit_code == 0
        assert 'Usage:' in result.output
