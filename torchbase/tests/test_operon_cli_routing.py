"""Tests for CLI routing of operon torches (docs/operon-strategy-plan.md §2, §6).

- `--strategy` is rejected for operon torches, mirroring the embedded-workflow guard.
- Without `--strategy`, an operon torch routes to the built-in operon_typing.wdl.
"""

import csv
from unittest.mock import patch

import pytest
import toml
from click.testing import CliRunner

from torchbase.cli import cli


OPERON_CFG = {
    "subunit_order": ["A", "B"],
    "intergenic_max": 36,
    "reference": {"file": "_resources/subunits.faa"},
    "identity_thresholds": {"default": 0.98},
}


@pytest.fixture
def operon_torch(tmp_path):
    torch_path = tmp_path / "test_ns" / "stx" / "1.0.0.torch"
    torch_path.mkdir(parents=True)

    metadata = {
        "namespace": "test_ns",
        "name": "stx",
        "version": "1.0.0",
        "version_meta": {"strategy": "semver", "timestamp": 1609459200},
        "manifest": {"profiles": "profiles.tsv"},
        "description": {"short": "stx operon torch"},
        "typing_model": "operon",
        "operon": OPERON_CFG,
    }
    with open(torch_path / "metadata.toml", "w") as f:
        toml.dump(metadata, f)

    with open(torch_path / "profiles.tsv", "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["subtype", "class"], delimiter="\t")
        writer.writeheader()
        writer.writerow({"subtype": "stx1a", "class": "1a"})

    (torch_path / "_resources").mkdir()
    (torch_path / "_resources" / "subunits.faa").write_text(">AAS07596.1|stxA2|stxA2\nMKC\n")

    return torch_path


@pytest.fixture
def contigs_file(tmp_path):
    f = tmp_path / "contigs.fasta"
    f.write_text(">NODE_1\nACGT\n")
    return f


class TestOperonStrategyRejection:
    def test_explicit_strategy_rejected(self, operon_torch, contigs_file):
        runner = CliRunner()
        result = runner.invoke(
            cli,
            ["run", str(operon_torch), "--strategy", "fast", "-c", str(contigs_file)],
        )
        assert result.exit_code != 0
        assert "operon" in result.output.lower()
        assert "--strategy" in result.output


class TestOperonWorkflowRouting:
    def test_routes_to_builtin_operon_workflow(self, operon_torch, contigs_file):
        runner = CliRunner()
        with patch("torchbase.cli.run") as mock_run:
            mock_run.return_value.returncode = 0
            result = runner.invoke(
                cli, ["run", str(operon_torch), "-c", str(contigs_file)]
            )
        assert result.exit_code == 0, result.output
        mock_run.assert_called_once()
        miniwdl_cmd = mock_run.call_args[0][0]
        assert any("operon_typing.wdl" in arg for arg in miniwdl_cmd)
        assert any(arg.startswith("subunit_reference=") for arg in miniwdl_cmd)
        assert any(arg.startswith("profiles_table=") for arg in miniwdl_cmd)
        assert any(arg.startswith("operon_config_json=") for arg in miniwdl_cmd)
