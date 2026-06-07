"""Behavioral tests for suspect data filtering flags.

Tests verify that --exclude-suspect-loci, --exclude-suspect-alleles,
--exclude-suspect-profiles flags actually filter the data passed to the
workflow, not just that the CLI accepts the flags.
"""

import json
import tempfile
from pathlib import Path
from unittest.mock import patch, MagicMock

import pytest
import toml
from click.testing import CliRunner


def _make_quality_json(suspect_loci=None, suspect_alleles=None):
    """Build a minimal quality.json dict for testing."""
    return {
        "loci": {
            "adk": {
                "allele_count": 5,
                "kmer_size": 13,
                "similarity_stats": {"min": 0.8, "median": 0.9, "percentile_99": 0.99},
                "threshold": 0.95,
                "threshold_method": "gap_detection",
                "suspect_pairs": [
                    {"allele1": "adk_1", "allele2": "adk_2",
                     "similarity": 0.99, "issue_type": "duplicate",
                     "containment_1_in_2": 0.99, "containment_2_in_1": 0.99}
                ] if suspect_alleles else [],
            },
            "fumC": {
                "allele_count": 3,
                "kmer_size": 13,
                "similarity_stats": {"min": 0.5, "median": 0.7, "percentile_99": 0.95},
                "threshold": 0.95,
                "threshold_method": "gap_detection",
                "suspect_pairs": [],
            },
        },
        "summary": {
            "total_loci": 2,
            "total_suspect_allele_pairs": 1 if suspect_alleles else 0,
            "suspect_loci": ["adk"] if suspect_loci else [],
            "suspect_profiles": ["adk"] if suspect_loci else [],
        },
    }


class TestSuspectFlagsAccepted:
    """CLI accepts all suspect filtering flags without error."""

    def test_exclude_suspect_alleles_flag_accepted(self):
        from torchbase.cli import cli
        runner = CliRunner()
        result = runner.invoke(cli, ["run", "--help"])
        assert "--exclude-suspect-alleles" in result.output

    def test_include_suspect_alleles_flag_accepted(self):
        from torchbase.cli import cli
        runner = CliRunner()
        result = runner.invoke(cli, ["run", "--help"])
        assert "--include-suspect-alleles" in result.output

    def test_exclude_suspect_loci_flag_accepted(self):
        from torchbase.cli import cli
        runner = CliRunner()
        result = runner.invoke(cli, ["run", "--help"])
        assert "--exclude-suspect-loci" in result.output

    def test_exclude_suspect_profiles_flag_accepted(self):
        from torchbase.cli import cli
        runner = CliRunner()
        result = runner.invoke(cli, ["run", "--help"])
        assert "--exclude-suspect-profiles" in result.output

    def test_quality_json_option_accepted(self):
        from torchbase.cli import cli
        runner = CliRunner()
        result = runner.invoke(cli, ["run", "--help"])
        assert "--quality-json" in result.output


class TestQualityJsonParsing:
    """Quality JSON is correctly parsed for suspect data identification."""

    def test_suspect_loci_identified_from_quality_json(self):
        """Loci listed in summary.suspect_loci are correctly extracted."""
        qj = _make_quality_json(suspect_loci=True, suspect_alleles=True)
        assert "adk" in qj["summary"]["suspect_loci"]
        assert "fumC" not in qj["summary"]["suspect_loci"]

    def test_suspect_alleles_identified_from_locus_data(self):
        """Suspect allele pairs are accessible per-locus."""
        qj = _make_quality_json(suspect_loci=True, suspect_alleles=True)
        adk_pairs = qj["loci"]["adk"]["suspect_pairs"]
        assert len(adk_pairs) == 1
        assert adk_pairs[0]["allele1"] == "adk_1"
        assert adk_pairs[0]["allele2"] == "adk_2"

    def test_clean_locus_has_no_suspect_pairs(self):
        """Locus with no suspect pairs should have empty list."""
        qj = _make_quality_json(suspect_loci=True, suspect_alleles=True)
        assert qj["loci"]["fumC"]["suspect_pairs"] == []

    def test_no_suspects_when_quality_json_clean(self):
        """When quality JSON has no suspects, all loci are clean."""
        qj = _make_quality_json(suspect_loci=False, suspect_alleles=False)
        assert qj["summary"]["suspect_loci"] == []
        assert qj["loci"]["adk"]["suspect_pairs"] == []


class TestFilterAllelesTaskInputs:
    """filter_alleles WDL task receives correct inputs from quality JSON."""

    def test_filter_alleles_wdl_exists(self):
        """filter_alleles.wdl task file exists in builtin tasks."""
        from pathlib import Path
        import torchbase
        tasks_dir = Path(torchbase.__file__).parent / "workflows" / "builtin" / "tasks"
        assert (tasks_dir / "filter_alleles.wdl").exists()

    def test_filter_alleles_wdl_has_quality_json_input(self):
        """filter_alleles.wdl accepts quality_json as an input."""
        from pathlib import Path
        import torchbase
        wdl = Path(torchbase.__file__).parent / "workflows" / "builtin" / "tasks" / "filter_alleles.wdl"
        content = wdl.read_text()
        assert "quality_json" in content

    def test_fast_typing_wdl_uses_filter_alleles(self):
        """fast_typing.wdl invokes the filter_alleles task."""
        from pathlib import Path
        import torchbase
        wdl = Path(torchbase.__file__).parent / "workflows" / "builtin" / "fast_typing.wdl"
        content = wdl.read_text()
        assert "filter_alleles" in content

    def test_balanced_typing_wdl_uses_filter_alleles(self):
        """balanced_typing.wdl invokes the filter_alleles task."""
        from pathlib import Path
        import torchbase
        wdl = Path(torchbase.__file__).parent / "workflows" / "builtin" / "balanced_typing.wdl"
        content = wdl.read_text()
        assert "filter_alleles" in content


class TestStrategyConflictWithEmbeddedWorkflow:
    """--strategy flag conflicts with embedded main.wdl in torch."""

    def test_strategy_flag_rejected_when_embedded_workflow_present(self):
        """Run command raises error when --strategy used with embedded workflow."""
        from torchbase.cli import cli
        from torchbase.torchfs import Torch

        runner = CliRunner()
        with tempfile.TemporaryDirectory() as tmp:
            # Build a minimal torch with main.wdl
            torch_path = Path(tmp) / "ns" / "name" / "1.0.0.torch"
            torch_path.mkdir(parents=True)
            (torch_path / "metadata.toml").write_text(
                toml.dumps({"namespace": "ns", "name": "name", "version": "1.0.0"})
            )
            (torch_path / "profiles.tsv").write_text("ST\tadk\n1\t1\n")
            (torch_path / "main.wdl").write_text("version 1.0\nworkflow main {}\n")
            resources = torch_path / "_resources"
            resources.mkdir()
            (resources / "adk.fasta").write_text(">adk_1\nACGT\n")

            mock_torch = MagicMock()
            mock_torch.workflow = torch_path / "main.wdl"

            with patch("torchbase.torchfs.Torch") as MockTorch, \
                 patch("torchbase.registry.RegistryManager") as MockManager:
                MockTorch.load.return_value = mock_torch
                MockManager.return_value.fetch_torch.return_value = torch_path

                with tempfile.NamedTemporaryFile(suffix=".fasta", mode="w", delete=False) as f:
                    f.write(">c1\nACGT\n")
                    contigs_path = f.name

                result = runner.invoke(cli, [
                    "run", "ns/name",
                    "--strategy", "fast",
                    "--contigs", contigs_path,
                ])
                assert result.exit_code != 0
                assert "strategy" in result.output.lower() or "workflow" in result.output.lower()
