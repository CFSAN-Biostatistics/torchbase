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


class TestSuspectFilteringReachesTheAlleleSet:
    """Quality filtering is applied before any compute is dispatched.

    This used to be a WDL task, and these tests used to assert that
    `filter_alleles.wdl` existed and was referenced by the typing workflows.
    Filtering is now a package-layer step (torchbase/quality_filters.py), so
    the contract worth defending is the observable one: with a quality.json and
    an exclusion flag, the excluded alleles are absent from the FASTA handed to
    the search, and the exclusions are reported in the result.
    """

    def _torch(self, tmp_path, alleles):
        from unittest.mock import MagicMock
        torch = MagicMock()
        allele_fasta = tmp_path / "alleles.fasta"
        allele_fasta.write_text(alleles)
        profiles = tmp_path / "profiles.tsv"
        profiles.write_text("ST\tadk\tfumC\n1\t1\t1\n")
        torch.get_unified_files.return_value = (allele_fasta, profiles)
        torch.path = tmp_path / "ns" / "name" / "1.0.0.torch"
        # New Torch fields (docs/adr/0003): a bare MagicMock auto-vivifies
        # any attribute access, so these must be pinned to their real
        # dataclass defaults or type_allelic misreads them as truthy.
        torch.calling_mode = "identity"
        torch.id_column = None
        return torch

    def _run_fast(self, torch, tmp_path, quality_path, **flags):
        """Run the fast strategy with compute stubbed, returning the FASTA it screened."""
        from unittest.mock import patch
        from torchbase import typing_run

        query = tmp_path / "query.fasta"
        query.write_text(">q1\nACGT\n")
        seen = {}

        def fake_screen(query_path, allele_fasta, strategy, threshold, engine):
            seen["alleles"] = Path(allele_fasta).read_text()
            return {}

        with patch.object(typing_run, "_screen", side_effect=fake_screen):
            result = typing_run.type_allelic(
                torch, str(query), strategy="fast",
                quality_json=str(quality_path) if quality_path else None,
                **flags
            )
        return seen.get("alleles", ""), result

    def test_suspect_locus_is_excluded_from_the_screened_alleles(self, tmp_path):
        # The schema quality_filters actually reads: a locus marked suspect.
        quality = tmp_path / "quality.json"
        quality.write_text(json.dumps({"loci": {"adk": {"suspect": True}}}))
        torch = self._torch(tmp_path, ">adk_1\nACGT\n>fumC_1\nACGT\n")

        screened, result = self._run_fast(
            torch, tmp_path, quality, exclude_suspect_loci=True
        )

        assert "fumC_1" in screened
        assert "adk_1" not in screened
        assert "adk" in json.dumps(result)

    def test_suspect_allele_is_excluded_from_the_screened_alleles(self, tmp_path):
        quality = tmp_path / "quality.json"
        quality.write_text(json.dumps(
            {"loci": {"adk": {"alleles": {"1": {"suspect": True}}}}}
        ))
        torch = self._torch(tmp_path, ">adk_1\nACGT\n>adk_2\nACGT\n")

        screened, _ = self._run_fast(
            torch, tmp_path, quality, exclude_suspect_alleles=True
        )

        assert "adk_1" not in screened
        assert "adk_2" in screened

    def test_nothing_is_excluded_without_a_flag(self, tmp_path):
        quality = tmp_path / "quality.json"
        quality.write_text(json.dumps({"loci": {"adk": {"suspect": True}}}))
        torch = self._torch(tmp_path, ">adk_1\nACGT\n>fumC_1\nACGT\n")

        screened, _ = self._run_fast(torch, tmp_path, quality)

        assert "adk_1" in screened and "fumC_1" in screened

    def test_no_quality_json_keeps_every_allele(self, tmp_path):
        torch = self._torch(tmp_path, ">adk_1\nACGT\n>fumC_1\nACGT\n")
        screened, _ = self._run_fast(
            torch, tmp_path, None, exclude_suspect_loci=True
        )
        assert "adk_1" in screened and "fumC_1" in screened

    def test_quality_report_schema_excludes_nothing(self, tmp_path):
        """Known defect, pinned: the producer and consumer schemas disagree.

        `torchbase/quality/report.py` writes `loci[x].suspect_pairs` plus a
        `summary.suspect_loci` list; the filtering code (inherited verbatim from
        the filter_alleles WDL task) reads `loci[x].suspect`,
        `loci[x].alleles[y].suspect`, `loci[x].similarities` and
        `profiles[p].suspect`. A real quality.json therefore excludes nothing,
        whatever flags are passed. Reconciling the two is a schema decision, not
        a refactor, so the current behaviour is asserted rather than changed.
        """
        quality = tmp_path / "quality.json"
        quality.write_text(json.dumps(_make_quality_json(
            suspect_loci=True, suspect_alleles=True
        )))
        torch = self._torch(tmp_path, ">adk_1\nACGT\n>adk_2\nACGT\n>fumC_1\nACGT\n")

        screened, _ = self._run_fast(
            torch, tmp_path, quality,
            exclude_suspect_alleles=True, exclude_suspect_loci=True,
            exclude_suspect_profiles=True,
        )

        assert "adk_1" in screened and "adk_2" in screened and "fumC_1" in screened


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
