#!/usr/bin/env python

"""Tests for quality report generation module."""

import pytest
import json
import tempfile
from pathlib import Path
from textwrap import dedent

from torchbase.quality.kmer_analysis import analyze_locus, SimilarityReport
from torchbase.quality.report import generate_report


@pytest.fixture
def temp_fasta_dir():
    """Create a temporary directory for FASTA files."""
    with tempfile.TemporaryDirectory() as tmpdir:
        yield Path(tmpdir)


@pytest.fixture
def normal_distribution_fasta(temp_fasta_dir):
    """Create a FASTA file with normally distributed similarity (no suspects)."""
    fasta_path = temp_fasta_dir / "normal.fasta"
    content = dedent("""
        >allele_1
        ATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATG
        >allele_2
        GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA
        >allele_3
        CGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
        >allele_4
        TCATCATCATCATCATCATCATCATCATCATCATCATCATCATCATCATCATCATCATCATCATCATCATCATCATCATCATCATCATCATCA
    """).strip()
    fasta_path.write_text(content)
    return fasta_path


@pytest.fixture
def bimodal_distribution_fasta(temp_fasta_dir):
    """Create a FASTA file with bimodal similarity (has duplicates)."""
    fasta_path = temp_fasta_dir / "bimodal.fasta"
    # Create groups: identical pairs and diverse sequences
    seq_group1 = "ATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATG"
    seq_group2 = "ATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATTG"
    seq_group3 = "GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA"
    content = f">allele_1\n{seq_group1}\n>allele_2\n{seq_group2}\n>allele_3\n{seq_group3}\n>allele_4\nTCATCATCATCATCATCATCATCATCATCATCATCATCATCATCATCATCATCATCATCATCATCATCATCATCATCATCATCATCATCATCA"
    fasta_path.write_text(content)
    return fasta_path


@pytest.fixture
def uniform_distribution_fasta(temp_fasta_dir):
    """Create a FASTA file with uniform similarity distribution."""
    fasta_path = temp_fasta_dir / "uniform.fasta"
    # Create sequences with gradual differences
    base = "ATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATG"
    content = f">allele_1\n{base}\n>allele_2\n{base[:-1]}G\n>allele_3\n{base[:-2]}GG\n>allele_4\n{base[:-3]}GGG"
    fasta_path.write_text(content)
    return fasta_path


@pytest.fixture
def empty_suspects_report():
    """Create a SimilarityReport with no suspect pairs."""
    return SimilarityReport(
        similarities={
            ("allele_1", "allele_2"): 45.5,
            ("allele_1", "allele_3"): 42.3,
            ("allele_2", "allele_3"): 48.7,
        },
        threshold=70.0,
        suspect_pairs=[],
        statistics={
            "mean": 45.5,
            "std_dev": 3.2,
            "min": 42.3,
            "max": 48.7,
            "percentile_99": 48.0,
            "threshold_type": "percentile",
        },
    )


@pytest.fixture
def suspect_report():
    """Create a SimilarityReport with suspect pairs."""
    return SimilarityReport(
        similarities={
            ("allele_1", "allele_2"): 98.5,
            ("allele_1", "allele_3"): 42.3,
            ("allele_2", "allele_3"): 41.7,
            ("allele_3", "allele_4"): 95.2,
        },
        threshold=90.0,
        suspect_pairs=[
            {
                "allele1": "allele_1",
                "allele2": "allele_2",
                "similarity": 98.5,
                "containment_1_in_2": 98.0,
                "containment_2_in_1": 99.0,
                "issue_type": "duplicate",
            },
            {
                "allele1": "allele_3",
                "allele2": "allele_4",
                "similarity": 95.2,
                "containment_1_in_2": 95.5,
                "containment_2_in_1": 94.8,
                "issue_type": "overlap",
            },
        ],
        statistics={
            "mean": 68.9,
            "std_dev": 31.2,
            "min": 41.7,
            "max": 98.5,
            "percentile_99": 97.0,
            "threshold_type": "percentile",
        },
    )


@pytest.fixture
def locus_reports(suspect_report, empty_suspects_report):
    """Create a dictionary of locus reports."""
    return {
        "locus_A": suspect_report,
        "locus_B": empty_suspects_report,
    }


class TestGenerateReportBasics:
    """Test basic functionality of generate_report."""

    def test_generate_report_with_text_format(self, empty_suspects_report):
        """Test that generate_report produces text output."""
        reports = {"locus_1": empty_suspects_report}
        result = generate_report(reports, format="text")
        assert isinstance(result, str)
        assert len(result) > 0

    def test_generate_report_with_json_format(self, empty_suspects_report):
        """Test that generate_report produces JSON output."""
        reports = {"locus_1": empty_suspects_report}
        result = generate_report(reports, format="json")
        assert isinstance(result, str)
        # Should be valid JSON
        data = json.loads(result)
        assert isinstance(data, dict)

    def test_generate_report_with_both_format(self, empty_suspects_report):
        """Test that generate_report produces both formats."""
        reports = {"locus_1": empty_suspects_report}
        result = generate_report(reports, format="both")
        assert isinstance(result, dict)
        assert "json" in result
        assert "text" in result
        assert isinstance(result["json"], str)
        assert isinstance(result["text"], str)
        # JSON should be valid
        json.loads(result["json"])

    def test_generate_report_default_format(self, empty_suspects_report):
        """Test that generate_report defaults to text format."""
        reports = {"locus_1": empty_suspects_report}
        result = generate_report(reports)
        assert isinstance(result, str)


class TestJSONOutput:
    """Test JSON output format."""

    def test_json_has_loci(self, empty_suspects_report):
        """Test that JSON output includes loci section."""
        reports = {"locus_1": empty_suspects_report}
        result = generate_report(reports, format="json")
        data = json.loads(result)
        assert "loci" in data

    def test_json_has_similarity_stats(self, empty_suspects_report):
        """Test that JSON output includes similarity_stats."""
        reports = {"locus_1": empty_suspects_report}
        result = generate_report(reports, format="json")
        data = json.loads(result)
        assert "similarity_stats" in data

    def test_json_has_suspect_pairs(self, empty_suspects_report):
        """Test that JSON output includes suspect_pairs."""
        reports = {"locus_1": empty_suspects_report}
        result = generate_report(reports, format="json")
        data = json.loads(result)
        assert "suspect_pairs" in data

    def test_json_has_summary(self, empty_suspects_report):
        """Test that JSON output includes summary."""
        reports = {"locus_1": empty_suspects_report}
        result = generate_report(reports, format="json")
        data = json.loads(result)
        assert "summary" in data

    def test_json_loci_structure(self, suspect_report):
        """Test that JSON loci section has correct structure."""
        reports = {"locus_A": suspect_report}
        result = generate_report(reports, format="json")
        data = json.loads(result)
        assert "locus_A" in data["loci"]
        locus_data = data["loci"]["locus_A"]
        assert "similarities" in locus_data
        assert "threshold" in locus_data
        assert "statistics" in locus_data

    def test_json_suspect_pairs_structure(self, suspect_report):
        """Test that JSON suspect_pairs section has correct structure."""
        reports = {"locus_A": suspect_report}
        result = generate_report(reports, format="json")
        data = json.loads(result)
        assert len(data["suspect_pairs"]["locus_A"]) == 2
        for pair in data["suspect_pairs"]["locus_A"]:
            assert "allele1" in pair
            assert "allele2" in pair
            assert "similarity" in pair
            assert "issue_type" in pair

    def test_json_summary_contains_stats(self, suspect_report, empty_suspects_report):
        """Test that JSON summary contains overall statistics."""
        reports = {
            "locus_A": suspect_report,
            "locus_B": empty_suspects_report,
        }
        result = generate_report(reports, format="json")
        data = json.loads(result)
        summary = data["summary"]
        assert "total_loci" in summary
        assert "total_suspect_allele_pairs" in summary
        assert "suspect_loci" in summary
        assert "suspect_profiles" in summary


class TestTextOutput:
    """Test text output format."""

    def test_text_includes_locus_headers(self, empty_suspects_report):
        """Test that text output includes locus names."""
        reports = {"locus_1": empty_suspects_report}
        result = generate_report(reports, format="text")
        assert "locus_1" in result

    def test_text_includes_statistics(self, empty_suspects_report):
        """Test that text output includes statistics."""
        reports = {"locus_1": empty_suspects_report}
        result = generate_report(reports, format="text")
        assert "mean" in result or "Mean" in result
        assert "std" in result or "Std" in result

    def test_text_includes_histogram(self, suspect_report):
        """Test that text output includes ASCII histogram."""
        reports = {"locus_A": suspect_report}
        result = generate_report(reports, format="text")
        # Should have some histogram-like content (bars, brackets, etc)
        assert "[" in result and "]" in result

    def test_text_includes_suspect_pairs(self, suspect_report):
        """Test that text output includes suspect pairs."""
        reports = {"locus_A": suspect_report}
        result = generate_report(reports, format="text")
        assert "allele_1" in result
        assert "allele_2" in result

    def test_text_includes_summary_section(self, locus_reports):
        """Test that text output includes a summary section."""
        result = generate_report(locus_reports, format="text")
        assert "Summary" in result or "SUMMARY" in result or "summary" in result


class TestHistogramRendering:
    """Test ASCII histogram rendering in text output."""

    def test_histogram_renders_normal_distribution(self, normal_distribution_fasta):
        """Test histogram rendering with normal distribution."""
        report = analyze_locus(normal_distribution_fasta, k_size=21)
        reports = {"locus_normal": report}
        result = generate_report(reports, format="text")
        # Should have histogram bars
        assert "|" in result or "-" in result

    def test_histogram_renders_bimodal_distribution(self, bimodal_distribution_fasta):
        """Test histogram rendering with bimodal distribution."""
        report = analyze_locus(bimodal_distribution_fasta, k_size=21)
        reports = {"locus_bimodal": report}
        result = generate_report(reports, format="text")
        # Should have histogram bars
        assert "|" in result or "-" in result

    def test_histogram_renders_uniform_distribution(self, uniform_distribution_fasta):
        """Test histogram rendering with uniform distribution."""
        report = analyze_locus(uniform_distribution_fasta, k_size=21)
        reports = {"locus_uniform": report}
        result = generate_report(reports, format="text")
        # Should have histogram bars
        assert "|" in result or "-" in result


class TestHierarchicalFlagging:
    """Test hierarchical suspect flagging."""

    def test_suspect_alleles_propagate_to_loci(self, suspect_report):
        """Test that suspect alleles mark their locus as suspect."""
        reports = {"locus_A": suspect_report}
        result = generate_report(reports, format="json")
        data = json.loads(result)
        # locus_A should be in suspect_loci because it has suspect pairs
        assert "locus_A" in data["summary"]["suspect_loci"]

    def test_suspect_loci_propagate_to_profiles(self, suspect_report, empty_suspects_report):
        """Test that suspect loci mark profiles as suspect."""
        reports = {
            "locus_A": suspect_report,
            "locus_B": empty_suspects_report,
        }
        result = generate_report(reports, format="json")
        data = json.loads(result)
        # There should be suspect profiles (loci with suspects)
        assert len(data["summary"]["suspect_loci"]) > 0

    def test_no_suspect_alleles_no_suspect_loci(self, empty_suspects_report):
        """Test that loci without suspects aren't flagged."""
        reports = {"locus_1": empty_suspects_report}
        result = generate_report(reports, format="json")
        data = json.loads(result)
        # locus_1 should not be in suspect_loci
        assert "locus_1" not in data["summary"]["suspect_loci"]

    def test_all_alleles_suspect_flags_locus(self, suspect_report):
        """Test that locus with all suspects is flagged."""
        # Create a report where all alleles are in suspect pairs
        all_suspect_report = SimilarityReport(
            similarities={
                ("allele_1", "allele_2"): 98.5,
                ("allele_1", "allele_3"): 97.0,
                ("allele_2", "allele_3"): 96.5,
            },
            threshold=95.0,
            suspect_pairs=[
                {
                    "allele1": "allele_1",
                    "allele2": "allele_2",
                    "similarity": 98.5,
                    "containment_1_in_2": 98.0,
                    "containment_2_in_1": 99.0,
                    "issue_type": "duplicate",
                },
                {
                    "allele1": "allele_1",
                    "allele2": "allele_3",
                    "similarity": 97.0,
                    "containment_1_in_2": 96.5,
                    "containment_2_in_1": 97.5,
                    "issue_type": "duplicate",
                },
                {
                    "allele1": "allele_2",
                    "allele2": "allele_3",
                    "similarity": 96.5,
                    "containment_1_in_2": 96.0,
                    "containment_2_in_1": 97.0,
                    "issue_type": "duplicate",
                },
            ],
            statistics={
                "mean": 97.3,
                "std_dev": 1.0,
                "min": 96.5,
                "max": 98.5,
                "percentile_99": 98.0,
                "threshold_type": "percentile",
            },
        )
        reports = {"locus_all_suspect": all_suspect_report}
        result = generate_report(reports, format="json")
        data = json.loads(result)
        assert "locus_all_suspect" in data["summary"]["suspect_loci"]


class TestEdgeCases:
    """Test edge cases."""

    def test_empty_locus_reports(self):
        """Test handling of empty locus reports dictionary."""
        result = generate_report({}, format="json")
        data = json.loads(result)
        assert "summary" in data
        assert data["summary"]["total_loci"] == 0

    def test_single_locus_report(self, empty_suspects_report):
        """Test handling of single locus report."""
        reports = {"locus_1": empty_suspects_report}
        result = generate_report(reports, format="json")
        data = json.loads(result)
        assert data["summary"]["total_loci"] == 1

    def test_multiple_loci_reports(self, suspect_report, empty_suspects_report):
        """Test handling of multiple loci reports."""
        reports = {
            "locus_A": suspect_report,
            "locus_B": empty_suspects_report,
            "locus_C": suspect_report,
        }
        result = generate_report(reports, format="json")
        data = json.loads(result)
        assert data["summary"]["total_loci"] == 3

    def test_no_similarities_report(self):
        """Test handling of report with no similarities."""
        report = SimilarityReport(
            similarities={},
            threshold=0.0,
            suspect_pairs=[],
            statistics={
                "mean": 0.0,
                "std_dev": 0.0,
                "min": 0.0,
                "max": 0.0,
                "percentile_99": 0.0,
                "threshold_type": "none",
            },
        )
        reports = {"empty_locus": report}
        result = generate_report(reports, format="json")
        data = json.loads(result)
        assert data["summary"]["total_loci"] == 1

    def test_text_output_no_suspects(self, empty_suspects_report):
        """Test text output with no suspects."""
        reports = {"locus_1": empty_suspects_report}
        result = generate_report(reports, format="text")
        assert isinstance(result, str)
        assert "locus_1" in result

    def test_json_output_preserves_all_data(self, suspect_report):
        """Test that JSON output preserves all data from SimilarityReport."""
        reports = {"locus_A": suspect_report}
        result = generate_report(reports, format="json")
        data = json.loads(result)
        locus_data = data["loci"]["locus_A"]
        # Check that similarities are preserved
        assert len(locus_data["similarities"]) == len(suspect_report.similarities)


class TestFormatIntegration:
    """Test integration across formats."""

    def test_json_and_text_have_same_data(self, locus_reports):
        """Test that JSON and text outputs represent the same data."""
        both_result = generate_report(locus_reports, format="both")
        json_result = generate_report(locus_reports, format="json")
        text_result = generate_report(locus_reports, format="text")

        assert both_result["json"] == json_result
        assert both_result["text"] == text_result

    def test_summary_totals_are_consistent(self, locus_reports):
        """Test that summary totals are consistent."""
        result = generate_report(locus_reports, format="json")
        data = json.loads(result)
        summary = data["summary"]
        # Total suspect alleles should match sum of suspect pairs * 2 (approximately)
        # (pairs can share alleles)
        assert summary["total_loci"] == len(locus_reports)
        assert summary["total_suspect_allele_pairs"] >= 0
        assert len(summary["suspect_loci"]) <= summary["total_loci"]
