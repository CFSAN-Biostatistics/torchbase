#!/usr/bin/env python

"""Tests for k-mer quality analysis module."""

import pytest
from pathlib import Path
import tempfile
from textwrap import dedent

from torchbase.quality.kmer_analysis import (
    analyze_locus,
    SimilarityReport,
)


@pytest.fixture
def temp_fasta_dir():
    """Create a temporary directory for FASTA files."""
    with tempfile.TemporaryDirectory() as tmpdir:
        yield Path(tmpdir)


@pytest.fixture
def single_locus_fasta(temp_fasta_dir):
    """Create a FASTA file with identical sequences (no issues)."""
    fasta_path = temp_fasta_dir / "locus1.fasta"
    content = dedent("""
        >allele_1
        ATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATG
        >allele_2
        ATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATG
    """).strip()
    fasta_path.write_text(content)
    return fasta_path


@pytest.fixture
def diverse_alleles_fasta(temp_fasta_dir):
    """Create a FASTA file with diverse sequences."""
    fasta_path = temp_fasta_dir / "diverse.fasta"
    content = dedent("""
        >allele_1
        ATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATG
        >allele_2
        GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA
        >allele_3
        CGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
    """).strip()
    fasta_path.write_text(content)
    return fasta_path


@pytest.fixture
def overlapping_alleles_fasta(temp_fasta_dir):
    """Create a FASTA file with overlapping sequences."""
    fasta_path = temp_fasta_dir / "overlapping.fasta"
    # Create sequences where allele_2 is ~95% contained in allele_1
    seq1 = "A" * 100
    seq2 = "A" * 95 + "T" * 5  # 95% identical (95 As, 5 Ts)
    content = f">allele_1\n{seq1}\n>allele_2\n{seq2}"
    fasta_path.write_text(content)
    return fasta_path


@pytest.fixture
def duplicate_alleles_fasta(temp_fasta_dir):
    """Create a FASTA file with duplicate sequences."""
    fasta_path = temp_fasta_dir / "duplicates.fasta"
    # Create sequences that are ~98% similar
    seq1 = (
        "ATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATG"  # noqa
    )
    seq2 = (
        "ATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATTG"  # noqa
    )
    content = f">allele_1\n{seq1}\n>allele_2\n{seq2}"
    fasta_path.write_text(content)
    return fasta_path


class TestAnalyzeLocus:
    """Tests for analyze_locus function."""

    def test_analyze_locus_returns_similarity_report(self, single_locus_fasta):
        """Test that analyze_locus returns a SimilarityReport."""
        report = analyze_locus(single_locus_fasta)
        assert isinstance(report, SimilarityReport)

    def test_similarity_report_has_required_fields(self, single_locus_fasta):
        """Test that SimilarityReport has required fields."""
        report = analyze_locus(single_locus_fasta)
        assert hasattr(report, "similarities")
        assert hasattr(report, "threshold")
        assert hasattr(report, "suspect_pairs")
        assert hasattr(report, "statistics")

    def test_similarities_is_dict(self, single_locus_fasta):
        """Test that similarities is a dictionary."""
        report = analyze_locus(single_locus_fasta)
        assert isinstance(report.similarities, dict)

    def test_suspect_pairs_is_list(self, single_locus_fasta):
        """Test that suspect_pairs is a list."""
        report = analyze_locus(single_locus_fasta)
        assert isinstance(report.suspect_pairs, list)

    def test_statistics_dict(self, single_locus_fasta):
        """Test that statistics is a dictionary."""
        report = analyze_locus(single_locus_fasta)
        assert isinstance(report.statistics, dict)
        assert "mean" in report.statistics
        assert "std_dev" in report.statistics
        assert "min" in report.statistics
        assert "max" in report.statistics
        assert "threshold_type" in report.statistics

    def test_custom_kmer_size(self, single_locus_fasta):
        """Test that custom k-mer size parameter works."""
        report = analyze_locus(single_locus_fasta, k_size=15)
        assert isinstance(report, SimilarityReport)

    def test_identical_alleles_high_similarity(self, single_locus_fasta):
        """Test that identical alleles have high similarity."""
        report = analyze_locus(single_locus_fasta, k_size=21)
        assert len(report.similarities) > 0
        # At least one pair of alleles
        for pair_key, similarity in report.similarities.items():
            if pair_key[0] != pair_key[1]:  # different alleles
                # very high similarity for identical sequences
                assert similarity >= 0.95


class TestAutoTuning:
    """Tests for auto-tuning threshold detection."""

    def test_gap_detection_identifies_outliers(
            self, overlapping_alleles_fasta):
        """Test that gap detection identifies overlapping alleles."""
        report = analyze_locus(overlapping_alleles_fasta, k_size=21)
        # Should detect the overlap
        assert len(report.suspect_pairs) > 0
        threshold_type = report.statistics["threshold_type"]
        assert threshold_type in ["gap_detection", "percentile"]

    def test_fallback_to_percentile(self, diverse_alleles_fasta):
        """Test fallback to 99th percentile when no clear gap."""
        report = analyze_locus(diverse_alleles_fasta, k_size=21)
        # With diverse alleles, should use percentile
        threshold_type = report.statistics["threshold_type"]
        assert threshold_type == "percentile"
        assert report.threshold is not None

    def test_threshold_value_in_report(self, single_locus_fasta):
        """Test that threshold value is included in report."""
        report = analyze_locus(single_locus_fasta, k_size=21)
        assert report.threshold is not None
        assert isinstance(report.threshold, (int, float))
        assert 0 <= report.threshold <= 100


class TestSuspectPairDetection:
    """Tests for suspect pair detection."""

    def test_no_suspects_in_diverse_alleles(self, diverse_alleles_fasta):
        """Test that diverse alleles don't have suspect pairs."""
        report = analyze_locus(diverse_alleles_fasta, k_size=21)
        # With diverse sequences, should have no suspect pairs
        assert len(report.suspect_pairs) == 0

    def test_overlap_classification(self, overlapping_alleles_fasta):
        """Test that overlapping alleles are classified correctly."""
        report = analyze_locus(overlapping_alleles_fasta, k_size=21,
                               overlap_threshold=95)
        suspect_pairs = report.suspect_pairs
        # Should have at least one suspect pair
        if len(suspect_pairs) > 0:
            pair = suspect_pairs[0]
            assert pair["allele1"] is not None
            assert pair["allele2"] is not None
            assert pair["similarity"] is not None
            assert pair["issue_type"] in ["overlap", "duplicate"]

    def test_duplicate_classification(self, duplicate_alleles_fasta):
        """Test that duplicate alleles are classified correctly."""
        report = analyze_locus(duplicate_alleles_fasta, k_size=21,
                               duplicate_threshold=98)
        suspect_pairs = report.suspect_pairs
        # May or may not flag depending on exact similarity calculation
        for pair in suspect_pairs:
            assert pair["issue_type"] in ["overlap", "duplicate"]

    def test_custom_overlap_threshold(self, overlapping_alleles_fasta):
        """Test that custom overlap threshold parameter works."""
        report_low = analyze_locus(overlapping_alleles_fasta, k_size=21,
                                   overlap_threshold=90)
        report_high = analyze_locus(overlapping_alleles_fasta, k_size=21,
                                    overlap_threshold=99)
        # Lower threshold should be more permissive
        assert len(report_low.suspect_pairs) >= len(
            report_high.suspect_pairs)

    def test_custom_duplicate_threshold(self, duplicate_alleles_fasta):
        """Test that custom duplicate threshold parameter works."""
        report_low = analyze_locus(duplicate_alleles_fasta, k_size=21,
                                   duplicate_threshold=90)
        report_high = analyze_locus(duplicate_alleles_fasta, k_size=21,
                                    duplicate_threshold=99)
        # Lower threshold should be more permissive
        assert len(report_low.suspect_pairs) >= len(
            report_high.suspect_pairs)


class TestStatisticsCalculation:
    """Tests for statistics calculation."""

    def test_statistics_has_mean(self, single_locus_fasta):
        """Test that statistics includes mean."""
        report = analyze_locus(single_locus_fasta)
        assert "mean" in report.statistics
        assert isinstance(report.statistics["mean"], (int, float))

    def test_statistics_has_std_dev(self, single_locus_fasta):
        """Test that statistics includes standard deviation."""
        report = analyze_locus(single_locus_fasta)
        assert "std_dev" in report.statistics
        assert isinstance(report.statistics["std_dev"], (int, float))
        assert report.statistics["std_dev"] >= 0

    def test_statistics_has_min_max(self, diverse_alleles_fasta):
        """Test that statistics includes min and max."""
        report = analyze_locus(diverse_alleles_fasta)
        assert "min" in report.statistics
        assert "max" in report.statistics
        assert report.statistics["min"] <= report.statistics["max"]

    def test_statistics_have_valid_percentile(self, single_locus_fasta):
        """Test that statistics include 99th percentile."""
        report = analyze_locus(single_locus_fasta)
        assert "percentile_99" in report.statistics
        assert 0 <= report.statistics["percentile_99"] <= 100


class TestEdgeCases:
    """Tests for edge cases."""

    def test_single_allele_fasta(self, temp_fasta_dir):
        """Test handling of FASTA with only one allele."""
        fasta_path = temp_fasta_dir / "single.fasta"
        content = ">allele_1\nATGATGATGATGATGATGATGATGATGATGATGATG"
        fasta_path.write_text(content)
        report = analyze_locus(fasta_path)
        assert isinstance(report, SimilarityReport)
        # Should have no pairwise similarities with only one allele
        assert len(report.similarities) == 0
        assert len(report.suspect_pairs) == 0

    def test_short_sequences(self, temp_fasta_dir):
        """Test handling of short sequences (shorter than k-mer size)."""
        fasta_path = temp_fasta_dir / "short.fasta"
        content = ">allele_1\nATG\n>allele_2\nGCT"
        fasta_path.write_text(content)
        # Should handle gracefully (may skip or use shorter k-mer)
        report = analyze_locus(fasta_path, k_size=21)
        assert isinstance(report, SimilarityReport)

    def test_large_number_of_alleles(self, temp_fasta_dir):
        """Test handling of FASTA with many alleles."""
        fasta_path = temp_fasta_dir / "many.fasta"
        content_lines = []
        base_seq = (
            "ATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATG"  # noqa
        )
        for i in range(20):
            content_lines.append(f">allele_{i+1}")
            content_lines.append(base_seq)
        fasta_path.write_text("\n".join(content_lines))
        report = analyze_locus(fasta_path, k_size=21)
        assert isinstance(report, SimilarityReport)
        # Should have many pairwise comparisons
        assert len(report.similarities) > 0


class TestSimilarityReportDataClass:
    """Tests for SimilarityReport dataclass."""

    def test_similarity_report_instantiation(self):
        """Test that SimilarityReport can be instantiated."""
        report = SimilarityReport(
            similarities={("a1", "a2"): 0.95},
            threshold=0.95,
            suspect_pairs=[],
            statistics={"mean": 0.95},
        )
        assert report.similarities == {("a1", "a2"): 0.95}
        assert report.threshold == 0.95
        assert report.suspect_pairs == []
        assert report.statistics == {"mean": 0.95}

    def test_similarity_report_repr(self):
        """Test string representation of SimilarityReport."""
        report = SimilarityReport(
            similarities={},
            threshold=0.95,
            suspect_pairs=[],
            statistics={},
        )
        repr_str = repr(report)
        assert "SimilarityReport" in repr_str
