#!/usr/bin/env python
"""RED tests verifying real allele sequence fetching (not stubs).

These tests should FAIL until BIGSdb client properly fetches allele sequences
and converter writes them to FASTA files.
"""

import pytest
from pathlib import Path
import toml

from torchbase.conversions.pubmlst import convert_scheme


@pytest.fixture
def mock_bigsdb_with_sequences(monkeypatch, tmp_path):
    """Mock BIGSdb client that returns actual allele sequences."""
    from torchbase.conversions.bigsdb_client import (
        BIGSdbClient, SchemeMetadata, LocusData, ProfileTable, SchemeData
    )
    from datetime import datetime, timezone

    # Define realistic allele sequences (different lengths, different sequences)
    allele_sequences = {
        "aroC": {
            "aroC_1": "ATGAGTCGCGTACTGGCAACGGAACAGCTGGATATGCGTTATGGCGCAGCGAACCCGAAAATCGACGACGTGGTGCTGGCAGAAGTGGTTGCCGAACTGGCAGATACCCTGCGCGAAGAGTTCAGCCAGGTGATGCTGGACCCACAGTTCGAGCTG",
            "aroC_2": "ATGAGTCGCGTACTGGCAACGGAACAGCTGGATATGCGTTATGGCGCAGCGAACCCGAAAATCGACGACGTGGTGCTGGCAGAAGTGGTTGCCGAACTGGCAGATACCCTGCGCGAAGAGTTCAGCCAGGTGATGCTGGACCCACAGTTCGAGCTA",
            "aroC_3": "ATGAGTCGCGTACTGGCAACGGAACAGCTGGATATGCGTTATGGCGCAGCGAACCCGAAAATCGACGACGTGGTGCTGGCAGAAGTGGTTGCCGAACTGGCAGATACCCTGCGCGAAGAGTTCAGCCAGGTGATGCTGGACCCACAGTTCGAGCTC",
        },
        "dnaN": {
            "dnaN_1": "ATGAAACCGGTAACTGTTGATATTGGTGATCGTGCCGTGCGCGACCGTCTGGAAACCGTCGTGAAAGATGGTGTGGAAGAAGGCGAACCGCTGCTGGCGGATGTGGGTATCCGCGCAATGCGTGAAGCGACCAAAGTGGTGGCAGAAGCGGCAA",
            "dnaN_2": "ATGAAACCGGTAACTGTTGATATTGGTGATCGTGCCGTGCGCGACCGTCTGGAAACCGTCGTGAAAGATGGTGTGGAAGAAGGCGAACCGCTGCTGGCGGATGTGGGTATCCGCGCAATGCGTGAAGCGACCAAAGTGGTGGCAGAAGCGGCAT",
        },
    }

    def mock_fetch_scheme(self, database, scheme_id):
        """Return mock scheme data."""
        metadata = SchemeMetadata(
            scheme_id=1,
            name="salmonella_mlst",
            description="MLST scheme",
            last_updated=datetime.now(timezone.utc)
        )

        loci = [
            LocusData(
                locus_id="aroC",
                locus_name="aroC",
                alleles_count=3,
                last_updated=datetime.now(timezone.utc)
            ),
            LocusData(
                locus_id="dnaN",
                locus_name="dnaN",
                alleles_count=2,
                last_updated=datetime.now(timezone.utc)
            ),
        ]

        profiles = ProfileTable(
            profiles=[
                {"ST": "1", "aroC": "1", "dnaN": "1"},
                {"ST": "2", "aroC": "2", "dnaN": "1"},
                {"ST": "3", "aroC": "1", "dnaN": "2"},
            ],
            row_count=3,
            last_updated=datetime.now(timezone.utc)
        )

        return SchemeData(metadata=metadata, loci=loci, profiles=profiles)

    def mock_fetch_alleles(self, database, locus_id):
        """Return allele sequences for a locus."""
        if locus_id not in allele_sequences:
            return {}
        return allele_sequences[locus_id]

    monkeypatch.setattr(BIGSdbClient, "fetch_scheme", mock_fetch_scheme)
    monkeypatch.setattr(BIGSdbClient, "fetch_alleles", mock_fetch_alleles)

    return allele_sequences


class TestRealAlleleSequences:
    """Test that real allele sequences are fetched and written."""

    def test_fasta_files_contain_real_sequences_not_stubs(self, mock_bigsdb_with_sequences, tmp_path):
        """FASTA files should contain actual sequences from BIGSdb, not stub sequences.

        This test will FAIL until the converter actually fetches allele sequences
        via BIGSdbClient.fetch_alleles() instead of using _write_stub_fasta().
        """
        output_dir = tmp_path / "output"

        torch_path = convert_scheme(
            database_url="http://test.pubmlst.org/api",
            scheme_id=1,
            output_path=str(output_dir),
        )

        torch_dir = Path(torch_path)
        schemes_dir = torch_dir / "schemes"
        organism_dir = list(schemes_dir.iterdir())[0]
        alleles_dir = organism_dir / "alleles"

        # Read aroC FASTA
        aroc_fasta = alleles_dir / "aroC.fasta"
        assert aroc_fasta.exists()

        content = aroc_fasta.read_text()

        # Should NOT contain the stub sequence pattern (repeated ATGATG...)
        stub_pattern = "ATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATG"
        assert stub_pattern not in content, "FASTA contains stub sequences, not real data from BIGSdb"

        # Should contain real allele headers
        assert ">aroC_1" in content
        assert ">aroC_2" in content
        assert ">aroC_3" in content

        # Should contain actual sequences (check for presence of real sequence start)
        assert "ATGAGTCGCGTACTGGCAACGGAACAGCTGGATATGCGTTATGGCGCAGCGAACCCGAAAAT" in content

    def test_different_loci_have_different_sequences(self, mock_bigsdb_with_sequences, tmp_path):
        """Different loci should have different sequences, not duplicates.

        This verifies that each locus is fetched independently from BIGSdb.
        """
        output_dir = tmp_path / "output"

        torch_path = convert_scheme(
            database_url="http://test.pubmlst.org/api",
            scheme_id=1,
            output_path=str(output_dir),
        )

        torch_dir = Path(torch_path)
        schemes_dir = torch_dir / "schemes"
        organism_dir = list(schemes_dir.iterdir())[0]
        alleles_dir = organism_dir / "alleles"

        # Read both FASTA files
        aroc_content = (alleles_dir / "aroC.fasta").read_text()
        dnaN_content = (alleles_dir / "dnaN.fasta").read_text()

        # Sequences should be different (not copy-paste)
        assert aroc_content != dnaN_content
        assert "aroC_1" in aroc_content
        assert "dnaN_1" in dnaN_content
        assert "aroC_1" not in dnaN_content
        assert "dnaN_1" not in aroc_content

    def test_allele_count_matches_sequences_in_fasta(self, mock_bigsdb_with_sequences, tmp_path):
        """Number of sequences in FASTA should match allele count from metadata.

        This verifies that all alleles are fetched and written.
        """
        output_dir = tmp_path / "output"

        torch_path = convert_scheme(
            database_url="http://test.pubmlst.org/api",
            scheme_id=1,
            output_path=str(output_dir),
        )

        torch_dir = Path(torch_path)
        schemes_dir = torch_dir / "schemes"
        organism_dir = list(schemes_dir.iterdir())[0]
        alleles_dir = organism_dir / "alleles"

        # aroC should have 3 alleles
        aroc_content = (alleles_dir / "aroC.fasta").read_text()
        aroc_count = aroc_content.count(">aroC_")
        assert aroc_count == 3, f"Expected 3 aroC alleles, found {aroc_count}"

        # dnaN should have 2 alleles
        dnaN_content = (alleles_dir / "dnaN.fasta").read_text()
        dnaN_count = dnaN_content.count(">dnaN_")
        assert dnaN_count == 2, f"Expected 2 dnaN alleles, found {dnaN_count}"

    def test_fasta_format_is_valid(self, mock_bigsdb_with_sequences, tmp_path):
        """FASTA files should be properly formatted (header line, sequence lines).

        This verifies the FASTA writer handles line wrapping correctly.
        """
        output_dir = tmp_path / "output"

        torch_path = convert_scheme(
            database_url="http://test.pubmlst.org/api",
            scheme_id=1,
            output_path=str(output_dir),
        )

        torch_dir = Path(torch_path)
        schemes_dir = torch_dir / "schemes"
        organism_dir = list(schemes_dir.iterdir())[0]
        alleles_dir = organism_dir / "alleles"

        aroc_content = (alleles_dir / "aroC.fasta").read_text()
        lines = aroc_content.strip().split("\n")

        # First line should be a header
        assert lines[0].startswith(">"), "First line should be FASTA header"

        # Should alternate between headers and sequences (roughly)
        header_count = sum(1 for line in lines if line.startswith(">"))
        assert header_count == 3, f"Expected 3 headers for 3 alleles, found {header_count}"

        # Sequence lines should not be empty
        for i, line in enumerate(lines):
            if not line.startswith(">"):
                assert len(line) > 0, f"Empty sequence line at line {i}"
                # Sequence lines should only contain valid nucleotides
                assert all(c in "ATGCNatgcn" for c in line), f"Invalid nucleotide characters at line {i}"
