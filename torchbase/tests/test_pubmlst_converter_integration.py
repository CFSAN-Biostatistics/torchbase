#!/usr/bin/env python
"""Real integration tests for PubMLST converter (Issue #8).

These tests use real file I/O and minimal mocking. They verify actual behavior
through public interfaces, not implementation details.

Tests should fail (RED) until proper implementation exists.
"""

import pytest
import json
import tempfile
from pathlib import Path
import toml

from torchbase import Torch
from torchbase.conversions.pubmlst import convert_scheme


@pytest.fixture
def mock_bigsdb_server(monkeypatch, tmp_path):
    """Mock HTTP server returning realistic PubMLST responses."""
    import requests

    # Create test data directory
    test_data_dir = tmp_path / "bigsdb_responses"
    test_data_dir.mkdir()

    # Write realistic scheme metadata
    scheme_metadata = {
        "scheme_id": 1,
        "name": "Salmonella enterica MLST",
        "description": "7-locus MLST scheme for Salmonella enterica",
        "loci": ["aroC", "dnaN", "hemD", "hisD", "purE", "sucA", "thrA"]
    }

    # Write realistic loci responses with actual allele sequences
    loci_data = {
        "aroC": {
            "alleles": [
                ("aroC_1", "ATGAGTCGCGTACTGGCAACGGAACAGCTGGATATGCGTTATGGCGCAGCGAACCCG"),
                ("aroC_2", "ATGAGTCGCGTACTGGCAACGGAACAGCTGGATATGCGTTATGGCGCAGCGAACCCG"),
                ("aroC_3", "ATGAGTCGCGTACTGGCAACGGAACAGCTGGATAT GCGTTATGGCGCAGCGAACCCG"),
            ],
            "count": 3
        },
        "dnaN": {
            "alleles": [
                ("dnaN_1", "ATGAAACCGGTAACTGTTGATATTGGTGATCGTGCCGTGCGCGACCGTCTGGAAACCG"),
                ("dnaN_2", "ATGAAACCGGTAACTGTTGATATTGGTGATCGTGCCGTGCGCGACCGTCTGGAAACCG"),
            ],
            "count": 2
        },
    }

    # Write realistic profiles
    profiles_data = [
        {"ST": "1", "aroC": "1", "dnaN": "1", "hemD": "1", "hisD": "1", "purE": "1", "sucA": "1", "thrA": "1"},
        {"ST": "2", "aroC": "2", "dnaN": "1", "hemD": "1", "hisD": "1", "purE": "1", "sucA": "1", "thrA": "1"},
        {"ST": "3", "aroC": "1", "dnaN": "2", "hemD": "1", "hisD": "1", "purE": "1", "sucA": "1", "thrA": "1"},
    ]

    # Mock the BIGSdbClient to return this data
    from torchbase.conversions.bigsdb_client import (
        BIGSdbClient, SchemeMetadata, LocusData, ProfileTable, SchemeData
    )
    from datetime import datetime, timezone

    original_fetch = BIGSdbClient.fetch_scheme

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
            profiles=profiles_data,
            row_count=len(profiles_data),
            last_updated=datetime.now(timezone.utc)
        )

        return SchemeData(metadata=metadata, loci=loci, profiles=profiles)

    def mock_fetch_alleles(self, database, locus_id):
        """Return stub allele sequences."""
        # Return simple stub sequences for mocked tests
        return {
            f"{locus_id}_1": "ATGCGTACGTAGCTAGCTAGCTAGCTAGCTAGCT",
            f"{locus_id}_2": "ATGCGTACGTAGCTAGCTAGCTAGCTAGCTAGCC",
        }

    monkeypatch.setattr(BIGSdbClient, "fetch_scheme", mock_fetch_scheme)
    monkeypatch.setattr(BIGSdbClient, "fetch_alleles", mock_fetch_alleles)

    return {
        "scheme_metadata": scheme_metadata,
        "loci_data": loci_data,
        "profiles_data": profiles_data,
    }


class TestRealConversion:
    """Test actual conversion with real file I/O."""

    def test_converts_scheme_to_torch_directory(self, mock_bigsdb_server, tmp_path):
        """Should create a torch directory with proper structure."""
        output_dir = tmp_path / "output"

        torch_path = convert_scheme(
            database_url="http://test.pubmlst.org/api",
            scheme_id=1,
            output_path=str(output_dir),
        )

        # Verify torch directory exists
        assert Path(torch_path).exists()
        assert Path(torch_path).is_dir()

        # Verify it follows namespace/name/version.torch pattern
        assert "pubmlst" in str(torch_path)
        assert ".torch" in str(torch_path)

    def test_creates_schemes_hierarchy(self, mock_bigsdb_server, tmp_path):
        """Should create schemes/<organism>/ subdirectory structure."""
        output_dir = tmp_path / "output"

        torch_path = convert_scheme(
            database_url="http://test.pubmlst.org/api",
            scheme_id=1,
            output_path=str(output_dir),
        )

        torch_dir = Path(torch_path)
        schemes_dir = torch_dir / "schemes"

        # schemes/ should exist
        assert schemes_dir.exists()
        assert schemes_dir.is_dir()

        # Should have at least one organism subdirectory
        organism_dirs = list(schemes_dir.iterdir())
        assert len(organism_dirs) > 0
        assert organism_dirs[0].is_dir()

    def test_creates_alleles_directory(self, mock_bigsdb_server, tmp_path):
        """Should create alleles/ subdirectory with FASTA files."""
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

        # alleles/ should exist
        assert alleles_dir.exists()
        assert alleles_dir.is_dir()

        # Should contain FASTA files
        fasta_files = list(alleles_dir.glob("*.fasta"))
        assert len(fasta_files) > 0

    def test_creates_profiles_tsv(self, mock_bigsdb_server, tmp_path):
        """Should create profiles.tsv in organism directory."""
        output_dir = tmp_path / "output"

        torch_path = convert_scheme(
            database_url="http://test.pubmlst.org/api",
            scheme_id=1,
            output_path=str(output_dir),
        )

        torch_dir = Path(torch_path)
        schemes_dir = torch_dir / "schemes"
        organism_dir = list(schemes_dir.iterdir())[0]
        profiles_path = organism_dir / "profiles.tsv"

        # profiles.tsv should exist
        assert profiles_path.exists()
        assert profiles_path.is_file()

        # Should be valid TSV with ST column
        content = profiles_path.read_text()
        assert "ST" in content
        assert "\t" in content  # Tab-separated

    def test_creates_metadata_toml(self, mock_bigsdb_server, tmp_path):
        """Should create metadata.toml at torch root."""
        output_dir = tmp_path / "output"

        torch_path = convert_scheme(
            database_url="http://test.pubmlst.org/api",
            scheme_id=1,
            output_path=str(output_dir),
        )

        torch_dir = Path(torch_path)
        metadata_path = torch_dir / "metadata.toml"

        # metadata.toml should exist
        assert metadata_path.exists()
        assert metadata_path.is_file()

        # Should be valid TOML
        metadata = toml.load(metadata_path)
        assert isinstance(metadata, dict)

    def test_metadata_contains_provenance(self, mock_bigsdb_server, tmp_path):
        """metadata.toml should have [provenance] section."""
        output_dir = tmp_path / "output"

        torch_path = convert_scheme(
            database_url="http://test.pubmlst.org/api",
            scheme_id=1,
            output_path=str(output_dir),
        )

        torch_dir = Path(torch_path)
        metadata = toml.load(torch_dir / "metadata.toml")

        assert "provenance" in metadata
        prov = metadata["provenance"]

        # Should capture source information
        assert "source" in prov or "database_url" in prov
        assert "scheme_id" in prov or "fetch_date" in prov

    def test_metadata_contains_data_quality(self, mock_bigsdb_server, tmp_path):
        """metadata.toml should have [data_quality] section."""
        output_dir = tmp_path / "output"

        torch_path = convert_scheme(
            database_url="http://test.pubmlst.org/api",
            scheme_id=1,
            output_path=str(output_dir),
        )

        torch_dir = Path(torch_path)
        metadata = toml.load(torch_dir / "metadata.toml")

        assert "data_quality" in metadata
        quality = metadata["data_quality"]

        # Should indicate if k-mer analysis was performed
        assert "kmer_analysis_performed" in quality or "kmer_size" in quality

    def test_metadata_contains_typing_section(self, mock_bigsdb_server, tmp_path):
        """metadata.toml should have [typing] section."""
        output_dir = tmp_path / "output"

        torch_path = convert_scheme(
            database_url="http://test.pubmlst.org/api",
            scheme_id=1,
            output_path=str(output_dir),
        )

        torch_dir = Path(torch_path)
        metadata = toml.load(torch_dir / "metadata.toml")

        assert "typing" in metadata
        typing = metadata["typing"]

        # Should describe the typing scheme
        assert "scheme_name" in typing or "loci_count" in typing

    def test_metadata_contains_schemes_section(self, mock_bigsdb_server, tmp_path):
        """metadata.toml should have [schemes] section."""
        output_dir = tmp_path / "output"

        torch_path = convert_scheme(
            database_url="http://test.pubmlst.org/api",
            scheme_id=1,
            output_path=str(output_dir),
        )

        torch_dir = Path(torch_path)
        metadata = toml.load(torch_dir / "metadata.toml")

        assert "schemes" in metadata
        schemes = metadata["schemes"]

        # Should list the organism scheme
        assert len(schemes) > 0

    def test_creates_quality_json(self, mock_bigsdb_server, tmp_path):
        """Should create quality.json at torch root."""
        output_dir = tmp_path / "output"

        torch_path = convert_scheme(
            database_url="http://test.pubmlst.org/api",
            scheme_id=1,
            output_path=str(output_dir),
        )

        torch_dir = Path(torch_path)
        quality_path = torch_dir / "quality.json"

        # quality.json should exist
        assert quality_path.exists()
        assert quality_path.is_file()

        # Should be valid JSON
        quality = json.loads(quality_path.read_text())
        assert isinstance(quality, dict)

    def test_quality_json_includes_kmer_analysis(self, mock_bigsdb_server, tmp_path):
        """quality.json should include k-mer analysis results."""
        output_dir = tmp_path / "output"

        torch_path = convert_scheme(
            database_url="http://test.pubmlst.org/api",
            scheme_id=1,
            output_path=str(output_dir),
        )

        torch_dir = Path(torch_path)
        quality = json.loads((torch_dir / "quality.json").read_text())

        # Should have k-mer analysis section
        assert "kmer_analysis" in quality
        kmer = quality["kmer_analysis"]

        assert "kmer_size" in kmer or "performed" in kmer

    def test_torch_is_loadable(self, mock_bigsdb_server, tmp_path):
        """Generated torch should be loadable via Torch.load()."""
        output_dir = tmp_path / "output"

        torch_path = convert_scheme(
            database_url="http://test.pubmlst.org/api",
            scheme_id=1,
            output_path=str(output_dir),
        )

        # Should load without error
        torch = Torch.load(Path(torch_path))

        # Should have schemes
        assert hasattr(torch, "schemes")
        assert len(torch.schemes) > 0

    def test_respects_kmer_size_parameter(self, mock_bigsdb_server, tmp_path):
        """Should respect --kmer-size parameter in quality analysis."""
        output_dir = tmp_path / "output"
        custom_kmer_size = 21

        torch_path = convert_scheme(
            database_url="http://test.pubmlst.org/api",
            scheme_id=1,
            output_path=str(output_dir),
            kmer_size=custom_kmer_size,
        )

        torch_dir = Path(torch_path)

        # Check metadata.toml
        metadata = toml.load(torch_dir / "metadata.toml")
        assert metadata["data_quality"]["kmer_size"] == custom_kmer_size

        # Check quality.json
        quality = json.loads((torch_dir / "quality.json").read_text())
        assert quality["kmer_analysis"]["kmer_size"] == custom_kmer_size

    def test_respects_overlap_threshold_parameter(self, mock_bigsdb_server, tmp_path):
        """Should respect --overlap-threshold parameter."""
        output_dir = tmp_path / "output"
        custom_threshold = 0.85

        torch_path = convert_scheme(
            database_url="http://test.pubmlst.org/api",
            scheme_id=1,
            output_path=str(output_dir),
            overlap_threshold=custom_threshold,
        )

        torch_dir = Path(torch_path)

        # Check metadata contains the threshold
        metadata = toml.load(torch_dir / "metadata.toml")
        assert metadata["data_quality"]["overlap_threshold"] == custom_threshold

    def test_respects_duplicate_threshold_parameter(self, mock_bigsdb_server, tmp_path):
        """Should respect --duplicate-threshold parameter."""
        output_dir = tmp_path / "output"
        custom_threshold = 0.98

        torch_path = convert_scheme(
            database_url="http://test.pubmlst.org/api",
            scheme_id=1,
            output_path=str(output_dir),
            duplicate_threshold=custom_threshold,
        )

        torch_dir = Path(torch_path)

        # Check metadata contains the threshold
        metadata = toml.load(torch_dir / "metadata.toml")
        assert metadata["data_quality"]["duplicate_threshold"] == custom_threshold
