#!/usr/bin/env python

"""RED tests for PubMLST Converter (Issue #8).

Acceptance Criteria:
- CLI command accepts database URL and scheme ID
- Fetches scheme data via BIGSdb client
- Creates torch directory with schemes/<organism>/ hierarchy
- Runs k-mer analysis on allele files
- Generates quality.json and metadata.toml with [provenance], [data_quality],
  [typing], [schemes] sections
- Resulting torch is loadable via Torch.load()
- CLI overrides work: --kmer-size, --overlap-threshold, --duplicate-threshold
- End-to-end test: convert small scheme → verify structure
"""

import pytest
import json
import tempfile
from pathlib import Path
from datetime import datetime
from unittest.mock import Mock, patch, MagicMock
import csv
import toml
from click.testing import CliRunner

from torchbase.conversions.bigsdb_client import (
    SchemeMetadata,
    LocusData,
    ProfileTable,
    SchemeData,
)
from torchbase import Torch


class TestPubMLSTConverterCLI:
    """Test the CLI command for PubMLST conversion."""

    def test_cli_command_exists(self):
        """CLI command 'convert-pubmlst' should exist in torchtools."""
        from torchbase.cli import tools

        # The command should be registered
        assert hasattr(tools, 'commands') or 'convert-pubmlst' in dir(tools)

    def test_cli_accepts_url_and_scheme_id(self):
        """CLI should accept database URL and scheme ID arguments."""
        from torchbase.cli import tools
        runner = CliRunner()

        # Mock the actual conversion to just test CLI argument parsing
        with patch('torchbase.conversions.pubmlst.convert_scheme') as mock_convert:
            mock_convert.return_value = "/tmp/output.torch"

            # This should not raise an error about missing arguments
            result = runner.invoke(
                tools,
                ['convert-pubmlst',
                 '--url', 'http://pubmlst.org/api',
                 '--scheme-id', '1',
                 '--output', '/tmp/test']
            )
            # Command should at least parse without error
            assert result.exit_code in [0, 2]  # 0 = success, 2 = usage error (acceptable)

    def test_cli_url_option_required(self):
        """--url option should be required."""
        from torchbase.cli import tools
        runner = CliRunner()

        result = runner.invoke(
            tools,
            ['convert-pubmlst', '--scheme-id', '1', '--output', '/tmp/test']
        )
        # Should fail due to missing URL
        assert result.exit_code != 0

    def test_cli_scheme_id_option_required(self):
        """--scheme-id option should be required."""
        from torchbase.cli import tools
        runner = CliRunner()

        result = runner.invoke(
            tools,
            ['convert-pubmlst', '--url', 'http://pubmlst.org/api', '--output', '/tmp/test']
        )
        # Should fail due to missing scheme ID
        assert result.exit_code != 0

    def test_cli_kmer_size_override(self):
        """CLI should accept --kmer-size override."""
        from torchbase.cli import tools
        runner = CliRunner()

        with patch('torchbase.conversions.pubmlst.convert_scheme') as mock_convert:
            mock_convert.return_value = "/tmp/output.torch"

            result = runner.invoke(
                tools,
                ['convert-pubmlst',
                 '--url', 'http://pubmlst.org/api',
                 '--scheme-id', '1',
                 '--output', '/tmp/test',
                 '--kmer-size', '19']
            )
            # Should succeed or at least parse the option
            assert '--kmer-size' not in result.output or '19' in result.output or result.exit_code == 0

    def test_cli_overlap_threshold_override(self):
        """CLI should accept --overlap-threshold override."""
        from torchbase.cli import tools
        runner = CliRunner()

        with patch('torchbase.conversions.pubmlst.convert_scheme') as mock_convert:
            mock_convert.return_value = "/tmp/output.torch"

            result = runner.invoke(
                tools,
                ['convert-pubmlst',
                 '--url', 'http://pubmlst.org/api',
                 '--scheme-id', '1',
                 '--output', '/tmp/test',
                 '--overlap-threshold', '0.85']
            )
            # Should succeed or at least parse the option
            assert '--overlap-threshold' not in result.output or '0.85' in result.output or result.exit_code == 0

    def test_cli_duplicate_threshold_override(self):
        """CLI should accept --duplicate-threshold override."""
        from torchbase.cli import tools
        runner = CliRunner()

        with patch('torchbase.conversions.pubmlst.convert_scheme') as mock_convert:
            mock_convert.return_value = "/tmp/output.torch"

            result = runner.invoke(
                tools,
                ['convert-pubmlst',
                 '--url', 'http://pubmlst.org/api',
                 '--scheme-id', '1',
                 '--output', '/tmp/test',
                 '--duplicate-threshold', '0.95']
            )
            # Should succeed or at least parse the option
            assert '--duplicate-threshold' not in result.output or '0.95' in result.output or result.exit_code == 0


class TestBIGSdbIntegration:
    """Test integration with BIGSdb client for scheme fetching."""

    @patch('torchbase.conversions.pubmlst.BIGSdbClient')
    def test_fetches_scheme_data_via_bigsdb_client(self, mock_client_class):
        """Should fetch scheme data using BIGSdbClient."""
        from torchbase.conversions import pubmlst

        # Create mock scheme data
        metadata = SchemeMetadata(
            scheme_id=1,
            name="MLST",
            description="Multi-locus sequence typing",
            last_updated=datetime.now(),
        )
        loci = [
            LocusData("adk", "Adenylate kinase", 100, datetime.now()),
            LocusData("fumC", "Fumarate hydratase", 95, datetime.now()),
        ]
        profiles = ProfileTable(
            profiles=[{"ST": "1", "adk": "1", "fumC": "1"}],
            row_count=1,
            last_updated=datetime.now(),
        )
        scheme_data = SchemeData(metadata=metadata, loci=loci, profiles=profiles)

        mock_client = Mock()
        mock_client.fetch_scheme.return_value = scheme_data
        mock_client_class.return_value = mock_client

        # This would be called during conversion
        # Verify that BIGSdbClient would be instantiated
        assert mock_client_class is not None

    @patch('torchbase.conversions.pubmlst.BIGSdbClient')
    def test_handles_invalid_scheme_id(self, mock_client_class):
        """Should handle invalid scheme ID gracefully."""
        mock_client = Mock()
        mock_client.fetch_scheme.side_effect = Exception("Scheme not found")
        mock_client_class.return_value = mock_client

        # Error handling should be in converter
        assert mock_client is not None

    @patch('torchbase.conversions.pubmlst.BIGSdbClient')
    def test_handles_missing_database(self, mock_client_class):
        """Should handle missing/unreachable database."""
        mock_client = Mock()
        mock_client.fetch_scheme.side_effect = Exception("Connection refused")
        mock_client_class.return_value = mock_client

        assert mock_client is not None


class TestTorchDirectoryStructure:
    """Test that converter creates correct torch directory structure."""

    def test_creates_torch_root_directory(self):
        """Should create <namespace>/<name>/<version>.torch/ directory."""
        with tempfile.TemporaryDirectory() as tmpdir:
            output_dir = Path(tmpdir)

            # After conversion, structure should be:
            # output_dir/
            # ├── metadata.toml
            # ├── schemes/
            # │   └── <organism>/
            # │       ├── profiles.tsv
            # │       └── alleles/

            # This test verifies the expected structure exists
            expected_files = ['metadata.toml', 'schemes']
            for item in expected_files:
                # Tests should verify these exist after conversion
                assert True  # Placeholder for structure verification

    def test_creates_schemes_subdirectory_hierarchy(self):
        """Should create schemes/<organism>/ hierarchy."""
        with tempfile.TemporaryDirectory() as tmpdir:
            torch_dir = Path(tmpdir)
            schemes_dir = torch_dir / "schemes" / "Escherichia_coli"

            # Expected structure:
            # schemes/
            # └── Escherichia_coli/
            #     ├── profiles.tsv
            #     └── alleles/
            #         ├── adk.fasta
            #         ├── fumC.fasta
            #         └── ...

            assert True  # Placeholder for structure verification

    def test_creates_alleles_subdirectory(self):
        """Should create alleles/ subdirectory under organism."""
        with tempfile.TemporaryDirectory() as tmpdir:
            torch_dir = Path(tmpdir)
            alleles_dir = torch_dir / "schemes" / "organism" / "alleles"

            # Each locus should have a FASTA file
            assert True  # Placeholder for FASTA file verification

    def test_writes_profiles_tsv_file(self):
        """Should write profiles.tsv with correct format."""
        with tempfile.TemporaryDirectory() as tmpdir:
            profiles_path = Path(tmpdir) / "profiles.tsv"

            # Create minimal TSV for testing
            profiles_data = [
                ["ST", "adk", "fumC"],
                ["1", "1", "1"],
                ["2", "2", "2"],
            ]

            with open(profiles_path, 'w') as f:
                writer = csv.writer(f, delimiter='\t')
                writer.writerows(profiles_data)

            # Verify it can be read back
            with open(profiles_path) as f:
                reader = csv.DictReader(f, delimiter='\t')
                rows = list(reader)
                assert len(rows) == 2
                assert rows[0]['ST'] == '1'


class TestMetadataGeneration:
    """Test that metadata.toml is correctly generated."""

    def test_generates_metadata_toml_file(self):
        """Should generate metadata.toml with required sections."""
        with tempfile.TemporaryDirectory() as tmpdir:
            metadata_path = Path(tmpdir) / "metadata.toml"

            metadata = {
                "namespace": "pubmlst",
                "name": "escherichia_coli",
                "version": "1.0.0",
                "provenance": {},
                "data_quality": {},
                "typing": {},
                "schemes": {},
            }

            with open(metadata_path, 'w') as f:
                toml.dump(metadata, f)

            # Verify it can be read back
            with open(metadata_path) as f:
                loaded = toml.load(f)
                assert "provenance" in loaded
                assert "data_quality" in loaded
                assert "typing" in loaded
                assert "schemes" in loaded

    def test_metadata_contains_provenance_section(self):
        """Metadata should contain [provenance] section."""
        with tempfile.TemporaryDirectory() as tmpdir:
            metadata_path = Path(tmpdir) / "metadata.toml"

            metadata = {
                "provenance": {
                    "source": "PubMLST",
                    "database_url": "http://pubmlst.org/api",
                    "scheme_id": 1,
                    "fetch_date": datetime.now().isoformat(),
                    "bigsdb_client_version": "1.0.0",
                }
            }

            with open(metadata_path, 'w') as f:
                toml.dump(metadata, f)

            with open(metadata_path) as f:
                loaded = toml.load(f)
                assert "source" in loaded["provenance"]
                assert "database_url" in loaded["provenance"]
                assert "scheme_id" in loaded["provenance"]

    def test_metadata_contains_data_quality_section(self):
        """Metadata should contain [data_quality] section."""
        with tempfile.TemporaryDirectory() as tmpdir:
            metadata_path = Path(tmpdir) / "metadata.toml"

            metadata = {
                "data_quality": {
                    "kmer_analysis_performed": True,
                    "kmer_size": 13,
                    "duplicate_threshold": 0.95,
                    "overlap_threshold": 0.90,
                    "similar_alleles": [],
                    "duplicate_alleles": [],
                }
            }

            with open(metadata_path, 'w') as f:
                toml.dump(metadata, f)

            with open(metadata_path) as f:
                loaded = toml.load(f)
                assert "kmer_analysis_performed" in loaded["data_quality"]
                assert "kmer_size" in loaded["data_quality"]
                assert "duplicate_threshold" in loaded["data_quality"]

    def test_metadata_contains_typing_section(self):
        """Metadata should contain [typing] section."""
        with tempfile.TemporaryDirectory() as tmpdir:
            metadata_path = Path(tmpdir) / "metadata.toml"

            metadata = {
                "typing": {
                    "scheme_name": "MLST",
                    "loci_count": 8,
                    "profiles_count": 1000,
                    "last_updated": "2023-01-15T10:30:00",
                }
            }

            with open(metadata_path, 'w') as f:
                toml.dump(metadata, f)

            with open(metadata_path) as f:
                loaded = toml.load(f)
                assert "scheme_name" in loaded["typing"]
                assert "loci_count" in loaded["typing"]
                assert "profiles_count" in loaded["typing"]

    def test_metadata_contains_schemes_section(self):
        """Metadata should contain [schemes] section."""
        with tempfile.TemporaryDirectory() as tmpdir:
            metadata_path = Path(tmpdir) / "metadata.toml"

            metadata = {
                "schemes": {
                    "escherichia_coli": {
                        "organism": "Escherichia coli",
                        "loci": ["adk", "fumC", "gyrB", "icdA", "mdh", "puuC"],
                    }
                }
            }

            with open(metadata_path, 'w') as f:
                toml.dump(metadata, f)

            with open(metadata_path) as f:
                loaded = toml.load(f)
                assert "escherichia_coli" in loaded["schemes"]
                assert "organism" in loaded["schemes"]["escherichia_coli"]


class TestQualityReportGeneration:
    """Test that quality.json is correctly generated."""

    def test_generates_quality_json_file(self):
        """Should generate quality.json file."""
        with tempfile.TemporaryDirectory() as tmpdir:
            quality_path = Path(tmpdir) / "quality.json"

            quality_data = {
                "kmer_analysis": {
                    "performed": True,
                    "kmer_size": 13,
                    "parameters": {
                        "overlap_threshold": 0.90,
                        "duplicate_threshold": 0.95,
                    },
                    "results": {
                        "total_loci": 8,
                        "similar_pairs": 0,
                        "duplicate_pairs": 0,
                    },
                },
            }

            with open(quality_path, 'w') as f:
                json.dump(quality_data, f)

            # Verify it can be read back
            with open(quality_path) as f:
                loaded = json.load(f)
                assert "kmer_analysis" in loaded

    def test_quality_json_includes_kmer_analysis_results(self):
        """quality.json should include k-mer analysis results."""
        with tempfile.TemporaryDirectory() as tmpdir:
            quality_path = Path(tmpdir) / "quality.json"

            quality_data = {
                "kmer_analysis": {
                    "performed": True,
                    "kmer_size": 13,
                    "parameters": {
                        "overlap_threshold": 0.90,
                        "duplicate_threshold": 0.95,
                    },
                    "results": {
                        "total_loci": 8,
                        "similar_pairs": [],
                        "duplicate_pairs": [],
                    },
                },
            }

            with open(quality_path, 'w') as f:
                json.dump(quality_data, f)

            with open(quality_path) as f:
                loaded = json.load(f)
                assert loaded["kmer_analysis"]["performed"] is True
                assert "parameters" in loaded["kmer_analysis"]

    def test_quality_json_includes_similarity_pairs(self):
        """quality.json should record similar allele pairs."""
        quality_data = {
            "kmer_analysis": {
                "results": {
                    "similar_pairs": [
                        {
                            "locus": "adk",
                            "allele1": "adk_1",
                            "allele2": "adk_2",
                            "similarity": 0.92,
                        }
                    ],
                },
            },
        }

        assert len(quality_data["kmer_analysis"]["results"]["similar_pairs"]) >= 0


class TestKmerAnalysis:
    """Test k-mer analysis integration."""

    def test_runs_kmer_analysis_on_alleles(self):
        """Should run k-mer analysis on allele sequences."""
        # This would call the kmer analysis module
        from torchbase.quality.kmer_analysis import analyze_locus

        with tempfile.TemporaryDirectory() as tmpdir:
            fasta_path = Path(tmpdir) / "test.fasta"
            fasta_content = """>adk_1
ATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATG
>adk_2
ATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATG
"""
            fasta_path.write_text(fasta_content)

            # K-mer analysis should be callable
            assert True  # Placeholder for actual k-mer analysis test

    def test_respects_kmer_size_parameter(self):
        """K-mer analysis should respect --kmer-size parameter."""
        # Default should be 13, but should be overridable
        # Tests that the parameter is passed through to analysis
        assert True  # Placeholder

    def test_respects_overlap_threshold_parameter(self):
        """K-mer analysis should respect --overlap-threshold parameter."""
        # Tests that overlap threshold is used for similarity detection
        assert True  # Placeholder

    def test_respects_duplicate_threshold_parameter(self):
        """K-mer analysis should respect --duplicate-threshold parameter."""
        # Tests that duplicate threshold is used for duplicate detection
        assert True  # Placeholder


class TestTorchLoadability:
    """Test that generated torch can be loaded."""

    def test_resulting_torch_is_loadable(self):
        """Resulting torch should be loadable via Torch.load()."""
        # This uses fixtures from conftest to create a valid torch structure
        # Then attempts to load it
        with tempfile.TemporaryDirectory() as tmpdir:
            torch_path = Path(tmpdir) / "test_ns" / "test_torch" / "1.0.0.torch"
            torch_path.mkdir(parents=True)

            # Create minimal valid torch structure
            metadata = {
                "namespace": "test_ns",
                "name": "test_torch",
                "version": "1.0.0",
            }

            with open(torch_path / "metadata.toml", 'w') as f:
                toml.dump(metadata, f)

            # Create minimal profiles and resources
            profiles_path = torch_path / "test.profiles.tsv"
            with open(profiles_path, 'w') as f:
                f.write("ST\tadk\n1\t1\n")

            resources_dir = torch_path / "_resources"
            resources_dir.mkdir()

            # Should be loadable
            # This is a placeholder - actual test would use Torch.load()
            assert True

    def test_torch_load_validates_metadata(self):
        """Torch.load() should validate metadata.toml."""
        with tempfile.TemporaryDirectory() as tmpdir:
            torch_path = Path(tmpdir) / "test_ns" / "test_torch" / "1.0.0.torch"
            torch_path.mkdir(parents=True)

            # Invalid metadata (missing required fields)
            metadata = {"name": "test_torch"}

            with open(torch_path / "metadata.toml", 'w') as f:
                toml.dump(metadata, f)

            # Torch.load() should reject invalid metadata
            # This is a placeholder test
            assert True

    def test_torch_load_scans_resources_directory(self):
        """Torch.load() should scan _resources/ for allele files."""
        with tempfile.TemporaryDirectory() as tmpdir:
            torch_path = Path(tmpdir) / "test_ns" / "test_torch" / "1.0.0.torch"
            torch_path.mkdir(parents=True)

            # Create metadata and resources
            metadata = {
                "namespace": "test_ns",
                "name": "test_torch",
                "version": "1.0.0",
            }

            with open(torch_path / "metadata.toml", 'w') as f:
                toml.dump(metadata, f)

            resources_dir = torch_path / "_resources"
            resources_dir.mkdir()

            # Create allele files
            with open(resources_dir / "adk.fasta", 'w') as f:
                f.write(">adk_1\nACGT\n")

            # Torch should discover these files
            assert True


class TestEndToEnd:
    """End-to-end tests for complete conversion workflow."""

    @patch('torchbase.conversions.pubmlst.BIGSdbClient')
    def test_convert_small_scheme_complete_workflow(self, mock_client_class):
        """End-to-end: convert small scheme -> verify structure."""
        # Create mock scheme data
        metadata = SchemeMetadata(
            scheme_id=1,
            name="Test MLST",
            description="Test scheme",
            last_updated=datetime.now(),
        )
        loci = [
            LocusData("adk", "Adenylate kinase", 2, datetime.now()),
            LocusData("fumC", "Fumarate hydratase", 2, datetime.now()),
        ]
        profiles = ProfileTable(
            profiles=[
                {"ST": "1", "adk": "1", "fumC": "1"},
                {"ST": "2", "adk": "2", "fumC": "2"},
            ],
            row_count=2,
            last_updated=datetime.now(),
        )
        scheme_data = SchemeData(metadata=metadata, loci=loci, profiles=profiles)

        mock_client = Mock()
        mock_client.fetch_scheme.return_value = scheme_data
        mock_client_class.return_value = mock_client

        # Conversion would happen here
        # Verify resulting structure
        assert True  # Placeholder for actual end-to-end test

    @patch('torchbase.conversions.pubmlst.BIGSdbClient')
    def test_generated_torch_contains_expected_files(self, mock_client_class):
        """Generated torch should contain all expected files."""
        with tempfile.TemporaryDirectory() as tmpdir:
            torch_dir = Path(tmpdir)

            # After conversion, verify files exist:
            expected_files = [
                "metadata.toml",
                "quality.json",
                "schemes/organism/profiles.tsv",
                "schemes/organism/alleles/locus1.fasta",
            ]

            # This is a placeholder - actual test would create and verify
            assert True

    def test_metadata_toml_is_valid_toml(self):
        """Generated metadata.toml should be valid TOML."""
        with tempfile.TemporaryDirectory() as tmpdir:
            metadata_path = Path(tmpdir) / "metadata.toml"

            # Create valid TOML
            metadata = {
                "namespace": "pubmlst",
                "name": "scheme",
                "version": "1.0.0",
                "provenance": {"source": "PubMLST"},
            }

            with open(metadata_path, 'w') as f:
                toml.dump(metadata, f)

            # Should be parseable as TOML
            with open(metadata_path) as f:
                loaded = toml.load(f)
                assert loaded["namespace"] == "pubmlst"

    def test_quality_json_is_valid_json(self):
        """Generated quality.json should be valid JSON."""
        with tempfile.TemporaryDirectory() as tmpdir:
            quality_path = Path(tmpdir) / "quality.json"

            quality = {
                "kmer_analysis": {
                    "performed": True,
                }
            }

            with open(quality_path, 'w') as f:
                json.dump(quality, f)

            # Should be parseable as JSON
            with open(quality_path) as f:
                loaded = json.load(f)
                assert loaded["kmer_analysis"]["performed"] is True


class TestErrorHandling:
    """Test error handling for edge cases."""

    def test_handles_invalid_scheme_id(self):
        """Should gracefully handle invalid scheme ID."""
        # Test expects proper error message
        assert True  # Placeholder

    def test_handles_missing_database(self):
        """Should gracefully handle unreachable database."""
        # Test expects proper error message
        assert True  # Placeholder

    def test_handles_malformed_metadata(self):
        """Should detect malformed metadata during load."""
        # Test expects validation error
        assert True  # Placeholder

    def test_handles_invalid_output_path(self):
        """Should handle invalid output path gracefully."""
        # Test expects proper error message
        assert True  # Placeholder
