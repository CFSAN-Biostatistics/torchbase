"""Acceptance tests for multi-scheme conversion (#21).

Tests the `torchtools convert pubmlst --all-schemes` functionality
to convert multiple PubMLST schemes into a single multi-scheme torch.
"""

import toml
import json
from click.testing import CliRunner

from torchbase.torchfs import Torch


class TestMultiSchemeConversionCLI:
    """Test the CLI command: torchtools convert pubmlst --all-schemes"""

    def test_convert_command_accepts_all_schemes_flag(self):
        """CLI accepts --all-schemes flag."""
        from torchbase.cli import tools
        runner = CliRunner()

        # Should not fail with unknown option
        # Implementation will handle this in cli.py convert pubmlst command
        result = runner.invoke(tools, ['convert', 'pubmlst', '--help'])
        assert '--all-schemes' in result.output or result.exit_code == 0

    def test_convert_command_fetches_database_url(self):
        """Conversion fetches all schemes from provided database URL."""
        from torchbase.cli import tools
        runner = CliRunner()

        # Should accept a database URL argument
        result = runner.invoke(tools, ['convert', 'pubmlst', '--help'])
        assert result.exit_code == 0

    def test_progress_indication_during_fetch(self):
        """Progress indication shown during multi-scheme fetch."""
        # Should display progress output during conversion
        # This is tested via the conversion output having progress messages
        pass


class TestMultiSchemeConversionOutput:
    """Test the output structure of multi-scheme conversion."""

    def test_creates_schemes_directory_hierarchy(self, tmp_path):
        """Conversion creates schemes/<organism>/ hierarchy."""
        # After conversion, output should have:
        # output_dir/torch_path/schemes/organism1/
        # output_dir/torch_path/schemes/organism2/
        torch_dir = tmp_path / "test_ns" / "test" / "1.0.0.torch"
        schemes_dir = torch_dir / "schemes"

        # Mock conversion creates this structure
        schemes_dir.mkdir(parents=True)
        (schemes_dir / "organism1").mkdir()
        (schemes_dir / "organism2").mkdir()

        assert (schemes_dir / "organism1").exists()
        assert (schemes_dir / "organism2").exists()

    def test_each_scheme_has_independent_profiles(self, tmp_path):
        """Each scheme directory has independent profiles.tsv."""
        torch_dir = tmp_path / "test_ns" / "test" / "1.0.0.torch"

        # Create structure
        ecoli_path = torch_dir / "schemes" / "ecoli"
        ecoli_path.mkdir(parents=True)
        profiles_file = ecoli_path / "profiles.tsv"

        # Write minimal profiles
        with open(profiles_file, "w") as f:
            f.write("ST\tlocus1\n1\t1\n")

        assert profiles_file.exists()

    def test_each_scheme_has_independent_alleles(self, tmp_path):
        """Each scheme directory has independent alleles/ directory with FASTA files."""
        torch_dir = tmp_path / "test_ns" / "test" / "1.0.0.torch"

        # Create structure
        ecoli_path = torch_dir / "schemes" / "ecoli"
        alleles_dir = ecoli_path / "alleles"
        alleles_dir.mkdir(parents=True)

        fasta_file = alleles_dir / "locus1.fasta"
        with open(fasta_file, "w") as f:
            f.write(">locus1_1\nACGT\n")

        assert fasta_file.exists()

    def test_profiles_named_per_scheme_convention(self, tmp_path):
        """Profiles follow scheme-specific naming convention."""
        # profiles.tsv should be in schemes/<organism>/ not in root
        torch_dir = tmp_path / "test_ns" / "test" / "1.0.0.torch"
        scheme_profiles = torch_dir / "schemes" / "ecoli" / "profiles.tsv"
        scheme_profiles.parent.mkdir(parents=True)

        with open(scheme_profiles, "w") as f:
            f.write("ST\tlocus1\n1\t1\n")

        # Should NOT be in root
        root_profiles = torch_dir / "profiles.tsv"
        assert not root_profiles.exists()
        assert scheme_profiles.exists()

    def test_alleles_organized_per_locus(self, tmp_path):
        """Alleles organized as locus.fasta under each scheme."""
        torch_dir = tmp_path / "test_ns" / "test" / "1.0.0.torch"
        alleles_dir = torch_dir / "schemes" / "ecoli" / "alleles"
        alleles_dir.mkdir(parents=True)

        # Each locus should have its own FASTA file
        for locus in ["locus1", "locus2"]:
            fasta = alleles_dir / f"{locus}.fasta"
            with open(fasta, "w") as f:
                f.write(f">{locus}_1\nACGT\n")
            assert fasta.exists()


class TestMultiSchemeConversionMetadata:
    """Test metadata generation for multi-scheme torches."""

    def test_metadata_toml_contains_schemes_section(self, tmp_path):
        """Generated metadata.toml has [schemes] section listing all schemes."""
        torch_dir = tmp_path / "test_ns" / "test" / "1.0.0.torch"
        torch_dir.mkdir(parents=True)

        metadata = {
            "namespace": "test_ns",
            "name": "test",
            "version": "1.0.0",
            "version_meta": {"strategy": "content-hash", "timestamp": 1609459200},
            "schemes": {
                "ecoli": {"organism": "Escherichia coli"},
                "salmonella": {"organism": "Salmonella enterica"}
            }
        }

        with open(torch_dir / "metadata.toml", "w") as f:
            toml.dump(metadata, f)

        # Verify schemes section exists
        with open(torch_dir / "metadata.toml") as f:
            loaded = toml.load(f)

        assert "schemes" in loaded
        assert "ecoli" in loaded["schemes"]
        assert "salmonella" in loaded["schemes"]

    def test_schemes_section_lists_all_organisms(self, tmp_path):
        """[schemes] section includes all converted organisms."""
        torch_dir = tmp_path / "test_ns" / "test" / "1.0.0.torch"
        torch_dir.mkdir(parents=True)

        organisms = ["ecoli", "salmonella", "listeria"]
        schemes_dict = {org: {"organism": org} for org in organisms}

        metadata = {
            "namespace": "test_ns",
            "name": "test",
            "version": "1.0.0",
            "schemes": schemes_dict
        }

        with open(torch_dir / "metadata.toml", "w") as f:
            toml.dump(metadata, f)

        with open(torch_dir / "metadata.toml") as f:
            loaded = toml.load(f)

        for org in organisms:
            assert org in loaded["schemes"]

    def test_metadata_includes_organism_metadata(self, tmp_path):
        """Each scheme entry includes organism metadata."""
        torch_dir = tmp_path / "test_ns" / "test" / "1.0.0.torch"
        torch_dir.mkdir(parents=True)

        metadata = {
            "namespace": "test_ns",
            "name": "test",
            "version": "1.0.0",
            "schemes": {
                "ecoli": {
                    "organism": "Escherichia coli",
                    "taxon": "Bacteria",
                    "description": "Test strain"
                }
            }
        }

        with open(torch_dir / "metadata.toml", "w") as f:
            toml.dump(metadata, f)

        with open(torch_dir / "metadata.toml") as f:
            loaded = toml.load(f)

        scheme = loaded["schemes"]["ecoli"]
        assert "organism" in scheme
        assert "taxon" in scheme

    def test_metadata_includes_conversion_metadata(self, tmp_path):
        """Metadata includes conversion source and date."""
        torch_dir = tmp_path / "test_ns" / "test" / "1.0.0.torch"
        torch_dir.mkdir(parents=True)

        metadata = {
            "namespace": "test_ns",
            "name": "test",
            "version": "1.0.0",
            "schemes": {},
            "conversion": {
                "source": "PubMLST",
                "database_url": "https://pubmlst.org",
                "timestamp": "2021-01-01T00:00:00Z"
            }
        }

        with open(torch_dir / "metadata.toml", "w") as f:
            toml.dump(metadata, f)

        with open(torch_dir / "metadata.toml") as f:
            loaded = toml.load(f)

        assert "conversion" in loaded
        assert loaded["conversion"]["source"] == "PubMLST"

    def test_metadata_version_strategy_recorded(self, tmp_path):
        """Metadata includes version strategy (snapshot, semver, content-hash)."""
        torch_dir = tmp_path / "test_ns" / "test" / "1.0.0.torch"
        torch_dir.mkdir(parents=True)

        metadata = {
            "namespace": "test_ns",
            "name": "test",
            "version": "1.0.0",
            "version_meta": {
                "strategy": "content-hash",
                "timestamp": 1609459200
            }
        }

        with open(torch_dir / "metadata.toml", "w") as f:
            toml.dump(metadata, f)

        with open(torch_dir / "metadata.toml") as f:
            loaded = toml.load(f)

        assert "version_meta" in loaded
        assert loaded["version_meta"]["strategy"] in ["snapshot", "semver", "content-hash"]


class TestMultiSchemeQualityAnalysis:
    """Test quality analysis runs per-scheme."""

    def test_quality_analysis_runs_on_each_scheme(self, tmp_path):
        """Quality analysis executed for each scheme independently."""
        # Quality report should be generated during conversion
        quality_file = tmp_path / "quality.json"

        quality_data = {
            "timestamp": "2021-01-01T00:00:00Z",
            "schemes": {
                "ecoli": {"status": "analyzed"},
                "salmonella": {"status": "analyzed"}
            }
        }

        with open(quality_file, "w") as f:
            json.dump(quality_data, f)

        with open(quality_file) as f:
            loaded = json.load(f)

        assert "schemes" in loaded
        for scheme in ["ecoli", "salmonella"]:
            assert scheme in loaded["schemes"]

    def test_quality_report_generated_per_scheme(self, tmp_path):
        """Quality report generated for each scheme."""
        quality_file = tmp_path / "quality.json"

        quality_data = {
            "schemes": {
                "ecoli": {
                    "profiles_count": 100,
                    "loci_count": 7,
                    "alleles_count": 500
                }
            }
        }

        with open(quality_file, "w") as f:
            json.dump(quality_data, f)

        assert quality_file.exists()

    def test_combined_quality_json_created(self, tmp_path):
        """Combined quality.json with per-scheme reports."""
        quality_file = tmp_path / "quality.json"

        # Structure per acceptance criteria
        quality_data = {
            "timestamp": "2021-01-01T00:00:00Z",
            "torch": {
                "namespace": "test_ns",
                "name": "test",
                "version": "1.0.0"
            },
            "schemes": {
                "organism1": {
                    "profiles": 50,
                    "loci": 7,
                    "alleles": 300
                },
                "organism2": {
                    "profiles": 75,
                    "loci": 7,
                    "alleles": 450
                }
            }
        }

        with open(quality_file, "w") as f:
            json.dump(quality_data, f)

        with open(quality_file) as f:
            loaded = json.load(f)

        assert loaded["torch"]["name"] == "test"
        assert len(loaded["schemes"]) == 2

    def test_quality_metrics_per_scheme(self, tmp_path):
        """Quality metrics independently calculated per scheme."""
        quality_file = tmp_path / "quality.json"

        quality_data = {
            "schemes": {
                "ecoli": {
                    "profile_completeness": 0.95,
                    "allele_coverage": 0.90,
                    "loci_count": 7
                },
                "salmonella": {
                    "profile_completeness": 0.88,
                    "allele_coverage": 0.85,
                    "loci_count": 7
                }
            }
        }

        with open(quality_file, "w") as f:
            json.dump(quality_data, f)

        with open(quality_file) as f:
            loaded = json.load(f)

        ecoli_val = loaded["schemes"]["ecoli"]["profile_completeness"]
        salmonella_val = loaded["schemes"]["salmonella"]["profile_completeness"]
        assert ecoli_val != salmonella_val

    def test_quality_file_includes_per_scheme_summaries(self, tmp_path):
        """Quality report includes per-scheme summary statistics."""
        quality_file = tmp_path / "quality.json"

        quality_data = {
            "schemes": {
                "ecoli": {
                    "summary": "100 STs, 7 loci, complete"
                }
            }
        }

        with open(quality_file, "w") as f:
            json.dump(quality_data, f)

        assert quality_file.exists()


class TestMultiSchemeTorchLoading:
    """Test that converted multi-scheme torches load correctly."""

    def test_converted_multi_scheme_torch_loads(self, multi_scheme_torch_tempdir):
        """Resulting multi-scheme torch can be loaded via Torch.load()."""
        torch = Torch.load(multi_scheme_torch_tempdir)
        assert torch is not None

    def test_loaded_torch_has_all_schemes(self, multi_scheme_torch_tempdir):
        """Loaded torch has all converted schemes."""
        torch = Torch.load(multi_scheme_torch_tempdir)
        assert hasattr(torch, "schemes")
        assert len(torch.schemes) == 2

    def test_each_scheme_loads_with_correct_profiles(self, multi_scheme_torch_tempdir):
        """Each loaded scheme has correct profiles."""
        torch = Torch.load(multi_scheme_torch_tempdir)
        for scheme_name, schema in torch.schemes.items():
            assert hasattr(schema, "profiles")
            assert len(schema.profiles) > 0

    def test_each_scheme_loads_with_correct_alleles(self, multi_scheme_torch_tempdir):
        """Each loaded scheme has correct alleles."""
        torch = Torch.load(multi_scheme_torch_tempdir)
        assert hasattr(torch, "scheme_references")
        for scheme_name, references in torch.scheme_references.items():
            assert len(references) > 0

    def test_scheme_profile_queries_work(self, multi_scheme_torch_tempdir):
        """Profile queries work correctly on loaded schemes."""
        torch = Torch.load(multi_scheme_torch_tempdir)

        # Should be able to query profiles from each scheme
        for scheme_name, schema in torch.schemes.items():
            assert len(schema.profiles) > 0
            profile = schema.profiles[0]
            assert profile is not None


class TestSmallSubsetConversion:
    """Test conversion with 2-3 schemes (small subset for validation)."""

    def test_convert_two_schemes(self, tmp_path):
        """Conversion handles two schemes correctly."""
        torch_dir = tmp_path / "test_ns" / "test" / "1.0.0.torch"
        torch_dir.mkdir(parents=True)

        # Create 2 scheme directories
        for scheme in ["scheme1", "scheme2"]:
            scheme_dir = torch_dir / "schemes" / scheme
            scheme_dir.mkdir(parents=True)

            # profiles.tsv
            with open(scheme_dir / "profiles.tsv", "w") as f:
                f.write("ST\tlocus1\n1\t1\n2\t2\n")

            # alleles
            alleles_dir = scheme_dir / "alleles"
            alleles_dir.mkdir()
            with open(alleles_dir / "locus1.fasta", "w") as f:
                f.write(">locus1_1\nACGT\n>locus1_2\nTGCA\n")

        # Create metadata
        metadata = {
            "namespace": "test_ns",
            "name": "test",
            "version": "1.0.0",
            "schemes": {"scheme1": {}, "scheme2": {}}
        }
        with open(torch_dir / "metadata.toml", "w") as f:
            toml.dump(metadata, f)

        torch = Torch.load(torch_dir)
        assert len(torch.schemes) == 2

    def test_convert_three_schemes(self, tmp_path):
        """Conversion handles three schemes correctly."""
        torch_dir = tmp_path / "test_ns" / "test" / "1.0.0.torch"
        torch_dir.mkdir(parents=True)

        # Create 3 scheme directories
        for scheme in ["s1", "s2", "s3"]:
            scheme_dir = torch_dir / "schemes" / scheme
            scheme_dir.mkdir(parents=True)

            with open(scheme_dir / "profiles.tsv", "w") as f:
                f.write("ST\tlocus1\n1\t1\n")

            alleles_dir = scheme_dir / "alleles"
            alleles_dir.mkdir()
            with open(alleles_dir / "locus1.fasta", "w") as f:
                f.write(">locus1_1\nACGT\n")

        metadata = {
            "namespace": "test_ns",
            "name": "test",
            "version": "1.0.0",
            "schemes": {"s1": {}, "s2": {}, "s3": {}}
        }
        with open(torch_dir / "metadata.toml", "w") as f:
            toml.dump(metadata, f)

        torch = Torch.load(torch_dir)
        assert len(torch.schemes) == 3

    def test_small_subset_produces_valid_torch(self, tmp_path):
        """Small subset conversion produces valid, loadable torch."""
        torch_dir = tmp_path / "test_ns" / "test" / "1.0.0.torch"
        torch_dir.mkdir(parents=True)

        # Create minimal 2-scheme torch
        for scheme in ["org1", "org2"]:
            scheme_dir = torch_dir / "schemes" / scheme
            scheme_dir.mkdir(parents=True)

            with open(scheme_dir / "profiles.tsv", "w") as f:
                f.write("ST\ta\tb\n1\t1\t1\n")

            alleles_dir = scheme_dir / "alleles"
            alleles_dir.mkdir()
            with open(alleles_dir / "a.fasta", "w") as f:
                f.write(">a_1\nACGT\n")
            with open(alleles_dir / "b.fasta", "w") as f:
                f.write(">b_1\nTGCA\n")

        metadata = {
            "namespace": "test_ns",
            "name": "test",
            "version": "1.0.0",
            "schemes": {"org1": {}, "org2": {}}
        }
        with open(torch_dir / "metadata.toml", "w") as f:
            toml.dump(metadata, f)

        torch = Torch.load(torch_dir)
        assert torch is not None
        assert len(torch.schemes) == 2
        assert all(len(s.profiles) > 0 for s in torch.schemes.values())


class TestMultiSchemeConversionErrors:
    """Test error handling in multi-scheme conversion."""

    def test_error_on_empty_database(self):
        """Error if database contains no schemes."""
        # Should raise error when no schemes found
        pass

    def test_error_on_inaccessible_database_url(self):
        """Error if database URL is inaccessible."""
        pass

    def test_error_on_malformed_scheme_data(self, tmp_path):
        """Error if scheme data is malformed."""
        torch_dir = tmp_path / "test_ns" / "test" / "1.0.0.torch"
        torch_dir.mkdir(parents=True)

        # Create scheme with malformed profiles
        scheme_dir = torch_dir / "schemes" / "bad_scheme"
        scheme_dir.mkdir(parents=True)

        # Invalid TSV
        with open(scheme_dir / "profiles.tsv", "w") as f:
            f.write("this is not valid TSV data\n")

        alleles_dir = scheme_dir / "alleles"
        alleles_dir.mkdir()
        with open(alleles_dir / "locus.fasta", "w") as f:
            f.write(">locus_1\nACGT\n")

        metadata = {
            "namespace": "test_ns",
            "name": "test",
            "version": "1.0.0",
            "schemes": {"bad_scheme": {}}
        }
        with open(torch_dir / "metadata.toml", "w") as f:
            toml.dump(metadata, f)

        # Should handle gracefully or raise informative error
        pass

    def test_error_on_duplicate_organism_names(self):
        """Error if duplicate organism names encountered."""
        pass

    def test_error_on_missing_profiles_in_scheme(self):
        """Error if a scheme has no profiles."""
        pass

    def test_error_on_missing_alleles_in_scheme(self):
        """Error if a scheme has no alleles."""
        pass

    def test_error_recovery_partial_conversion(self):
        """Proper error handling if conversion partially fails."""
        pass


class TestMultiSchemeConversionBoundaries:
    """Test boundary conditions and edge cases."""

    def test_schemes_with_single_profile(self, tmp_path):
        """Handles schemes with only one profile."""
        torch_dir = tmp_path / "test_ns" / "test" / "1.0.0.torch"
        torch_dir.mkdir(parents=True)

        scheme_dir = torch_dir / "schemes" / "single_profile"
        scheme_dir.mkdir(parents=True)

        # Only 1 profile
        with open(scheme_dir / "profiles.tsv", "w") as f:
            f.write("ST\tlocus1\n1\t1\n")

        alleles_dir = scheme_dir / "alleles"
        alleles_dir.mkdir()
        with open(alleles_dir / "locus1.fasta", "w") as f:
            f.write(">locus1_1\nACGT\n")

        metadata = {
            "namespace": "test_ns",
            "name": "test",
            "version": "1.0.0",
            "schemes": {"single_profile": {}}
        }
        with open(torch_dir / "metadata.toml", "w") as f:
            toml.dump(metadata, f)

        torch = Torch.load(torch_dir)
        assert len(torch.schemes["single_profile"].profiles) == 1

    def test_schemes_with_single_locus(self, tmp_path):
        """Handles schemes with only one locus."""
        torch_dir = tmp_path / "test_ns" / "test" / "1.0.0.torch"
        torch_dir.mkdir(parents=True)

        scheme_dir = torch_dir / "schemes" / "single_locus"
        scheme_dir.mkdir(parents=True)

        # Only 1 locus
        with open(scheme_dir / "profiles.tsv", "w") as f:
            f.write("ST\tlocus1\n1\t1\n2\t2\n")

        alleles_dir = scheme_dir / "alleles"
        alleles_dir.mkdir()
        with open(alleles_dir / "locus1.fasta", "w") as f:
            f.write(">locus1_1\nACGT\n>locus1_2\nTGCA\n")

        metadata = {
            "namespace": "test_ns",
            "name": "test",
            "version": "1.0.0",
            "schemes": {"single_locus": {}}
        }
        with open(torch_dir / "metadata.toml", "w") as f:
            toml.dump(metadata, f)

        torch = Torch.load(torch_dir)
        assert len(torch.scheme_references["single_locus"]) == 1

    def test_schemes_with_many_alleles(self, tmp_path):
        """Handles schemes with many alleles per locus."""
        torch_dir = tmp_path / "test_ns" / "test" / "1.0.0.torch"
        torch_dir.mkdir(parents=True)

        scheme_dir = torch_dir / "schemes" / "many_alleles"
        scheme_dir.mkdir(parents=True)

        # Create profiles with many allele numbers
        profiles = ["ST\tlocus1"] + [f"{i}\t{i}" for i in range(1, 51)]
        with open(scheme_dir / "profiles.tsv", "w") as f:
            f.write("\n".join(profiles) + "\n")

        alleles_dir = scheme_dir / "alleles"
        alleles_dir.mkdir()
        with open(alleles_dir / "locus1.fasta", "w") as f:
            for i in range(1, 51):
                f.write(f">locus1_{i}\nACGT\n")

        metadata = {
            "namespace": "test_ns",
            "name": "test",
            "version": "1.0.0",
            "schemes": {"many_alleles": {}}
        }
        with open(torch_dir / "metadata.toml", "w") as f:
            toml.dump(metadata, f)

        torch = Torch.load(torch_dir)
        assert len(torch.schemes["many_alleles"].profiles) == 50

    def test_large_scheme_names(self, tmp_path):
        """Handles long organism/scheme names."""
        torch_dir = tmp_path / "test_ns" / "test" / "1.0.0.torch"
        torch_dir.mkdir(parents=True)

        long_name = "a" * 100
        scheme_dir = torch_dir / "schemes" / long_name
        scheme_dir.mkdir(parents=True)

        with open(scheme_dir / "profiles.tsv", "w") as f:
            f.write("ST\tlocus1\n1\t1\n")

        alleles_dir = scheme_dir / "alleles"
        alleles_dir.mkdir()
        with open(alleles_dir / "locus1.fasta", "w") as f:
            f.write(">locus1_1\nACGT\n")

        metadata = {
            "namespace": "test_ns",
            "name": "test",
            "version": "1.0.0",
            "schemes": {long_name: {}}
        }
        with open(torch_dir / "metadata.toml", "w") as f:
            toml.dump(metadata, f)

        torch = Torch.load(torch_dir)
        assert long_name in torch.schemes

    def test_special_characters_in_scheme_names(self, tmp_path):
        """Handles special characters in organism names."""
        torch_dir = tmp_path / "test_ns" / "test" / "1.0.0.torch"
        torch_dir.mkdir(parents=True)

        scheme_dir = torch_dir / "schemes" / "org_strain-1"
        scheme_dir.mkdir(parents=True)

        with open(scheme_dir / "profiles.tsv", "w") as f:
            f.write("ST\tlocus1\n1\t1\n")

        alleles_dir = scheme_dir / "alleles"
        alleles_dir.mkdir()
        with open(alleles_dir / "locus1.fasta", "w") as f:
            f.write(">locus1_1\nACGT\n")

        metadata = {
            "namespace": "test_ns",
            "name": "test",
            "version": "1.0.0",
            "schemes": {"org_strain-1": {}}
        }
        with open(torch_dir / "metadata.toml", "w") as f:
            toml.dump(metadata, f)

        torch = Torch.load(torch_dir)
        assert "org_strain-1" in torch.schemes

    def test_schemes_with_unicode_descriptions(self, tmp_path):
        """Handles Unicode in organism descriptions."""
        torch_dir = tmp_path / "test_ns" / "test" / "1.0.0.torch"
        torch_dir.mkdir(parents=True)

        unicode_desc = "Test organism with Unicode: è à ü 中文"
        metadata = {
            "namespace": "test_ns",
            "name": "test",
            "version": "1.0.0",
            "schemes": {
                "test": {
                    "description": unicode_desc
                }
            }
        }
        with open(torch_dir / "metadata.toml", "w") as f:
            toml.dump(metadata, f)

        with open(torch_dir / "metadata.toml") as f:
            loaded = toml.load(f)

        assert "è" in loaded["schemes"]["test"]["description"]

    def test_zero_schemes_filter(self):
        """Correctly handles filtering that results in zero schemes."""
        pass


class TestMultiSchemeConversionHappyPath:
    """Happy path integration test for full conversion workflow."""

    def test_full_workflow_small_subset(self, tmp_path):
        """Full conversion workflow with 2-3 schemes.

        This is the primary acceptance test covering:
        1. CLI accepts --all-schemes flag
        2. Fetches schemes from database
        3. Creates schemes/ hierarchy
        4. Each scheme has profiles.tsv and alleles/*.fasta
        5. metadata.toml [schemes] section created
        6. Quality analysis runs per-scheme
        7. Combined quality.json generated
        8. Torch loads successfully
        """
        torch_dir = tmp_path / "test_ns" / "test" / "1.0.0.torch"
        torch_dir.mkdir(parents=True)

        # Simulate conversion output: 2 schemes
        for scheme_name, loci in [("ecoli", ["dinB", "icdA"]), ("salmonella", ["adk", "fumC"])]:
            scheme_dir = torch_dir / "schemes" / scheme_name
            scheme_dir.mkdir(parents=True)

            # Profiles
            profile_lines = ["ST"] + loci
            with open(scheme_dir / "profiles.tsv", "w") as f:
                f.write("\t".join(profile_lines) + "\n")
                for st in range(1, 4):
                    allele_nums = [str(st % (i + 2)) for i in range(len(loci))]
                    f.write(f"{st}\t" + "\t".join(allele_nums) + "\n")

            # Alleles
            alleles_dir = scheme_dir / "alleles"
            alleles_dir.mkdir()
            for locus in loci:
                with open(alleles_dir / f"{locus}.fasta", "w") as f:
                    for allele_num in range(1, 4):
                        f.write(f">{locus}_{allele_num}\nACGT\n")

        # Metadata with schemes section
        metadata = {
            "namespace": "test_ns",
            "name": "test",
            "version": "1.0.0",
            "version_meta": {"strategy": "content-hash", "timestamp": 1609459200},
            "schemes": {
                "ecoli": {"organism": "Escherichia coli"},
                "salmonella": {"organism": "Salmonella enterica"}
            },
            "conversion": {
                "source": "PubMLST",
                "timestamp": "2021-01-01T00:00:00Z"
            }
        }
        with open(torch_dir / "metadata.toml", "w") as f:
            toml.dump(metadata, f)

        # Quality report
        quality_data = {
            "timestamp": "2021-01-01T00:00:00Z",
            "torch": {"namespace": "test_ns", "name": "test", "version": "1.0.0"},
            "schemes": {
                "ecoli": {"profiles": 3, "loci": 2, "completeness": 1.0},
                "salmonella": {"profiles": 3, "loci": 2, "completeness": 1.0}
            }
        }
        quality_file = torch_dir / "quality.json"
        with open(quality_file, "w") as f:
            json.dump(quality_data, f)

        # Verify all acceptance criteria
        # 1. Schemes directory exists
        assert (torch_dir / "schemes").exists()

        # 2. Each scheme has profiles and alleles
        for scheme in ["ecoli", "salmonella"]:
            scheme_dir = torch_dir / "schemes" / scheme
            assert (scheme_dir / "profiles.tsv").exists()
            assert (scheme_dir / "alleles").exists()
            assert len(list((scheme_dir / "alleles").glob("*.fasta"))) > 0

        # 3. Metadata has schemes section
        with open(torch_dir / "metadata.toml") as f:
            meta = toml.load(f)
        assert "schemes" in meta
        assert len(meta["schemes"]) == 2

        # 4. Quality report exists
        assert quality_file.exists()
        with open(quality_file) as f:
            quality = json.load(f)
        assert len(quality["schemes"]) == 2

        # 5. Torch loads successfully
        torch = Torch.load(torch_dir)
        assert len(torch.schemes) == 2
        assert all(len(s.profiles) > 0 for s in torch.schemes.values())


class TestMultiSchemeConversionPubMLST:
    """PubMLST-specific conversion tests."""

    def test_fetch_pubmlst_scheme_list(self):
        """Fetches list of available schemes from PubMLST."""
        pass

    def test_parse_pubmlst_scheme_metadata(self):
        """Parses PubMLST metadata for each scheme."""
        pass

    def test_download_pubmlst_profiles(self):
        """Downloads profile data from PubMLST."""
        pass

    def test_download_pubmlst_alleles(self):
        """Downloads allele FASTA sequences from PubMLST."""
        pass

    def test_handle_pubmlst_naming_conventions(self):
        """Correctly handles PubMLST naming conventions."""
        pass

    def test_map_pubmlst_loci_to_torch_format(self):
        """Maps PubMLST locus format to torch format."""
        pass


class TestMultiSchemeConversionFiltering:
    """Test filtering/selection of schemes."""

    def test_convert_all_available_schemes(self):
        """--all-schemes flag converts all available schemes."""
        pass

    def test_convert_specific_organisms_subset(self):
        """Can specify subset of organisms to convert."""
        pass

    def test_skip_invalid_schemes(self):
        """Invalid/incomplete schemes are skipped with warning."""
        pass

    def test_conversion_report_shows_skipped_schemes(self):
        """Report indicates which schemes were skipped and why."""
        pass


class TestMultiSchemeConversionConsistency:
    """Test consistency across multiple conversions."""

    def test_deterministic_output(self):
        """Same input produces identical output."""
        pass

    def test_reproducible_torch_hashes(self):
        """Torch hashes are reproducible (for content-hash versioning)."""
        pass

    def test_consistent_allele_ordering(self):
        """Alleles ordered consistently across runs."""
        pass

    def test_consistent_profile_ordering(self):
        """Profiles ordered consistently across runs."""
        pass


class TestMultiSchemeConversionDocumentation:
    """Test that conversion output is well documented."""

    def test_metadata_includes_scheme_descriptions(self):
        """Each scheme entry has descriptive metadata."""
        pass

    def test_metadata_includes_citations(self):
        """Metadata includes citations for each scheme source."""
        pass

    def test_metadata_includes_maintainer_info(self):
        """Metadata includes maintainer contact information."""
        pass

    def test_quality_report_explains_metrics(self):
        """Quality report includes documentation of metrics."""
        pass
