"""Tests for synthetic example torches (#61).

Tests that example torches exist, are correctly structured, and can be loaded
and used for testing and documentation purposes.
"""

from pathlib import Path
from torchbase.torchfs import Torch
from torchbase.torchbase import Schema
import toml
import csv


class TestSimpleMLSTExampleStructure:
    """Test the simple_mlst example torch exists and has correct structure."""

    def test_simple_mlst_directory_exists(self):
        """examples/simple_mlst directory exists."""
        examples_dir = Path(__file__).parent.parent.parent / "examples"
        simple_mlst = examples_dir / "simple_mlst"
        assert simple_mlst.exists(), "examples/simple_mlst directory not found"
        assert simple_mlst.is_dir(), "examples/simple_mlst must be a directory"

    def test_simple_mlst_has_torch_structure(self):
        """examples/simple_mlst has proper torch directory structure."""
        examples_dir = Path(__file__).parent.parent.parent / "examples"
        # Should be: examples/simple_mlst/<version>.torch/
        simple_mlst = examples_dir / "simple_mlst"

        # Find the torch directory (should end with .torch)
        torch_dirs = list(simple_mlst.glob("*.torch"))
        assert len(torch_dirs) == 1, (
            "examples/simple_mlst should have exactly one .torch directory"
        )

        torch_dir = torch_dirs[0]
        assert torch_dir.is_dir(), "Torch path must be a directory"

    def test_simple_mlst_has_metadata_toml(self):
        """examples/simple_mlst torch has metadata.toml."""
        examples_dir = Path(__file__).parent.parent.parent / "examples"
        simple_mlst = examples_dir / "simple_mlst"
        torch_dir = list(simple_mlst.glob("*.torch"))[0]

        metadata_file = torch_dir / "metadata.toml"
        assert metadata_file.exists(), "metadata.toml not found"
        assert metadata_file.is_file(), "metadata.toml must be a file"

    def test_simple_mlst_has_profiles_tsv(self):
        """examples/simple_mlst torch has profiles.tsv."""
        examples_dir = Path(__file__).parent.parent.parent / "examples"
        simple_mlst = examples_dir / "simple_mlst"
        torch_dir = list(simple_mlst.glob("*.torch"))[0]

        profiles_file = torch_dir / "profiles.tsv"
        assert profiles_file.exists(), "profiles.tsv not found"
        assert profiles_file.is_file(), "profiles.tsv must be a file"

    def test_simple_mlst_has_resources_directory(self):
        """examples/simple_mlst torch has _resources/ directory."""
        examples_dir = Path(__file__).parent.parent.parent / "examples"
        simple_mlst = examples_dir / "simple_mlst"
        torch_dir = list(simple_mlst.glob("*.torch"))[0]

        resources_dir = torch_dir / "_resources"
        assert resources_dir.exists(), "_resources directory not found"
        assert resources_dir.is_dir(), "_resources must be a directory"

    def test_simple_mlst_not_in_schemes_format(self):
        """examples/simple_mlst uses flat structure, not schemes/."""
        examples_dir = Path(__file__).parent.parent.parent / "examples"
        simple_mlst = examples_dir / "simple_mlst"
        torch_dir = list(simple_mlst.glob("*.torch"))[0]

        schemes_dir = torch_dir / "schemes"
        assert not schemes_dir.exists(), (
            "simple_mlst should use flat structure, not schemes/"
        )


class TestSimpleMLSTExampleContent:
    """Test the content of the simple_mlst example torch."""

    def test_simple_mlst_has_three_loci(self):
        """simple_mlst has exactly 3 loci (adk, fumC, gyrB)."""
        examples_dir = Path(__file__).parent.parent.parent / "examples"
        simple_mlst = examples_dir / "simple_mlst"
        torch_dir = list(simple_mlst.glob("*.torch"))[0]

        with open(torch_dir / "profiles.tsv") as f:
            reader = csv.reader(f, delimiter="\t")
            header = next(reader)
            # First column is ST, rest are loci
            loci = header[1:]

        assert len(loci) == 3, "simple_mlst should have exactly 3 loci"
        assert set(loci) == {"adk", "fumC", "gyrB"}, (
            "Loci should be adk, fumC, gyrB"
        )

    def test_simple_mlst_has_three_alleles_per_locus(self):
        """Each locus in simple_mlst has 3 alleles."""
        examples_dir = Path(__file__).parent.parent.parent / "examples"
        simple_mlst = examples_dir / "simple_mlst"
        torch_dir = list(simple_mlst.glob("*.torch"))[0]
        resources_dir = torch_dir / "_resources"

        for locus in ["adk", "fumC", "gyrB"]:
            fasta_file = resources_dir / f"{locus}.fasta"
            assert fasta_file.exists(), f"{locus}.fasta not found"

            # Count alleles (>header lines)
            with open(fasta_file) as f:
                allele_count = sum(1 for line in f if line.startswith(">"))

            assert allele_count == 3, (
                f"{locus} should have exactly 3 alleles, found {allele_count}"
            )

    def test_simple_mlst_has_five_profiles(self):
        """simple_mlst has exactly 5 profiles."""
        examples_dir = Path(__file__).parent.parent.parent / "examples"
        simple_mlst = examples_dir / "simple_mlst"
        torch_dir = list(simple_mlst.glob("*.torch"))[0]

        with open(torch_dir / "profiles.tsv") as f:
            reader = csv.reader(f, delimiter="\t")
            next(reader)  # skip header
            profiles = list(reader)

        assert len(profiles) == 5, (
            f"simple_mlst should have exactly 5 profiles, found {len(profiles)}"
        )

    def test_simple_mlst_metadata_valid(self):
        """simple_mlst metadata.toml is valid and minimal."""
        examples_dir = Path(__file__).parent.parent.parent / "examples"
        simple_mlst = examples_dir / "simple_mlst"
        torch_dir = list(simple_mlst.glob("*.torch"))[0]

        with open(torch_dir / "metadata.toml") as f:
            metadata = toml.load(f)

        # Required fields
        assert "namespace" in metadata
        assert "name" in metadata
        assert "version" in metadata
        assert "manifest" in metadata
        assert metadata["namespace"] == "examples"
        assert metadata["name"] == "simple_mlst"


class TestSimpleMLSTExampleLoadable:
    """Test that simple_mlst torch can be loaded via Torch.load()."""

    def test_simple_mlst_loads_via_torch_load(self):
        """Torch.load() successfully loads simple_mlst example."""
        examples_dir = Path(__file__).parent.parent.parent / "examples"
        simple_mlst = examples_dir / "simple_mlst"
        torch_dir = list(simple_mlst.glob("*.torch"))[0]

        # Should not raise
        torch = Torch.load(torch_dir)
        assert torch is not None

    def test_simple_mlst_loads_as_single_scheme_format(self):
        """simple_mlst loads as single-scheme (legacy) format."""
        examples_dir = Path(__file__).parent.parent.parent / "examples"
        simple_mlst = examples_dir / "simple_mlst"
        torch_dir = list(simple_mlst.glob("*.torch"))[0]

        torch = Torch.load(torch_dir)

        # Single-scheme format has profile attribute
        assert hasattr(torch, "profile")
        assert torch.profile is not None
        assert isinstance(torch.profile, Schema)

    def test_simple_mlst_has_correct_profile_count(self):
        """Loaded simple_mlst has 5 profiles."""
        examples_dir = Path(__file__).parent.parent.parent / "examples"
        simple_mlst = examples_dir / "simple_mlst"
        torch_dir = list(simple_mlst.glob("*.torch"))[0]

        torch = Torch.load(torch_dir)
        assert len(torch.profile.profiles) == 5

    def test_simple_mlst_references_loaded(self):
        """simple_mlst allele references are loaded."""
        examples_dir = Path(__file__).parent.parent.parent / "examples"
        simple_mlst = examples_dir / "simple_mlst"
        torch_dir = list(simple_mlst.glob("*.torch"))[0]

        torch = Torch.load(torch_dir)

        assert hasattr(torch, "references")
        ref_names = {ref.name for ref in torch.references}
        assert "adk.fasta" in ref_names
        assert "fumC.fasta" in ref_names
        assert "gyrB.fasta" in ref_names


class TestMultiOrganismExampleStructure:
    """Test the multi_organism example torch exists and has correct structure."""

    def test_multi_organism_directory_exists(self):
        """examples/multi_organism directory exists."""
        examples_dir = Path(__file__).parent.parent.parent / "examples"
        multi_organism = examples_dir / "multi_organism"
        assert multi_organism.exists(), (
            "examples/multi_organism directory not found"
        )
        assert multi_organism.is_dir(), (
            "examples/multi_organism must be a directory"
        )

    def test_multi_organism_has_torch_structure(self):
        """examples/multi_organism has proper torch directory structure."""
        examples_dir = Path(__file__).parent.parent.parent / "examples"
        multi_organism = examples_dir / "multi_organism"

        torch_dirs = list(multi_organism.glob("*.torch"))
        assert len(torch_dirs) == 1, (
            "examples/multi_organism should have exactly one .torch directory"
        )

        torch_dir = torch_dirs[0]
        assert torch_dir.is_dir(), "Torch path must be a directory"

    def test_multi_organism_has_metadata_toml(self):
        """examples/multi_organism torch has metadata.toml."""
        examples_dir = Path(__file__).parent.parent.parent / "examples"
        multi_organism = examples_dir / "multi_organism"
        torch_dir = list(multi_organism.glob("*.torch"))[0]

        metadata_file = torch_dir / "metadata.toml"
        assert metadata_file.exists(), "metadata.toml not found"
        assert metadata_file.is_file(), "metadata.toml must be a file"

    def test_multi_organism_has_schemes_directory(self):
        """examples/multi_organism uses schemes/ format."""
        examples_dir = Path(__file__).parent.parent.parent / "examples"
        multi_organism = examples_dir / "multi_organism"
        torch_dir = list(multi_organism.glob("*.torch"))[0]

        schemes_dir = torch_dir / "schemes"
        assert schemes_dir.exists(), "schemes/ directory not found"
        assert schemes_dir.is_dir(), "schemes/ must be a directory"

    def test_multi_organism_has_two_schemes(self):
        """multi_organism has exactly 2 scheme subdirectories."""
        examples_dir = Path(__file__).parent.parent.parent / "examples"
        multi_organism = examples_dir / "multi_organism"
        torch_dir = list(multi_organism.glob("*.torch"))[0]
        schemes_dir = torch_dir / "schemes"

        # Count non-hidden subdirectories
        scheme_dirs = [
            d for d in schemes_dir.iterdir()
            if d.is_dir() and not d.name.startswith(".")
        ]

        assert len(scheme_dirs) == 2, (
            f"Should have exactly 2 schemes, found {len(scheme_dirs)}"
        )

    def test_multi_organism_has_salmonella_and_ecoli_schemes(self):
        """multi_organism has salmonella and ecoli schemes."""
        examples_dir = Path(__file__).parent.parent.parent / "examples"
        multi_organism = examples_dir / "multi_organism"
        torch_dir = list(multi_organism.glob("*.torch"))[0]
        schemes_dir = torch_dir / "schemes"

        scheme_names = {
            d.name for d in schemes_dir.iterdir()
            if d.is_dir() and not d.name.startswith(".")
        }

        assert scheme_names == {"salmonella", "ecoli"}, (
            f"Expected salmonella and ecoli, found {scheme_names}"
        )


class TestMultiOrganismSchemeContent:
    """Test the content of each scheme in multi_organism torch."""

    def test_salmonella_scheme_has_profiles_and_alleles(self):
        """salmonella scheme has profiles.tsv and alleles/ directory."""
        examples_dir = Path(__file__).parent.parent.parent / "examples"
        multi_organism = examples_dir / "multi_organism"
        torch_dir = list(multi_organism.glob("*.torch"))[0]
        salmonella = torch_dir / "schemes" / "salmonella"

        assert (salmonella / "profiles.tsv").exists()
        assert (salmonella / "alleles").is_dir()

    def test_ecoli_scheme_has_profiles_and_alleles(self):
        """ecoli scheme has profiles.tsv and alleles/ directory."""
        examples_dir = Path(__file__).parent.parent.parent / "examples"
        multi_organism = examples_dir / "multi_organism"
        torch_dir = list(multi_organism.glob("*.torch"))[0]
        ecoli = torch_dir / "schemes" / "ecoli"

        assert (ecoli / "profiles.tsv").exists()
        assert (ecoli / "alleles").is_dir()

    def test_each_scheme_has_three_loci(self):
        """Both salmonella and ecoli schemes have exactly 3 loci."""
        examples_dir = Path(__file__).parent.parent.parent / "examples"
        multi_organism = examples_dir / "multi_organism"
        torch_dir = list(multi_organism.glob("*.torch"))[0]

        for scheme_name in ["salmonella", "ecoli"]:
            scheme_path = torch_dir / "schemes" / scheme_name

            with open(scheme_path / "profiles.tsv") as f:
                reader = csv.reader(f, delimiter="\t")
                header = next(reader)
                loci = header[1:]  # First column is ST

            assert len(loci) == 3, (
                f"{scheme_name} should have 3 loci, found {len(loci)}"
            )

    def test_each_scheme_has_three_alleles_per_locus(self):
        """Each locus in each scheme has 3 alleles."""
        examples_dir = Path(__file__).parent.parent.parent / "examples"
        multi_organism = examples_dir / "multi_organism"
        torch_dir = list(multi_organism.glob("*.torch"))[0]

        for scheme_name in ["salmonella", "ecoli"]:
            scheme_path = torch_dir / "schemes" / scheme_name
            alleles_dir = scheme_path / "alleles"

            # Each scheme should have 3 FASTA files (one per locus)
            fasta_files = list(alleles_dir.glob("*.fasta"))
            assert len(fasta_files) == 3, (
                f"{scheme_name} should have 3 FASTA files"
            )

            for fasta_file in fasta_files:
                with open(fasta_file) as f:
                    allele_count = sum(1 for line in f if line.startswith(">"))

                assert allele_count == 3, (
                    f"{fasta_file.name} should have 3 alleles"
                )

    def test_each_scheme_has_three_profiles(self):
        """Both schemes have exactly 3 profiles."""
        examples_dir = Path(__file__).parent.parent.parent / "examples"
        multi_organism = examples_dir / "multi_organism"
        torch_dir = list(multi_organism.glob("*.torch"))[0]

        for scheme_name in ["salmonella", "ecoli"]:
            scheme_path = torch_dir / "schemes" / scheme_name

            with open(scheme_path / "profiles.tsv") as f:
                reader = csv.reader(f, delimiter="\t")
                next(reader)  # skip header
                profiles = list(reader)

            assert len(profiles) == 3, (
                f"{scheme_name} should have 3 profiles, found {len(profiles)}"
            )


class TestMultiOrganismExampleLoadable:
    """Test that multi_organism torch can be loaded via Torch.load()."""

    def test_multi_organism_loads_via_torch_load(self):
        """Torch.load() successfully loads multi_organism example."""
        examples_dir = Path(__file__).parent.parent.parent / "examples"
        multi_organism = examples_dir / "multi_organism"
        torch_dir = list(multi_organism.glob("*.torch"))[0]

        # Should not raise
        torch = Torch.load(torch_dir)
        assert torch is not None

    def test_multi_organism_loads_as_multi_scheme_format(self):
        """multi_organism loads as multi-scheme format."""
        examples_dir = Path(__file__).parent.parent.parent / "examples"
        multi_organism = examples_dir / "multi_organism"
        torch_dir = list(multi_organism.glob("*.torch"))[0]

        torch = Torch.load(torch_dir)

        # Multi-scheme format has schemes attribute
        assert hasattr(torch, "schemes")
        assert isinstance(torch.schemes, dict)

    def test_multi_organism_has_two_loaded_schemes(self):
        """Loaded multi_organism has 2 schemes."""
        examples_dir = Path(__file__).parent.parent.parent / "examples"
        multi_organism = examples_dir / "multi_organism"
        torch_dir = list(multi_organism.glob("*.torch"))[0]

        torch = Torch.load(torch_dir)
        assert len(torch.schemes) == 2
        assert "salmonella" in torch.schemes
        assert "ecoli" in torch.schemes

    def test_each_loaded_scheme_is_a_schema(self):
        """Each loaded scheme is a Schema object."""
        examples_dir = Path(__file__).parent.parent.parent / "examples"
        multi_organism = examples_dir / "multi_organism"
        torch_dir = list(multi_organism.glob("*.torch"))[0]

        torch = Torch.load(torch_dir)

        assert isinstance(torch.schemes["salmonella"], Schema)
        assert isinstance(torch.schemes["ecoli"], Schema)

    def test_each_loaded_scheme_has_correct_profile_count(self):
        """Each loaded scheme has 3 profiles."""
        examples_dir = Path(__file__).parent.parent.parent / "examples"
        multi_organism = examples_dir / "multi_organism"
        torch_dir = list(multi_organism.glob("*.torch"))[0]

        torch = Torch.load(torch_dir)

        assert len(torch.schemes["salmonella"].profiles) == 3
        assert len(torch.schemes["ecoli"].profiles) == 3

    def test_scheme_references_loaded(self):
        """Both schemes have allele references loaded."""
        examples_dir = Path(__file__).parent.parent.parent / "examples"
        multi_organism = examples_dir / "multi_organism"
        torch_dir = list(multi_organism.glob("*.torch"))[0]

        torch = Torch.load(torch_dir)

        assert hasattr(torch, "scheme_references")
        assert "salmonella" in torch.scheme_references
        assert "ecoli" in torch.scheme_references

        # Each scheme should have 3 reference files
        assert len(torch.scheme_references["salmonella"]) == 3
        assert len(torch.scheme_references["ecoli"]) == 3


class TestOldMLSTTorchRemoved:
    """Test that old MLST example torch is removed."""

    def test_old_mlst_workflow_directory_removed(self):
        """torchbase/workflows/mlst/ directory is removed."""
        workflows_dir = Path(__file__).parent.parent / "workflows"
        mlst_dir = workflows_dir / "mlst"

        # The implementation should remove this directory
        # For now, it should exist (this test should fail in RED phase)
        # After implementation, it should not exist
        assert not mlst_dir.exists(), (
            "Old torchbase/workflows/mlst/ directory should be removed. "
            "Found at: " + str(mlst_dir)
        )


class TestExampleTorchesInTests:
    """Test that tests can use the new example torches."""

    def test_tests_can_import_and_use_simple_mlst(self):
        """Tests can load and use simple_mlst example."""
        examples_dir = Path(__file__).parent.parent.parent / "examples"
        simple_mlst = examples_dir / "simple_mlst"
        torch_dir = list(simple_mlst.glob("*.torch"))[0]

        torch = Torch.load(torch_dir)

        # Should be usable for testing
        assert torch.profile is not None
        assert len(torch.profile.profiles) > 0
        assert len(torch.references) > 0

    def test_tests_can_import_and_use_multi_organism(self):
        """Tests can load and use multi_organism example."""
        examples_dir = Path(__file__).parent.parent.parent / "examples"
        multi_organism = examples_dir / "multi_organism"
        torch_dir = list(multi_organism.glob("*.torch"))[0]

        torch = Torch.load(torch_dir)

        # Should be usable for testing
        assert len(torch.schemes) > 0
        assert len(torch.scheme_references) > 0


class TestExampleDocumentation:
    """Test that examples are documented or referenced."""

    def test_readme_or_docs_mention_examples(self):
        """README or docs reference the example torches."""
        repo_root = Path(__file__).parent.parent.parent

        # Check README files
        readme_candidates = [
            repo_root / "README.md",
            repo_root / "README.rst",
            repo_root / "docs" / "examples.md",
            repo_root / "examples" / "README.md",
        ]

        found_docs = False
        for readme in readme_candidates:
            if readme.exists():
                with open(readme) as f:
                    content = f.read().lower()
                    if "example" in content or "simple_mlst" in content:
                        found_docs = True
                        break

        assert found_docs, (
            "Examples should be documented in README or docs/"
        )
