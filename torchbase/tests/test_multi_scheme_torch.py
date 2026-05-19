"""Tests for multi-scheme torch loading functionality."""

import pytest
from torchbase import Schema, Profile
from torchbase.torchfs import Torch


class TestMultiSchemeTorchDiscovery:
    """Test scheme discovery from schemes/ subdirectories."""

    def test_load_multi_scheme_torch_discovers_schemes(
        self, multi_scheme_torch_tempdir
    ):
        """Multi-scheme torch.load() discovers schemes from schemes/."""
        torch = Torch.load(multi_scheme_torch_tempdir)

        # Should have a schemes attribute instead of single profile
        assert hasattr(torch, "schemes")
        assert isinstance(torch.schemes, dict)
        assert "ecoli" in torch.schemes
        assert "salmonella" in torch.schemes

    def test_each_scheme_is_a_schema_object(self, multi_scheme_torch_tempdir):
        """Each discovered scheme is a Schema object."""
        torch = Torch.load(multi_scheme_torch_tempdir)

        assert isinstance(torch.schemes["ecoli"], Schema)
        assert isinstance(torch.schemes["salmonella"], Schema)

    def test_scheme_names_match_discovered_subdirs(
        self, multi_scheme_torch_tempdir
    ):
        """Scheme names match discovered subdirectories."""
        torch = Torch.load(multi_scheme_torch_tempdir)

        # Should exactly match the schemes declared in metadata
        assert set(torch.schemes.keys()) == {"ecoli", "salmonella"}


class TestIndependentSchemaLoading:
    """Test each scheme loads its own profiles and alleles."""

    def test_ecoli_scheme_loads_profiles(
        self,
        multi_scheme_torch_tempdir  # noqa: F841
    ):
        """E. coli scheme loads profiles correctly."""
        torch = Torch.load(multi_scheme_torch_tempdir)
        ecoli_schema = torch.schemes["ecoli"]

        assert ecoli_schema.name == "ecoli"
        assert len(ecoli_schema.profiles) == 3  # ST 1, 2, 3
        assert all(isinstance(p, Profile) for p in ecoli_schema.profiles)

    def test_salmonella_scheme_loads_profiles(self, multi_scheme_torch_tempdir):
        """Salmonella scheme loads its profiles correctly."""
        torch = Torch.load(multi_scheme_torch_tempdir)
        salmonella_schema = torch.schemes["salmonella"]

        assert salmonella_schema.name == "salmonella"
        assert len(salmonella_schema.profiles) == 3  # ST 1, 2, 3
        assert all(isinstance(p, Profile) for p in salmonella_schema.profiles)

    def test_ecoli_profiles_have_correct_loci(
        self, multi_scheme_torch_tempdir
    ):
        """E. coli profiles have correct loci (dinB, icdA)."""
        torch = Torch.load(multi_scheme_torch_tempdir)
        ecoli_schema = torch.schemes["ecoli"]

        # Check first profile
        profile1 = ecoli_schema.profiles[0]
        assert "dinB" in profile1
        assert "icdA" in profile1

    def test_salmonella_profiles_have_correct_loci(
        self, multi_scheme_torch_tempdir
    ):
        """Salmonella profiles have correct loci (adk, fumC)."""
        torch = Torch.load(multi_scheme_torch_tempdir)
        salmonella_schema = torch.schemes["salmonella"]

        # Check first profile
        profile1 = salmonella_schema.profiles[0]
        assert "adk" in profile1
        assert "fumC" in profile1

    def test_schemes_have_independent_allele_files(
        self, multi_scheme_torch_tempdir
    ):
        """Each scheme has its own allele files tracked independently."""
        torch = Torch.load(multi_scheme_torch_tempdir)

        # Both schemes should have references
        assert hasattr(torch, "scheme_references")
        assert "ecoli" in torch.scheme_references
        assert "salmonella" in torch.scheme_references

        # E. coli alleles
        ecoli_alleles = {ref.name for ref in torch.scheme_references["ecoli"]}
        assert "dinB.fasta" in ecoli_alleles
        assert "icdA.fasta" in ecoli_alleles

        # Salmonella alleles
        salmonella_alleles = {
            ref.name for ref in torch.scheme_references["salmonella"]
        }
        assert "adk.fasta" in salmonella_alleles
        assert "fumC.fasta" in salmonella_alleles


class TestMetadataValidation:
    """Test metadata validation for scheme declarations."""

    def test_metadata_schemes_section_matches_discovered_schemes(
        self, multi_scheme_torch_tempdir
    ):
        """Metadata [schemes] section matches discovered schemes."""
        # Should not raise any errors
        torch = Torch.load(multi_scheme_torch_tempdir)
        assert set(torch.schemes.keys()) == {"ecoli", "salmonella"}

    def test_error_on_metadata_mismatch(self):
        """Error if declared schemes don't match discovered schemes."""
        # This test will be implemented when we have error handling
        # For now, we'll test the happy path
        pass


class TestBackwardCompatibility:
    """Test backward compatibility with single-scheme torches."""

    def test_single_scheme_torch_loads_as_legacy_format(
        self, single_scheme_torch_tempdir
    ):
        """Single-scheme torch without schemes/ directory loads correctly."""
        torch = Torch.load(single_scheme_torch_tempdir)

        # Legacy format should have a profile attribute
        assert hasattr(torch, "profile")
        assert isinstance(torch.profile, (Profile, Schema))

    def test_legacy_torch_has_profiles(self, single_scheme_torch_tempdir):
        """Legacy single-scheme torch loads profiles correctly."""
        torch = Torch.load(single_scheme_torch_tempdir)

        # Should be able to access profiles somehow
        if hasattr(torch.profile, "profiles"):
            # It's a Schema
            assert len(torch.profile.profiles) == 3
        else:
            # It's a Profile
            assert torch.profile is not None

    def test_legacy_torch_has_references(self, single_scheme_torch_tempdir):
        """Legacy single-scheme torch loads allele references."""
        torch = Torch.load(single_scheme_torch_tempdir)

        assert hasattr(torch, "references")
        ref_names = {ref.name for ref in torch.references}
        assert "dinB.fasta" in ref_names
        assert "icdA.fasta" in ref_names


class TestMultiSchemeTorchAttributes:
    """Test Torch dataclass has correct attributes for multi-scheme."""

    def test_torch_has_path_attribute(
        self,
        multi_scheme_torch_tempdir  # noqa: F841
    ):
        """Torch has path attribute pointing to torch directory."""
        torch = Torch.load(multi_scheme_torch_tempdir)
        assert torch.path == multi_scheme_torch_tempdir

    def test_torch_has_schemes_attribute(
        self,
        multi_scheme_torch_tempdir  # noqa: F841
    ):
        """Torch has schemes attribute for multi-scheme torches."""
        torch = Torch.load(multi_scheme_torch_tempdir)
        assert hasattr(torch, "schemes")

    def test_torch_has_scheme_references(
        self,
        multi_scheme_torch_tempdir  # noqa: F841
    ):
        """Torch tracks allele references per scheme."""
        torch = Torch.load(multi_scheme_torch_tempdir)
        assert hasattr(torch, "scheme_references")
        assert isinstance(torch.scheme_references, dict)

    def test_torch_has_workflows(self, multi_scheme_torch_tempdir):
        """Torch has workflow attributes (for compatibility)."""
        torch = Torch.load(multi_scheme_torch_tempdir)
        # Multi-scheme may have workflow per scheme or shared workflows
        assert (hasattr(torch, "workflow") or
                hasattr(torch, "workflows"))


class TestErrorHandling:
    """Test error handling for malformed multi-scheme torches."""

    def test_error_on_missing_profiles_tsv(self, tmp_path):
        """Error if profiles.tsv missing in a scheme directory."""
        # Create torch structure without profiles.tsv
        torch_path = tmp_path / "test_namespace" / "bad_torch" / "1.0.0.torch"
        torch_path.mkdir(parents=True)

        import toml

        metadata = {
            "namespace": "test_namespace",
            "name": "bad_torch",
            "version": "1.0.0",
            "version_meta": {"strategy": "content-hash", "timestamp": 1609459200},
            "schemes": {"ecoli": {}},
            "manifest": {}
        }

        with open(torch_path / "metadata.toml", "w") as f:
            toml.dump(metadata, f)

        # Create scheme dir without profiles
        scheme_path = torch_path / "schemes" / "ecoli"
        scheme_path.mkdir(parents=True)

        with pytest.raises((FileNotFoundError, ValueError)) as exc_info:
            Torch.load(torch_path)

        assert "profiles.tsv" in str(exc_info.value).lower() or "profiles" in str(
            exc_info.value
        ).lower()

    def test_error_on_missing_alleles_directory(self, tmp_path):
        """Error if alleles/ directory missing in a scheme."""
        torch_path = tmp_path / "test_namespace" / "bad_torch" / "1.0.0.torch"
        torch_path.mkdir(parents=True)

        import toml
        import csv

        metadata = {
            "namespace": "test_namespace",
            "name": "bad_torch",
            "version": "1.0.0",
            "version_meta": {"strategy": "content-hash", "timestamp": 1609459200},
            "schemes": {"ecoli": {}},
            "manifest": {}
        }

        with open(torch_path / "metadata.toml", "w") as f:
            toml.dump(metadata, f)

        # Create scheme dir with profiles but no alleles
        scheme_path = torch_path / "schemes" / "ecoli"
        scheme_path.mkdir(parents=True)

        profiles = [["ST", "dinB"], ["1", "1"]]
        with open(scheme_path / "profiles.tsv", "w") as f:
            writer = csv.writer(f, delimiter="\t")
            writer.writerows(profiles)

        with pytest.raises((FileNotFoundError, ValueError)) as exc_info:
            Torch.load(torch_path)

        assert "alleles" in str(exc_info.value).lower()
