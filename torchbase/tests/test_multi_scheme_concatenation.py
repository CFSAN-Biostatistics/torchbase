"""Tests for multi-scheme torch allele concatenation and profile transformation.

This test suite covers the acceptance criteria for issue #53:
Multi-scheme torch support requires concatenating allele databases with
scheme-prefixed locus names and transforming profile tables to include
scheme identification.
"""

import pytest
from pathlib import Path
from torchbase.torchfs import Torch


class TestAlleleConcatenation:
    """Test allele file concatenation with scheme-prefixed headers."""

    def test_concatenate_alleles_returns_unified_fasta(
        self, multi_scheme_torch_tempdir
    ):
        """Concatenate allele files from all schemes into single FASTA."""
        torch = Torch.load(multi_scheme_torch_tempdir)

        # Call concatenation method - should return Path to unified file
        unified_alleles = torch.concatenate_alleles()

        assert unified_alleles is not None
        assert isinstance(unified_alleles, Path)
        assert unified_alleles.exists()
        assert unified_alleles.suffix == ".fasta"

    def test_concatenated_alleles_have_scheme_prefixes(
        self, multi_scheme_torch_tempdir
    ):
        """Allele headers in concatenated file have scheme prefixes."""
        torch = Torch.load(multi_scheme_torch_tempdir)
        unified_alleles = torch.concatenate_alleles()

        # Read concatenated file
        with open(unified_alleles) as f:
            content = f.read()

        # Check for scheme-prefixed headers
        assert ">ecoli_dinB_1" in content
        assert ">ecoli_dinB_2" in content
        assert ">ecoli_icdA_1" in content
        assert ">salmonella_adk_1" in content
        assert ">salmonella_adk_2" in content
        assert ">salmonella_fumC_1" in content
        assert ">salmonella_fumC_2" in content

    def test_concatenated_alleles_preserve_sequences(
        self, multi_scheme_torch_tempdir
    ):
        """Sequences in concatenated file match original alleles."""
        torch = Torch.load(multi_scheme_torch_tempdir)
        unified_alleles = torch.concatenate_alleles()

        with open(unified_alleles) as f:
            content = f.read()

        # Verify sequences from fixtures are preserved
        assert "ACGT" in content  # ecoli dinB_1
        assert "TGCA" in content  # ecoli dinB_2
        assert "GATC" in content  # ecoli icdA_1
        assert "CCCC" in content  # salmonella adk_1
        assert "GGGG" in content  # salmonella adk_2
        assert "AAAA" in content  # salmonella fumC_1
        assert "TTTT" in content  # salmonella fumC_2

    def test_concatenated_alleles_no_name_collisions(
        self, multi_scheme_torch_tempdir
    ):
        """Scheme prefixes prevent locus name collisions between schemes."""
        torch = Torch.load(multi_scheme_torch_tempdir)
        unified_alleles = torch.concatenate_alleles()

        with open(unified_alleles) as f:
            lines = f.readlines()

        headers = [line for line in lines if line.startswith(">")]

        # Count unique headers - should be 7 total
        unique_headers = set(h.strip() for h in headers)
        assert len(unique_headers) == 7

        # Verify no unprefixed headers exist
        for header in headers:
            header = header.strip()
            assert "_" in header  # All headers should have underscore
            # Headers should start with scheme name
            assert (
                header.startswith(">ecoli_") or
                header.startswith(">salmonella_")
            )

    def test_concatenate_alleles_empty_scheme_raises_error(self, tmp_path):
        """Error if a scheme has no allele files."""
        # Create torch with empty alleles directory
        torch_path = tmp_path / "test_namespace" / "bad_torch" / "1.0.0.torch"
        torch_path.mkdir(parents=True)

        import toml
        import csv

        metadata = {
            "namespace": "test_namespace",
            "name": "bad_torch",
            "version": "1.0.0",
            "version_meta": {
                "strategy": "content-hash",
                "timestamp": 1609459200
            },
            "schemes": {"ecoli": {}},
            "manifest": {}
        }

        with open(torch_path / "metadata.toml", "w") as f:
            toml.dump(metadata, f)

        # Create scheme with profiles but empty alleles
        scheme_path = torch_path / "schemes" / "ecoli"
        scheme_path.mkdir(parents=True)

        profiles = [["ST", "dinB"], ["1", "1"]]
        with open(scheme_path / "profiles.tsv", "w") as f:
            writer = csv.writer(f, delimiter="\t")
            writer.writerows(profiles)

        alleles_path = scheme_path / "alleles"
        alleles_path.mkdir(parents=True)
        # Leave empty

        torch = Torch.load(torch_path)

        with pytest.raises((ValueError, RuntimeError)) as exc_info:
            torch.concatenate_alleles()

        assert "allele" in str(exc_info.value).lower()


class TestProfileTransformation:
    """Test profile table transformation with scheme identification."""

    def test_transform_profiles_returns_unified_tsv(
        self, multi_scheme_torch_tempdir
    ):
        """Transform profiles from all schemes into single TSV."""
        torch = Torch.load(multi_scheme_torch_tempdir)

        unified_profiles = torch.transform_profiles()

        assert unified_profiles is not None
        assert isinstance(unified_profiles, Path)
        assert unified_profiles.exists()
        assert unified_profiles.suffix in (".tsv", ".txt")

    def test_transformed_profiles_have_scheme_column(
        self, multi_scheme_torch_tempdir
    ):
        """Transformed profile table includes 'scheme' column."""
        torch = Torch.load(multi_scheme_torch_tempdir)
        unified_profiles = torch.transform_profiles()

        import csv
        with open(unified_profiles) as f:
            reader = csv.DictReader(f, delimiter="\t")
            rows = list(reader)

        # Check header has 'scheme' column
        assert "scheme" in rows[0].keys()

        # Check scheme values are populated
        scheme_values = {row["scheme"] for row in rows}
        assert "ecoli" in scheme_values
        assert "salmonella" in scheme_values

    def test_transformed_profiles_have_prefixed_locus_names(
        self, multi_scheme_torch_tempdir
    ):
        """Locus column names in transformed table have scheme prefixes."""
        torch = Torch.load(multi_scheme_torch_tempdir)
        unified_profiles = torch.transform_profiles()

        import csv
        with open(unified_profiles) as f:
            reader = csv.reader(f, delimiter="\t")
            header = next(reader)

        # Check for scheme-prefixed locus columns
        assert "ecoli_dinB" in header
        assert "ecoli_icdA" in header
        assert "salmonella_adk" in header
        assert "salmonella_fumC" in header

        # Verify original unprefixed names are NOT in header
        # (would cause ambiguity)
        locus_columns = [col for col in header if col not in ("ST", "scheme")]
        for col in locus_columns:
            assert "_" in col  # All locus names should be prefixed

    def test_transformed_profiles_preserve_allele_values(
        self, multi_scheme_torch_tempdir
    ):
        """Allele values in transformed table match original profiles."""
        torch = Torch.load(multi_scheme_torch_tempdir)
        unified_profiles = torch.transform_profiles()

        import csv
        with open(unified_profiles) as f:
            reader = csv.DictReader(f, delimiter="\t")
            rows = list(reader)

        # Find E. coli ST1 row
        ecoli_st1 = next(
            r for r in rows if r["scheme"] == "ecoli" and r["ST"] == "1"
        )
        assert ecoli_st1["ecoli_dinB"] == "1"
        assert ecoli_st1["ecoli_icdA"] == "1"

        # Find Salmonella ST2 row
        salmonella_st2 = next(
            r for r in rows if r["scheme"] == "salmonella" and r["ST"] == "2"
        )
        assert salmonella_st2["salmonella_adk"] == "1"
        assert salmonella_st2["salmonella_fumC"] == "2"

    def test_transformed_profiles_use_empty_for_cross_scheme_loci(
        self, multi_scheme_torch_tempdir
    ):
        """Cross-scheme loci (not in scheme) have empty or null values."""
        torch = Torch.load(multi_scheme_torch_tempdir)
        unified_profiles = torch.transform_profiles()

        import csv
        with open(unified_profiles) as f:
            reader = csv.DictReader(f, delimiter="\t")
            rows = list(reader)

        # E. coli rows should have empty salmonella loci
        ecoli_st1 = next(
            r for r in rows if r["scheme"] == "ecoli" and r["ST"] == "1"
        )
        # Salmonella loci should be empty/null for E. coli profiles
        assert (
            ecoli_st1.get("salmonella_adk", "") == "" or
            ecoli_st1.get("salmonella_adk") is None
        )
        assert (
            ecoli_st1.get("salmonella_fumC", "") == "" or
            ecoli_st1.get("salmonella_fumC") is None
        )

        # Salmonella rows should have empty E. coli loci
        salmonella_st1 = next(
            r for r in rows if r["scheme"] == "salmonella" and r["ST"] == "1"
        )
        assert (
            salmonella_st1.get("ecoli_dinB", "") == "" or
            salmonella_st1.get("ecoli_dinB") is None
        )
        assert (
            salmonella_st1.get("ecoli_icdA", "") == "" or
            salmonella_st1.get("ecoli_icdA") is None
        )

    def test_transformed_profiles_maintain_st_uniqueness_per_scheme(
        self, multi_scheme_torch_tempdir
    ):
        """ST numbers are unique within each scheme (not globally)."""
        torch = Torch.load(multi_scheme_torch_tempdir)
        unified_profiles = torch.transform_profiles()

        import csv
        with open(unified_profiles) as f:
            reader = csv.DictReader(f, delimiter="\t")
            rows = list(reader)

        # Both schemes can have ST1, ST2, ST3
        ecoli_sts = {r["ST"] for r in rows if r["scheme"] == "ecoli"}
        salmonella_sts = {r["ST"] for r in rows if r["scheme"] == "salmonella"}

        assert "1" in ecoli_sts
        assert "2" in ecoli_sts
        assert "3" in ecoli_sts

        assert "1" in salmonella_sts
        assert "2" in salmonella_sts
        assert "3" in salmonella_sts


class TestUnifiedFilesForWorkflows:
    """Test unified files can be passed to workflows."""

    def test_unified_files_both_exist(self, multi_scheme_torch_tempdir):
        """Both unified alleles and profiles are created."""
        torch = Torch.load(multi_scheme_torch_tempdir)

        alleles = torch.concatenate_alleles()
        profiles = torch.transform_profiles()

        assert alleles.exists()
        assert profiles.exists()

    def test_unified_alleles_and_profiles_are_consistent(
        self, multi_scheme_torch_tempdir
    ):
        """Locus names in unified profiles match headers in unified alleles."""
        torch = Torch.load(multi_scheme_torch_tempdir)

        alleles = torch.concatenate_alleles()
        profiles = torch.transform_profiles()

        # Extract locus names from profile header
        import csv
        with open(profiles) as f:
            reader = csv.reader(f, delimiter="\t")
            profile_header = next(reader)

        locus_columns = [
            col for col in profile_header if col not in ("ST", "scheme")
        ]

        # Extract locus names from allele headers
        with open(alleles) as f:
            allele_lines = f.readlines()

        allele_headers = [
            line.strip()[1:]  # Remove '>'
            for line in allele_lines
            if line.startswith(">")
        ]

        # Extract unique loci from allele headers (remove allele numbers)
        # e.g., "ecoli_dinB_1" -> "ecoli_dinB"
        allele_loci = set()
        for header in allele_headers:
            parts = header.rsplit("_", 1)  # Split on last underscore
            if len(parts) == 2:
                locus = parts[0]
                allele_loci.add(locus)

        # Every locus in profiles should have alleles in FASTA
        for locus in locus_columns:
            assert locus in allele_loci, (
                f"Locus {locus} in profiles but not in alleles"
            )

    def test_get_unified_files_returns_both_paths(
        self, multi_scheme_torch_tempdir
    ):
        """Convenience method returns both unified files."""
        torch = Torch.load(multi_scheme_torch_tempdir)

        # Method that returns both files at once
        alleles_path, profiles_path = torch.get_unified_files()

        assert alleles_path.exists()
        assert profiles_path.exists()
        assert alleles_path.suffix == ".fasta"
        assert profiles_path.suffix in (".tsv", ".txt")


class TestBackwardCompatibilitySingleScheme:
    """Test single-scheme torches still work with new API."""

    def test_single_scheme_concatenate_alleles_works(
        self, single_scheme_torch_tempdir
    ):
        """Single-scheme torch concatenation doesn't break."""
        torch = Torch.load(single_scheme_torch_tempdir)

        # Should work without errors
        unified_alleles = torch.concatenate_alleles()

        assert unified_alleles.exists()

    def test_single_scheme_alleles_no_scheme_prefix(
        self, single_scheme_torch_tempdir
    ):
        """Single-scheme alleles don't get redundant scheme prefix.

        For backward compatibility, single-scheme torches should either:
        1. Not add a prefix (preferred), OR
        2. Use the torch name as prefix (acceptable)
        """
        torch = Torch.load(single_scheme_torch_tempdir)
        unified_alleles = torch.concatenate_alleles()

        with open(unified_alleles) as f:
            content = f.read()

        # Should have allele headers - either unprefixed or torch-name-prefixed
        # But NOT scheme-prefixed (since there's no schemes/ directory)
        assert (
            ">dinB_1" in content or  # No prefix (preferred)
            ">legacy_torch_dinB_1" in content  # Torch name prefix (acceptable)
        )

    def test_single_scheme_transform_profiles_works(
        self, single_scheme_torch_tempdir
    ):
        """Single-scheme torch profile transformation doesn't break."""
        torch = Torch.load(single_scheme_torch_tempdir)

        unified_profiles = torch.transform_profiles()

        assert unified_profiles.exists()

    def test_single_scheme_profiles_optional_scheme_column(
        self, single_scheme_torch_tempdir
    ):
        """Single-scheme profiles may or may not have scheme column.

        For backward compatibility, single-scheme can either:
        1. Omit the scheme column (preferred), OR
        2. Add it with torch name as value (acceptable)
        """
        torch = Torch.load(single_scheme_torch_tempdir)
        unified_profiles = torch.transform_profiles()

        import csv
        with open(unified_profiles) as f:
            reader = csv.DictReader(f, delimiter="\t")
            rows = list(reader)

        # Should have profiles
        assert len(rows) > 0

        # Either has no scheme column OR has it with consistent value
        if "scheme" in rows[0].keys():
            # If present, all should be same
            schemes = {row["scheme"] for row in rows}
            assert len(schemes) == 1


class TestConcatenationEdgeCases:
    """Test edge cases in concatenation logic."""

    def test_concatenate_preserves_allele_order(
        self, multi_scheme_torch_tempdir
    ):
        """Alleles within each scheme maintain consistent ordering."""
        torch = Torch.load(multi_scheme_torch_tempdir)
        unified_alleles = torch.concatenate_alleles()

        with open(unified_alleles) as f:
            lines = f.readlines()

        headers = [
            line.strip() for line in lines if line.startswith(">")
        ]

        # ecoli alleles should come before salmonella (alphabetical scheme order)
        ecoli_indices = [
            i for i, h in enumerate(headers) if h.startswith(">ecoli_")
        ]
        salmonella_indices = [
            i for i, h in enumerate(headers) if h.startswith(">salmonella_")
        ]

        if ecoli_indices and salmonella_indices:
            # All ecoli should come before salmonella
            assert max(ecoli_indices) < min(salmonella_indices)

    def test_concatenate_handles_multiline_fasta(self, tmp_path):
        """Concatenation preserves multiline FASTA sequences."""
        # Create torch with multiline sequences
        torch_path = tmp_path / "test_namespace" / "multiline_torch" / "1.0.0.torch"
        torch_path.mkdir(parents=True)

        import toml
        import csv

        metadata = {
            "namespace": "test_namespace",
            "name": "multiline_torch",
            "version": "1.0.0",
            "version_meta": {
                "strategy": "content-hash",
                "timestamp": 1609459200
            },
            "schemes": {"test": {}},
            "manifest": {}
        }

        with open(torch_path / "metadata.toml", "w") as f:
            toml.dump(metadata, f)

        scheme_path = torch_path / "schemes" / "test"
        scheme_path.mkdir(parents=True)

        profiles = [["ST", "geneA"], ["1", "1"]]
        with open(scheme_path / "profiles.tsv", "w") as f:
            writer = csv.writer(f, delimiter="\t")
            writer.writerows(profiles)

        alleles_path = scheme_path / "alleles"
        alleles_path.mkdir(parents=True)

        # Write multiline FASTA
        with open(alleles_path / "geneA.fasta", "w") as f:
            f.write(">geneA_1\n")
            f.write("ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n")
            f.write("ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n")

        torch = Torch.load(torch_path)
        unified_alleles = torch.concatenate_alleles()

        with open(unified_alleles) as f:
            content = f.read()

        # Should have prefixed header
        assert ">test_geneA_1" in content

        # Should preserve the full sequence
        sequence_lines = [
            line for line in content.split("\n")
            if not line.startswith(">") and line.strip()
        ]
        total_sequence = "".join(sequence_lines)
        assert len(total_sequence) >= 104  # 52 * 2 = 104 bases

    def test_transform_handles_special_allele_values(
        self, tmp_path
    ):
        """Profile transformation preserves special values (?, X)."""
        torch_path = tmp_path / "test_namespace" / "special_torch" / "1.0.0.torch"
        torch_path.mkdir(parents=True)

        import toml
        import csv

        metadata = {
            "namespace": "test_namespace",
            "name": "special_torch",
            "version": "1.0.0",
            "version_meta": {
                "strategy": "content-hash",
                "timestamp": 1609459200
            },
            "schemes": {"test": {}},
            "manifest": {}
        }

        with open(torch_path / "metadata.toml", "w") as f:
            toml.dump(metadata, f)

        scheme_path = torch_path / "schemes" / "test"
        scheme_path.mkdir(parents=True)

        # Profiles with special values
        profiles = [
            ["ST", "locus1", "locus2"],
            ["1", "1", "?"],  # IGNORE wildcard
            ["2", "X", "2"],  # EXCLUDE marker
        ]
        with open(scheme_path / "profiles.tsv", "w") as f:
            writer = csv.writer(f, delimiter="\t")
            writer.writerows(profiles)

        alleles_path = scheme_path / "alleles"
        alleles_path.mkdir(parents=True)
        with open(alleles_path / "locus1.fasta", "w") as f:
            f.write(">locus1_1\nACGT\n")
        with open(alleles_path / "locus2.fasta", "w") as f:
            f.write(">locus2_2\nTGCA\n")

        torch = Torch.load(torch_path)
        unified_profiles = torch.transform_profiles()

        with open(unified_profiles) as f:
            reader = csv.DictReader(f, delimiter="\t")
            rows = list(reader)

        # Find profiles and check special values preserved
        st1 = next(r for r in rows if r["ST"] == "1")
        assert st1["test_locus2"] == "?"

        st2 = next(r for r in rows if r["ST"] == "2")
        assert st2["test_locus1"] == "X"
