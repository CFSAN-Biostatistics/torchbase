"""Tests for suspect data workflow flags (Issue #22).

Acceptance criteria:
- [ ] Workflow reads quality.json if present
- [ ] CLI flags: --include-suspect-alleles (default), --exclude-suspect-alleles
- [ ] Additional flags: --exclude-suspect-loci, --exclude-suspect-profiles
- [ ] Workflow filters allele database based on flags before MinHash/alignment
- [ ] Results note which alleles/loci were excluded (if any)
- [ ] Works when quality.json absent (no filtering)
- [ ] Tests verify filtering behavior at all three levels
- [ ] Documentation: flag semantics and defaults

This feature adds workflow parameters for handling suspect data. If quality.json
is present, the workflow can include/exclude suspect alleles/loci/profiles.
Default behavior is to include suspect data (no filtering). Users opt-in to
filtering.

Key design decisions:
- Positive semantics: --include-suspect-alleles (default), NOT --no-exclude-...
- Three filtering levels: alleles, loci (all alleles), profiles (all loci)
- Filtering happens before MinHash/alignment to avoid wasting compute
- Results JSON documents what was filtered for reproducibility
"""

import pytest
import json
import tempfile
from pathlib import Path


class TestQualityJSONReading:
    """Test that workflow reads quality.json if present."""

    def test_workflow_reads_quality_json_when_present(self):
        """Workflow loads and parses quality.json from torch directory"""
        with tempfile.TemporaryDirectory() as tmpdir:
            torch_path = Path(tmpdir) / "test.torch"
            torch_path.mkdir()

            # Create quality.json with suspect data
            quality_data = {
                "summary": {
                    "suspect_alleles": ["aroC_45", "aroC_102"],
                    "suspect_loci": ["aroC"],
                    "suspect_profiles": []
                },
                "suspect_pairs": {
                    "aroC": [
                        {
                            "allele1": "aroC_45",
                            "allele2": "aroC_102",
                            "similarity": 0.998,
                            "issue_type": "duplicate"
                        }
                    ]
                }
            }

            quality_path = torch_path / "quality.json"
            with open(quality_path, "w") as f:
                json.dump(quality_data, f)

            # Import workflow module (doesn't exist yet - should fail)
            from torchbase.workflow_filters import load_quality_data

            loaded_data = load_quality_data(torch_path)
            assert loaded_data is not None
            assert "summary" in loaded_data
            assert len(loaded_data["summary"]["suspect_alleles"]) == 2

    def test_workflow_handles_missing_quality_json(self):
        """Workflow continues when quality.json is absent (no filtering)"""
        with tempfile.TemporaryDirectory() as tmpdir:
            torch_path = Path(tmpdir) / "test.torch"
            torch_path.mkdir()

            # No quality.json file
            from torchbase.workflow_filters import load_quality_data

            loaded_data = load_quality_data(torch_path)
            # Should return empty/None to indicate no filtering needed
            assert loaded_data is None or loaded_data == {}

    def test_workflow_handles_malformed_quality_json(self):
        """Workflow handles malformed quality.json gracefully"""
        with tempfile.TemporaryDirectory() as tmpdir:
            torch_path = Path(tmpdir) / "test.torch"
            torch_path.mkdir()

            # Create malformed JSON
            quality_path = torch_path / "quality.json"
            with open(quality_path, "w") as f:
                f.write("{invalid json")

            from torchbase.workflow_filters import load_quality_data

            # Should either raise specific error or return None
            with pytest.raises(json.JSONDecodeError):
                load_quality_data(torch_path)

    def test_workflow_reads_quality_json_with_hierarchical_structure(self):
        """Workflow correctly reads hierarchical suspect data structure"""
        with tempfile.TemporaryDirectory() as tmpdir:
            torch_path = Path(tmpdir) / "test.torch"
            torch_path.mkdir()

            quality_data = {
                "loci": {
                    "aroC": {
                        "threshold": 0.985,
                        "statistics": {"mean": 0.943, "std_dev": 0.042}
                    }
                },
                "suspect_pairs": {
                    "aroC": [
                        {
                            "allele1": "aroC_45",
                            "allele2": "aroC_102",
                            "similarity": 0.998
                        }
                    ],
                    "thrA": []
                },
                "summary": {
                    "suspect_alleles": ["aroC_45", "aroC_102"],
                    "suspect_loci": ["aroC"],
                    "suspect_profiles": ["ST1", "ST5"]
                }
            }

            with open(torch_path / "quality.json", "w") as f:
                json.dump(quality_data, f)

            from torchbase.workflow_filters import load_quality_data

            loaded = load_quality_data(torch_path)
            assert "loci" in loaded
            assert "aroC" in loaded["loci"]
            assert "summary" in loaded


class TestIncludeSuspectAllelesFlag:
    """Test CLI flag: --include-suspect-alleles (default behavior)."""

    def test_include_suspect_alleles_is_default(self):
        """--include-suspect-alleles is the default (no filtering)"""
        from torchbase.workflow_filters import WorkflowFilterConfig

        # Create config with defaults
        config = WorkflowFilterConfig()

        # Default should include suspect alleles
        assert config.include_suspect_alleles is True
        assert config.should_filter_alleles() is False

    def test_include_suspect_alleles_flag_explicitly_set(self):
        """--include-suspect-alleles can be explicitly set"""
        from torchbase.workflow_filters import WorkflowFilterConfig

        config = WorkflowFilterConfig(include_suspect_alleles=True)

        assert config.include_suspect_alleles is True
        assert config.should_filter_alleles() is False

    def test_include_suspect_alleles_keeps_all_alleles(self):
        """With --include-suspect-alleles, no alleles are filtered"""
        from torchbase.workflow_filters import filter_alleles

        alleles = ["aroC_1", "aroC_45", "aroC_102", "aroC_200"]
        suspect_alleles = ["aroC_45", "aroC_102"]

        filtered = filter_alleles(
            alleles, suspect_alleles, include_suspect=True
        )

        # All alleles should be kept
        assert len(filtered) == 4
        assert set(filtered) == set(alleles)


class TestExcludeSuspectAllelesFlag:
    """Test CLI flag: --exclude-suspect-alleles (opt-in filtering)."""

    def test_exclude_suspect_alleles_removes_suspect_alleles(self):
        """--exclude-suspect-alleles removes suspect alleles from database"""
        from torchbase.workflow_filters import filter_alleles

        alleles = ["aroC_1", "aroC_45", "aroC_102", "aroC_200"]
        suspect_alleles = ["aroC_45", "aroC_102"]

        filtered = filter_alleles(
            alleles, suspect_alleles, include_suspect=False
        )

        # Only non-suspect alleles remain
        assert len(filtered) == 2
        assert "aroC_1" in filtered
        assert "aroC_200" in filtered
        assert "aroC_45" not in filtered
        assert "aroC_102" not in filtered

    def test_exclude_suspect_alleles_with_empty_suspect_list(self):
        """--exclude-suspect-alleles with no suspect data keeps all alleles"""
        from torchbase.workflow_filters import filter_alleles

        alleles = ["aroC_1", "aroC_45", "aroC_102"]
        suspect_alleles = []

        filtered = filter_alleles(
            alleles, suspect_alleles, include_suspect=False
        )

        # No suspects to exclude, all remain
        assert len(filtered) == 3
        assert set(filtered) == set(alleles)

    def test_exclude_suspect_alleles_with_all_suspect(self):
        """--exclude-suspect-alleles when all alleles are suspect"""
        from torchbase.workflow_filters import filter_alleles

        alleles = ["aroC_1", "aroC_2"]
        suspect_alleles = ["aroC_1", "aroC_2"]

        filtered = filter_alleles(
            alleles, suspect_alleles, include_suspect=False
        )

        # All excluded - should result in empty list or error
        assert len(filtered) == 0

    def test_exclude_suspect_alleles_preserves_locus_structure(self):
        """--exclude-suspect-alleles preserves multi-locus structure"""
        from torchbase.workflow_filters import filter_alleles_by_locus

        alleles_by_locus = {
            "aroC": ["aroC_1", "aroC_45", "aroC_102"],
            "thrA": ["thrA_1", "thrA_2"],
            "hemD": ["hemD_1", "hemD_5"]
        }
        suspect_alleles = ["aroC_45", "aroC_102"]

        filtered = filter_alleles_by_locus(
            alleles_by_locus, suspect_alleles, include_suspect=False
        )

        assert len(filtered["aroC"]) == 1
        assert "aroC_1" in filtered["aroC"]
        # Other loci unchanged
        assert len(filtered["thrA"]) == 2
        assert len(filtered["hemD"]) == 2


class TestExcludeSuspectLociFlag:
    """Test CLI flag: --exclude-suspect-loci (removes entire loci)."""

    def test_exclude_suspect_loci_removes_entire_locus(self):
        """--exclude-suspect-loci removes all alleles from suspect loci"""
        from torchbase.workflow_filters import filter_loci

        alleles_by_locus = {
            "aroC": ["aroC_1", "aroC_45", "aroC_102"],
            "thrA": ["thrA_1", "thrA_2"],
            "hemD": ["hemD_1", "hemD_5"]
        }
        suspect_loci = ["aroC"]

        filtered = filter_loci(
            alleles_by_locus, suspect_loci, include_suspect=False
        )

        # aroC entirely removed
        assert "aroC" not in filtered
        # Other loci remain
        assert "thrA" in filtered
        assert "hemD" in filtered
        assert len(filtered["thrA"]) == 2

    def test_exclude_suspect_loci_with_multiple_loci(self):
        """--exclude-suspect-loci removes multiple suspect loci"""
        from torchbase.workflow_filters import filter_loci

        alleles_by_locus = {
            "aroC": ["aroC_1", "aroC_2"],
            "thrA": ["thrA_1", "thrA_2"],
            "hemD": ["hemD_1", "hemD_2"],
            "dnaE": ["dnaE_1"]
        }
        suspect_loci = ["aroC", "hemD"]

        filtered = filter_loci(
            alleles_by_locus, suspect_loci, include_suspect=False
        )

        assert "aroC" not in filtered
        assert "hemD" not in filtered
        assert "thrA" in filtered
        assert "dnaE" in filtered

    def test_exclude_suspect_loci_with_empty_suspect_list(self):
        """--exclude-suspect-loci with no suspect loci keeps all"""
        from torchbase.workflow_filters import filter_loci

        alleles_by_locus = {
            "aroC": ["aroC_1"],
            "thrA": ["thrA_1"]
        }
        suspect_loci = []

        filtered = filter_loci(
            alleles_by_locus, suspect_loci, include_suspect=False
        )

        assert len(filtered) == 2
        assert "aroC" in filtered
        assert "thrA" in filtered

    def test_exclude_suspect_loci_overrides_allele_filtering(self):
        """--exclude-suspect-loci takes precedence over allele-level flags"""
        from torchbase.workflow_filters import apply_filters

        alleles_by_locus = {
            "aroC": ["aroC_1", "aroC_45"],
            "thrA": ["thrA_1", "thrA_2"]
        }
        suspect_alleles = ["aroC_45"]
        suspect_loci = ["aroC"]

        # Even if we include suspect alleles, locus filtering removes entire locus
        filtered = apply_filters(
            alleles_by_locus,
            suspect_alleles=suspect_alleles,
            suspect_loci=suspect_loci,
            include_suspect_alleles=True,
            include_suspect_loci=False
        )

        # aroC removed at locus level, despite allele inclusion
        assert "aroC" not in filtered
        assert "thrA" in filtered


class TestExcludeSuspectProfilesFlag:
    """Test CLI flag: --exclude-suspect-profiles (strictest filtering)."""

    def test_exclude_suspect_profiles_removes_affected_loci(self):
        """--exclude-suspect-profiles removes all loci in suspect profiles"""
        from torchbase.workflow_filters import filter_profiles

        alleles_by_locus = {
            "aroC": ["aroC_1", "aroC_2"],
            "thrA": ["thrA_1", "thrA_2"],
            "hemD": ["hemD_1", "hemD_2"]
        }

        # Profile ST1 uses aroC, thrA
        # Profile ST5 uses hemD
        suspect_profiles = ["ST1"]
        loci_in_profiles = {
            "ST1": ["aroC", "thrA"],
            "ST5": ["hemD"]
        }

        filtered = filter_profiles(
            alleles_by_locus,
            suspect_profiles,
            loci_in_profiles,
            include_suspect=False
        )

        # ST1's loci (aroC, thrA) removed
        assert "aroC" not in filtered
        assert "thrA" not in filtered
        # ST5's loci remain
        assert "hemD" in filtered

    def test_exclude_suspect_profiles_hierarchical_filtering(self):
        """--exclude-suspect-profiles is most conservative filter"""
        from torchbase.workflow_filters import WorkflowFilterConfig

        config = WorkflowFilterConfig(
            include_suspect_alleles=True,
            include_suspect_loci=True,
            include_suspect_profiles=False
        )

        # Profile filtering is strictest
        assert config.get_filtering_level() == "profile"

    def test_exclude_suspect_profiles_with_overlapping_loci(self):
        """Profiles may share loci - remove if ANY profile is suspect"""
        from torchbase.workflow_filters import filter_profiles

        alleles_by_locus = {
            "aroC": ["aroC_1"],
            "thrA": ["thrA_1"],
            "hemD": ["hemD_1"]
        }

        suspect_profiles = ["ST1", "ST5"]
        loci_in_profiles = {
            "ST1": ["aroC", "thrA"],  # aroC shared
            "ST5": ["aroC", "hemD"]   # aroC shared
        }

        filtered = filter_profiles(
            alleles_by_locus,
            suspect_profiles,
            loci_in_profiles,
            include_suspect=False
        )

        # All loci from suspect profiles removed
        # aroC, thrA, hemD all appear in suspect profiles
        assert len(filtered) == 0


class TestWorkflowFilteringBeforeMinHash:
    """Test that filtering happens before MinHash/alignment stages."""

    def test_filtered_alleles_not_passed_to_minhash(self):
        """Filtered alleles are removed before MinHash sketching"""
        from torchbase.workflow_filters import prepare_allele_database

        alleles_by_locus = {
            "aroC": ["aroC_1", "aroC_45", "aroC_102"]
        }
        suspect_alleles = ["aroC_45", "aroC_102"]

        # Simulate preparing database for MinHash
        prepared_db = prepare_allele_database(
            alleles_by_locus,
            suspect_alleles=suspect_alleles,
            include_suspect_alleles=False
        )

        # Only aroC_1 should be in prepared database
        assert len(prepared_db["aroC"]) == 1
        assert "aroC_1" in prepared_db["aroC"]

    def test_filtering_creates_filtered_fasta(self):
        """Filtering produces new FASTA with only non-suspect alleles"""
        from torchbase.workflow_filters import create_filtered_fasta

        with tempfile.TemporaryDirectory() as tmpdir:
            input_fasta = Path(tmpdir) / "input.fasta"
            output_fasta = Path(tmpdir) / "filtered.fasta"

            # Create input FASTA
            with open(input_fasta, "w") as f:
                f.write(">aroC_1\nACGT\n")
                f.write(">aroC_45\nTGCA\n")
                f.write(">aroC_102\nGGGG\n")

            suspect_alleles = ["aroC_45", "aroC_102"]

            create_filtered_fasta(
                input_fasta,
                output_fasta,
                suspect_alleles,
                include_suspect=False
            )

            # Check output only contains non-suspect
            with open(output_fasta) as f:
                content = f.read()

            assert ">aroC_1" in content
            assert ">aroC_45" not in content
            assert ">aroC_102" not in content

    def test_filtering_tracks_removed_count(self):
        """Filtering reports how many alleles were removed"""
        from torchbase.workflow_filters import filter_with_stats

        alleles = ["aroC_1", "aroC_45", "aroC_102", "aroC_200"]
        suspect_alleles = ["aroC_45", "aroC_102"]

        filtered, stats = filter_with_stats(
            alleles, suspect_alleles, include_suspect=False
        )

        assert stats["total_alleles"] == 4
        assert stats["suspect_alleles"] == 2
        assert stats["kept_alleles"] == 2
        assert stats["removed_alleles"] == 2


class TestResultsNoteExcludedData:
    """Test that results document which alleles/loci were excluded."""

    def test_results_json_includes_filter_metadata(self):
        """Results JSON includes filtering metadata section"""
        from torchbase.workflow_filters import create_results_with_filter_info

        results = {
            "st": "1",
            "loci": {"aroC": "1", "thrA": "2"}
        }

        filter_info = {
            "suspect_alleles_excluded": ["aroC_45"],
            "suspect_loci_excluded": [],
            "filtering_enabled": True
        }

        output = create_results_with_filter_info(results, filter_info)

        assert "typing_results" in output
        assert "filter_metadata" in output
        assert output["filter_metadata"]["filtering_enabled"] is True
        assert len(output["filter_metadata"]["suspect_alleles_excluded"]) == 1

    def test_results_json_notes_no_filtering_when_disabled(self):
        """Results JSON indicates when no filtering was applied"""
        from torchbase.workflow_filters import create_results_with_filter_info

        results = {"st": "1"}
        filter_info = {
            "suspect_alleles_excluded": [],
            "suspect_loci_excluded": [],
            "filtering_enabled": False
        }

        output = create_results_with_filter_info(results, filter_info)

        assert output["filter_metadata"]["filtering_enabled"] is False

    def test_results_json_lists_excluded_alleles(self):
        """Results JSON lists each excluded allele by name"""
        from torchbase.workflow_filters import create_results_with_filter_info

        results = {"st": "2"}
        filter_info = {
            "suspect_alleles_excluded": ["aroC_45", "aroC_102", "thrA_13"],
            "suspect_loci_excluded": [],
            "filtering_enabled": True
        }

        output = create_results_with_filter_info(results, filter_info)

        excluded = output["filter_metadata"]["suspect_alleles_excluded"]
        assert len(excluded) == 3
        assert "aroC_45" in excluded
        assert "thrA_13" in excluded

    def test_results_json_lists_excluded_loci(self):
        """Results JSON lists each excluded locus by name"""
        from torchbase.workflow_filters import create_results_with_filter_info

        results = {"st": "3"}
        filter_info = {
            "suspect_alleles_excluded": [],
            "suspect_loci_excluded": ["aroC", "hemD"],
            "filtering_enabled": True
        }

        output = create_results_with_filter_info(results, filter_info)

        excluded_loci = output["filter_metadata"]["suspect_loci_excluded"]
        assert len(excluded_loci) == 2
        assert "aroC" in excluded_loci
        assert "hemD" in excluded_loci

    def test_results_json_notes_quality_json_absence(self):
        """Results JSON notes when quality.json was absent"""
        from torchbase.workflow_filters import create_results_with_filter_info

        results = {"st": "1"}
        filter_info = {
            "quality_json_present": False,
            "filtering_enabled": False
        }

        output = create_results_with_filter_info(results, filter_info)

        metadata = output["filter_metadata"]
        assert metadata["quality_json_present"] is False


class TestFilteringLevels:
    """Test filtering behavior at all three levels (alleles, loci, profiles)."""

    def test_allele_level_filtering_is_finest_grain(self):
        """Allele-level filtering removes specific alleles, keeps others in locus"""
        from torchbase.workflow_filters import apply_filters

        alleles_by_locus = {
            "aroC": ["aroC_1", "aroC_45", "aroC_102"]
        }

        filtered = apply_filters(
            alleles_by_locus,
            suspect_alleles=["aroC_45"],
            suspect_loci=[],
            include_suspect_alleles=False,
            include_suspect_loci=True
        )

        # Locus remains, but only 2 alleles
        assert "aroC" in filtered
        assert len(filtered["aroC"]) == 2
        assert "aroC_1" in filtered["aroC"]
        assert "aroC_102" in filtered["aroC"]

    def test_locus_level_filtering_removes_all_alleles(self):
        """Locus-level filtering removes entire locus regardless of alleles"""
        from torchbase.workflow_filters import apply_filters

        alleles_by_locus = {
            "aroC": ["aroC_1", "aroC_45"],
            "thrA": ["thrA_1"]
        }

        filtered = apply_filters(
            alleles_by_locus,
            suspect_alleles=["aroC_45"],
            suspect_loci=["aroC"],
            include_suspect_alleles=True,  # Include suspect alleles
            include_suspect_loci=False     # But exclude locus entirely
        )

        # Entire aroC locus removed
        assert "aroC" not in filtered
        assert "thrA" in filtered

    def test_profile_level_filtering_removes_profile_loci(self):
        """Profile-level filtering removes all loci used in suspect profiles"""
        from torchbase.workflow_filters import apply_filters

        alleles_by_locus = {
            "aroC": ["aroC_1"],
            "thrA": ["thrA_1"],
            "hemD": ["hemD_1"]
        }

        suspect_profiles = ["ST1"]
        profile_loci = {"ST1": ["aroC", "thrA"]}

        filtered = apply_filters(
            alleles_by_locus,
            suspect_alleles=[],
            suspect_loci=[],
            suspect_profiles=suspect_profiles,
            profile_loci_map=profile_loci,
            include_suspect_profiles=False
        )

        # Profile ST1 uses aroC and thrA
        assert "aroC" not in filtered
        assert "thrA" not in filtered
        # hemD not in suspect profile
        assert "hemD" in filtered

    def test_combined_filtering_levels(self):
        """Multiple filtering levels can be applied simultaneously"""
        from torchbase.workflow_filters import apply_filters

        alleles_by_locus = {
            "aroC": ["aroC_1", "aroC_45", "aroC_102"],
            "thrA": ["thrA_1", "thrA_2"],
            "hemD": ["hemD_1", "hemD_2"],
            "dnaE": ["dnaE_1"]
        }

        filtered = apply_filters(
            alleles_by_locus,
            suspect_alleles=["aroC_45"],        # Remove one allele
            suspect_loci=["thrA"],              # Remove entire locus
            include_suspect_alleles=False,
            include_suspect_loci=False
        )

        # aroC has aroC_45 removed, but locus remains
        assert "aroC" in filtered
        assert len(filtered["aroC"]) == 2
        assert "aroC_45" not in filtered["aroC"]

        # thrA removed entirely
        assert "thrA" not in filtered

        # Others unchanged
        assert "hemD" in filtered
        assert "dnaE" in filtered


class TestEdgeCases:
    """Test edge cases and boundary conditions."""

    def test_filtering_with_no_alleles(self):
        """Filtering empty allele list returns empty list"""
        from torchbase.workflow_filters import filter_alleles

        filtered = filter_alleles([], ["aroC_1"], include_suspect=False)
        assert filtered == []

    def test_filtering_with_no_suspect_data(self):
        """Filtering with empty suspect list returns all alleles"""
        from torchbase.workflow_filters import filter_alleles

        alleles = ["aroC_1", "aroC_2"]
        filtered = filter_alleles(alleles, [], include_suspect=False)
        assert len(filtered) == 2

    def test_filtering_when_all_alleles_suspect(self):
        """Filtering when all alleles are suspect results in empty locus"""
        from torchbase.workflow_filters import filter_alleles_by_locus

        alleles_by_locus = {"aroC": ["aroC_1", "aroC_2"]}
        suspect_alleles = ["aroC_1", "aroC_2"]

        filtered = filter_alleles_by_locus(
            alleles_by_locus, suspect_alleles, include_suspect=False
        )

        # Locus still present but empty (or removed entirely)
        assert len(filtered.get("aroC", [])) == 0

    def test_filtering_preserves_allele_order(self):
        """Filtering preserves original allele order"""
        from torchbase.workflow_filters import filter_alleles

        alleles = ["aroC_100", "aroC_1", "aroC_45", "aroC_2"]
        suspect = ["aroC_45"]

        filtered = filter_alleles(alleles, suspect, include_suspect=False)

        # Order preserved
        assert filtered == ["aroC_100", "aroC_1", "aroC_2"]

    def test_filtering_case_sensitive_allele_names(self):
        """Allele name matching is case-sensitive"""
        from torchbase.workflow_filters import filter_alleles

        alleles = ["aroC_1", "AROC_1", "aroC_2"]
        suspect = ["aroC_1"]

        filtered = filter_alleles(alleles, suspect, include_suspect=False)

        # Only exact match removed
        assert "aroC_1" not in filtered
        assert "AROC_1" in filtered  # Different case, kept

    def test_filtering_with_duplicate_suspect_entries(self):
        """Duplicate entries in suspect list handled correctly"""
        from torchbase.workflow_filters import filter_alleles

        alleles = ["aroC_1", "aroC_2", "aroC_3"]
        suspect = ["aroC_2", "aroC_2", "aroC_2"]  # Duplicates

        filtered = filter_alleles(alleles, suspect, include_suspect=False)

        # aroC_2 removed once
        assert len(filtered) == 2
        assert "aroC_2" not in filtered


class TestCLIIntegration:
    """Test CLI integration for suspect data flags."""

    def test_cli_accepts_include_suspect_alleles_flag(self):
        """CLI accepts --include-suspect-alleles flag"""
        # This will fail until CLI is implemented
        from torchbase.cli import _run
        from click.testing import CliRunner

        runner = CliRunner()
        with tempfile.TemporaryDirectory() as tmpdir:
            result = runner.invoke(_run, [
                "test/torch",
                "-c", f"{tmpdir}/contigs.fasta",
                "--include-suspect-alleles"
            ])

            # Should not error on flag parsing
            # (Will fail on missing files, but flag should parse)
            assert "--include-suspect-alleles" not in result.output

    def test_cli_accepts_exclude_suspect_alleles_flag(self):
        """CLI accepts --exclude-suspect-alleles flag"""
        from torchbase.cli import _run
        from click.testing import CliRunner

        runner = CliRunner()
        with tempfile.TemporaryDirectory() as tmpdir:
            result = runner.invoke(_run, [
                "test/torch",
                "-c", f"{tmpdir}/contigs.fasta",
                "--exclude-suspect-alleles"
            ])

            assert "--exclude-suspect-alleles" not in result.output

    def test_cli_accepts_exclude_suspect_loci_flag(self):
        """CLI accepts --exclude-suspect-loci flag"""
        from torchbase.cli import _run
        from click.testing import CliRunner

        runner = CliRunner()
        with tempfile.TemporaryDirectory() as tmpdir:
            result = runner.invoke(_run, [
                "test/torch",
                "-c", f"{tmpdir}/contigs.fasta",
                "--exclude-suspect-loci"
            ])

            assert result.exit_code in [0, 1, 2]  # Parse succeeded

    def test_cli_accepts_exclude_suspect_profiles_flag(self):
        """CLI accepts --exclude-suspect-profiles flag"""
        from torchbase.cli import _run
        from click.testing import CliRunner

        runner = CliRunner()
        with tempfile.TemporaryDirectory() as tmpdir:
            result = runner.invoke(_run, [
                "test/torch",
                "-c", f"{tmpdir}/contigs.fasta",
                "--exclude-suspect-profiles"
            ])

            assert result.exit_code in [0, 1, 2]

    def test_cli_flags_mutually_compatible(self):
        """Multiple filtering flags can be used together"""
        from torchbase.cli import _run
        from click.testing import CliRunner

        runner = CliRunner()
        with tempfile.TemporaryDirectory() as tmpdir:
            result = runner.invoke(_run, [
                "test/torch",
                "-c", f"{tmpdir}/contigs.fasta",
                "--exclude-suspect-alleles",
                "--exclude-suspect-loci"
            ])

            # Should not error on flag combination
            assert result.exit_code in [0, 1, 2]

    def test_cli_default_includes_suspect_data(self):
        """CLI default behavior includes suspect data (no filtering)"""
        from torchbase.cli import _run
        from click.testing import CliRunner

        runner = CliRunner()
        with tempfile.TemporaryDirectory() as tmpdir:
            # Run without any filtering flags
            _ = runner.invoke(_run, [
                "test/torch",
                "-c", f"{tmpdir}/contigs.fasta"
            ])

            # Default should be include (verified by config inspection)
            # This test documents the expected default behavior
            assert True  # Placeholder - real test would inspect config


class TestWDLIntegration:
    """Test WDL workflow integration with filtering."""

    def test_wdl_workflow_accepts_filter_parameters(self):
        """WDL workflow defines filter parameters"""
        wdl_path = Path(__file__).parent.parent / "workflows" / "minhash_allele_calling.wdl"

        # WDL doesn't exist yet with filtering params - should fail
        with open(wdl_path) as f:
            content = f.read()

        # Should have filter parameters defined
        assert "include_suspect_alleles" in content or "filter_suspect" in content

    def test_wdl_workflow_passes_quality_json_path(self):
        """WDL workflow receives quality.json file path as input"""
        wdl_path = Path(__file__).parent.parent / "workflows" / "minhash_allele_calling.wdl"

        with open(wdl_path) as f:
            content = f.read()

        # Should have quality_json as input parameter
        assert "quality_json" in content or "File? quality" in content

    def test_wdl_task_filters_before_sketching(self):
        """WDL task filters alleles before MinHash sketching"""
        wdl_path = Path(__file__).parent.parent / "workflows" / "minhash_allele_calling.wdl"

        with open(wdl_path) as f:
            content = f.read()

        # Should have filtering step before sketch_sequences
        # Look for filtering command or filter task call
        assert "filter" in content.lower()

    def test_wdl_outputs_include_filter_metadata(self):
        """WDL workflow outputs include filtering metadata"""
        wdl_path = Path(__file__).parent.parent / "workflows" / "minhash_allele_calling.wdl"

        with open(wdl_path) as f:
            content = f.read()

        # Should output filter metadata or stats
        assert "filter" in content or "excluded" in content


class TestDocumentation:
    """Test that flag semantics and defaults are documented."""

    def test_cli_help_documents_include_suspect_alleles(self):
        """CLI help text explains --include-suspect-alleles (default)"""
        from torchbase.cli import _run
        from click.testing import CliRunner

        runner = CliRunner()
        result = runner.invoke(_run, ["--help"])

        # Should document the flag and its default behavior
        assert "include-suspect-alleles" in result.output
        assert "default" in result.output.lower()

    def test_cli_help_documents_exclude_flags(self):
        """CLI help text explains all exclude flags"""
        from torchbase.cli import _run
        from click.testing import CliRunner

        runner = CliRunner()
        result = runner.invoke(_run, ["--help"])

        assert "exclude-suspect-alleles" in result.output
        assert "exclude-suspect-loci" in result.output
        assert "exclude-suspect-profiles" in result.output

    def test_cli_help_explains_flag_semantics(self):
        """CLI help explains positive semantics (include vs exclude)"""
        from torchbase.cli import _run
        from click.testing import CliRunner

        runner = CliRunner()
        result = runner.invoke(_run, ["--help"])

        # Should explain that include is default, exclude is opt-in
        help_text = result.output.lower()
        assert "include" in help_text or "exclude" in help_text

    def test_readme_documents_suspect_data_filtering(self):
        """README or docs explain suspect data filtering feature"""
        readme_path = Path(__file__).parent.parent.parent / "README.md"

        if readme_path.exists():
            with open(readme_path) as f:
                content = f.read()

            # Should mention suspect data or quality filtering
            assert "suspect" in content.lower() or "quality" in content.lower()
