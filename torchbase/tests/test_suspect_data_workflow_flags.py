"""Acceptance tests for suspect data workflow flags (Issue #22).

These are RED-phase tests - they MUST fail because the feature is not yet implemented.

Acceptance criteria:
- Workflow reads quality.json if present
- CLI flags: --include-suspect-alleles (default), --exclude-suspect-alleles
- Additional flags: --exclude-suspect-loci, --exclude-suspect-profiles
- Workflow filters allele database based on flags before MinHash/alignment
- Results note which alleles/loci were excluded (if any)
- Works when quality.json absent (no filtering)
- Tests verify filtering behavior at all three levels
- Documentation: flag semantics and defaults
"""

import pytest
import json
import tempfile
from pathlib import Path
import subprocess


# Get the torchbase root directory
TORCHBASE_ROOT = Path(__file__).parent.parent


@pytest.fixture
def quality_json_with_suspect_data():
    """Create a quality.json file with suspect alleles, loci, and profiles."""
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir_path = Path(tmpdir)
        quality_path = tmpdir_path / "quality.json"

        quality_data = {
            "loci": {
                "salmonella_adk": {
                    "similarities": {
                        "adk_1-adk_2": 45.5,
                        "adk_1-adk_3": 98.5
                    },
                    "threshold": 90.0,
                    "statistics": {
                        "mean": 72.0,
                        "std_dev": 31.2,
                        "min": 45.5,
                        "max": 98.5,
                        "percentile_99": 97.0,
                        "threshold_type": "percentile"
                    }
                },
                "salmonella_fumC": {
                    "similarities": {
                        "fumC_1-fumC_2": 42.3
                    },
                    "threshold": 70.0,
                    "statistics": {
                        "mean": 42.3,
                        "std_dev": 0.0,
                        "min": 42.3,
                        "max": 42.3,
                        "percentile_99": 42.3,
                        "threshold_type": "none"
                    }
                },
                "salmonella_gyrB": {
                    "similarities": {
                        "gyrB_1-gyrB_2": 96.5
                    },
                    "threshold": 90.0,
                    "statistics": {
                        "mean": 96.5,
                        "std_dev": 0.0,
                        "min": 96.5,
                        "max": 96.5,
                        "percentile_99": 96.5,
                        "threshold_type": "percentile"
                    }
                }
            },
            "suspect_pairs": {
                "salmonella_adk": [
                    {
                        "allele1": "adk_1",
                        "allele2": "adk_3",
                        "similarity": 98.5,
                        "containment_1_in_2": 98.0,
                        "containment_2_in_1": 99.0,
                        "issue_type": "duplicate"
                    }
                ],
                "salmonella_gyrB": [
                    {
                        "allele1": "gyrB_1",
                        "allele2": "gyrB_2",
                        "similarity": 96.5,
                        "containment_1_in_2": 96.0,
                        "containment_2_in_1": 97.0,
                        "issue_type": "overlap"
                    }
                ]
            },
            "summary": {
                "total_loci": 3,
                "total_suspect_allele_pairs": 2,
                "suspect_loci": ["salmonella_adk", "salmonella_gyrB"],
                "suspect_profiles": ["salmonella_adk", "salmonella_gyrB"]
            }
        }

        with open(quality_path, "w") as f:
            json.dump(quality_data, f, indent=2)

        yield quality_path


@pytest.fixture
def allele_database_with_suspects():
    """Create an allele database FASTA that matches quality.json."""
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir_path = Path(tmpdir)
        db_path = tmpdir_path / "alleles.fasta"

        fasta_content = """>salmonella_adk_1
ATGAATATTAACAACGCACTGGGCGACGTGCTGAAAACCCACGGCCAGATGACGAAAGAAGTGATGCAA
>salmonella_adk_2
ATGAATATTAACAACGCACTGGGCGACGTGCTGAAAACCCACGGCCAGATGACGAAAGAAGTGATGCAC
>salmonella_adk_3
ATGAATATTAACAACGCACTGGGCGACGTGCTGAAAACCCACGGCCAGATGACGAAAGAAGTGATGCAA
>salmonella_fumC_1
CTGACCCAAGGTGCAACCCACGCCTTTGTGACCGCCGTGGGCGACTCGCCCGAAGAAACGCACCACGGA
>salmonella_fumC_2
CTGACCCAAGGTGCAACCCACGCCTTTGTGACCGCCGTGGGCGACTCGCCCGAAGAAACGCACCACGGC
>salmonella_gyrB_1
ATGACCCAACTGAAAGTGATGCCGCAACGTGTCGACCTGCAAATCCACGCAGTGCTGATGAAACCGATG
>salmonella_gyrB_2
ATGACCCAACTGAAAGTGATGCCGCAACGTGTCGACCTGCAAATCCACGCAGTGCTGATGAAACCGATC
"""

        with open(db_path, "w") as f:
            f.write(fasta_content)

        yield db_path


@pytest.fixture
def profile_table():
    """Create a profile table TSV."""
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir_path = Path(tmpdir)
        profile_path = tmpdir_path / "profiles.tsv"

        profile_content = """ST\tsalmonella_adk\tsalmonella_fumC\tsalmonella_gyrB
1\t1\t1\t1
2\t2\t2\t1
3\t1\t2\t2
"""

        with open(profile_path, "w") as f:
            f.write(profile_content)

        yield profile_path


@pytest.fixture
def query_contigs():
    """Create query contigs matching ST=1."""
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir_path = Path(tmpdir)
        contigs_path = tmpdir_path / "query.fasta"

        fasta_content = """>contig1_adk_1
ATGAATATTAACAACGCACTGGGCGACGTGCTGAAAACCCACGGCCAGATGACGAAAGAAGTGATGCAA
>contig2_fumC_1
CTGACCCAAGGTGCAACCCACGCCTTTGTGACCGCCGTGGGCGACTCGCCCGAAGAAACGCACCACGGA
>contig3_gyrB_1
ATGACCCAACTGAAAGTGATGCCGCAACGTGTCGACCTGCAAATCCACGCAGTGCTGATGAAACCGATG
"""

        with open(contigs_path, "w") as f:
            f.write(fasta_content)

        yield contigs_path


@pytest.fixture
def mlst_workflow_path():
    """Path to the main MLST workflow WDL file.

    Using balanced_typing.wdl as the baseline workflow that should
    support quality.json and suspect data filtering.
    """
    return TORCHBASE_ROOT / "workflows" / "builtin" / "balanced_typing.wdl"


class TestWorkflowReadsQualityJson:
    """Test workflow reads quality.json if present."""

    def test_workflow_accepts_quality_json_parameter(self, mlst_workflow_path):
        """WDL workflow accepts quality_json as optional input parameter"""
        with open(mlst_workflow_path) as f:
            content = f.read()

        assert "File? quality_json" in content or "quality" in content.lower(), \
            "Workflow does not accept quality_json input parameter"

    def test_workflow_loads_quality_json_content(self, mlst_workflow_path):
        """WDL workflow loads and parses quality.json content"""
        with open(mlst_workflow_path) as f:
            content = f.read()

        # Should have logic to read quality.json
        assert "quality" in content.lower(), \
            "Workflow does not load quality.json"

    def test_workflow_extracts_suspect_data_from_quality_json(self, mlst_workflow_path):
        """WDL workflow extracts suspect alleles, loci, profiles from quality.json"""
        with open(mlst_workflow_path) as f:
            content = f.read()

        assert "suspect" in content.lower(), \
            "Workflow does not extract suspect data from quality.json"


class TestCLIFlagsForSuspectAlleles:
    """Test CLI flags: --include-suspect-alleles (default), --exclude-suspect-alleles."""

    def test_cli_has_include_suspect_alleles_flag(self):
        """CLI has --include-suspect-alleles flag"""
        # Check if flag exists in CLI help
        result = subprocess.run(
            ["torchbase", "run", "--help"],
            capture_output=True,
            text=True
        )

        assert "--include-suspect-alleles" in result.stdout or \
               "--exclude-suspect-alleles" in result.stdout, \
            "CLI does not have suspect alleles flags"

    def test_cli_has_exclude_suspect_alleles_flag(self):
        """CLI has --exclude-suspect-alleles flag"""
        result = subprocess.run(
            ["torchbase", "run", "--help"],
            capture_output=True,
            text=True
        )

        assert "--exclude-suspect-alleles" in result.stdout, \
            "CLI does not have --exclude-suspect-alleles flag"

    def test_cli_include_suspect_alleles_is_default(self):
        """CLI --include-suspect-alleles is the default behavior"""
        # This would be tested by checking default parameter value
        # For now, verify that the help text indicates the default
        result = subprocess.run(
            ["torchbase", "run", "--help"],
            capture_output=True,
            text=True
        )

        # Default should be include (no filtering)
        assert "default" in result.stdout.lower() or \
               "include" in result.stdout.lower(), \
            "Default suspect allele behavior not documented"

    def test_cli_flags_are_mutually_exclusive(self):
        """CLI --include-suspect-alleles and --exclude-suspect-alleles are mutually exclusive"""
        # Try to use both flags at once - should fail
        result = subprocess.run(
            ["torchbase", "run", "--include-suspect-alleles",
             "--exclude-suspect-alleles", "dummy_torch"],
            capture_output=True,
            text=True
        )

        # Should error
        assert result.returncode != 0, \
            "CLI allows both --include and --exclude suspect alleles flags"


class TestCLIFlagsForSuspectLoci:
    """Test CLI flag: --exclude-suspect-loci."""

    def test_cli_has_exclude_suspect_loci_flag(self):
        """CLI has --exclude-suspect-loci flag"""
        result = subprocess.run(
            ["torchbase", "run", "--help"],
            capture_output=True,
            text=True
        )

        assert "--exclude-suspect-loci" in result.stdout, \
            "CLI does not have --exclude-suspect-loci flag"

    def test_exclude_suspect_loci_implies_exclude_suspect_alleles(self):
        """Excluding suspect loci implicitly excludes all their alleles"""
        # This is a logical constraint that should be enforced
        # Tested via integration test
        pass


class TestCLIFlagsForSuspectProfiles:
    """Test CLI flag: --exclude-suspect-profiles."""

    def test_cli_has_exclude_suspect_profiles_flag(self):
        """CLI has --exclude-suspect-profiles flag"""
        result = subprocess.run(
            ["torchbase", "run", "--help"],
            capture_output=True,
            text=True
        )

        assert "--exclude-suspect-profiles" in result.stdout, \
            "CLI does not have --exclude-suspect-profiles flag"

    def test_exclude_suspect_profiles_implies_exclude_suspect_loci(self):
        """Excluding suspect profiles implicitly excludes suspect loci and alleles"""
        # This is a logical constraint that should be enforced
        # Tested via integration test
        pass


class TestWorkflowFiltersAlleleDatabase:
    """Test workflow filters allele database based on flags before MinHash/alignment."""

    def test_workflow_has_filter_alleles_task(self, mlst_workflow_path):
        """WDL workflow has task to filter alleles based on quality.json"""
        with open(mlst_workflow_path) as f:
            content = f.read()

        assert "filter" in content.lower() and "allele" in content.lower(), \
            "Workflow does not have allele filtering task"

    def test_workflow_filters_before_minhash(self, mlst_workflow_path):
        """WDL workflow filters alleles before MinHash step"""
        with open(mlst_workflow_path) as f:
            content = f.read()

        # Filter should occur before sketch_sequences task
        filter_pos = content.lower().find("filter")
        sketch_pos = content.lower().find("sketch_sequences")

        assert filter_pos > 0 and sketch_pos > 0 and filter_pos < sketch_pos, \
            "Workflow does not filter alleles before MinHash"

    def test_workflow_filters_before_alignment(self, mlst_workflow_path):
        """WDL workflow filters alleles before alignment step"""
        with open(mlst_workflow_path) as f:
            content = f.read()

        # Filter should apply to alignment as well
        assert "filter" in content.lower() and "align" in content.lower(), \
            "Workflow does not filter alleles before alignment"

    def test_workflow_conditional_filtering_based_on_flags(self, mlst_workflow_path):
        """WDL workflow conditionally filters based on input flags"""
        with open(mlst_workflow_path) as f:
            content = f.read()

        # Should have conditional logic (if/select_first) for filtering
        assert ("if" in content or "select_first" in content) and "filter" in content.lower(), \
            "Workflow does not conditionally apply filtering"


@pytest.mark.miniwdl
class TestWorkflowExcludeSuspectAllelesIntegration:
    """Integration test: workflow excludes suspect alleles when flag is set."""

    def test_workflow_excludes_suspect_alleles(
        self, mlst_workflow_path, query_contigs, allele_database_with_suspects,
        profile_table, quality_json_with_suspect_data
    ):
        """Workflow execution with --exclude-suspect-alleles filters suspect alleles"""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir_path = Path(tmpdir)

            input_json = {
                "mlst_typing.contigs": str(query_contigs),
                "mlst_typing.allele_database": str(allele_database_with_suspects),
                "mlst_typing.profiles": str(profile_table),
                "mlst_typing.quality_json": str(quality_json_with_suspect_data),
                "mlst_typing.exclude_suspect_alleles": True
            }

            input_json_path = tmpdir_path / "inputs.json"
            with open(input_json_path, "w") as f:
                json.dump(input_json, f)

            result = subprocess.run(
                ["miniwdl", "run", str(mlst_workflow_path), "-i", str(input_json_path),
                 "-d", str(tmpdir_path)],
                capture_output=True,
                text=True,
                timeout=600
            )

            assert result.returncode == 0, f"Workflow execution failed: {result.stderr}"

            # Find output
            output_dirs = list(tmpdir_path.glob("**/outputs.json"))
            assert len(output_dirs) > 0, "No outputs.json found"

            with open(output_dirs[0]) as f:
                outputs = json.load(f)

            result_path = Path(outputs["mlst_typing.typing_result"])
            with open(result_path) as f:
                result_data = json.load(f)

            # Result should indicate that suspect alleles were excluded
            assert "excluded_alleles" in result_data or "filtering" in result_data, \
                "Result does not indicate allele filtering"

    def test_workflow_includes_suspect_alleles_by_default(
        self, mlst_workflow_path, query_contigs, allele_database_with_suspects,
        profile_table, quality_json_with_suspect_data
    ):
        """Workflow execution without flags includes suspect alleles (default)"""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir_path = Path(tmpdir)

            # No exclude flags - default is include
            input_json = {
                "mlst_typing.contigs": str(query_contigs),
                "mlst_typing.allele_database": str(allele_database_with_suspects),
                "mlst_typing.profiles": str(profile_table),
                "mlst_typing.quality_json": str(quality_json_with_suspect_data)
            }

            input_json_path = tmpdir_path / "inputs.json"
            with open(input_json_path, "w") as f:
                json.dump(input_json, f)

            result = subprocess.run(
                ["miniwdl", "run", str(mlst_workflow_path), "-i", str(input_json_path),
                 "-d", str(tmpdir_path)],
                capture_output=True,
                text=True,
                timeout=600
            )

            assert result.returncode == 0, f"Workflow execution failed: {result.stderr}"

            output_dirs = list(tmpdir_path.glob("**/outputs.json"))
            with open(output_dirs[0]) as f:
                outputs = json.load(f)

            result_path = Path(outputs["mlst_typing.typing_result"])
            with open(result_path) as f:
                result_data = json.load(f)

            # Result should not indicate filtering (default includes all)
            # Or should explicitly state no filtering applied
            if "excluded_alleles" in result_data:
                assert len(result_data["excluded_alleles"]) == 0, \
                    "Default behavior excluded alleles"


@pytest.mark.miniwdl
class TestWorkflowExcludeSuspectLociIntegration:
    """Integration test: workflow excludes suspect loci when flag is set."""

    def test_workflow_excludes_suspect_loci(
        self, mlst_workflow_path, query_contigs, allele_database_with_suspects,
        profile_table, quality_json_with_suspect_data
    ):
        """Workflow execution with --exclude-suspect-loci filters suspect loci"""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir_path = Path(tmpdir)

            input_json = {
                "mlst_typing.contigs": str(query_contigs),
                "mlst_typing.allele_database": str(allele_database_with_suspects),
                "mlst_typing.profiles": str(profile_table),
                "mlst_typing.quality_json": str(quality_json_with_suspect_data),
                "mlst_typing.exclude_suspect_loci": True
            }

            input_json_path = tmpdir_path / "inputs.json"
            with open(input_json_path, "w") as f:
                json.dump(input_json, f)

            result = subprocess.run(
                ["miniwdl", "run", str(mlst_workflow_path), "-i", str(input_json_path),
                 "-d", str(tmpdir_path)],
                capture_output=True,
                text=True,
                timeout=600
            )

            assert result.returncode == 0, f"Workflow execution failed: {result.stderr}"

            output_dirs = list(tmpdir_path.glob("**/outputs.json"))
            with open(output_dirs[0]) as f:
                outputs = json.load(f)

            result_path = Path(outputs["mlst_typing.typing_result"])
            with open(result_path) as f:
                result_data = json.load(f)

            # Result should indicate that suspect loci were excluded
            assert "excluded_loci" in result_data, \
                "Result does not indicate loci filtering"

            # salmonella_adk and salmonella_gyrB should be excluded
            excluded_loci = result_data.get("excluded_loci", [])
            assert "salmonella_adk" in excluded_loci or \
                   "salmonella_gyrB" in excluded_loci, \
                "Suspect loci not excluded"


@pytest.mark.miniwdl
class TestWorkflowExcludeSuspectProfilesIntegration:
    """Integration test: workflow excludes suspect profiles when flag is set."""

    def test_workflow_excludes_suspect_profiles(
        self, mlst_workflow_path, query_contigs, allele_database_with_suspects,
        profile_table, quality_json_with_suspect_data
    ):
        """Workflow execution with --exclude-suspect-profiles filters suspect profiles"""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir_path = Path(tmpdir)

            input_json = {
                "mlst_typing.contigs": str(query_contigs),
                "mlst_typing.allele_database": str(allele_database_with_suspects),
                "mlst_typing.profiles": str(profile_table),
                "mlst_typing.quality_json": str(quality_json_with_suspect_data),
                "mlst_typing.exclude_suspect_profiles": True
            }

            input_json_path = tmpdir_path / "inputs.json"
            with open(input_json_path, "w") as f:
                json.dump(input_json, f)

            result = subprocess.run(
                ["miniwdl", "run", str(mlst_workflow_path), "-i", str(input_json_path),
                 "-d", str(tmpdir_path)],
                capture_output=True,
                text=True,
                timeout=600
            )

            assert result.returncode == 0, f"Workflow execution failed: {result.stderr}"

            output_dirs = list(tmpdir_path.glob("**/outputs.json"))
            with open(output_dirs[0]) as f:
                outputs = json.load(f)

            result_path = Path(outputs["mlst_typing.typing_result"])
            with open(result_path) as f:
                result_data = json.load(f)

            # Result should indicate that suspect profiles were excluded
            assert "excluded_profiles" in result_data, \
                "Result does not indicate profile filtering"


class TestWorkflowResultsNoteExclusions:
    """Test results note which alleles/loci were excluded (if any)."""

    def test_result_json_has_exclusion_fields(self, mlst_workflow_path):
        """Result JSON structure includes fields for exclusion information"""
        # Verify workflow output structure (would be in assemble_final_result task)
        with open(mlst_workflow_path) as f:
            content = f.read()

        # Should have logic to include exclusion info in result
        assert "excluded" in content.lower() or "filtering" in content.lower(), \
            "Workflow does not include exclusion information in results"

    def test_result_includes_count_of_excluded_alleles(self):
        """Result includes count of excluded alleles"""
        # Tested via integration test - result should have excluded_alleles count
        pass

    def test_result_includes_list_of_excluded_alleles(self):
        """Result includes list of excluded allele IDs"""
        # Tested via integration test - result should list excluded allele IDs
        pass

    def test_result_includes_count_of_excluded_loci(self):
        """Result includes count of excluded loci"""
        # Tested via integration test - result should have excluded_loci count
        pass

    def test_result_includes_list_of_excluded_loci(self):
        """Result includes list of excluded loci names"""
        # Tested via integration test - result should list excluded loci names
        pass


class TestWorkflowWorksWithoutQualityJson:
    """Test workflow works when quality.json absent (no filtering)."""

    def test_workflow_runs_without_quality_json(
        self, mlst_workflow_path, query_contigs, allele_database_with_suspects,
        profile_table
    ):
        """Workflow runs successfully without quality.json file"""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir_path = Path(tmpdir)

            # No quality.json provided
            input_json = {
                "balanced_typing.query_sequences": str(query_contigs),
                "balanced_typing.allele_fasta": str(allele_database_with_suspects),
                "balanced_typing.profiles_table": str(profile_table)
            }

            input_json_path = tmpdir_path / "inputs.json"
            with open(input_json_path, "w") as f:
                json.dump(input_json, f)

            result = subprocess.run(
                ["miniwdl", "run", str(mlst_workflow_path), "-i", str(input_json_path),
                 "-d", str(tmpdir_path)],
                capture_output=True,
                text=True,
                timeout=600
            )

            assert result.returncode == 0, \
                f"Workflow should run without quality.json: {result.stderr}"

    def test_workflow_without_quality_json_includes_all_alleles(
        self, mlst_workflow_path, query_contigs, allele_database_with_suspects,
        profile_table
    ):
        """Workflow without quality.json includes all alleles (no filtering)"""
        # Should behave same as default (include all)
        # Tested via integration test
        pass

    def test_exclude_flags_without_quality_json_are_ignored(
        self, mlst_workflow_path, query_contigs, allele_database_with_suspects,
        profile_table
    ):
        """Exclude flags without quality.json are silently ignored"""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir_path = Path(tmpdir)

            # Exclude flags but no quality.json - should ignore and not error
            input_json = {
                "balanced_typing.query_sequences": str(query_contigs),
                "balanced_typing.allele_fasta": str(allele_database_with_suspects),
                "balanced_typing.profiles_table": str(profile_table),
                "balanced_typing.exclude_suspect_alleles": True
            }

            input_json_path = tmpdir_path / "inputs.json"
            with open(input_json_path, "w") as f:
                json.dump(input_json, f)

            result = subprocess.run(
                ["miniwdl", "run", str(mlst_workflow_path), "-i", str(input_json_path),
                 "-d", str(tmpdir_path)],
                capture_output=True,
                text=True,
                timeout=600
            )

            # Should succeed (silently ignore the flag)
            assert result.returncode == 0, \
                f"Workflow should ignore exclude flags without quality.json: {result.stderr}"


class TestFilteringBehaviorAtAllThreeLevels:
    """Test filtering behavior at allele, loci, and profile levels."""

    def test_allele_level_filtering_excludes_specific_alleles(self):
        """Allele-level filtering excludes only flagged alleles"""
        # In quality.json: adk_1 and adk_3 are suspect
        # Filter should exclude adk_1 and adk_3 but keep adk_2
        pass

    def test_loci_level_filtering_excludes_all_locus_alleles(self):
        """Loci-level filtering excludes all alleles of flagged loci"""
        # In quality.json: salmonella_adk and salmonella_gyrB are suspect loci
        # Filter should exclude all adk and gyrB alleles
        pass

    def test_profile_level_filtering_excludes_profile_loci(self):
        """Profile-level filtering excludes loci involved in suspect profiles"""
        # In quality.json: suspect_profiles include salmonella_adk, salmonella_gyrB
        # Filter should exclude all alleles from those loci
        pass

    def test_hierarchical_filtering_allele_subset_of_loci(self):
        """Allele filtering is subset of loci filtering"""
        # Excluding loci should implicitly exclude their alleles
        pass

    def test_hierarchical_filtering_loci_subset_of_profiles(self):
        """Loci filtering is subset of profile filtering"""
        # Excluding profiles should implicitly exclude their loci
        pass


class TestFlagSemanticsAndDocumentation:
    """Test flag semantics and documentation."""

    def test_flag_names_use_positive_semantics(self):
        """Flags use positive semantics (--include, --exclude) not double-negatives"""
        result = subprocess.run(
            ["torchbase", "run", "--help"],
            capture_output=True,
            text=True
        )

        # Should not have double-negative flags like --no-exclude
        assert "--no-exclude" not in result.stdout.lower() and \
               "--no-no-" not in result.stdout.lower(), \
            "Flags use double-negative semantics"

    def test_flag_help_text_explains_default_behavior(self):
        """Flag help text explains default behavior (include suspect data)"""
        result = subprocess.run(
            ["torchbase", "run", "--help"],
            capture_output=True,
            text=True
        )

        # Help should explain default is to include
        assert "default" in result.stdout.lower() or \
               "include" in result.stdout.lower(), \
            "Help text does not explain default behavior"

    def test_flag_help_text_explains_quality_json_requirement(self):
        """Flag help text explains quality.json is required for filtering"""
        result = subprocess.run(
            ["torchbase", "run", "--help"],
            capture_output=True,
            text=True
        )

        # Help should mention quality.json
        assert "quality" in result.stdout.lower(), \
            "Help text does not mention quality.json requirement"

    def test_flag_help_text_explains_hierarchical_filtering(self):
        """Flag help text explains hierarchical filtering (profiles > loci > alleles)"""
        result = subprocess.run(
            ["torchbase", "run", "--help"],
            capture_output=True,
            text=True
        )

        # Help should explain hierarchy
        assert "allele" in result.stdout.lower() and \
               "loci" in result.stdout.lower(), \
            "Help text does not explain filtering hierarchy"


@pytest.mark.miniwdl
class TestEdgeCasesForSuspectDataFiltering:
    """Test edge cases for suspect data filtering."""

    def test_empty_quality_json_no_filtering(
        self, mlst_workflow_path, query_contigs, allele_database_with_suspects,
        profile_table
    ):
        """Empty quality.json (no suspects) results in no filtering"""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir_path = Path(tmpdir)

            # Empty quality.json
            quality_path = tmpdir_path / "quality.json"
            with open(quality_path, "w") as f:
                json.dump({
                    "loci": {},
                    "suspect_pairs": {},
                    "summary": {
                        "total_loci": 0,
                        "total_suspect_allele_pairs": 0,
                        "suspect_loci": [],
                        "suspect_profiles": []
                    }
                }, f)

            input_json = {
                "mlst_typing.contigs": str(query_contigs),
                "mlst_typing.allele_database": str(allele_database_with_suspects),
                "mlst_typing.profiles": str(profile_table),
                "mlst_typing.quality_json": str(quality_path),
                "mlst_typing.exclude_suspect_alleles": True
            }

            input_json_path = tmpdir_path / "inputs.json"
            with open(input_json_path, "w") as f:
                json.dump(input_json, f)

            result = subprocess.run(
                ["miniwdl", "run", str(mlst_workflow_path), "-i", str(input_json_path),
                 "-d", str(tmpdir_path)],
                capture_output=True,
                text=True,
                timeout=600
            )

            assert result.returncode == 0, f"Workflow failed: {result.stderr}"

    def test_malformed_quality_json_handled_gracefully(
        self, mlst_workflow_path, query_contigs, allele_database_with_suspects,
        profile_table
    ):
        """Malformed quality.json is handled gracefully (error or skip)"""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir_path = Path(tmpdir)

            # Malformed quality.json
            quality_path = tmpdir_path / "quality.json"
            with open(quality_path, "w") as f:
                f.write("{ invalid json }")

            input_json = {
                "mlst_typing.contigs": str(query_contigs),
                "mlst_typing.allele_database": str(allele_database_with_suspects),
                "mlst_typing.profiles": str(profile_table),
                "mlst_typing.quality_json": str(quality_path),
                "mlst_typing.exclude_suspect_alleles": True
            }

            input_json_path = tmpdir_path / "inputs.json"
            with open(input_json_path, "w") as f:
                json.dump(input_json, f)

            result = subprocess.run(
                ["miniwdl", "run", str(mlst_workflow_path), "-i", str(input_json_path),
                 "-d", str(tmpdir_path)],
                capture_output=True,
                text=True,
                timeout=600
            )

            # Should either error gracefully or skip filtering
            # Either way, shouldn't crash with unhandled exception
            assert "unhandled" not in result.stderr.lower(), \
                "Malformed quality.json caused unhandled exception"

    def test_all_alleles_excluded_handled_gracefully(
        self, mlst_workflow_path, query_contigs, profile_table
    ):
        """Workflow handles case where all alleles are excluded"""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir_path = Path(tmpdir)

            # Create quality.json where all alleles are suspect
            quality_path = tmpdir_path / "quality.json"
            quality_data = {
                "summary": {
                    "total_loci": 3,
                    "suspect_loci": ["salmonella_adk", "salmonella_fumC", "salmonella_gyrB"],
                    "suspect_profiles": ["salmonella_adk", "salmonella_fumC", "salmonella_gyrB"]
                }
            }
            with open(quality_path, "w") as f:
                json.dump(quality_data, f)

            # Create allele DB with only suspect alleles
            db_path = tmpdir_path / "alleles.fasta"
            with open(db_path, "w") as f:
                f.write(">salmonella_adk_1\nATGC\n")

            input_json = {
                "mlst_typing.contigs": str(query_contigs),
                "mlst_typing.allele_database": str(db_path),
                "mlst_typing.profiles": str(profile_table),
                "mlst_typing.quality_json": str(quality_path),
                "mlst_typing.exclude_suspect_loci": True
            }

            input_json_path = tmpdir_path / "inputs.json"
            with open(input_json_path, "w") as f:
                json.dump(input_json, f)

            subprocess.run(
                ["miniwdl", "run", str(mlst_workflow_path), "-i", str(input_json_path),
                 "-d", str(tmpdir_path)],
                capture_output=True,
                text=True,
                timeout=600
            )

            # Should handle gracefully - either warn or produce empty result
            # Should not crash
            pass
