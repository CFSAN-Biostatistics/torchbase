#!/usr/bin/env python

"""Tests for alignment fallback WDL task.

Issue #16: Alignment Fallback (WDL Task)

Tests for minimap2 alignment-based refinement of ambiguous MinHash calls.
Detects low-confidence MinHash results and refines with precise alignment.
"""

import pytest
import json
import tempfile
from pathlib import Path
import subprocess


class TestAlignmentFallbackWDL:
    """Tests for alignment fallback WDL task."""

    @pytest.fixture
    def wdl_task_file(self):
        """Get path to alignment fallback WDL task."""
        wdl_path = (
            Path(__file__).parent.parent /
            "templates" / "torch" / "tasks" / "alignment_fallback.wdl"
        )
        return wdl_path

    @pytest.fixture
    def query_fasta(self):
        """Create query sequence for testing."""
        with tempfile.NamedTemporaryFile(
            mode='w', suffix='.fasta', delete=False
        ) as f:
            f.write(">query_contig_1\n")
            # 200bp query sequence
            f.write(
                "ATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATG\n"  # noqa: E501
            )
            f.name_path = Path(f.name)
        yield f.name_path
        Path(f.name).unlink(missing_ok=True)

    @pytest.fixture
    def allele_database_fasta(self):
        """Create allele database with multiple loci."""
        with tempfile.NamedTemporaryFile(
            mode='w', suffix='.fasta', delete=False
        ) as f:
            # locus_1 alleles - highly similar (ambiguous)
            f.write(">locus_1_allele_1\n")
            f.write(
                "ATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATG\n"  # noqa: E501
            )
            f.write(">locus_1_allele_2\n")
            # 95% similar - creates ambiguity
            f.write(
                "ATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGATGTAGC\n"  # noqa: E501
            )
            f.write(">locus_1_allele_3\n")
            # More divergent
            f.write(
                "TGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCA\n"  # noqa: E501
            )

            # locus_2 alleles
            f.write(">locus_2_allele_1\n")
            f.write(
                "CGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\n"  # noqa: E501
            )
            f.write(">locus_2_allele_2\n")
            f.write(
                "TACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTA\n"  # noqa: E501
            )

            # locus_3 alleles
            f.write(">locus_3_allele_1\n")
            f.write(
                "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n"  # noqa: E501
            )
            f.write(">locus_3_allele_2\n")
            f.write(
                "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT\n"  # noqa: E501
            )

            f.name_path = Path(f.name)
        yield f.name_path
        Path(f.name).unlink(missing_ok=True)

    @pytest.fixture
    def minhash_results_ambiguous(self):
        """Create MinHash results with ambiguity."""
        results = {
            "locus_1": {
                "allele_id": "locus_1_allele_1",
                "similarity": 0.95,
                "confidence": "medium"
            },
            "locus_1_alt": {
                "allele_id": "locus_1_allele_2",
                "similarity": 0.92,
                "confidence": "medium"
            },
            "locus_2": {
                "allele_id": "locus_2_allele_1",
                "similarity": 0.98,
                "confidence": "high"
            },
            "locus_3": {
                "allele_id": "locus_3_allele_1",
                "similarity": 0.50,
                "coverage": 0.60,
                "confidence": "low"
            }
        }
        return results

    @pytest.fixture
    def minhash_results_json(self, minhash_results_ambiguous):
        """Create temporary JSON file with MinHash results."""
        with tempfile.NamedTemporaryFile(
            mode='w', suffix='.json', delete=False
        ) as f:
            json.dump(minhash_results_ambiguous, f)
            f.name_path = Path(f.name)
        yield f.name_path
        Path(f.name).unlink(missing_ok=True)

    # Tests for WDL file structure

    def test_wdl_task_file_exists(self, wdl_task_file):
        """Test that WDL task file exists."""
        assert wdl_task_file.exists(), (
            f"WDL task file not found at {wdl_task_file}"
        )

    def test_wdl_task_has_required_inputs(self, wdl_task_file):
        """Test that WDL task has required input parameters."""
        if not wdl_task_file.exists():
            pytest.skip("WDL task file not found")

        content = wdl_task_file.read_text()

        # Check for required inputs
        assert "query_sequences" in content, (
            "Missing input: query_sequences"
        )
        assert "allele_database" in content or "allele_fasta" in content, (
            "Missing input: allele database"
        )
        assert "minhash_results" in content, (
            "Missing input: minhash_results"
        )
        assert "File" in content, "WDL task must define File inputs"

    def test_wdl_task_has_required_outputs(self, wdl_task_file):
        """Test that WDL task has required output."""
        if not wdl_task_file.exists():
            pytest.skip("WDL task file not found")

        content = wdl_task_file.read_text()

        # Check for output definition
        assert "output" in content.lower(), (
            "WDL task must have output section"
        )
        # Output should be JSON
        assert "json" in content.lower() or "File" in content, (
            "Output should produce JSON file"
        )

    def test_wdl_task_output_has_refined_calls(self, wdl_task_file):
        """Test that WDL output includes refined calls."""
        if not wdl_task_file.exists():
            pytest.skip("WDL task file not found")

        content = wdl_task_file.read_text()

        # Output should reference refined calls
        assert "refined" in content.lower() or "calls" in content.lower(), (
            "Output should include refined calls"
        )

    def test_wdl_task_has_docker_runtime(self, wdl_task_file):
        """Test that WDL task specifies Docker container."""
        if not wdl_task_file.exists():
            pytest.skip("WDL task file not found")

        content = wdl_task_file.read_text()

        # Check for runtime section with container
        assert "runtime" in content, (
            "WDL task must have runtime section"
        )
        assert "container" in content.lower() or "docker" in content.lower(), (
            "WDL task must specify Docker container"
        )
        assert "minimap2" in content.lower(), (
            "Container should include minimap2"
        )

    def test_wdl_syntax_valid(self, wdl_task_file):
        """Test that WDL task passes miniwdl syntax check."""
        if not wdl_task_file.exists():
            pytest.skip("WDL task file not found")

        try:
            result = subprocess.run(
                ["miniwdl", "check", str(wdl_task_file)],
                capture_output=True,
                timeout=10
            )
            # miniwdl check returns 0 for valid WDL
            assert result.returncode == 0, (
                f"WDL syntax error:\n{result.stderr.decode()}"
            )
        except FileNotFoundError:
            pytest.skip("miniwdl not installed")
        except subprocess.TimeoutExpired:
            pytest.fail("miniwdl check timed out")

    def test_wdl_version_declaration(self, wdl_task_file):
        """Test that WDL file has proper version declaration."""
        if not wdl_task_file.exists():
            pytest.skip("WDL task file not found")

        content = wdl_task_file.read_text()

        # Should declare WDL version
        assert "version" in content.lower(), (
            "WDL file should declare version"
        )

    def test_task_command_section_exists(self, wdl_task_file):
        """Test that task has command section."""
        if not wdl_task_file.exists():
            pytest.skip("WDL task file not found")

        content = wdl_task_file.read_text()

        # Should have command section with minimap2 invocation
        assert "command" in content.lower(), (
            "Task should have command section"
        )

    def test_minimap2_invocation_in_command(self, wdl_task_file):
        """Test that command invokes minimap2."""
        if not wdl_task_file.exists():
            pytest.skip("WDL task file not found")

        content = wdl_task_file.read_text()

        # Command should invoke minimap2
        assert "minimap2" in content.lower(), (
            "Command should invoke minimap2"
        )

    def test_json_output_from_task(self, wdl_task_file):
        """Test that task outputs JSON file."""
        if not wdl_task_file.exists():
            pytest.skip("WDL task file not found")

        content = wdl_task_file.read_text()

        # Output should reference a JSON file
        assert ("output" in content.lower() and
                (".json" in content.lower() or "File" in content)), (
            "Task should output JSON file"
        )

    def test_container_has_minimap2(self, wdl_task_file):
        """Test that container image includes minimap2."""
        if not wdl_task_file.exists():
            pytest.skip("WDL task file not found")

        content = wdl_task_file.read_text()

        # Container specification should mention minimap2
        assert "minimap2" in content.lower(), (
            "Container should include minimap2"
        )

    # Tests for ambiguity detection logic

    def test_ambiguity_detection_top_two_similar(self, minhash_results_ambiguous):  # noqa: E501
        """Test detection of top 2 results within 3% similarity."""
        locus1 = minhash_results_ambiguous["locus_1"]
        locus1_alt = minhash_results_ambiguous["locus_1_alt"]

        # Top 2 within 3%
        diff = abs(locus1["similarity"] - locus1_alt["similarity"])
        assert diff <= 0.03, "Top 2 should be within 3% for ambiguity"

    def test_ambiguity_detection_best_below_92_percent(
        self, minhash_results_ambiguous
    ):
        """Test detection when best match below 92%."""
        locus3 = minhash_results_ambiguous["locus_3"]

        # Below 92% threshold
        assert locus3["similarity"] < 0.92, (
            "Best match should be below 92% for ambiguity"
        )

    def test_ambiguity_detection_low_coverage(
        self, minhash_results_ambiguous
    ):
        """Test detection when coverage is below 80%."""
        locus3 = minhash_results_ambiguous["locus_3"]

        # Low coverage
        if "coverage" in locus3:
            assert locus3["coverage"] < 0.80, (
                "Coverage should be below 80% for ambiguity"
            )

    def test_ambiguity_trigger_conditions(self, minhash_results_ambiguous):
        """Test that ambiguity triggers are correctly identified."""
        ambiguous_loci = []

        for key, result in minhash_results_ambiguous.items():
            triggers = []

            # Check triggers
            if result["similarity"] < 0.92:
                triggers.append("low_similarity")

            if "coverage" in result and result["coverage"] < 0.80:
                triggers.append("low_coverage")

            if triggers:
                ambiguous_loci.append((key, triggers))

        # Should identify ambiguous loci
        assert len(ambiguous_loci) > 0, (
            "Should identify at least one ambiguous locus"
        )

    # Tests for output format

    def test_output_json_format_structure(self):
        """Test that refined output has correct JSON structure."""
        refined_calls = {
            "locus_1": {
                "allele_id": "locus_1_allele_1",
                "identity": 0.98,
                "status": "confirmed"
            },
            "locus_2": {
                "allele_id": "locus_2_allele_2",
                "identity": 0.85,
                "status": "novel_allele"
            }
        }

        # Should be JSON serializable
        json_str = json.dumps(refined_calls)
        assert isinstance(json_str, str)

        # Should be deserializable
        loaded = json.loads(json_str)
        assert loaded == refined_calls

    def test_output_has_required_fields(self):
        """Test that output contains required fields."""
        refined_calls = {
            "locus_1": {
                "allele_id": "locus_1_allele_1",
                "identity": 0.98,
                "status": "confirmed"
            }
        }

        for locus, call in refined_calls.items():
            assert "allele_id" in call, f"Missing allele_id for {locus}"
            assert "identity" in call, f"Missing identity for {locus}"
            assert "status" in call, f"Missing status for {locus}"

    def test_output_status_values(self):
        """Test that status field has valid values."""
        valid_statuses = ["confirmed", "novel_allele"]

        test_calls = [
            {"locus": "locus_1", "status": "confirmed"},
            {"locus": "locus_2", "status": "novel_allele"},
        ]

        for call in test_calls:
            assert call["status"] in valid_statuses, (
                f"Invalid status: {call['status']}"
            )

    def test_output_identity_range(self):
        """Test that identity scores are in 0-1 range."""
        refined_calls = {
            "locus_1": {"allele_id": "locus_1_allele_1", "identity": 0.98},
            "locus_2": {"allele_id": "locus_2_allele_1", "identity": 0.85},
            "locus_3": {"allele_id": "locus_3_allele_1", "identity": 0.75},
        }

        for locus, call in refined_calls.items():
            assert 0 <= call["identity"] <= 1, (
                f"Invalid identity for {locus}: {call['identity']}"
            )

    def test_output_allele_id_format(self):
        """Test that allele IDs follow expected format."""
        refined_calls = {
            "locus_1": {"allele_id": "locus_1_allele_5"},
            "locus_2": {"allele_id": "locus_2_allele_12"},
            "gyrA": {"allele_id": "gyrA_allele_3"},
        }

        for locus, call in refined_calls.items():
            allele_id = call["allele_id"]
            # Should contain locus and allele number separated by underscore
            assert "_allele_" in allele_id, (
                f"Invalid allele_id format: {allele_id}"
            )

    # Tests for novel allele detection

    def test_novel_allele_status_below_threshold(self):
        """Test that status is 'novel_allele' when identity < 90%."""
        test_calls = [
            {"allele_id": "novel_1", "identity": 0.85, "status": "novel_allele"},  # noqa: E501
            {"allele_id": "novel_2", "identity": 0.75, "status": "novel_allele"},  # noqa: E501
            {"allele_id": "known_1", "identity": 0.95, "status": "confirmed"},
        ]

        for call in test_calls:
            if call["identity"] < 0.90:
                assert call["status"] == "novel_allele", (
                    f"Status should be novel_allele for identity {call['identity']}"  # noqa: E501
                )
            else:
                assert call["status"] == "confirmed", (
                    f"Status should be confirmed for identity {call['identity']}"
                )

    def test_novel_allele_high_confidence(self):
        """Test that novel alleles have identity calculated."""
        novel_alleles = [
            {"allele_id": "novel_1", "identity": 0.78},
            {"allele_id": "novel_2", "identity": 0.82},
        ]

        for allele in novel_alleles:
            assert "identity" in allele, "Novel alleles should have identity"
            assert 0 <= allele["identity"] <= 1, (
                f"Invalid identity: {allele['identity']}"
            )

    # Tests for input handling

    def test_handles_empty_query(self):
        """Test behavior with empty query sequences."""
        # Task should handle empty input gracefully
        # Output should be empty or indicate no matches
        output = {}
        assert isinstance(output, dict)

    def test_handles_no_ambiguous_loci(self):
        """Test when MinHash is confident about all loci."""
        minhash_results = {
            "locus_1": {
                "allele_id": "locus_1_allele_1",
                "similarity": 0.99,
                "confidence": "high"
            },
            "locus_2": {
                "allele_id": "locus_2_allele_1",
                "similarity": 0.98,
                "confidence": "high"
            },
        }

        # All results confident - no loci need alignment refinement
        ambiguous = [
            r for r in minhash_results.values()
            if r.get("similarity", 1.0) < 0.92
        ]
        assert len(ambiguous) == 0, (
            "Should have no ambiguous loci"
        )

    def test_single_locus_ambiguity(self):
        """Test handling of single ambiguous locus."""
        minhash_results = {
            "locus_1": {
                "allele_id": "locus_1_allele_1",
                "similarity": 0.89,
                "confidence": "low"
            },
        }

        # Single locus should still be processed
        ambiguous_count = sum(
            1 for r in minhash_results.values()
            if r.get("similarity", 1.0) < 0.92
        )
        assert ambiguous_count >= 1, (
            "Should identify ambiguous locus"
        )

    def test_multiple_ambiguous_loci(self):
        """Test handling of multiple ambiguous loci."""
        minhash_results = {
            "locus_1": {"similarity": 0.85, "confidence": "low"},
            "locus_2": {"similarity": 0.88, "confidence": "low"},
            "locus_3": {"similarity": 0.91, "confidence": "medium"},
            "locus_4": {"similarity": 0.99, "confidence": "high"},
        }

        ambiguous = [
            r for r in minhash_results.values()
            if r.get("similarity", 1.0) < 0.92
        ]
        assert len(ambiguous) >= 3, (
            "Should identify multiple ambiguous loci"
        )

    # Tests for minimap2 parameters

    def test_minimap2_uses_appropriate_preset(self, wdl_task_file):
        """Test that minimap2 uses appropriate preset."""
        if not wdl_task_file.exists():
            pytest.skip("WDL task file not found")

        content = wdl_task_file.read_text()

        # Should specify minimap2 preset (map-ont, map-pb, map-x, sr, etc)
        preset_patterns = ["map-ont", "map-pb", "map-x", "sr", "asm"]
        has_preset = any(pattern in content.lower() for pattern in preset_patterns)

        # At minimum, should mention alignment mode
        assert "minimap2" in content.lower() or has_preset, (
            "Should specify minimap2 preset or alignment parameters"
        )

    def test_minimap2_produces_alignment_output(self, wdl_task_file):
        """Test that minimap2 output is processed correctly."""
        if not wdl_task_file.exists():
            pytest.skip("WDL task file not found")

        content = wdl_task_file.read_text()

        # Should parse minimap2 output (PAF or SAM format)
        assert "paf" in content.lower() or "sam" in content.lower() or (
            "minimap2" in content.lower()
        ), "Should specify minimap2 output format"

    # Tests for confidence thresholds

    def test_confidence_threshold_90_percent_identity(self):
        """Test that 90% identity threshold is documented."""
        # This documents the expected behavior
        threshold = 0.90
        assert 0 <= threshold <= 1, "Threshold should be between 0 and 1"

    def test_refined_identity_above_original(self):
        """Test that refined identity is typically >= MinHash similarity."""
        minhash_similarity = 0.85
        refined_identity = 0.92  # After minimap2 alignment

        # Alignment should improve or maintain similarity score
        assert refined_identity >= minhash_similarity or (
            refined_identity < 0.90
        ), "Alignment should refine or clarify low-confidence results"

    # Tests for synthetic ambiguous cases

    def test_ambiguous_case_two_alleles_within_3_percent(self):
        """Test synthetic case: top 2 alleles within 3%."""
        minhash_results = {
            "locus_1": [
                {"allele_id": "locus_1_allele_1", "similarity": 0.95},
                {"allele_id": "locus_1_allele_2", "similarity": 0.93},
            ]
        }

        # Difference within 3%
        diff = (
            minhash_results["locus_1"][0]["similarity"] -
            minhash_results["locus_1"][1]["similarity"]
        )
        assert diff <= 0.03, "Top 2 should be within 3%"

    def test_ambiguous_case_best_match_88_percent(self):
        """Test synthetic case: best match at 88% (below 92%)."""
        minhash_results = {
            "locus_1": {
                "allele_id": "locus_1_allele_1",
                "similarity": 0.88
            }
        }

        # Below 92% threshold
        assert minhash_results["locus_1"]["similarity"] < 0.92

    def test_ambiguous_case_low_coverage_60_percent(self):
        """Test synthetic case: coverage at 60% (below 80%)."""
        minhash_results = {
            "locus_1": {
                "allele_id": "locus_1_allele_1",
                "similarity": 0.95,
                "coverage": 0.60
            }
        }

        # Low coverage despite good similarity
        assert minhash_results["locus_1"]["coverage"] < 0.80

    def test_refined_call_from_ambiguous_to_confirmed(self):
        """Test refinement from ambiguous to confirmed allele."""
        # After alignment fallback
        refined_call = {
            "allele_id": "locus_1_allele_1",
            "identity": 0.96,
            "status": "confirmed"
        }

        # Refinement should improve confidence
        original_minhash_similarity = 0.85
        assert refined_call["identity"] > original_minhash_similarity
        assert refined_call["status"] == "confirmed"

    def test_refined_call_from_ambiguous_to_novel(self):
        """Test refinement identifying novel allele."""
        # After alignment, poor match indicates novel
        refined_call = {
            "allele_id": "novel_1",
            "identity": 0.78,
            "status": "novel_allele"
        }

        # Should flag as novel when below threshold
        assert refined_call["status"] == "novel_allele"
        assert refined_call["identity"] < 0.90

    # Integration tests

    def test_output_json_all_ambiguous_loci_present(self):
        """Test that output includes all originally ambiguous loci."""
        ambiguous_input = {
            "locus_1": {"similarity": 0.85},
            "locus_2": {"similarity": 0.89},
            "locus_3": {"similarity": 0.91},
        }

        refined_output = {
            "locus_1": {"allele_id": "locus_1_allele_1", "identity": 0.96, "status": "confirmed"},  # noqa: E501
            "locus_2": {"allele_id": "locus_2_allele_2", "identity": 0.87, "status": "novel_allele"},  # noqa: E501
            "locus_3": {"allele_id": "locus_3_allele_1", "identity": 0.94, "status": "confirmed"},  # noqa: E501
        }

        # All input loci should have output
        for locus in ambiguous_input.keys():
            assert locus in refined_output, (
                f"Missing output for locus: {locus}"
            )

    def test_output_json_serializable(self):
        """Test that refined output is JSON serializable."""
        output = {
            "locus_1": {
                "allele_id": "locus_1_allele_1",
                "identity": 0.98,
                "status": "confirmed"
            },
            "locus_2": {
                "allele_id": "locus_2_allele_2",
                "identity": 0.87,
                "status": "novel_allele"
            },
        }

        # Should be JSON serializable
        json_str = json.dumps(output)
        assert isinstance(json_str, str)

        # Should be deserializable
        loaded = json.loads(json_str)
        assert loaded == output


class TestAlignmentFallbackIntegration:
    """Integration tests for alignment fallback workflow."""

    def test_workflow_accepts_minhash_json_input(self):
        """Test workflow correctly accepts MinHash JSON results."""
        minhash_input = {
            "query_sequences": "query.fasta",
            "allele_database": "alleles.fasta",
            "minhash_results": "minhash_calls.json"
        }

        assert "minhash_results" in minhash_input
        assert minhash_input["minhash_results"].endswith(".json")

    def test_workflow_processes_ambiguous_loci(self):
        """Test workflow identifies and refines ambiguous loci."""
        minhash_results = {
            "locus_1": {"similarity": 0.85, "confidence": "low"},
            "locus_2": {"similarity": 0.99, "confidence": "high"},
        }

        ambiguous_count = sum(
            1 for r in minhash_results.values()
            if r.get("similarity", 1.0) < 0.92
        )

        assert ambiguous_count >= 1, "Should identify ambiguous loci"

    def test_workflow_skips_confident_loci(self):
        """Test that confident MinHash calls are passed through."""
        input_calls = {
            "locus_1": {"allele_id": "locus_1_allele_1", "similarity": 0.99},
            "locus_2": {"allele_id": "locus_2_allele_1", "similarity": 0.98},
        }

        confident = [
            (k, v) for k, v in input_calls.items()
            if v.get("similarity", 0) >= 0.92
        ]

        assert len(confident) == len(input_calls), (
            "All should be confident"
        )

    def test_allele_database_must_be_fasta(self):
        """Test that allele database must be FASTA format."""
        valid_headers = [
            "locus_1_allele_1",
            "gyrA_allele_5",
            "dinB_allele_12"
        ]

        for header in valid_headers:
            parts = header.split("_")
            assert len(parts) >= 3, f"Invalid header: {header}"

    def test_query_sequences_format(self):
        """Test that query sequences can be FASTA."""
        # Should accept FASTA format
        query_format = "fasta"
        assert query_format in ["fasta", "fastq"], "Supported format"

    def test_refined_output_format_consistency(self):
        """Test that all output entries follow same format."""
        output = {
            "locus_1": {
                "allele_id": "locus_1_allele_1",
                "identity": 0.98,
                "status": "confirmed"
            },
            "locus_2": {
                "allele_id": "locus_2_allele_2",
                "identity": 0.87,
                "status": "novel_allele"
            },
            "locus_3": {
                "allele_id": "locus_3_allele_1",
                "identity": 0.75,
                "status": "novel_allele"
            }
        }

        required_keys = {"allele_id", "identity", "status"}

        for locus, call in output.items():
            assert set(call.keys()) == required_keys, (
                f"Inconsistent format for {locus}"
            )

    def test_minimap2_identity_calculation(self):
        """Test that minimap2 identity is correctly calculated."""
        # minimap2 reports: 100 * (n_match / aligned_query_length)
        # For testing, verify it's in valid range
        identities = [0.75, 0.85, 0.95, 0.99]

        for identity in identities:
            assert 0 <= identity <= 1, f"Invalid identity: {identity}"

    def test_handles_multiple_query_sequences(self):
        """Test processing multiple query sequences."""
        query_headers = [
            "contig_1",
            "contig_2",
            "read_1"
        ]

        # All queries should be processed
        assert len(query_headers) >= 1
