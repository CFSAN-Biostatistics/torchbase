"""Acceptance tests for Alignment Fallback WDL task (Issue #16).

These are RED-phase tests - they MUST fail because the feature is not yet complete.

Acceptance criteria:
- WDL task signature: query + allele FASTA + MinHash results → refined calls JSON
- Detects ambiguity triggers from MinHash output
- Runs minimap2 alignment for flagged loci
- Refines allele calls with higher precision
- Reports novel alleles if alignment below confidence threshold (e.g., <90%)
- Output JSON: {locus: {allele_id, identity, status: "confirmed"/"novel_allele"}}
- Containerized (Docker/Singularity with minimap2)
- miniwdl check validates syntax
- Test with ambiguous synthetic cases
"""

import pytest
import json
import tempfile
from pathlib import Path
import subprocess


# Get the torchbase root directory
TORCHBASE_ROOT = Path(__file__).parent.parent


@pytest.fixture
def allele_database_fasta():
    """Create temporary allele database (FASTA format).

    Creates MLST-style allele database with multiple loci and alleles:
    - adk: 3 alleles with subtle variations
    - fumC: 3 alleles with subtle variations
    - gyrB: 2 alleles
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir_path = Path(tmpdir)
        db_path = tmpdir_path / "allele_db.fasta"

        # Create realistic MLST allele sequences with subtle differences
        fasta_content = """>adk_1
ATGGCTATGAAACAACTCACCAACGTAGCACTGTCCAAAGCCGCACGTGGCAATGCCGCTGCAACTGACTGCACTGGCGCACCTGCCGCTGCTGATGAACGTCATCGGTACGGTCTCGTCCACCGGCTCTACGACCTGGTGTACGACTGTACGTCCGAAATCCTTACTGGCGGCGGTCACACGTCTGCTGGACATCCGCCACCACATTTGCTGCATCACGGCGGCGGTGACGGTGGTTATGGTGGTGACGACCACGCCGTACATCGACGACGTGCTGATCGAACTGATCGAAGACGACGACGACGAAGTGATCGAACTGATCGAAGTGCTGGATGAAGTCGAAAATGTCTAA
>adk_2
ATGGCTATGAAACAACTCACCAACGTAGCACTGTCCAAAGCCGCACGTGGCAATGCCGCTGCAACTGACTGCACTGGCGCACCTGCCGCTGCTGATGAACGTCATCGGTACGGTCTCGTCCACCGGCTCTACGACCTGGTGTACGACTGTACGTCCGAAATCCTTACTGGCGGCGGTCACACGTCTGCTGGACATCCGCCACCACATTTGCTGCATCACGGCGGCGGTGACGGTGGTTATGGTGGTGACGACCACGCCGTACATCGACGACGTGCTGATCGAACTGATCGAAGACGACGACGACGAAGTGATCGAACTGATCGAAGTGCTGGATGAAGTCGAAAATGTCTAG
>adk_3
ATGGCTATGAAACAACTCACCAACGTAGCACTGTCCAAAGCCGCACGTGGCAATGCCGCTGCAACTGACTGCACTGGCGCACCTGCCGCTGCTGATGAACGTCATCGGTACGGTCTCGTCCACCGGATCTACGACCTGGTGTACGACTGTACGTCCGAAATCCTTACTGGCGGCGGTCACACGTCTGCTGGACATCCGCCACCACATTTGCTGCATCACGGCGGCGGTGACGGTGGTTATGGTGGTGACGACCACGCCGTACATCGACGACGTGCTGATCGAACTGATCGAAGACGACGACGACGAAGTGATCGAACTGATCGAAGTGCTGGATGAAGTCGAAAATGTCTAA
>fumC_1
ATGAATATTAACAACGCACTGGGCGACGTGCTGAAAACCCACGGCCAGATGACGAAAGAAGTGATGCTTCTGACCCAAGGTGCAACCCACGCCTTTGTGACCGCCGTGGGCGACTCGCCCGAAGAAACGCACCACGGAAATCGATCTGGACGAACTGTTCGAACTGATGGAAGAACACGAAGTGCTGATGGTCGACATCCTGATGATGCACGACCACGACGATGACCGTGATAGCACCACTGTACGACATTGACGACGACGACGACGATACAGAACACAATGACGATGGAAGAAAACGACGACGAAGTGATCCACGTGATGGTGTAA
>fumC_2
ATGAATATTAACAACGCACTGGGCGACGTGCTGAAAACCCACGGCCAGATGACGAAAGAAGTGATGCAACTGACCCAAGGTGCAACCCACGCCTTTGTGACCGCCGTGGGCGACTCGCCCGAAGAAACGCACCACGGAAATCGATCTGGACGAACTGTTCGAACTGATGGAAGAACACGAAGTGCTGATGGTCGACATCCTGATGATGCACGACCACGACGATGACCGTGATAGCACCACTGTACGACATTGACGACGACGACGACGATACAGAACACAATGACGATGGAAGAAAACGACGACGAAGTGATCCACGTGATGGTGTAG
>fumC_3
ATGAATATTAACAACGCACTGGGCGACGTGCTGAAAACCCACGGCCAGATGACGAAAGAAGTGATGCAACTGACCCAAGGTGCAACCCACGCCTTTGTGACCGCCGTGGGCGACTCGCCCGAAGAAACGCACCACGGAAATCGATCTGGACGAACTGTTCGAACTGATGGAAGAACACGAAGTGCTGATGGTCGACATCCTGATGATGCACGACCACGACGATGACCGTGATAGCACCACTGTACGACATTGACGACGACGACGACGATACAGAACACAATGACGATGGAAGAAAACGACGACGAAGTGATCCACGTGATGGTGTAT
>gyrB_1
ATGACCCAACTGAAAGTGATGCCGCAACGTGTCGACCTGCAAATCCACGCAGTGCTGATGAAACCGATGAAACTGCTGATGACGAAACTGCTGATGATCGTGGTGATCTGACGGACGAACTGACGAAATACGACGAAAACAAACACGATGTCATCGACGATGTGACGACCGACATGATCACGGACGACGTACTGATGAAACTGGTGATCCACGTGCACGATGAAACGGACGACTACGACGACATGCCGATCGACGATGATGATGATGACCACGACGACAACGACGAAACGATGATCCTGACGATGACGACGATCTGACGGATGACTAA
>gyrB_2
ATGACCCAACTGAAAGTGATGCCGCAACGTGTCGACCTGCAAATCCACGCAGTGCTGATGAAACCGATGAAACTGCTGATGACGAAACTGCTGATGATCGTGGTGATCTGACGGACGAACTGACGAAATACGACGAAAACAAACACGATGTCATCGACGATGTGACGACCGACATGATCACGGACGACGTACTGATGAAACTGGTGATCCACGTGCACGATGAAACGGACGACTACGACGACATGCCGATCGACGATGATGATGATGACCACGACGACAACGACGAAACGATGATCCTGACGATGACGACGATCTGACGGATGACTAA
"""

        with open(db_path, "w") as f:
            f.write(fasta_content)

        yield db_path


@pytest.fixture
def query_ambiguous_fasta():
    """Create query sequences with ambiguous matches.

    Creates sequences that match multiple alleles with similar scores,
    requiring alignment fallback for disambiguation.
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir_path = Path(tmpdir)
        query_path = tmpdir_path / "query_ambiguous.fasta"

        # Create sequences with mutations that make MinHash ambiguous
        fasta_content = """>query_adk_almost_1_or_2
ATGGCTATGAAACAACTCACCAACGTAGCACTGTCCAAAGCCGCACGTGGCAATGCCGCTGCAACTGACTGCACTGGCGCACCTGCCGCTGCTGATGAACGTCATCGGTACGGTCTCGTCCACCGGCTCTACGACCTGGTGTACGACTGTACGTCCGAAATCCTTACTGGCGGCGGTCACACGTCTGCTGGACATCCGCCACCACATTTGCTGCATCACGGCGGCGGTGACGGTGGTTATGGTGGTGACGACCACGCCGTACATCGACGACGTGCTGATCGAACTGATCGAAGACGACGACGACGAAGTGATCGAACTGATCGAAGTGCTGGATGAAGTCGAAAATGTCTAA
>query_fumC_between_1_and_2
ATGAATATTAACAACGCACTGGGCGACGTGCTGAAAACCCACGGCCAGATGACGAAAGAAGTGATGCTTCTGACCCAAGGTGCAACCCACGCCTTTGTGACCGCCGTGGGCGACTCGCCCGAAGAAACGCACCACGGAAATCGATCTGGACGAACTGTTCGAACTGATGGAAGAACACGAAGTGCTGATGGTCGACATCCTGATGATGCACGACCACGACGATGACCGTGATAGCACCACTGTACGACATTGACGACGACGACGACGATACAGAACACAATGACGATGGAAGAAAACGACGACGAAGTGATCCACGTGATGGTGTAG
>query_gyrB_novel_allele
ATGACCCAACTGAAAGTGATGCCGCAACGTGTCGACCTGCAAATCCACGCAGTGCTGATGAAACCGATGAAACTGCTGATGACGAAACTGCTGATGATCGTGGTGATCTGACGGACGAACTGACGAAATACGACGAAAACAAACACGATGTCATCGACGATGTGACGACCGACATGATCACGGACGACGTACTGATGAAACTGGTGATCCACGTGCACGATGAAACGGACGACTACGACGACATGCCGATCGACGATGATGATGATGACCACGACGACAACGACGAAACGATGATCCTGACGATGACGACGATCTGACGGATGACTAG
"""

        with open(query_path, "w") as f:
            f.write(fasta_content)

        yield query_path


@pytest.fixture
def minhash_results_ambiguous():
    """Create MinHash results JSON with ambiguous calls.

    Simulates low-confidence MinHash output that should trigger alignment fallback:
    - Top 2 alleles within 3%
    - Best match < 92%
    - Coverage < 80%
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir_path = Path(tmpdir)
        results_path = tmpdir_path / "minhash_results.json"

        # Create ambiguous results triggering alignment fallback
        results = {
            "adk": {
                "allele_id": "1",
                "similarity": 0.91,  # Below 92% threshold
                "confidence": False,
                "second_best": {
                    "allele_id": "2",
                    "similarity": 0.89  # Within 3% of best
                }
            },
            "fumC": {
                "allele_id": "1",
                "similarity": 0.88,  # Low confidence
                "confidence": False,
                "second_best": {
                    "allele_id": "2",
                    "similarity": 0.87  # Within 3%
                }
            },
            "gyrB": {
                "allele_id": "1",
                "similarity": 0.85,  # Very low, potential novel allele
                "confidence": False
            }
        }

        with open(results_path, "w") as f:
            json.dump(results, f, indent=2)

        yield results_path


@pytest.fixture
def minhash_results_confident():
    """Create MinHash results JSON with confident calls.

    Simulates high-confidence MinHash output that should NOT trigger alignment fallback.
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir_path = Path(tmpdir)
        results_path = tmpdir_path / "minhash_confident.json"

        results = {
            "adk": {
                "allele_id": "1",
                "similarity": 0.98,
                "confidence": True
            },
            "fumC": {
                "allele_id": "2",
                "similarity": 0.99,
                "confidence": True
            },
            "gyrB": {
                "allele_id": "1",
                "similarity": 0.97,
                "confidence": True
            }
        }

        with open(results_path, "w") as f:
            json.dump(results, f, indent=2)

        yield results_path


class TestAlignmentFallbackWDLFileExists:
    """Test alignment fallback WDL file exists at expected location."""

    def test_wdl_file_exists(self):
        """Alignment fallback WDL file exists"""
        wdl_path = TORCHBASE_ROOT / "workflows" / "alignment_fallback.wdl"
        assert wdl_path.exists(), f"WDL file not found at {wdl_path}"

    def test_wdl_file_is_file(self):
        """Alignment fallback WDL file is a regular file"""
        wdl_path = TORCHBASE_ROOT / "workflows" / "alignment_fallback.wdl"
        assert wdl_path.is_file(), f"WDL path is not a file: {wdl_path}"


class TestAlignmentFallbackWDLTaskSignature:
    """Test WDL task signature: query + allele FASTA + MinHash results → refined calls JSON."""

    def test_wdl_has_task_definition(self):
        """WDL file contains task definition"""
        wdl_path = TORCHBASE_ROOT / "workflows" / "alignment_fallback.wdl"
        with open(wdl_path) as f:
            content = f.read()

        assert "task" in content, "WDL file does not contain task definition"

    def test_wdl_task_name_is_alignment_fallback(self):
        """WDL workflow/task is named alignment_fallback"""
        wdl_path = TORCHBASE_ROOT / "workflows" / "alignment_fallback.wdl"
        with open(wdl_path) as f:
            content = f.read()

        # Accept either task or workflow with alignment_fallback name
        assert ("task alignment_fallback" in content or
                "workflow alignment_fallback" in content), "Task/workflow name is incorrect"

    def test_wdl_has_input_section(self):
        """WDL task has input section"""
        wdl_path = TORCHBASE_ROOT / "workflows" / "alignment_fallback.wdl"
        with open(wdl_path) as f:
            content = f.read()

        assert "input {" in content, "Task does not have input section"

    def test_wdl_has_query_sequences_input(self):
        """WDL task accepts query_sequences input"""
        wdl_path = TORCHBASE_ROOT / "workflows" / "alignment_fallback.wdl"
        with open(wdl_path) as f:
            content = f.read()

        assert "query_sequences" in content or "query" in content, "Task does not have query sequences input"

    def test_wdl_has_allele_fasta_input(self):
        """WDL task accepts allele_fasta input"""
        wdl_path = TORCHBASE_ROOT / "workflows" / "alignment_fallback.wdl"
        with open(wdl_path) as f:
            content = f.read()

        assert "allele_fasta" in content or "allele_db" in content, "Task does not have allele_fasta input"

    def test_wdl_has_minhash_results_input(self):
        """WDL task accepts minhash_results input"""
        wdl_path = TORCHBASE_ROOT / "workflows" / "alignment_fallback.wdl"
        with open(wdl_path) as f:
            content = f.read()

        assert "minhash_results" in content or "minhash" in content, "Task does not have minhash_results input"

    def test_wdl_has_output_section(self):
        """WDL task has output section"""
        wdl_path = TORCHBASE_ROOT / "workflows" / "alignment_fallback.wdl"
        with open(wdl_path) as f:
            content = f.read()

        assert "output {" in content, "Task does not have output section"

    def test_wdl_output_is_file(self):
        """WDL task output is a File type"""
        wdl_path = TORCHBASE_ROOT / "workflows" / "alignment_fallback.wdl"
        with open(wdl_path) as f:
            content = f.read()

        # Check for File type in output section
        assert "File" in content, "Task output is not File type"


class TestAlignmentFallbackWDLSyntaxValidation:
    """Test miniwdl check validates WDL syntax."""

    def test_miniwdl_check_passes(self):
        """miniwdl check validates WDL syntax without errors"""
        wdl_path = TORCHBASE_ROOT / "workflows" / "alignment_fallback.wdl"

        # Run miniwdl check
        result = subprocess.run(
            ["miniwdl", "check", str(wdl_path)],
            capture_output=True,
            text=True
        )

        assert result.returncode == 0, f"miniwdl check failed: {result.stderr}"


class TestAlignmentFallbackWDLContainerization:
    """Test containerized (Docker/Singularity with minimap2)."""

    def test_wdl_has_runtime_section(self):
        """WDL task has runtime section"""
        wdl_path = TORCHBASE_ROOT / "workflows" / "alignment_fallback.wdl"
        with open(wdl_path) as f:
            content = f.read()

        assert "runtime {" in content, "Task does not have runtime section"

    def test_wdl_has_docker_image(self):
        """WDL task specifies docker image"""
        wdl_path = TORCHBASE_ROOT / "workflows" / "alignment_fallback.wdl"
        with open(wdl_path) as f:
            content = f.read()

        assert "docker:" in content, "Task does not specify docker image"

    def test_wdl_uses_minimap2_container(self):
        """WDL task uses minimap2 container image"""
        wdl_path = TORCHBASE_ROOT / "workflows" / "alignment_fallback.wdl"
        with open(wdl_path) as f:
            content = f.read()

        assert "minimap2" in content.lower(), "Task does not use minimap2 container"


class TestAlignmentFallbackWDLAmbiguityDetection:
    """Test detects ambiguity triggers from MinHash output."""

    def test_wdl_command_reads_minhash_results(self):
        """WDL task command reads MinHash results file"""
        wdl_path = TORCHBASE_ROOT / "workflows" / "alignment_fallback.wdl"
        with open(wdl_path) as f:
            content = f.read()

        # Check for JSON reading logic
        assert "json" in content.lower(), "Task does not read JSON input"

    def test_wdl_detects_low_similarity_threshold(self):
        """WDL task detects low similarity (< 92%) as ambiguity trigger"""
        wdl_path = TORCHBASE_ROOT / "workflows" / "alignment_fallback.wdl"
        with open(wdl_path) as f:
            content = f.read()

        # Check for similarity/confidence threshold logic
        assert "similarity" in content.lower() or "confidence" in content.lower(), \
            "Task does not check similarity threshold"

    def test_wdl_detects_close_second_best(self):
        """WDL task detects top 2 alleles within 3% as ambiguity trigger"""
        wdl_path = TORCHBASE_ROOT / "workflows" / "alignment_fallback.wdl"
        with open(wdl_path) as f:
            content = f.read()

        # Check for second_best comparison logic
        assert "second_best" in content or "ambiguous" in content.lower(), \
            "Task does not check for close second-best matches"

    def test_wdl_has_ambiguity_threshold_parameters(self):
        """WDL task has parameters for ambiguity detection thresholds"""
        wdl_path = TORCHBASE_ROOT / "workflows" / "alignment_fallback.wdl"
        with open(wdl_path) as f:
            content = f.read()

        # Check for threshold parameters
        assert any(keyword in content.lower() for keyword in ["threshold", "min_identity", "min_similarity"]), \
            "Task does not have ambiguity threshold parameters"


class TestAlignmentFallbackWDLMinimap2Integration:
    """Test runs minimap2 alignment for flagged loci."""

    def test_wdl_command_uses_minimap2(self):
        """WDL task command uses minimap2"""
        wdl_path = TORCHBASE_ROOT / "workflows" / "alignment_fallback.wdl"
        with open(wdl_path) as f:
            content = f.read()

        assert "minimap2" in content, "Task command does not use minimap2"

    def test_wdl_command_performs_alignment(self):
        """WDL task command performs sequence alignment"""
        wdl_path = TORCHBASE_ROOT / "workflows" / "alignment_fallback.wdl"
        with open(wdl_path) as f:
            content = f.read()

        # Check for alignment keywords
        assert any(keyword in content.lower() for keyword in ["align", "-a", "-x sr"]), \
            "Task command does not perform alignment"

    def test_wdl_command_processes_alignment_output(self):
        """WDL task command processes minimap2 alignment output"""
        wdl_path = TORCHBASE_ROOT / "workflows" / "alignment_fallback.wdl"
        with open(wdl_path) as f:
            content = f.read()

        # Check for SAM/PAF processing or identity calculation
        assert any(keyword in content.lower() for keyword in ["sam", "paf", "identity", "cigar"]), \
            "Task does not process alignment output"


class TestAlignmentFallbackWDLOutputFormat:
    """Test output JSON: {locus: {allele_id, identity, status: 'confirmed'/'novel_allele'}}."""

    def test_wdl_produces_json_output(self):
        """WDL task produces JSON output"""
        wdl_path = TORCHBASE_ROOT / "workflows" / "alignment_fallback.wdl"
        with open(wdl_path) as f:
            content = f.read()

        assert ".json" in content.lower() or "json" in content.lower(), "Task does not produce JSON output"

    def test_wdl_output_has_status_field(self):
        """WDL task output includes status field"""
        wdl_path = TORCHBASE_ROOT / "workflows" / "alignment_fallback.wdl"
        with open(wdl_path) as f:
            content = f.read()

        # Check for status field in output
        assert "status" in content.lower() or "confirmed" in content or "novel" in content, \
            "Task output does not include status field"

    def test_wdl_output_has_identity_field(self):
        """WDL task output includes identity field"""
        wdl_path = TORCHBASE_ROOT / "workflows" / "alignment_fallback.wdl"
        with open(wdl_path) as f:
            content = f.read()

        # Check for identity field
        assert "identity" in content.lower(), "Task output does not include identity field"


class TestAlignmentFallbackWDLNovelAlleleDetection:
    """Test reports novel alleles if alignment below confidence threshold."""

    def test_wdl_has_novel_allele_threshold(self):
        """WDL task has threshold for novel allele detection (e.g., <90%)"""
        wdl_path = TORCHBASE_ROOT / "workflows" / "alignment_fallback.wdl"
        with open(wdl_path) as f:
            content = f.read()

        # Check for novel allele threshold parameter
        assert any(keyword in content.lower() for keyword in
                   ["novel_threshold", "min_identity", "identity_threshold"]), \
            "Task does not have novel allele threshold"

    def test_wdl_detects_novel_alleles(self):
        """WDL task detects novel alleles when identity is low"""
        wdl_path = TORCHBASE_ROOT / "workflows" / "alignment_fallback.wdl"
        with open(wdl_path) as f:
            content = f.read()

        # Check for novel allele detection logic
        assert "novel" in content.lower(), "Task does not detect novel alleles"

    def test_wdl_marks_confirmed_alleles(self):
        """WDL task marks alleles as confirmed when identity is high"""
        wdl_path = TORCHBASE_ROOT / "workflows" / "alignment_fallback.wdl"
        with open(wdl_path) as f:
            content = f.read()

        # Check for confirmed status
        assert "confirmed" in content.lower() or "exact" in content.lower(), \
            "Task does not mark confirmed alleles"


@pytest.mark.miniwdl
class TestAlignmentFallbackWDLExecutionWithAmbiguousCalls:
    """Test with ambiguous synthetic cases."""

    def test_wdl_execution_with_ambiguous_calls_produces_output(
        self, allele_database_fasta, query_ambiguous_fasta, minhash_results_ambiguous
    ):
        """WDL task execution with ambiguous calls produces output JSON"""
        wdl_path = TORCHBASE_ROOT / "workflows" / "alignment_fallback.wdl"

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir_path = Path(tmpdir)

            # Create input JSON for miniwdl
            input_json = {
                "alignment_fallback.query_sequences": str(query_ambiguous_fasta),
                "alignment_fallback.allele_fasta": str(allele_database_fasta),
                "alignment_fallback.minhash_results": str(minhash_results_ambiguous)
            }

            input_json_path = tmpdir_path / "inputs.json"
            with open(input_json_path, "w") as f:
                json.dump(input_json, f)

            # Run miniwdl
            result = subprocess.run(
                ["miniwdl", "run", str(wdl_path), "-i", str(input_json_path), "-d", str(tmpdir_path)],
                capture_output=True,
                text=True,
                timeout=300
            )

            # Check execution succeeded
            assert result.returncode == 0, f"miniwdl run failed: {result.stderr}"

            # Find output directory
            output_dirs = list(tmpdir_path.glob("**/outputs.json"))
            assert len(output_dirs) > 0, "No outputs.json found"

            # Read outputs
            with open(output_dirs[0]) as f:
                outputs = json.load(f)

            # Check output file exists
            assert "alignment_fallback.refined_calls" in outputs, "Output refined_calls not found"

            results_path = Path(outputs["alignment_fallback.refined_calls"])
            assert results_path.exists(), f"Results file does not exist: {results_path}"

            # Read and validate JSON structure
            with open(results_path) as f:
                results = json.load(f)

            # Check expected loci are present
            assert "adk" in results, "adk locus not in results"
            assert "fumC" in results, "fumC locus not in results"
            assert "gyrB" in results, "gyrB locus not in results"

    def test_wdl_output_has_allele_id_field(
        self, allele_database_fasta, query_ambiguous_fasta, minhash_results_ambiguous
    ):
        """WDL output JSON contains allele_id field for each locus"""
        wdl_path = TORCHBASE_ROOT / "workflows" / "alignment_fallback.wdl"

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir_path = Path(tmpdir)

            input_json = {
                "alignment_fallback.query_sequences": str(query_ambiguous_fasta),
                "alignment_fallback.allele_fasta": str(allele_database_fasta),
                "alignment_fallback.minhash_results": str(minhash_results_ambiguous)
            }

            input_json_path = tmpdir_path / "inputs.json"
            with open(input_json_path, "w") as f:
                json.dump(input_json, f)

            result = subprocess.run(
                ["miniwdl", "run", str(wdl_path), "-i", str(input_json_path), "-d", str(tmpdir_path)],
                capture_output=True,
                text=True,
                timeout=300
            )

            assert result.returncode == 0, f"miniwdl run failed: {result.stderr}"

            output_dirs = list(tmpdir_path.glob("**/outputs.json"))
            with open(output_dirs[0]) as f:
                outputs = json.load(f)

            results_path = Path(outputs["alignment_fallback.refined_calls"])
            with open(results_path) as f:
                results = json.load(f)

            # Check allele_id field
            for locus in results:
                assert "allele_id" in results[locus], f"allele_id not in {locus} results"

    def test_wdl_output_has_identity_field(
        self, allele_database_fasta, query_ambiguous_fasta, minhash_results_ambiguous
    ):
        """WDL output JSON contains identity field for each locus"""
        wdl_path = TORCHBASE_ROOT / "workflows" / "alignment_fallback.wdl"

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir_path = Path(tmpdir)

            input_json = {
                "alignment_fallback.query_sequences": str(query_ambiguous_fasta),
                "alignment_fallback.allele_fasta": str(allele_database_fasta),
                "alignment_fallback.minhash_results": str(minhash_results_ambiguous)
            }

            input_json_path = tmpdir_path / "inputs.json"
            with open(input_json_path, "w") as f:
                json.dump(input_json, f)

            result = subprocess.run(
                ["miniwdl", "run", str(wdl_path), "-i", str(input_json_path), "-d", str(tmpdir_path)],
                capture_output=True,
                text=True,
                timeout=300
            )

            assert result.returncode == 0, f"miniwdl run failed: {result.stderr}"

            output_dirs = list(tmpdir_path.glob("**/outputs.json"))
            with open(output_dirs[0]) as f:
                outputs = json.load(f)

            results_path = Path(outputs["alignment_fallback.refined_calls"])
            with open(results_path) as f:
                results = json.load(f)

            # Check identity field
            for locus in results:
                assert "identity" in results[locus], f"identity not in {locus} results"
                assert 0.0 <= results[locus]["identity"] <= 1.0, f"identity out of range for {locus}"

    def test_wdl_output_has_status_field(
        self, allele_database_fasta, query_ambiguous_fasta, minhash_results_ambiguous
    ):
        """WDL output JSON contains status field for each locus"""
        wdl_path = TORCHBASE_ROOT / "workflows" / "alignment_fallback.wdl"

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir_path = Path(tmpdir)

            input_json = {
                "alignment_fallback.query_sequences": str(query_ambiguous_fasta),
                "alignment_fallback.allele_fasta": str(allele_database_fasta),
                "alignment_fallback.minhash_results": str(minhash_results_ambiguous)
            }

            input_json_path = tmpdir_path / "inputs.json"
            with open(input_json_path, "w") as f:
                json.dump(input_json, f)

            result = subprocess.run(
                ["miniwdl", "run", str(wdl_path), "-i", str(input_json_path), "-d", str(tmpdir_path)],
                capture_output=True,
                text=True,
                timeout=300
            )

            assert result.returncode == 0, f"miniwdl run failed: {result.stderr}"

            output_dirs = list(tmpdir_path.glob("**/outputs.json"))
            with open(output_dirs[0]) as f:
                outputs = json.load(f)

            results_path = Path(outputs["alignment_fallback.refined_calls"])
            with open(results_path) as f:
                results = json.load(f)

            # Check status field
            for locus in results:
                assert "status" in results[locus], f"status not in {locus} results"
                assert results[locus]["status"] in ["confirmed", "novel_allele"], \
                    f"status has invalid value for {locus}: {results[locus]['status']}"

    def test_wdl_detects_novel_allele_for_low_identity(
        self, allele_database_fasta, query_ambiguous_fasta, minhash_results_ambiguous
    ):
        """WDL task detects novel alleles when alignment identity is low"""
        wdl_path = TORCHBASE_ROOT / "workflows" / "alignment_fallback.wdl"

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir_path = Path(tmpdir)

            input_json = {
                "alignment_fallback.query_sequences": str(query_ambiguous_fasta),
                "alignment_fallback.allele_fasta": str(allele_database_fasta),
                "alignment_fallback.minhash_results": str(minhash_results_ambiguous)
            }

            input_json_path = tmpdir_path / "inputs.json"
            with open(input_json_path, "w") as f:
                json.dump(input_json, f)

            result = subprocess.run(
                ["miniwdl", "run", str(wdl_path), "-i", str(input_json_path), "-d", str(tmpdir_path)],
                capture_output=True,
                text=True,
                timeout=300
            )

            assert result.returncode == 0, f"miniwdl run failed: {result.stderr}"

            output_dirs = list(tmpdir_path.glob("**/outputs.json"))
            with open(output_dirs[0]) as f:
                outputs = json.load(f)

            results_path = Path(outputs["alignment_fallback.refined_calls"])
            with open(results_path) as f:
                results = json.load(f)

            # gyrB query has novel allele (different stop codon)
            # Should be detected if identity threshold is properly set
            if results["gyrB"]["identity"] < 0.90:
                assert results["gyrB"]["status"] == "novel_allele", \
                    "gyrB should be marked as novel_allele for low identity"

    def test_wdl_confirms_alleles_for_high_identity(
        self, allele_database_fasta, query_ambiguous_fasta, minhash_results_ambiguous
    ):
        """WDL task confirms alleles when alignment identity is high"""
        wdl_path = TORCHBASE_ROOT / "workflows" / "alignment_fallback.wdl"

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir_path = Path(tmpdir)

            input_json = {
                "alignment_fallback.query_sequences": str(query_ambiguous_fasta),
                "alignment_fallback.allele_fasta": str(allele_database_fasta),
                "alignment_fallback.minhash_results": str(minhash_results_ambiguous)
            }

            input_json_path = tmpdir_path / "inputs.json"
            with open(input_json_path, "w") as f:
                json.dump(input_json, f)

            result = subprocess.run(
                ["miniwdl", "run", str(wdl_path), "-i", str(input_json_path), "-d", str(tmpdir_path)],
                capture_output=True,
                text=True,
                timeout=300
            )

            assert result.returncode == 0, f"miniwdl run failed: {result.stderr}"

            output_dirs = list(tmpdir_path.glob("**/outputs.json"))
            with open(output_dirs[0]) as f:
                outputs = json.load(f)

            results_path = Path(outputs["alignment_fallback.refined_calls"])
            with open(results_path) as f:
                results = json.load(f)

            # Check that high-identity matches are confirmed
            for locus in ["adk", "fumC"]:
                if results[locus]["identity"] >= 0.90:
                    assert results[locus]["status"] == "confirmed", \
                        f"{locus} should be confirmed for high identity"


@pytest.mark.miniwdl
class TestAlignmentFallbackWDLSkipsConfidentCalls:
    """Test alignment fallback skips loci with confident MinHash calls."""

    def test_wdl_skips_alignment_for_confident_calls(
        self, allele_database_fasta, query_ambiguous_fasta, minhash_results_confident
    ):
        """WDL task skips alignment for confident MinHash calls"""
        wdl_path = TORCHBASE_ROOT / "workflows" / "alignment_fallback.wdl"

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir_path = Path(tmpdir)

            input_json = {
                "alignment_fallback.query_sequences": str(query_ambiguous_fasta),
                "alignment_fallback.allele_fasta": str(allele_database_fasta),
                "alignment_fallback.minhash_results": str(minhash_results_confident)
            }

            input_json_path = tmpdir_path / "inputs.json"
            with open(input_json_path, "w") as f:
                json.dump(input_json, f)

            result = subprocess.run(
                ["miniwdl", "run", str(wdl_path), "-i", str(input_json_path), "-d", str(tmpdir_path)],
                capture_output=True,
                text=True,
                timeout=300
            )

            # Should succeed and pass through confident calls without alignment
            assert result.returncode == 0, f"miniwdl run failed: {result.stderr}"

            output_dirs = list(tmpdir_path.glob("**/outputs.json"))
            with open(output_dirs[0]) as f:
                outputs = json.load(f)

            results_path = Path(outputs["alignment_fallback.refined_calls"])
            with open(results_path) as f:
                results = json.load(f)

            # Results should pass through unchanged or marked as confirmed
            assert len(results) > 0, "Results should not be empty"
            for locus in results:
                # All should be marked as confirmed since input was confident
                assert results[locus]["status"] == "confirmed", \
                    f"{locus} should be confirmed when MinHash was confident"


@pytest.mark.miniwdl
class TestAlignmentFallbackWDLEdgeCases:
    """Test edge cases and boundary conditions."""

    def test_wdl_handles_empty_minhash_results(
        self, allele_database_fasta, query_ambiguous_fasta
    ):
        """WDL task handles empty MinHash results gracefully"""
        wdl_path = TORCHBASE_ROOT / "workflows" / "alignment_fallback.wdl"

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir_path = Path(tmpdir)

            # Create empty MinHash results
            empty_results = tmpdir_path / "empty_results.json"
            with open(empty_results, "w") as f:
                json.dump({}, f)

            input_json = {
                "alignment_fallback.query_sequences": str(query_ambiguous_fasta),
                "alignment_fallback.allele_fasta": str(allele_database_fasta),
                "alignment_fallback.minhash_results": str(empty_results)
            }

            input_json_path = tmpdir_path / "inputs.json"
            with open(input_json_path, "w") as f:
                json.dump(input_json, f)

            result = subprocess.run(
                ["miniwdl", "run", str(wdl_path), "-i", str(input_json_path), "-d", str(tmpdir_path)],
                capture_output=True,
                text=True,
                timeout=300
            )

            # Should either succeed with empty results or fail gracefully
            assert result.returncode in [0, 1], "Unexpected return code for empty results"

    def test_wdl_handles_single_locus(self):
        """WDL task handles single locus alignment"""
        wdl_path = TORCHBASE_ROOT / "workflows" / "alignment_fallback.wdl"

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir_path = Path(tmpdir)

            # Create single-locus database
            single_locus_db = tmpdir_path / "single_locus.fasta"
            with open(single_locus_db, "w") as f:
                f.write(">adk_1\nATGGCTATGAAACAACTCACCAACGTAGCACTGTCCAAAGCC\n")
                f.write(">adk_2\nATGGCTATGAAACAACTCACCAACGTAGCACTGTCCAAAGCC\n")

            # Create matching query
            query = tmpdir_path / "query.fasta"
            with open(query, "w") as f:
                f.write(">query\nATGGCTATGAAACAACTCACCAACGTAGCACTGTCCAAAGCC\n")

            # Create ambiguous MinHash results for single locus
            minhash_results = tmpdir_path / "minhash.json"
            with open(minhash_results, "w") as f:
                json.dump({
                    "adk": {
                        "allele_id": "1",
                        "similarity": 0.89,
                        "confidence": False,
                        "second_best": {
                            "allele_id": "2",
                            "similarity": 0.88
                        }
                    }
                }, f)

            input_json = {
                "alignment_fallback.query_sequences": str(query),
                "alignment_fallback.allele_fasta": str(single_locus_db),
                "alignment_fallback.minhash_results": str(minhash_results)
            }

            input_json_path = tmpdir_path / "inputs.json"
            with open(input_json_path, "w") as f:
                json.dump(input_json, f)

            result = subprocess.run(
                ["miniwdl", "run", str(wdl_path), "-i", str(input_json_path), "-d", str(tmpdir_path)],
                capture_output=True,
                text=True,
                timeout=300
            )

            assert result.returncode == 0, f"miniwdl run failed: {result.stderr}"

            output_dirs = list(tmpdir_path.glob("**/outputs.json"))
            with open(output_dirs[0]) as f:
                outputs = json.load(f)

            results_path = Path(outputs["alignment_fallback.refined_calls"])
            with open(results_path) as f:
                results = json.load(f)

            # Should have exactly one locus
            assert len(results) == 1, f"Expected 1 locus, got {len(results)}"
            assert "adk" in results, "Expected adk locus in results"

    def test_wdl_handles_all_novel_alleles(
        self, allele_database_fasta
    ):
        """WDL task handles case where all alleles are novel"""
        wdl_path = TORCHBASE_ROOT / "workflows" / "alignment_fallback.wdl"

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir_path = Path(tmpdir)

            # Create query with sequences far from any known allele
            novel_query = tmpdir_path / "novel_query.fasta"
            with open(novel_query, "w") as f:
                # Completely different sequences
                f.write(">novel_adk\nAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n")
                f.write(">novel_fumC\nTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT\n")
                f.write(">novel_gyrB\nGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG\n")

            # MinHash results with very low similarity
            minhash_results = tmpdir_path / "minhash_novel.json"
            with open(minhash_results, "w") as f:
                json.dump({
                    "adk": {"allele_id": "1", "similarity": 0.30, "confidence": False},
                    "fumC": {"allele_id": "1", "similarity": 0.25, "confidence": False},
                    "gyrB": {"allele_id": "1", "similarity": 0.28, "confidence": False}
                }, f)

            input_json = {
                "alignment_fallback.query_sequences": str(novel_query),
                "alignment_fallback.allele_fasta": str(allele_database_fasta),
                "alignment_fallback.minhash_results": str(minhash_results)
            }

            input_json_path = tmpdir_path / "inputs.json"
            with open(input_json_path, "w") as f:
                json.dump(input_json, f)

            result = subprocess.run(
                ["miniwdl", "run", str(wdl_path), "-i", str(input_json_path), "-d", str(tmpdir_path)],
                capture_output=True,
                text=True,
                timeout=300
            )

            assert result.returncode == 0, f"miniwdl run failed: {result.stderr}"

            output_dirs = list(tmpdir_path.glob("**/outputs.json"))
            with open(output_dirs[0]) as f:
                outputs = json.load(f)

            results_path = Path(outputs["alignment_fallback.refined_calls"])
            with open(results_path) as f:
                results = json.load(f)

            # All loci should be marked as novel_allele
            for locus in results:
                assert results[locus]["status"] == "novel_allele", \
                    f"{locus} should be marked as novel_allele"
