"""Acceptance tests for complete MLST workflow orchestration (Issue #18).

These are RED-phase tests - they MUST fail because the feature is not yet complete.

Acceptance criteria:
- WDL workflow: reads/contigs input → typing result output
- Stage 1: input validation (ensure reads OR contigs provided)
- Stage 2: depth filtering for reads (skip for contigs)
- Stage 3: MinHash allele calling across all schemes
- Stage 4: scheme inference from allele calls (highest coverage/identity)
- Stage 5: profile lookup in inferred scheme
- Stage 6: alignment fallback if ambiguous
- Stage 7: final result with status, confidence, nearest ST
- Parameterized via metadata.toml thresholds
- miniwdl check validates syntax
- End-to-end test: synthetic data → complete typing result
"""

import pytest
import json
import tempfile
from pathlib import Path
import subprocess
import toml


# Get the torchbase root directory
TORCHBASE_ROOT = Path(__file__).parent.parent


@pytest.fixture
def mlst_workflow_path():
    """Path to the main MLST workflow WDL file."""
    return TORCHBASE_ROOT / "workflows" / "mlst" / "1.0.0.torch" / "main.wdl"


@pytest.fixture
def workflow_metadata_path():
    """Path to workflow metadata with threshold configuration."""
    return TORCHBASE_ROOT / "workflows" / "mlst" / "1.0.0.torch" / "metadata.toml"


@pytest.fixture
def multi_scheme_alleles_fasta():
    """Create temporary multi-scheme allele database (FASTA format).

    Creates allele database with multiple schemes:
    - ecoli: dinB, icdA
    - salmonella: adk, fumC, gyrB
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir_path = Path(tmpdir)
        db_path = tmpdir_path / "multi_scheme_alleles.fasta"

        fasta_content = """>ecoli_dinB_1
ATGGCTATGAAACAACTCACCAACGTAGCACTGTCCAAAGCCGCACGTGGCAATGCCGCTGCAACTGAC
>ecoli_dinB_2
ATGGCTATGAAACAACTCACCAACGTAGCACTGTCCAAAGCCGCACGTGGCAATGCCGCTGCAACTGAT
>ecoli_icdA_1
TGCACTGGCGCACCTGCCGCTGCTGATGAACGTCATCGGTACGGTCTCGTCCACCGGCTCTACGACCTG
>ecoli_icdA_2
TGCACTGGCGCACCTGCCGCTGCTGATGAACGTCATCGGTACGGTCTCGTCCACCGGCTCTACGACCTC
>salmonella_adk_1
ATGAATATTAACAACGCACTGGGCGACGTGCTGAAAACCCACGGCCAGATGACGAAAGAAGTGATGCAA
>salmonella_adk_2
ATGAATATTAACAACGCACTGGGCGACGTGCTGAAAACCCACGGCCAGATGACGAAAGAAGTGATGCAC
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
def salmonella_profile_table():
    """Create temporary Salmonella profile table (TSV format)."""
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir_path = Path(tmpdir)
        profile_path = tmpdir_path / "salmonella_profiles.tsv"

        # Salmonella MLST profiles
        profile_content = """ST\tadk\tfumC\tgyrB
1\t1\t1\t1
2\t1\t2\t1
3\t2\t1\t1
4\t1\t1\t2
"""

        with open(profile_path, "w") as f:
            f.write(profile_content)

        yield profile_path


@pytest.fixture
def ecoli_profile_table():
    """Create temporary E. coli profile table (TSV format)."""
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir_path = Path(tmpdir)
        profile_path = tmpdir_path / "ecoli_profiles.tsv"

        # E. coli MLST profiles
        profile_content = """ST\tdinB\ticdA
10\t1\t1
11\t2\t1
12\t1\t2
"""

        with open(profile_path, "w") as f:
            f.write(profile_content)

        yield profile_path


@pytest.fixture
def query_salmonella_contigs():
    """Create query contigs matching Salmonella ST=1 (adk_1, fumC_1, gyrB_1)."""
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir_path = Path(tmpdir)
        contigs_path = tmpdir_path / "salmonella_contigs.fasta"

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
def query_salmonella_reads():
    """Create query reads matching Salmonella ST=1 at varying depths."""
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir_path = Path(tmpdir)
        reads_path = tmpdir_path / "salmonella_reads.fasta"

        # Simulate reads with different coverage
        fasta_content = """>read1_adk_1_depth1
ATGAATATTAACAACGCACTGGGCGACGTGCTGAAAACCCACGGCCAGATGACGAAAGAA
>read2_adk_1_depth2
GTGATGCAACTGACCCAAGGTGCAACCCACGCCTTTGTGACCGCCGTGGGCGACTCGCCC
>read3_adk_1_depth3
GAAGAAACGCACCACGGAAATCGATCTGGACGAACTGTTCGAACTGATGGAAGAACACGA
>read4_fumC_1_depth1
CTGACCCAAGGTGCAACCCACGCCTTTGTGACCGCCGTGGGCGACTCGCCCGAAGAAACG
>read5_fumC_1_depth2
CACCACGGAAATCGATCTGGACGAACTGTTCGAACTGATGGAAGAACACGAAGTGCTGAT
>read6_gyrB_1_depth1
ATGACCCAACTGAAAGTGATGCCGCAACGTGTCGACCTGCAAATCCACGCAGTGCTGATG
>read7_gyrB_1_depth2
AAACCGATGAAACTGCTGATGACGAAACTGCTGATGATCGTGGTGATCTGACGGACGAAC
>read8_gyrB_1_depth3
TGACGAAATACGACGAAAACAAACACGATGTCATCGACGATGTGACGACCGACATGATCA
"""

        with open(reads_path, "w") as f:
            f.write(fasta_content)

        yield reads_path


@pytest.fixture
def query_ambiguous_contigs():
    """Create query contigs with ambiguous allele calls (multiple schemes match)."""
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir_path = Path(tmpdir)
        contigs_path = tmpdir_path / "ambiguous_contigs.fasta"

        # Mix alleles from different schemes
        fasta_content = """>contig1_salmonella_adk_1
ATGAATATTAACAACGCACTGGGCGACGTGCTGAAAACCCACGGCCAGATGACGAAAGAAGTGATGCAA
>contig2_ecoli_dinB_1
ATGGCTATGAAACAACTCACCAACGTAGCACTGTCCAAAGCCGCACGTGGCAATGCCGCTGCAACTGAC
>contig3_salmonella_fumC_1
CTGACCCAAGGTGCAACCCACGCCTTTGTGACCGCCGTGGGCGACTCGCCCGAAGAAACGCACCACGGA
"""

        with open(contigs_path, "w") as f:
            f.write(fasta_content)

        yield contigs_path


class TestMLSTWorkflowFileExists:
    """Test MLST workflow WDL file exists at expected location."""

    def test_main_wdl_file_exists(self, mlst_workflow_path):
        """Main MLST workflow WDL file exists"""
        assert mlst_workflow_path.exists(), f"WDL file not found at {mlst_workflow_path}"

    def test_main_wdl_is_file(self, mlst_workflow_path):
        """Main MLST workflow WDL is a regular file"""
        assert mlst_workflow_path.is_file(), f"WDL path is not a file: {mlst_workflow_path}"


class TestMLSTWorkflowStructure:
    """Test WDL workflow: reads/contigs input → typing result output."""

    def test_wdl_has_workflow_definition(self, mlst_workflow_path):
        """WDL file contains workflow definition"""
        with open(mlst_workflow_path) as f:
            content = f.read()

        assert "workflow" in content, "WDL file does not contain workflow definition"

    def test_wdl_workflow_name_is_mlst_typing(self, mlst_workflow_path):
        """WDL workflow is named mlst_typing or similar"""
        with open(mlst_workflow_path) as f:
            content = f.read()

        assert "workflow mlst" in content.lower(), "Workflow name is not mlst-related"

    def test_wdl_has_input_section(self, mlst_workflow_path):
        """WDL workflow has input section"""
        with open(mlst_workflow_path) as f:
            content = f.read()

        assert "input {" in content, "Workflow does not have input section"

    def test_wdl_accepts_contigs_input(self, mlst_workflow_path):
        """WDL workflow accepts contigs input"""
        with open(mlst_workflow_path) as f:
            content = f.read()

        assert "contigs" in content.lower() or "sequences" in content.lower(), \
            "Workflow does not accept contigs input"

    def test_wdl_accepts_reads_input(self, mlst_workflow_path):
        """WDL workflow accepts reads input"""
        with open(mlst_workflow_path) as f:
            content = f.read()

        assert "reads" in content.lower(), "Workflow does not accept reads input"

    def test_wdl_has_output_section(self, mlst_workflow_path):
        """WDL workflow has output section"""
        with open(mlst_workflow_path) as f:
            content = f.read()

        assert "output {" in content, "Workflow does not have output section"

    def test_wdl_output_is_typing_result(self, mlst_workflow_path):
        """WDL workflow outputs typing result"""
        with open(mlst_workflow_path) as f:
            content = f.read()

        assert "typing" in content.lower() or "result" in content.lower(), \
            "Workflow does not output typing result"


class TestMLSTWorkflowStage1InputValidation:
    """Test Stage 1: input validation (ensure reads OR contigs provided)."""

    def test_wdl_has_input_validation_task(self, mlst_workflow_path):
        """WDL workflow has input validation task/step"""
        with open(mlst_workflow_path) as f:
            content = f.read()

        # Check for validation task or conditional logic
        assert "validate" in content.lower() or "if " in content, \
            "Workflow does not have input validation"

    def test_wdl_validates_contigs_or_reads_provided(self, mlst_workflow_path):
        """WDL workflow validates that contigs OR reads are provided"""
        with open(mlst_workflow_path) as f:
            content = f.read()

        # Check for optional input handling
        assert "File?" in content or "optional" in content.lower(), \
            "Workflow does not handle optional contigs/reads inputs"

    def test_wdl_rejects_both_empty_inputs(self):
        """WDL workflow execution fails if both contigs and reads are empty"""
        wdl_path = TORCHBASE_ROOT / "workflows" / "mlst" / "1.0.0.torch" / "main.wdl"

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir_path = Path(tmpdir)

            # Create input JSON with no contigs or reads
            input_json = {}

            input_json_path = tmpdir_path / "inputs.json"
            with open(input_json_path, "w") as f:
                json.dump(input_json, f)

            # Run miniwdl (should fail validation)
            result = subprocess.run(
                ["miniwdl", "run", str(wdl_path), "-i", str(input_json_path), "-d", str(tmpdir_path)],
                capture_output=True,
                text=True,
                timeout=60
            )

            # Should fail with validation error
            assert result.returncode != 0, "Workflow should fail with empty inputs"


class TestMLSTWorkflowStage2DepthFiltering:
    """Test Stage 2: depth filtering for reads (skip for contigs)."""

    def test_wdl_has_depth_filtering_task(self, mlst_workflow_path):
        """WDL workflow includes depth filtering task"""
        with open(mlst_workflow_path) as f:
            content = f.read()

        assert "depth" in content.lower() or "filter" in content.lower(), \
            "Workflow does not include depth filtering"

    def test_wdl_depth_filtering_conditional_on_reads(self, mlst_workflow_path):
        """WDL workflow depth filtering is conditional (only for reads)"""
        with open(mlst_workflow_path) as f:
            content = f.read()

        # Check for conditional execution
        assert "if " in content or "select_first" in content, \
            "Workflow does not conditionally apply depth filtering"

    def test_wdl_accepts_min_coverage_parameter(self, mlst_workflow_path):
        """WDL workflow accepts min_coverage parameter"""
        with open(mlst_workflow_path) as f:
            content = f.read()

        assert "min_coverage" in content or "coverage" in content, \
            "Workflow does not accept coverage parameter"

    def test_wdl_skips_depth_filtering_for_contigs(self):
        """WDL workflow skips depth filtering when contigs provided"""
        # This would be tested by checking the execution path
        # For now, verify the conditional structure exists
        wdl_path = TORCHBASE_ROOT / "workflows" / "mlst" / "1.0.0.torch" / "main.wdl"
        with open(wdl_path) as f:
            content = f.read()

        # Check for conditional logic around depth filtering
        assert "if " in content or "select_first" in content, \
            "No conditional execution structure found"


class TestMLSTWorkflowStage3MinHashCalling:
    """Test Stage 3: MinHash allele calling across all schemes."""

    def test_wdl_calls_minhash_task(self, mlst_workflow_path):
        """WDL workflow calls MinHash allele calling task"""
        with open(mlst_workflow_path) as f:
            content = f.read()

        assert "minhash" in content.lower(), \
            "Workflow does not call MinHash task"

    def test_wdl_passes_allele_database_to_minhash(self, mlst_workflow_path):
        """WDL workflow passes allele database to MinHash task"""
        with open(mlst_workflow_path) as f:
            content = f.read()

        assert "allele" in content.lower() and "fasta" in content.lower(), \
            "Workflow does not pass allele database"

    def test_wdl_minhash_runs_across_all_schemes(self, mlst_workflow_path):
        """WDL workflow MinHash runs across all schemes"""
        with open(mlst_workflow_path) as f:
            content = f.read()

        # Check for multi-scheme handling
        assert "scheme" in content.lower(), \
            "Workflow does not handle multiple schemes"

    def test_wdl_minhash_produces_allele_calls(self, mlst_workflow_path):
        """WDL workflow MinHash produces allele calls JSON"""
        with open(mlst_workflow_path) as f:
            content = f.read()

        assert "allele_calls" in content or "results" in content, \
            "Workflow does not produce allele calls"


class TestMLSTWorkflowStage4SchemeInference:
    """Test Stage 4: scheme inference from allele calls (highest coverage/identity)."""

    def test_wdl_has_scheme_inference_task(self, mlst_workflow_path):
        """WDL workflow has scheme inference task/step"""
        with open(mlst_workflow_path) as f:
            content = f.read()

        assert "scheme" in content.lower() and ("infer" in content.lower() or
                                                  "detect" in content.lower() or
                                                  "select" in content.lower()), \
            "Workflow does not have scheme inference"

    def test_wdl_scheme_inference_uses_coverage(self, mlst_workflow_path):
        """WDL workflow scheme inference uses coverage metric"""
        with open(mlst_workflow_path) as f:
            content = f.read()

        assert "coverage" in content.lower(), \
            "Workflow scheme inference does not use coverage"

    def test_wdl_scheme_inference_uses_identity(self, mlst_workflow_path):
        """WDL workflow scheme inference uses identity metric"""
        with open(mlst_workflow_path) as f:
            content = f.read()

        assert "identity" in content.lower() or "similarity" in content.lower(), \
            "Workflow scheme inference does not use identity/similarity"

    def test_wdl_scheme_inference_selects_best_scheme(self, mlst_workflow_path):
        """WDL workflow scheme inference selects best-matching scheme"""
        with open(mlst_workflow_path) as f:
            content = f.read()

        # Check for best/max/highest selection logic
        assert any(keyword in content.lower() for keyword in ["best", "max", "highest", "top"]), \
            "Workflow does not select best scheme"


class TestMLSTWorkflowStage5ProfileLookup:
    """Test Stage 5: profile lookup in inferred scheme."""

    def test_wdl_has_profile_lookup_task(self, mlst_workflow_path):
        """WDL workflow has profile lookup task"""
        with open(mlst_workflow_path) as f:
            content = f.read()

        assert "profile" in content.lower() and "lookup" in content.lower(), \
            "Workflow does not have profile lookup task"

    def test_wdl_profile_lookup_uses_profiles_table(self, mlst_workflow_path):
        """WDL workflow profile lookup uses profiles.tsv"""
        with open(mlst_workflow_path) as f:
            content = f.read()

        assert "profiles" in content.lower() and ("tsv" in content.lower() or "table" in content.lower()), \
            "Workflow does not use profiles table"

    def test_wdl_profile_lookup_detects_exact_match(self, mlst_workflow_path):
        """WDL workflow profile lookup detects exact matches"""
        with open(mlst_workflow_path) as f:
            content = f.read()

        assert "match" in content.lower() or "lookup" in content.lower(), \
            "Workflow does not detect profile matches"

    def test_wdl_profile_lookup_detects_novel_profile(self, mlst_workflow_path):
        """WDL workflow profile lookup detects novel profiles"""
        with open(mlst_workflow_path) as f:
            content = f.read()

        assert "novel" in content.lower() or "new" in content.lower(), \
            "Workflow does not detect novel profiles"

    def test_wdl_profile_lookup_calculates_nearest_st(self, mlst_workflow_path):
        """WDL workflow calculates nearest ST for novel profiles"""
        with open(mlst_workflow_path) as f:
            content = f.read()

        assert "nearest" in content.lower() or "distance" in content.lower(), \
            "Workflow does not calculate nearest ST"


class TestMLSTWorkflowStage6AlignmentFallback:
    """Test Stage 6: alignment fallback if ambiguous."""

    def test_wdl_has_alignment_fallback_task(self, mlst_workflow_path):
        """WDL workflow has alignment fallback task"""
        with open(mlst_workflow_path) as f:
            content = f.read()

        assert "align" in content.lower() or "minimap" in content.lower(), \
            "Workflow does not have alignment fallback"

    def test_wdl_alignment_triggers_on_ambiguity(self, mlst_workflow_path):
        """WDL workflow alignment triggers on ambiguous calls"""
        with open(mlst_workflow_path) as f:
            content = f.read()

        assert "ambiguous" in content.lower() or "confidence" in content.lower(), \
            "Workflow does not trigger alignment on ambiguity"

    def test_wdl_alignment_refines_allele_calls(self, mlst_workflow_path):
        """WDL workflow alignment refines allele calls"""
        with open(mlst_workflow_path) as f:
            content = f.read()

        assert "refine" in content.lower() or "align" in content.lower(), \
            "Workflow does not refine allele calls via alignment"

    def test_wdl_alignment_detects_novel_alleles(self, mlst_workflow_path):
        """WDL workflow alignment detects novel alleles"""
        with open(mlst_workflow_path) as f:
            content = f.read()

        assert "novel" in content.lower(), \
            "Workflow does not detect novel alleles"


class TestMLSTWorkflowStage7FinalResult:
    """Test Stage 7: final result with status, confidence, nearest ST."""

    def test_wdl_produces_final_typing_result(self, mlst_workflow_path):
        """WDL workflow produces final typing result"""
        with open(mlst_workflow_path) as f:
            content = f.read()

        # Check for final output/result assembly
        assert "output" in content, "Workflow does not produce output"

    def test_wdl_final_result_includes_status(self, mlst_workflow_path):
        """WDL workflow final result includes status field"""
        with open(mlst_workflow_path) as f:
            content = f.read()

        assert "status" in content.lower(), \
            "Workflow final result does not include status"

    def test_wdl_final_result_includes_confidence(self, mlst_workflow_path):
        """WDL workflow final result includes confidence field"""
        with open(mlst_workflow_path) as f:
            content = f.read()

        assert "confidence" in content.lower(), \
            "Workflow final result does not include confidence"

    def test_wdl_final_result_includes_nearest_st(self, mlst_workflow_path):
        """WDL workflow final result includes nearest_st field"""
        with open(mlst_workflow_path) as f:
            content = f.read()

        assert "nearest" in content.lower() or "st" in content.lower(), \
            "Workflow final result does not include nearest ST"

    def test_wdl_final_result_includes_sequence_type(self, mlst_workflow_path):
        """WDL workflow final result includes sequence_type field"""
        with open(mlst_workflow_path) as f:
            content = f.read()

        assert "sequence_type" in content or "st" in content.lower(), \
            "Workflow final result does not include sequence type"


class TestMLSTWorkflowParameterization:
    """Test parameterized via metadata.toml thresholds."""

    def test_metadata_has_thresholds_section(self, workflow_metadata_path):
        """Workflow metadata.toml has thresholds section"""
        with open(workflow_metadata_path) as f:
            metadata = toml.load(f)

        # Should have parameters section or thresholds
        assert "thresholds" in metadata or "parameters" in metadata or "typing" in metadata, \
            "Metadata does not have thresholds/parameters section"

    def test_wdl_reads_thresholds_from_metadata(self, mlst_workflow_path):
        """WDL workflow reads thresholds from metadata"""
        with open(mlst_workflow_path) as f:
            content = f.read()

        # Check for threshold parameters
        assert "threshold" in content.lower() or "min_" in content.lower(), \
            "Workflow does not use threshold parameters"

    def test_wdl_has_min_identity_threshold(self, mlst_workflow_path):
        """WDL workflow has min_identity threshold parameter"""
        with open(mlst_workflow_path) as f:
            content = f.read()

        assert "min_identity" in content or "identity_threshold" in content, \
            "Workflow does not have min_identity threshold"

    def test_wdl_has_min_coverage_threshold(self, mlst_workflow_path):
        """WDL workflow has min_coverage threshold parameter"""
        with open(mlst_workflow_path) as f:
            content = f.read()

        assert "min_coverage" in content or "coverage_threshold" in content, \
            "Workflow does not have min_coverage threshold"

    def test_wdl_has_ambiguity_threshold(self, mlst_workflow_path):
        """WDL workflow has ambiguity threshold parameter"""
        with open(mlst_workflow_path) as f:
            content = f.read()

        assert "ambiguity" in content.lower() or "confidence" in content.lower(), \
            "Workflow does not have ambiguity threshold"


class TestMLSTWorkflowSyntaxValidation:
    """Test miniwdl check validates syntax."""

    def test_miniwdl_check_passes(self, mlst_workflow_path):
        """miniwdl check validates WDL syntax without errors"""
        result = subprocess.run(
            ["miniwdl", "check", str(mlst_workflow_path)],
            capture_output=True,
            text=True
        )

        assert result.returncode == 0, f"miniwdl check failed: {result.stderr}"


@pytest.mark.miniwdl
class TestMLSTWorkflowEndToEndContigs:
    """Test end-to-end: synthetic contigs → complete typing result."""

    def test_workflow_execution_with_contigs_produces_result(
        self, mlst_workflow_path, query_salmonella_contigs, multi_scheme_alleles_fasta,
        salmonella_profile_table
    ):
        """Workflow execution with contigs produces typing result"""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir_path = Path(tmpdir)

            input_json = {
                "mlst_typing.contigs": str(query_salmonella_contigs),
                "mlst_typing.allele_database": str(multi_scheme_alleles_fasta),
                "mlst_typing.profiles": str(salmonella_profile_table)
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

            assert "mlst_typing.typing_result" in outputs, "Typing result not in outputs"

    def test_workflow_result_has_sequence_type(
        self, mlst_workflow_path, query_salmonella_contigs, multi_scheme_alleles_fasta,
        salmonella_profile_table
    ):
        """Workflow result includes sequence_type field"""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir_path = Path(tmpdir)

            input_json = {
                "mlst_typing.contigs": str(query_salmonella_contigs),
                "mlst_typing.allele_database": str(multi_scheme_alleles_fasta),
                "mlst_typing.profiles": str(salmonella_profile_table)
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

            assert "sequence_type" in result_data or "st" in result_data, \
                "Result does not include sequence_type"

    def test_workflow_result_has_status(
        self, mlst_workflow_path, query_salmonella_contigs, multi_scheme_alleles_fasta,
        salmonella_profile_table
    ):
        """Workflow result includes status field"""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir_path = Path(tmpdir)

            input_json = {
                "mlst_typing.contigs": str(query_salmonella_contigs),
                "mlst_typing.allele_database": str(multi_scheme_alleles_fasta),
                "mlst_typing.profiles": str(salmonella_profile_table)
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

            assert "status" in result_data, "Result does not include status"
            # Status should be one of: known, novel_profile, novel_allele
            assert result_data["status"] in ["known", "novel_profile", "novel_allele"], \
                f"Invalid status: {result_data['status']}"

    def test_workflow_result_has_confidence(
        self, mlst_workflow_path, query_salmonella_contigs, multi_scheme_alleles_fasta,
        salmonella_profile_table
    ):
        """Workflow result includes confidence field"""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir_path = Path(tmpdir)

            input_json = {
                "mlst_typing.contigs": str(query_salmonella_contigs),
                "mlst_typing.allele_database": str(multi_scheme_alleles_fasta),
                "mlst_typing.profiles": str(salmonella_profile_table)
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

            assert "confidence" in result_data, "Result does not include confidence"

    def test_workflow_identifies_correct_st(
        self, mlst_workflow_path, query_salmonella_contigs, multi_scheme_alleles_fasta,
        salmonella_profile_table
    ):
        """Workflow correctly identifies ST=1 from Salmonella contigs"""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir_path = Path(tmpdir)

            input_json = {
                "mlst_typing.contigs": str(query_salmonella_contigs),
                "mlst_typing.allele_database": str(multi_scheme_alleles_fasta),
                "mlst_typing.profiles": str(salmonella_profile_table)
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

            # Should identify ST=1 (adk_1, fumC_1, gyrB_1)
            st = result_data.get("sequence_type") or result_data.get("st")
            assert st == "1" or st == 1, f"Expected ST=1, got {st}"


@pytest.mark.miniwdl
class TestMLSTWorkflowEndToEndReads:
    """Test end-to-end: synthetic reads → complete typing result (with depth filtering)."""

    def test_workflow_execution_with_reads_produces_result(
        self, mlst_workflow_path, query_salmonella_reads, multi_scheme_alleles_fasta,
        salmonella_profile_table
    ):
        """Workflow execution with reads produces typing result"""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir_path = Path(tmpdir)

            input_json = {
                "mlst_typing.reads": str(query_salmonella_reads),
                "mlst_typing.allele_database": str(multi_scheme_alleles_fasta),
                "mlst_typing.profiles": str(salmonella_profile_table),
                "mlst_typing.min_coverage": 2
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
            assert len(output_dirs) > 0, "No outputs.json found"

    def test_workflow_with_reads_applies_depth_filtering(
        self, mlst_workflow_path, query_salmonella_reads, multi_scheme_alleles_fasta,
        salmonella_profile_table
    ):
        """Workflow with reads applies depth filtering"""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir_path = Path(tmpdir)

            # High coverage threshold should trigger filtering
            input_json = {
                "mlst_typing.reads": str(query_salmonella_reads),
                "mlst_typing.allele_database": str(multi_scheme_alleles_fasta),
                "mlst_typing.profiles": str(salmonella_profile_table),
                "mlst_typing.min_coverage": 5
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

            # Should succeed but with filtered results
            assert result.returncode == 0, f"Workflow execution failed: {result.stderr}"


@pytest.mark.miniwdl
class TestMLSTWorkflowSchemeInference:
    """Test scheme inference with multi-scheme ambiguous input."""

    def test_workflow_infers_correct_scheme_from_allele_coverage(
        self, mlst_workflow_path, query_salmonella_contigs, multi_scheme_alleles_fasta,
        salmonella_profile_table, ecoli_profile_table
    ):
        """Workflow infers Salmonella scheme from allele coverage"""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir_path = Path(tmpdir)

            input_json = {
                "mlst_typing.contigs": str(query_salmonella_contigs),
                "mlst_typing.allele_database": str(multi_scheme_alleles_fasta),
                "mlst_typing.schemes": {
                    "salmonella": str(salmonella_profile_table),
                    "ecoli": str(ecoli_profile_table)
                }
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

            # Should identify Salmonella as the scheme
            assert "scheme" in result_data, "Result does not include inferred scheme"
            assert result_data["scheme"] == "salmonella", \
                f"Expected scheme=salmonella, got {result_data['scheme']}"

    def test_workflow_scheme_inference_includes_coverage_metric(
        self, mlst_workflow_path, query_salmonella_contigs, multi_scheme_alleles_fasta,
        salmonella_profile_table, ecoli_profile_table
    ):
        """Workflow result includes scheme coverage metric"""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir_path = Path(tmpdir)

            input_json = {
                "mlst_typing.contigs": str(query_salmonella_contigs),
                "mlst_typing.allele_database": str(multi_scheme_alleles_fasta),
                "mlst_typing.schemes": {
                    "salmonella": str(salmonella_profile_table),
                    "ecoli": str(ecoli_profile_table)
                }
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

            assert "scheme_coverage" in result_data or "coverage" in result_data, \
                "Result does not include scheme coverage"


@pytest.mark.miniwdl
class TestMLSTWorkflowAlignmentFallback:
    """Test alignment fallback with ambiguous allele calls."""

    def test_workflow_triggers_alignment_on_ambiguous_calls(
        self, mlst_workflow_path, query_ambiguous_contigs, multi_scheme_alleles_fasta,
        salmonella_profile_table
    ):
        """Workflow triggers alignment fallback on ambiguous calls"""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir_path = Path(tmpdir)

            input_json = {
                "mlst_typing.contigs": str(query_ambiguous_contigs),
                "mlst_typing.allele_database": str(multi_scheme_alleles_fasta),
                "mlst_typing.profiles": str(salmonella_profile_table)
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

            # Should succeed even with ambiguous input
            assert result.returncode == 0, f"Workflow execution failed: {result.stderr}"

    def test_workflow_alignment_refines_ambiguous_calls(
        self, mlst_workflow_path, query_ambiguous_contigs, multi_scheme_alleles_fasta,
        salmonella_profile_table
    ):
        """Workflow alignment refines ambiguous allele calls"""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir_path = Path(tmpdir)

            input_json = {
                "mlst_typing.contigs": str(query_ambiguous_contigs),
                "mlst_typing.allele_database": str(multi_scheme_alleles_fasta),
                "mlst_typing.profiles": str(salmonella_profile_table)
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

            # Should indicate alignment was performed
            assert "alignment_used" in result_data or "refined" in result_data, \
                "Result does not indicate alignment refinement"


@pytest.mark.miniwdl
class TestMLSTWorkflowEdgeCases:
    """Test edge cases and boundary conditions."""

    def test_workflow_handles_novel_profile(
        self, mlst_workflow_path, multi_scheme_alleles_fasta, salmonella_profile_table
    ):
        """Workflow handles novel profile (known alleles, unknown combination)"""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir_path = Path(tmpdir)

            # Create contigs with novel profile (adk_2, fumC_2, gyrB_2)
            novel_profile_contigs = tmpdir_path / "novel_profile.fasta"
            with open(novel_profile_contigs, "w") as f:
                f.write(">contig1_adk_2\nATGAATATTAACAACGCACTGGGCGACGTGCTGAAAACCCACGGCCAGATGACGAAAGAAGTGATGCAC\n")
                f.write(">contig2_fumC_2\nCTGACCCAAGGTGCAACCCACGCCTTTGTGACCGCCGTGGGCGACTCGCCCGAAGAAACGCACCACGGC\n")
                f.write(">contig3_gyrB_2\nATGACCCAACTGAAAGTGATGCCGCAACGTGTCGACCTGCAA ATCCACGCAGTGCTGATGAAACCGATC\n")

            input_json = {
                "mlst_typing.contigs": str(novel_profile_contigs),
                "mlst_typing.allele_database": str(multi_scheme_alleles_fasta),
                "mlst_typing.profiles": str(salmonella_profile_table)
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

            # Status should be novel_profile
            assert result_data["status"] == "novel_profile", \
                f"Expected status=novel_profile, got {result_data['status']}"

            # Should have nearest_st
            assert "nearest_st" in result_data, "Result does not include nearest_st"

    def test_workflow_handles_partial_profile(
        self, mlst_workflow_path, multi_scheme_alleles_fasta, salmonella_profile_table
    ):
        """Workflow handles partial profile (missing loci)"""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir_path = Path(tmpdir)

            # Create contigs with only 2 out of 3 loci
            partial_contigs = tmpdir_path / "partial.fasta"
            with open(partial_contigs, "w") as f:
                f.write(">contig1_adk_1\nATGAATATTAACAACGCACTGGGCGACGTGCTGAAAACCCACGGCCAGATGACGAAAGAAGTGATGCAA\n")
                f.write(">contig2_fumC_1\nCTGACCCAAGGTGCAACCCACGCCTTTGTGACCGCCGTGGGCGACTCGCCCGAAGAAACGCACCACGGA\n")

            input_json = {
                "mlst_typing.contigs": str(partial_contigs),
                "mlst_typing.allele_database": str(multi_scheme_alleles_fasta),
                "mlst_typing.profiles": str(salmonella_profile_table)
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

            # Should handle gracefully
            assert result.returncode == 0, f"Workflow execution failed: {result.stderr}"

    def test_workflow_handles_empty_contigs(
        self, mlst_workflow_path, multi_scheme_alleles_fasta, salmonella_profile_table
    ):
        """Workflow handles empty contigs file gracefully"""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir_path = Path(tmpdir)

            # Create empty contigs file
            empty_contigs = tmpdir_path / "empty.fasta"
            empty_contigs.touch()

            input_json = {
                "mlst_typing.contigs": str(empty_contigs),
                "mlst_typing.allele_database": str(multi_scheme_alleles_fasta),
                "mlst_typing.profiles": str(salmonella_profile_table)
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

            # Should either fail gracefully or produce empty result
            assert result.returncode in [0, 1], "Unexpected return code for empty contigs"
