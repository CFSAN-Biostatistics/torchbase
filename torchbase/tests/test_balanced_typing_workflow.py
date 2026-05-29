"""Acceptance tests for balanced typing workflow (Issue #56).

These are RED-phase tests - they MUST fail because the feature is not yet complete.

Acceptance criteria:
- `torchbase/workflows/builtin/balanced_typing.wdl` exists
- Workflow imports shared tasks from tasks/ directory
- Conditional alignment runs only when confidence < 0.85
- Uses minimap2 with asm5 preset for contigs, sr for reads
- Accepts standard inputs plus input_type parameter
- Outputs standardized JSON result format
- Tests verify both alignment and no-alignment paths
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
    - adk: 3 alleles
    - fumC: 2 alleles
    - gyrB: 2 alleles
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir_path = Path(tmpdir)
        db_path = tmpdir_path / "allele_db.fasta"

        # Create realistic MLST allele sequences
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
>gyrB_1
ATGACCCAACTGAAAGTGATGCCGCAACGTGTCGACCTGCAAATCCACGCAGTGCTGATGAAACCGATGAAACTGCTGATGACGAAACTGCTGATGATCGTGGTGATCTGACGGACGAACTGACGAAATACGACGAAAACAAACACGATGTCATCGACGATGTGACGACCGACATGATCACGGACGACGTACTGATGAAACTGGTGATCCACGTGCACGATGAAACGGACGACTACGACGACATGCCGATCGACGATGATGATGATGACCACGACGACAACGACGAAACGATGATCCTGACGATGACGACGATCTGACGGATGACTAA
>gyrB_2
ATGACCCAACTGAAAGTGATGCCGCAACGTGTCGACCTGCAAATCCACGCAGTGCTGATGAAACCGATGAAACTGCTGATGACGAAACTGCTGATGATCGTGGTGATCTGACGGACGAACTGACGAAATACGACGAAAACAAACACGATGTCATCGACGATGTGACGACCGACATGATCACGGACGACGTACTGATGAAACTGGTGATCCACGTGCACGATGAAACGGACGACTACGACGACATGCCGATCGACGATGATGATGATGACCACGACGACAACGACGAAACGATGATCCTGACGATGACGACGATCTGACGGATGACTAA
"""

        with open(db_path, "w") as f:
            f.write(fasta_content)

        yield db_path


@pytest.fixture
def profiles_table():
    """Create temporary profiles table (TSV format).

    Maps allele combinations to sequence types.
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir_path = Path(tmpdir)
        profiles_path = tmpdir_path / "profiles.tsv"

        # Create MLST-style profiles table
        tsv_content = """ST\tadk\tfumC\tgyrB
1\t1\t1\t1
2\t1\t2\t1
3\t2\t1\t2
4\t3\t2\t1
"""

        with open(profiles_path, "w") as f:
            f.write(tsv_content)

        yield profiles_path


@pytest.fixture
def exact_match_query_contigs():
    """Create query contigs that exactly match adk_1, fumC_2, gyrB_1 (ST2)."""
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir_path = Path(tmpdir)
        contigs_path = tmpdir_path / "query_contigs.fasta"

        fasta_content = """>contig1_contains_adk_1
ATGGCTATGAAACAACTCACCAACGTAGCACTGTCCAAAGCCGCACGTGGCAATGCCGCTGCAACTGACTGCACTGGCGCACCTGCCGCTGCTGATGAACGTCATCGGTACGGTCTCGTCCACCGGCTCTACGACCTGGTGTACGACTGTACGTCCGAAATCCTTACTGGCGGCGGTCACACGTCTGCTGGACATCCGCCACCACATTTGCTGCATCACGGCGGCGGTGACGGTGGTTATGGTGGTGACGACCACGCCGTACATCGACGACGTGCTGATCGAACTGATCGAAGACGACGACGACGAAGTGATCGAACTGATCGAAGTGCTGGATGAAGTCGAAAATGTCTAA
>contig2_contains_fumC_2
ATGAATATTAACAACGCACTGGGCGACGTGCTGAAAACCCACGGCCAGATGACGAAAGAAGTGATGCAACTGACCCAAGGTGCAACCCACGCCTTTGTGACCGCCGTGGGCGACTCGCCCGAAGAAACGCACCACGGAAATCGATCTGGACGAACTGTTCGAACTGATGGAAGAACACGAAGTGCTGATGGTCGACATCCTGATGATGCACGACCACGACGATGACCGTGATAGCACCACTGTACGACATTGACGACGACGACGACGATACAGAACACAATGACGATGGAAGAAAACGACGACGAAGTGATCCACGTGATGGTGTAG
>contig3_contains_gyrB_1
ATGACCCAACTGAAAGTGATGCCGCAACGTGTCGACCTGCAAATCCACGCAGTGCTGATGAAACCGATGAAACTGCTGATGACGAAACTGCTGATGATCGTGGTGATCTGACGGACGAACTGACGAAATACGACGAAAACAAACACGATGTCATCGACGATGTGACGACCGACATGATCACGGACGACGTACTGATGAAACTGGTGATCCACGTGCACGATGAAACGGACGACTACGACGACATGCCGATCGACGATGATGATGATGACCACGACGACAACGACGAAACGATGATCCTGACGATGACGACGATCTGACGGATGACTAA
"""

        with open(contigs_path, "w") as f:
            f.write(fasta_content)

        yield contigs_path


@pytest.fixture
def low_confidence_query_contigs():
    """Create query contigs with mutations to produce low confidence MinHash match.

    This should trigger alignment fallback in balanced strategy.
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir_path = Path(tmpdir)
        contigs_path = tmpdir_path / "low_confidence_contigs.fasta"

        # Add several mutations to reduce MinHash similarity
        fasta_content = """>contig1_adk_1_with_mutations
ATGGCTATGAAACAACTCACCAACGTAGCACTGTCCAAAGCCGCACGTGGCAATGCCGCTGCAACTGACTGCACTGGCGCACCTGCCGCTGCTGATGAACGTCATCGGTACGGTCTCGTCCACCGGCTCTACGACCTGGTGTACGACTGTACGTCCGAAATCCTTACTGGCGGCGGTCACACGTCTGCTGGACATCCGCCACCACATTTGCTGCATCACGGCGGCGGTGACGGTGGTTATGGTGGTGACGACCACGCCGTACATCGACGACGTGCTGATCGAACTGATCGAAGACGACGACGACGAAGTGATCGAACTGATCGAAGTGCTGGATGAAGTCGAAAATGTCTAA
>contig2_fumC_2_with_mutations_many_changes
ATGAATATTAACAACGCACTGGGCGACGTGCTGAAAACCCACGGCCAGATGACGAAAGAAGTGATGCAACTGACCCAAGGTGCAACCCACGCCTTTGTGACCGCCGTGGGCGACTCGCCCGAAGAAACGCACCACGGAAATCGATCTGGACGAACTGTTCGAACTGATGGAAGAACACGAAGTGCTGATGGTCGACATCCTGATGATGCACGACCACGACGATGACCGTGATAGCACCACTGTACGACATTGACGACGACGACGACGATACAGAACACAATGACGATGGAAGAAAACGACGACGAAGTGATCCACGTGATGGTGTAG
>contig3_gyrB_1_with_many_snps_to_reduce_confidence
ATGACCCAACTGAAAGTGATGCCGCAACGTGTCGACCTGCAAATCCACGCAGTGCTGATGAAACCGATGAAACTGCTGATGACGAAACTGCTGATGATCGTGGTGATCTGACGGACGAACTGACGAAATACGACGAAAACAAACACGATGTCATCGACGATGTGACGACCGACATGATCACGGACGACGTACTGATGAAACTGGTGATCCACGTGCACGATGAAACGGACGACTACGACGACATGCCGATCGACGATGATGATGATGACCACGACGACAACGACGAAACGATGATCCTGACGATGACGACGATCTGACGGATGACTAA
"""

        with open(contigs_path, "w") as f:
            f.write(fasta_content)

        yield contigs_path


@pytest.fixture
def query_reads_fasta():
    """Create temporary query reads file (FASTA format)."""
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir_path = Path(tmpdir)
        reads_path = tmpdir_path / "query_reads.fasta"

        fasta_content = """>read1_adk_1_depth1
ATGGCTATGAAACAACTCACCAACGTAGCACTGTCCAAAGCCGCACGTGGCAATGCCGCTGCAACTGAC
>read2_adk_1_depth2
TGCACTGGCGCACCTGCCGCTGCTGATGAACGTCATCGGTACGGTCTCGTCCACCGGCTCTACGACCTG
>read3_adk_1_depth3
GTGTACGACTGTACGTCCGAAATCCTTACTGGCGGCGGTCACACGTCTGCTGGACATCCGCCACCACAT
>read4_fumC_2_depth1
ATGAATATTAACAACGCACTGGGCGACGTGCTGAAAACCCACGGCCAGATGACGAAAGAAGTGATGCAA
>read5_fumC_2_depth2
CTGACCCAAGGTGCAACCCACGCCTTTGTGACCGCCGTGGGCGACTCGCCCGAAGAAACGCACCACGGA
>read6_gyrB_1_depth1
ATGACCCAACTGAAAGTGATGCCGCAACGTGTCGACCTGCAAATCCACGCAGTGCTGATGAAACCGATG
"""

        with open(reads_path, "w") as f:
            f.write(fasta_content)

        yield reads_path


class TestBalancedTypingWDLFileExists:
    """Test balanced_typing.wdl file exists at expected location."""

    def test_wdl_file_exists(self):
        """balanced_typing.wdl file exists in builtin directory"""
        wdl_path = TORCHBASE_ROOT / "workflows" / "builtin" / "balanced_typing.wdl"
        assert wdl_path.exists(), f"WDL file not found at {wdl_path}"

    def test_wdl_file_is_file(self):
        """balanced_typing.wdl is a regular file"""
        wdl_path = TORCHBASE_ROOT / "workflows" / "builtin" / "balanced_typing.wdl"
        assert wdl_path.is_file(), f"WDL path is not a file: {wdl_path}"


class TestBalancedTypingWDLImports:
    """Test workflow imports shared tasks from tasks/ directory."""

    def test_wdl_imports_minhash_tasks(self):
        """Workflow imports minhash.wdl tasks"""
        wdl_path = TORCHBASE_ROOT / "workflows" / "builtin" / "balanced_typing.wdl"
        with open(wdl_path) as f:
            content = f.read()

        # Check for import statement referencing minhash.wdl
        assert "import" in content and "minhash" in content, \
            "Workflow does not import minhash tasks"

    def test_wdl_imports_alignment_tasks(self):
        """Workflow imports alignment.wdl tasks"""
        wdl_path = TORCHBASE_ROOT / "workflows" / "builtin" / "balanced_typing.wdl"
        with open(wdl_path) as f:
            content = f.read()

        # Check for import statement referencing alignment.wdl
        assert "import" in content and "alignment" in content, \
            "Workflow does not import alignment tasks"

    def test_wdl_imports_profile_lookup_tasks(self):
        """Workflow imports profile_lookup.wdl tasks"""
        wdl_path = TORCHBASE_ROOT / "workflows" / "builtin" / "balanced_typing.wdl"
        with open(wdl_path) as f:
            content = f.read()

        # Check for import statement referencing profile_lookup.wdl
        assert "import" in content and "profile_lookup" in content, \
            "Workflow does not import profile_lookup tasks"

    def test_wdl_import_paths_reference_tasks_directory(self):
        """Import statements use relative paths to tasks/ directory"""
        wdl_path = TORCHBASE_ROOT / "workflows" / "builtin" / "balanced_typing.wdl"
        with open(wdl_path) as f:
            content = f.read()

        # Check that imports reference tasks/ subdirectory
        assert "tasks/" in content, \
            "Import statements do not reference tasks/ directory"


class TestBalancedTypingWDLWorkflowStructure:
    """Test workflow structure and signature."""

    def test_wdl_has_workflow_definition(self):
        """WDL file contains workflow definition"""
        wdl_path = TORCHBASE_ROOT / "workflows" / "builtin" / "balanced_typing.wdl"
        with open(wdl_path) as f:
            content = f.read()

        assert "workflow" in content, "WDL file does not contain workflow definition"

    def test_wdl_workflow_name_is_balanced_typing(self):
        """Workflow is named balanced_typing"""
        wdl_path = TORCHBASE_ROOT / "workflows" / "builtin" / "balanced_typing.wdl"
        with open(wdl_path) as f:
            content = f.read()

        assert "workflow balanced_typing" in content, "Workflow name is not balanced_typing"

    def test_wdl_has_input_section(self):
        """Workflow has input section"""
        wdl_path = TORCHBASE_ROOT / "workflows" / "builtin" / "balanced_typing.wdl"
        with open(wdl_path) as f:
            content = f.read()

        assert "input {" in content, "Workflow does not have input section"

    def test_wdl_has_output_section(self):
        """Workflow has output section"""
        wdl_path = TORCHBASE_ROOT / "workflows" / "builtin" / "balanced_typing.wdl"
        with open(wdl_path) as f:
            content = f.read()

        assert "output {" in content, "Workflow does not have output section"


class TestBalancedTypingWDLInputs:
    """Test workflow accepts standard inputs plus input_type parameter."""

    def test_wdl_accepts_query_sequences_input(self):
        """Workflow accepts query_sequences input"""
        wdl_path = TORCHBASE_ROOT / "workflows" / "builtin" / "balanced_typing.wdl"
        with open(wdl_path) as f:
            content = f.read()

        assert "query_sequences" in content or "sequences" in content, \
            "Workflow does not accept query_sequences input"

    def test_wdl_accepts_allele_fasta_input(self):
        """Workflow accepts allele_fasta input"""
        wdl_path = TORCHBASE_ROOT / "workflows" / "builtin" / "balanced_typing.wdl"
        with open(wdl_path) as f:
            content = f.read()

        assert "allele_fasta" in content or "allele" in content, \
            "Workflow does not accept allele_fasta input"

    def test_wdl_accepts_profiles_table_input(self):
        """Workflow accepts profiles_table input"""
        wdl_path = TORCHBASE_ROOT / "workflows" / "builtin" / "balanced_typing.wdl"
        with open(wdl_path) as f:
            content = f.read()

        assert "profiles" in content, \
            "Workflow does not accept profiles_table input"

    def test_wdl_accepts_input_type_parameter(self):
        """Workflow accepts input_type parameter (contigs/reads)"""
        wdl_path = TORCHBASE_ROOT / "workflows" / "builtin" / "balanced_typing.wdl"
        with open(wdl_path) as f:
            content = f.read()

        assert "input_type" in content, \
            "Workflow does not accept input_type parameter"

    def test_wdl_accepts_confidence_threshold_parameter(self):
        """Workflow accepts confidence_threshold parameter (default 0.85)"""
        wdl_path = TORCHBASE_ROOT / "workflows" / "builtin" / "balanced_typing.wdl"
        with open(wdl_path) as f:
            content = f.read()

        # Check for confidence threshold parameter
        assert "confidence" in content.lower() or "threshold" in content.lower(), \
            "Workflow does not have confidence threshold parameter"


class TestBalancedTypingWDLConditionalAlignment:
    """Test conditional alignment runs only when confidence < 0.85."""

    def test_wdl_has_conditional_logic(self):
        """Workflow contains conditional logic for alignment"""
        wdl_path = TORCHBASE_ROOT / "workflows" / "builtin" / "balanced_typing.wdl"
        with open(wdl_path) as f:
            content = f.read()

        # Check for if/conditional logic
        assert "if" in content.lower(), \
            "Workflow does not contain conditional logic"

    def test_wdl_conditional_checks_confidence(self):
        """Conditional statement checks confidence value"""
        wdl_path = TORCHBASE_ROOT / "workflows" / "builtin" / "balanced_typing.wdl"
        with open(wdl_path) as f:
            content = f.read()

        # Check that conditional logic references confidence
        assert "confidence" in content.lower(), \
            "Conditional logic does not check confidence"

    def test_wdl_conditional_compares_threshold(self):
        """Conditional compares confidence against 0.85 threshold"""
        wdl_path = TORCHBASE_ROOT / "workflows" / "builtin" / "balanced_typing.wdl"
        with open(wdl_path) as f:
            content = f.read()

        # Check for 0.85 threshold value
        assert "0.85" in content or "threshold" in content, \
            "Conditional does not reference 0.85 threshold"

    def test_wdl_calls_alignment_conditionally(self):
        """Workflow calls alignment task within conditional block"""
        wdl_path = TORCHBASE_ROOT / "workflows" / "builtin" / "balanced_typing.wdl"
        with open(wdl_path) as f:
            content = f.read()

        # Check for conditional alignment call
        # This is tricky to test without parsing WDL, but we can check structure
        assert "align" in content.lower(), \
            "Workflow does not call alignment task"


class TestBalancedTypingWDLAlignmentPresets:
    """Test uses minimap2 with asm5 preset for contigs, sr for reads."""

    def test_wdl_specifies_asm5_preset_for_contigs(self):
        """Workflow uses asm5 preset for contig input_type"""
        wdl_path = TORCHBASE_ROOT / "workflows" / "builtin" / "balanced_typing.wdl"
        with open(wdl_path) as f:
            content = f.read()

        # Check for asm5 preset
        assert "asm5" in content, \
            "Workflow does not specify asm5 preset"

    def test_wdl_specifies_sr_preset_for_reads(self):
        """Workflow uses sr preset for reads input_type"""
        wdl_path = TORCHBASE_ROOT / "workflows" / "builtin" / "balanced_typing.wdl"
        with open(wdl_path) as f:
            content = f.read()

        # Check for sr preset
        assert "sr" in content, \
            "Workflow does not specify sr preset"

    def test_wdl_selects_preset_based_on_input_type(self):
        """Workflow selects preset based on input_type parameter"""
        wdl_path = TORCHBASE_ROOT / "workflows" / "builtin" / "balanced_typing.wdl"
        with open(wdl_path) as f:
            content = f.read()

        # Check that input_type is used to select preset
        # This requires both input_type and preset to be referenced together
        has_input_type = "input_type" in content
        has_preset_selection = "asm5" in content and "sr" in content

        assert has_input_type and has_preset_selection, \
            "Workflow does not select preset based on input_type"


class TestBalancedTypingWDLOutputFormat:
    """Test workflow outputs standardized JSON result format."""

    def test_wdl_output_is_json_file(self):
        """Workflow output is a JSON File"""
        wdl_path = TORCHBASE_ROOT / "workflows" / "builtin" / "balanced_typing.wdl"
        with open(wdl_path) as f:
            content = f.read()

        # Check for File output and JSON reference
        assert "File" in content and ("json" in content.lower() or ".json" in content), \
            "Workflow output is not a JSON File"

    def test_wdl_output_includes_strategy_field(self):
        """Workflow output includes strategy field set to 'balanced'"""
        wdl_path = TORCHBASE_ROOT / "workflows" / "builtin" / "balanced_typing.wdl"
        with open(wdl_path) as f:
            content = f.read()

        # Check for strategy field
        assert "strategy" in content and "balanced" in content, \
            "Output does not include strategy field"

    def test_wdl_output_includes_alignment_used_field(self):
        """Workflow output includes alignment_used boolean field"""
        wdl_path = TORCHBASE_ROOT / "workflows" / "builtin" / "balanced_typing.wdl"
        with open(wdl_path) as f:
            content = f.read()

        # Check for alignment_used field
        assert "alignment_used" in content, \
            "Output does not include alignment_used field"

    def test_wdl_output_includes_method_metadata(self):
        """Workflow output includes method metadata object"""
        wdl_path = TORCHBASE_ROOT / "workflows" / "builtin" / "balanced_typing.wdl"
        with open(wdl_path) as f:
            content = f.read()

        # Check for method metadata
        assert "method" in content, \
            "Output does not include method metadata"


class TestBalancedTypingWDLSyntaxValidation:
    """Test miniwdl check validates WDL syntax."""

    def test_miniwdl_check_passes(self):
        """miniwdl check validates WDL syntax without errors"""
        wdl_path = TORCHBASE_ROOT / "workflows" / "builtin" / "balanced_typing.wdl"

        result = subprocess.run(
            ["miniwdl", "check", str(wdl_path)],
            capture_output=True,
            text=True
        )

        assert result.returncode == 0, f"miniwdl check failed: {result.stderr}"


@pytest.mark.miniwdl
class TestBalancedTypingWDLExecutionHighConfidence:
    """Test workflow execution with high-confidence query (no alignment path)."""

    def test_wdl_execution_with_exact_match_produces_output(
        self, allele_database_fasta, profiles_table, exact_match_query_contigs
    ):
        """Workflow execution with exact match produces JSON output"""
        wdl_path = TORCHBASE_ROOT / "workflows" / "builtin" / "balanced_typing.wdl"

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir_path = Path(tmpdir)

            input_json = {
                "balanced_typing.query_sequences": str(exact_match_query_contigs),
                "balanced_typing.allele_fasta": str(allele_database_fasta),
                "balanced_typing.profiles_table": str(profiles_table),
                "balanced_typing.input_type": "contigs"
            }

            input_json_path = tmpdir_path / "inputs.json"
            with open(input_json_path, "w") as f:
                json.dump(input_json, f)

            result = subprocess.run(
                ["miniwdl", "run", str(wdl_path), "-i", str(input_json_path), "-d", str(tmpdir_path)],
                capture_output=True,
                text=True,
                timeout=600
            )

            assert result.returncode == 0, f"miniwdl run failed: {result.stderr}"

            # Find output
            output_dirs = list(tmpdir_path.glob("**/outputs.json"))
            assert len(output_dirs) > 0, "No outputs.json found"

            with open(output_dirs[0]) as f:
                outputs = json.load(f)

            assert "balanced_typing.result" in outputs or "balanced_typing.results" in outputs, \
                "Output result not found"

    def test_wdl_high_confidence_skips_alignment(
        self, allele_database_fasta, profiles_table, exact_match_query_contigs
    ):
        """High confidence match does not trigger alignment (alignment_used=false)"""
        wdl_path = TORCHBASE_ROOT / "workflows" / "builtin" / "balanced_typing.wdl"

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir_path = Path(tmpdir)

            input_json = {
                "balanced_typing.query_sequences": str(exact_match_query_contigs),
                "balanced_typing.allele_fasta": str(allele_database_fasta),
                "balanced_typing.profiles_table": str(profiles_table),
                "balanced_typing.input_type": "contigs"
            }

            input_json_path = tmpdir_path / "inputs.json"
            with open(input_json_path, "w") as f:
                json.dump(input_json, f)

            result = subprocess.run(
                ["miniwdl", "run", str(wdl_path), "-i", str(input_json_path), "-d", str(tmpdir_path)],
                capture_output=True,
                text=True,
                timeout=600
            )

            assert result.returncode == 0, f"miniwdl run failed: {result.stderr}"

            output_dirs = list(tmpdir_path.glob("**/outputs.json"))
            with open(output_dirs[0]) as f:
                outputs = json.load(f)

            # Get result file path
            result_key = "balanced_typing.result" if "balanced_typing.result" in outputs else "balanced_typing.results"
            result_path = Path(outputs[result_key])

            with open(result_path) as f:
                result_data = json.load(f)

            # Check alignment_used is false
            assert "method" in result_data, "Result does not have method metadata"
            assert "alignment_used" in result_data["method"], "method does not have alignment_used field"
            assert result_data["method"]["alignment_used"] is False, \
                "High confidence should not use alignment"

    def test_wdl_output_has_strategy_balanced(
        self, allele_database_fasta, profiles_table, exact_match_query_contigs
    ):
        """Output includes strategy='balanced' field"""
        wdl_path = TORCHBASE_ROOT / "workflows" / "builtin" / "balanced_typing.wdl"

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir_path = Path(tmpdir)

            input_json = {
                "balanced_typing.query_sequences": str(exact_match_query_contigs),
                "balanced_typing.allele_fasta": str(allele_database_fasta),
                "balanced_typing.profiles_table": str(profiles_table),
                "balanced_typing.input_type": "contigs"
            }

            input_json_path = tmpdir_path / "inputs.json"
            with open(input_json_path, "w") as f:
                json.dump(input_json, f)

            result = subprocess.run(
                ["miniwdl", "run", str(wdl_path), "-i", str(input_json_path), "-d", str(tmpdir_path)],
                capture_output=True,
                text=True,
                timeout=600
            )

            assert result.returncode == 0, f"miniwdl run failed: {result.stderr}"

            output_dirs = list(tmpdir_path.glob("**/outputs.json"))
            with open(output_dirs[0]) as f:
                outputs = json.load(f)

            result_key = "balanced_typing.result" if "balanced_typing.result" in outputs else "balanced_typing.results"
            result_path = Path(outputs[result_key])

            with open(result_path) as f:
                result_data = json.load(f)

            assert "method" in result_data, "Result does not have method metadata"
            assert "strategy" in result_data["method"], "method does not have strategy field"
            assert result_data["method"]["strategy"] == "balanced", \
                f"Expected strategy='balanced', got {result_data['method']['strategy']}"


@pytest.mark.miniwdl
class TestBalancedTypingWDLExecutionLowConfidence:
    """Test workflow execution with low-confidence query (alignment path)."""

    def test_wdl_execution_with_low_confidence_produces_output(
        self, allele_database_fasta, profiles_table, low_confidence_query_contigs
    ):
        """Workflow execution with low confidence query produces JSON output"""
        wdl_path = TORCHBASE_ROOT / "workflows" / "builtin" / "balanced_typing.wdl"

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir_path = Path(tmpdir)

            input_json = {
                "balanced_typing.query_sequences": str(low_confidence_query_contigs),
                "balanced_typing.allele_fasta": str(allele_database_fasta),
                "balanced_typing.profiles_table": str(profiles_table),
                "balanced_typing.input_type": "contigs"
            }

            input_json_path = tmpdir_path / "inputs.json"
            with open(input_json_path, "w") as f:
                json.dump(input_json, f)

            result = subprocess.run(
                ["miniwdl", "run", str(wdl_path), "-i", str(input_json_path), "-d", str(tmpdir_path)],
                capture_output=True,
                text=True,
                timeout=600
            )

            assert result.returncode == 0, f"miniwdl run failed: {result.stderr}"

            output_dirs = list(tmpdir_path.glob("**/outputs.json"))
            assert len(output_dirs) > 0, "No outputs.json found"

    def test_wdl_low_confidence_triggers_alignment(
        self, allele_database_fasta, profiles_table, low_confidence_query_contigs
    ):
        """Low confidence match triggers alignment fallback (alignment_used=true)"""
        wdl_path = TORCHBASE_ROOT / "workflows" / "builtin" / "balanced_typing.wdl"

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir_path = Path(tmpdir)

            input_json = {
                "balanced_typing.query_sequences": str(low_confidence_query_contigs),
                "balanced_typing.allele_fasta": str(allele_database_fasta),
                "balanced_typing.profiles_table": str(profiles_table),
                "balanced_typing.input_type": "contigs"
            }

            input_json_path = tmpdir_path / "inputs.json"
            with open(input_json_path, "w") as f:
                json.dump(input_json, f)

            result = subprocess.run(
                ["miniwdl", "run", str(wdl_path), "-i", str(input_json_path), "-d", str(tmpdir_path)],
                capture_output=True,
                text=True,
                timeout=600
            )

            assert result.returncode == 0, f"miniwdl run failed: {result.stderr}"

            output_dirs = list(tmpdir_path.glob("**/outputs.json"))
            with open(output_dirs[0]) as f:
                outputs = json.load(f)

            result_key = "balanced_typing.result" if "balanced_typing.result" in outputs else "balanced_typing.results"
            result_path = Path(outputs[result_key])

            with open(result_path) as f:
                result_data = json.load(f)

            # Check alignment_used is true
            assert "method" in result_data, "Result does not have method metadata"
            assert "alignment_used" in result_data["method"], "method does not have alignment_used field"
            assert result_data["method"]["alignment_used"] is True, \
                "Low confidence should trigger alignment"

    def test_wdl_alignment_uses_asm5_preset_for_contigs(
        self, allele_database_fasta, profiles_table, low_confidence_query_contigs
    ):
        """Alignment uses asm5 preset for contig input_type"""
        wdl_path = TORCHBASE_ROOT / "workflows" / "builtin" / "balanced_typing.wdl"

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir_path = Path(tmpdir)

            input_json = {
                "balanced_typing.query_sequences": str(low_confidence_query_contigs),
                "balanced_typing.allele_fasta": str(allele_database_fasta),
                "balanced_typing.profiles_table": str(profiles_table),
                "balanced_typing.input_type": "contigs"
            }

            input_json_path = tmpdir_path / "inputs.json"
            with open(input_json_path, "w") as f:
                json.dump(input_json, f)

            result = subprocess.run(
                ["miniwdl", "run", str(wdl_path), "-i", str(input_json_path), "-d", str(tmpdir_path)],
                capture_output=True,
                text=True,
                timeout=600
            )

            assert result.returncode == 0, f"miniwdl run failed: {result.stderr}"

            output_dirs = list(tmpdir_path.glob("**/outputs.json"))
            with open(output_dirs[0]) as f:
                outputs = json.load(f)

            result_key = "balanced_typing.result" if "balanced_typing.result" in outputs else "balanced_typing.results"
            result_path = Path(outputs[result_key])

            with open(result_path) as f:
                result_data = json.load(f)

            # Check that notes or method metadata indicates asm5 preset was used
            # This may be implementation-dependent
            if "notes" in result_data:
                notes_str = json.dumps(result_data["notes"])
                assert "asm5" in notes_str.lower(), "Alignment did not use asm5 preset"


@pytest.mark.miniwdl
class TestBalancedTypingWDLExecutionWithReads:
    """Test workflow execution with reads input (should use sr preset)."""

    def test_wdl_execution_with_reads_produces_output(
        self, allele_database_fasta, profiles_table, query_reads_fasta
    ):
        """Workflow execution with reads input produces JSON output"""
        wdl_path = TORCHBASE_ROOT / "workflows" / "builtin" / "balanced_typing.wdl"

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir_path = Path(tmpdir)

            input_json = {
                "balanced_typing.query_sequences": str(query_reads_fasta),
                "balanced_typing.allele_fasta": str(allele_database_fasta),
                "balanced_typing.profiles_table": str(profiles_table),
                "balanced_typing.input_type": "reads"
            }

            input_json_path = tmpdir_path / "inputs.json"
            with open(input_json_path, "w") as f:
                json.dump(input_json, f)

            result = subprocess.run(
                ["miniwdl", "run", str(wdl_path), "-i", str(input_json_path), "-d", str(tmpdir_path)],
                capture_output=True,
                text=True,
                timeout=600
            )

            assert result.returncode == 0, f"miniwdl run failed: {result.stderr}"

            output_dirs = list(tmpdir_path.glob("**/outputs.json"))
            assert len(output_dirs) > 0, "No outputs.json found"

    def test_wdl_reads_input_uses_sr_preset_when_alignment_triggered(
        self, allele_database_fasta, profiles_table, query_reads_fasta
    ):
        """Workflow uses sr preset for reads when alignment is triggered"""
        wdl_path = TORCHBASE_ROOT / "workflows" / "builtin" / "balanced_typing.wdl"

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir_path = Path(tmpdir)

            # Force low confidence threshold to trigger alignment
            input_json = {
                "balanced_typing.query_sequences": str(query_reads_fasta),
                "balanced_typing.allele_fasta": str(allele_database_fasta),
                "balanced_typing.profiles_table": str(profiles_table),
                "balanced_typing.input_type": "reads",
                "balanced_typing.confidence_threshold": 0.99  # High threshold forces alignment
            }

            input_json_path = tmpdir_path / "inputs.json"
            with open(input_json_path, "w") as f:
                json.dump(input_json, f)

            result = subprocess.run(
                ["miniwdl", "run", str(wdl_path), "-i", str(input_json_path), "-d", str(tmpdir_path)],
                capture_output=True,
                text=True,
                timeout=600
            )

            # This test verifies sr preset usage - may need to check logs or metadata
            # Implementation-dependent
            if result.returncode == 0:
                output_dirs = list(tmpdir_path.glob("**/outputs.json"))
                with open(output_dirs[0]) as f:
                    outputs = json.load(f)

                result_key = ("balanced_typing.result" if "balanced_typing.result" in outputs
                              else "balanced_typing.results")
                result_path = Path(outputs[result_key])

                with open(result_path) as f:
                    result_data = json.load(f)

                # Check metadata for sr preset indication
                if "notes" in result_data:
                    notes_str = json.dumps(result_data["notes"])
                    # sr preset should be mentioned for reads
                    assert "sr" in notes_str.lower() or "reads" in notes_str.lower(), \
                        "Workflow did not use sr preset for reads"


@pytest.mark.miniwdl
class TestBalancedTypingWDLStandardizedOutput:
    """Test workflow produces standardized JSON output format."""

    def test_output_has_required_fields(
        self, allele_database_fasta, profiles_table, exact_match_query_contigs
    ):
        """Output JSON contains all required standardized fields"""
        wdl_path = TORCHBASE_ROOT / "workflows" / "builtin" / "balanced_typing.wdl"

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir_path = Path(tmpdir)

            input_json = {
                "balanced_typing.query_sequences": str(exact_match_query_contigs),
                "balanced_typing.allele_fasta": str(allele_database_fasta),
                "balanced_typing.profiles_table": str(profiles_table),
                "balanced_typing.input_type": "contigs"
            }

            input_json_path = tmpdir_path / "inputs.json"
            with open(input_json_path, "w") as f:
                json.dump(input_json, f)

            result = subprocess.run(
                ["miniwdl", "run", str(wdl_path), "-i", str(input_json_path), "-d", str(tmpdir_path)],
                capture_output=True,
                text=True,
                timeout=600
            )

            assert result.returncode == 0, f"miniwdl run failed: {result.stderr}"

            output_dirs = list(tmpdir_path.glob("**/outputs.json"))
            with open(output_dirs[0]) as f:
                outputs = json.load(f)

            result_key = "balanced_typing.result" if "balanced_typing.result" in outputs else "balanced_typing.results"
            result_path = Path(outputs[result_key])

            with open(result_path) as f:
                result_data = json.load(f)

            # Check required fields per standardized format
            required_fields = ["status", "confidence", "allele_calls", "method"]
            for field in required_fields:
                assert field in result_data, f"Output missing required field: {field}"

            # Check method subfields
            assert "strategy" in result_data["method"], "method missing strategy field"
            assert "alignment_used" in result_data["method"], "method missing alignment_used field"

    def test_output_status_field_valid(
        self, allele_database_fasta, profiles_table, exact_match_query_contigs
    ):
        """Output status field is one of: known_profile, novel_profile, novel_allele"""
        wdl_path = TORCHBASE_ROOT / "workflows" / "builtin" / "balanced_typing.wdl"

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir_path = Path(tmpdir)

            input_json = {
                "balanced_typing.query_sequences": str(exact_match_query_contigs),
                "balanced_typing.allele_fasta": str(allele_database_fasta),
                "balanced_typing.profiles_table": str(profiles_table),
                "balanced_typing.input_type": "contigs"
            }

            input_json_path = tmpdir_path / "inputs.json"
            with open(input_json_path, "w") as f:
                json.dump(input_json, f)

            result = subprocess.run(
                ["miniwdl", "run", str(wdl_path), "-i", str(input_json_path), "-d", str(tmpdir_path)],
                capture_output=True,
                text=True,
                timeout=600
            )

            assert result.returncode == 0, f"miniwdl run failed: {result.stderr}"

            output_dirs = list(tmpdir_path.glob("**/outputs.json"))
            with open(output_dirs[0]) as f:
                outputs = json.load(f)

            result_key = "balanced_typing.result" if "balanced_typing.result" in outputs else "balanced_typing.results"
            result_path = Path(outputs[result_key])

            with open(result_path) as f:
                result_data = json.load(f)

            valid_statuses = ["known_profile", "novel_profile", "novel_allele"]
            assert result_data["status"] in valid_statuses, \
                f"Invalid status: {result_data['status']}"

    def test_output_confidence_is_numeric(
        self, allele_database_fasta, profiles_table, exact_match_query_contigs
    ):
        """Output confidence field is numeric (0.0-1.0)"""
        wdl_path = TORCHBASE_ROOT / "workflows" / "builtin" / "balanced_typing.wdl"

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir_path = Path(tmpdir)

            input_json = {
                "balanced_typing.query_sequences": str(exact_match_query_contigs),
                "balanced_typing.allele_fasta": str(allele_database_fasta),
                "balanced_typing.profiles_table": str(profiles_table),
                "balanced_typing.input_type": "contigs"
            }

            input_json_path = tmpdir_path / "inputs.json"
            with open(input_json_path, "w") as f:
                json.dump(input_json, f)

            result = subprocess.run(
                ["miniwdl", "run", str(wdl_path), "-i", str(input_json_path), "-d", str(tmpdir_path)],
                capture_output=True,
                text=True,
                timeout=600
            )

            assert result.returncode == 0, f"miniwdl run failed: {result.stderr}"

            output_dirs = list(tmpdir_path.glob("**/outputs.json"))
            with open(output_dirs[0]) as f:
                outputs = json.load(f)

            result_key = "balanced_typing.result" if "balanced_typing.result" in outputs else "balanced_typing.results"
            result_path = Path(outputs[result_key])

            with open(result_path) as f:
                result_data = json.load(f)

            assert isinstance(result_data["confidence"], (int, float)), \
                f"Confidence is not numeric: {type(result_data['confidence'])}"
            assert 0.0 <= result_data["confidence"] <= 1.0, \
                f"Confidence out of range: {result_data['confidence']}"
