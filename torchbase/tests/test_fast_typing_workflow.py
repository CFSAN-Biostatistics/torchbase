"""Acceptance tests for fast typing workflow (Issue #55).

These are RED-phase tests - they MUST fail because the feature is not yet complete.

Acceptance criteria:
- torchbase/workflows/builtin/fast_typing.wdl exists
- Workflow imports shared tasks from tasks/ directory
- Pipeline completes without alignment stage
- Accepts standard inputs (query_sequences, allele_database, profiles_table)
- Outputs standardized JSON result format
- Can be executed via miniwdl independently
- Tests verify fast workflow completes successfully
- Result format includes method: {strategy: "fast", alignment_used: false}
"""

import pytest
import json
import tempfile
from pathlib import Path
import subprocess


# Get the torchbase root directory
TORCHBASE_ROOT = Path(__file__).parent.parent


@pytest.fixture
def fast_workflow_path():
    """Path to the fast typing workflow WDL file."""
    return TORCHBASE_ROOT / "workflows" / "builtin" / "fast_typing.wdl"


@pytest.fixture
def tasks_directory():
    """Path to the shared tasks directory."""
    return TORCHBASE_ROOT / "workflows" / "builtin" / "tasks"


@pytest.fixture
def allele_database_fasta():
    """Create temporary allele database (FASTA format).

    Creates multi-scheme allele database:
    - ecoli: dinB, icdA
    - salmonella: adk, fumC, gyrB
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir_path = Path(tmpdir)
        db_path = tmpdir_path / "alleles.fasta"

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
def profiles_table_tsv():
    """Create temporary profiles table (TSV format)."""
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir_path = Path(tmpdir)
        profile_path = tmpdir_path / "profiles.tsv"

        # Salmonella MLST profiles
        profile_content = """ST\tsalmonella_adk\tsalmonella_fumC\tsalmonella_gyrB
1\t1\t1\t1
2\t1\t2\t1
3\t2\t1\t1
4\t1\t1\t2
"""

        with open(profile_path, "w") as f:
            f.write(profile_content)

        yield profile_path


@pytest.fixture
def query_contigs_salmonella_st1():
    """Create query contigs matching Salmonella ST=1 (adk_1, fumC_1, gyrB_1)."""
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
def query_contigs_novel_profile():
    """Create query contigs with novel profile (adk_2, fumC_2, gyrB_2)."""
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir_path = Path(tmpdir)
        contigs_path = tmpdir_path / "novel_query.fasta"

        fasta_content = """>contig1_adk_2
ATGAATATTAACAACGCACTGGGCGACGTGCTGAAAACCCACGGCCAGATGACGAAAGAAGTGATGCAC
>contig2_fumC_2
CTGACCCAAGGTGCAACCCACGCCTTTGTGACCGCCGTGGGCGACTCGCCCGAAGAAACGCACCACGGC
>contig3_gyrB_2
ATGACCCAACTGAAAGTGATGCCGCAACGTGTCGACCTGCAAATCCACGCAGTGCTGATGAAACCGATC
"""

        with open(contigs_path, "w") as f:
            f.write(fasta_content)

        yield contigs_path


class TestFastWorkflowFileExists:
    """Test fast_typing.wdl exists at expected location."""

    def test_fast_workflow_file_exists(self, fast_workflow_path):
        """fast_typing.wdl file exists in builtin workflows directory"""
        assert fast_workflow_path.exists(), \
            f"Fast workflow WDL not found at {fast_workflow_path}"

    def test_fast_workflow_is_file(self, fast_workflow_path):
        """fast_typing.wdl is a regular file"""
        assert fast_workflow_path.is_file(), \
            f"Fast workflow path is not a file: {fast_workflow_path}"

    def test_builtin_workflows_directory_exists(self):
        """workflows/builtin directory exists"""
        builtin_dir = TORCHBASE_ROOT / "workflows" / "builtin"
        assert builtin_dir.exists(), f"Builtin workflows directory not found at {builtin_dir}"

    def test_tasks_directory_exists(self, tasks_directory):
        """workflows/builtin/tasks directory exists for shared tasks"""
        assert tasks_directory.exists(), \
            f"Tasks directory not found at {tasks_directory}"


class TestFastWorkflowStructure:
    """Test workflow structure and naming."""

    def test_wdl_has_workflow_definition(self, fast_workflow_path):
        """WDL file contains workflow definition"""
        with open(fast_workflow_path) as f:
            content = f.read()

        assert "workflow" in content, "WDL file does not contain workflow definition"

    def test_wdl_workflow_name_is_fast_typing(self, fast_workflow_path):
        """WDL workflow is named fast_typing"""
        with open(fast_workflow_path) as f:
            content = f.read()

        assert "workflow fast_typing" in content, \
            "Workflow is not named fast_typing"

    def test_wdl_has_version_declaration(self, fast_workflow_path):
        """WDL file declares version (should be 1.0 or higher)"""
        with open(fast_workflow_path) as f:
            content = f.read()

        assert "version" in content, "WDL file does not declare version"

    def test_wdl_version_is_1_0_or_higher(self, fast_workflow_path):
        """WDL file declares version 1.0 or higher"""
        with open(fast_workflow_path) as f:
            first_line = f.readline().strip()

        assert "version 1." in first_line or "version 2." in first_line, \
            "WDL version should be 1.0 or higher"


class TestFastWorkflowImports:
    """Test workflow imports shared tasks from tasks/ directory."""

    def test_wdl_imports_minhash_tasks(self, fast_workflow_path):
        """WDL imports minhash tasks"""
        with open(fast_workflow_path) as f:
            content = f.read()

        assert "import" in content and "minhash" in content, \
            "Workflow does not import minhash tasks"

    def test_wdl_imports_from_tasks_directory(self, fast_workflow_path):
        """WDL imports use tasks/ directory"""
        with open(fast_workflow_path) as f:
            content = f.read()

        assert "tasks/" in content, \
            "Workflow does not import from tasks/ directory"

    def test_wdl_imports_profile_lookup(self, fast_workflow_path):
        """WDL imports profile_lookup task"""
        with open(fast_workflow_path) as f:
            content = f.read()

        assert "profile_lookup" in content or "profile" in content.lower(), \
            "Workflow does not import or reference profile_lookup"

    def test_wdl_does_not_import_alignment(self, fast_workflow_path):
        """WDL does not import alignment tasks (fast strategy skips alignment)"""
        with open(fast_workflow_path) as f:
            content = f.read()

        # Should not import alignment.wdl
        assert "import" not in content or "alignment" not in content, \
            "Fast workflow should not import alignment tasks"


class TestFastWorkflowInputs:
    """Test workflow accepts standard inputs."""

    def test_wdl_has_input_section(self, fast_workflow_path):
        """WDL workflow has input section"""
        with open(fast_workflow_path) as f:
            content = f.read()

        assert "input {" in content, "Workflow does not have input section"

    def test_wdl_accepts_query_sequences_input(self, fast_workflow_path):
        """WDL workflow accepts query_sequences input"""
        with open(fast_workflow_path) as f:
            content = f.read()

        assert "query_sequences" in content, \
            "Workflow does not accept query_sequences input"

    def test_wdl_accepts_allele_database_input(self, fast_workflow_path):
        """WDL workflow accepts allele_database input"""
        with open(fast_workflow_path) as f:
            content = f.read()

        assert "allele_database" in content or "allele_fasta" in content, \
            "Workflow does not accept allele_database input"

    def test_wdl_accepts_profiles_table_input(self, fast_workflow_path):
        """WDL workflow accepts profiles_table input"""
        with open(fast_workflow_path) as f:
            content = f.read()

        assert "profiles" in content, \
            "Workflow does not accept profiles_table input"

    def test_wdl_query_sequences_is_file_type(self, fast_workflow_path):
        """WDL query_sequences input is File type"""
        with open(fast_workflow_path) as f:
            content = f.read()

        assert "File" in content, "Workflow inputs are not File type"


class TestFastWorkflowPipeline:
    """Test pipeline structure: MinHash → Allele calling → Profile lookup → Result."""

    def test_wdl_calls_minhash_sketch_task(self, fast_workflow_path):
        """WDL workflow calls MinHash sketch task"""
        with open(fast_workflow_path) as f:
            content = f.read()

        assert "sketch" in content.lower(), \
            "Workflow does not call MinHash sketch task"

    def test_wdl_calls_allele_calling_task(self, fast_workflow_path):
        """WDL workflow calls allele calling task"""
        with open(fast_workflow_path) as f:
            content = f.read()

        has_allele = "allele" in content.lower()
        has_call = "call" in content.lower() or "calling" in content.lower()
        assert has_allele and has_call, "Workflow does not call allele calling task"

    def test_wdl_calls_profile_lookup_task(self, fast_workflow_path):
        """WDL workflow calls profile lookup task"""
        with open(fast_workflow_path) as f:
            content = f.read()

        assert "profile" in content.lower() and "lookup" in content.lower(), \
            "Workflow does not call profile lookup task"

    def test_wdl_does_not_call_alignment_task(self, fast_workflow_path):
        """WDL workflow does not call alignment task (fast strategy)"""
        with open(fast_workflow_path) as f:
            content = f.read()

        # Should not have minimap2 or alignment calls
        assert "minimap" not in content.lower() and "align" not in content.lower(), \
            "Fast workflow should not call alignment tasks"

    def test_wdl_pipeline_is_linear(self, fast_workflow_path):
        """WDL workflow has linear pipeline (no conditionals for alignment fallback)"""
        with open(fast_workflow_path) as f:
            content = f.read()

        # Fast strategy should not have complex branching
        # (No "if" statements for alignment fallback)
        if "if " in content:
            # If present, should not be for alignment fallback
            assert "alignment" not in content.lower() or "minimap" not in content.lower(), \
                "Fast workflow should not have conditional alignment logic"


class TestFastWorkflowOutputs:
    """Test workflow outputs standardized JSON result format."""

    def test_wdl_has_output_section(self, fast_workflow_path):
        """WDL workflow has output section"""
        with open(fast_workflow_path) as f:
            content = f.read()

        assert "output {" in content, "Workflow does not have output section"

    def test_wdl_output_is_typing_result(self, fast_workflow_path):
        """WDL workflow outputs typing result"""
        with open(fast_workflow_path) as f:
            content = f.read()

        assert "typing_result" in content or "result" in content, \
            "Workflow does not output typing result"

    def test_wdl_output_is_file_type(self, fast_workflow_path):
        """WDL workflow output is File type (JSON result)"""
        with open(fast_workflow_path) as f:
            content = f.read()

        # Check output section has File type
        lines = content.split('\n')
        in_output = False
        has_file_output = False

        for line in lines:
            if "output {" in line:
                in_output = True
            elif in_output and "File" in line:
                has_file_output = True
                break
            elif in_output and "}" in line:
                break

        assert has_file_output, "Workflow output is not File type"

    def test_wdl_result_includes_strategy_metadata(self, fast_workflow_path):
        """WDL workflow result includes strategy metadata"""
        with open(fast_workflow_path) as f:
            content = f.read()

        # Should have logic to add strategy: "fast" to result
        assert "fast" in content or "strategy" in content, \
            "Workflow does not include strategy metadata"


class TestFastWorkflowSyntaxValidation:
    """Test miniwdl check validates syntax."""

    def test_miniwdl_check_passes(self, fast_workflow_path):
        """miniwdl check validates WDL syntax without errors"""
        result = subprocess.run(
            ["miniwdl", "check", str(fast_workflow_path)],
            capture_output=True,
            text=True
        )

        assert result.returncode == 0, \
            f"miniwdl check failed: {result.stderr}"


@pytest.mark.miniwdl
class TestFastWorkflowExecution:
    """Test workflow can be executed via miniwdl independently."""

    def test_workflow_executes_successfully(
        self, fast_workflow_path, query_contigs_salmonella_st1,
        allele_database_fasta, profiles_table_tsv
    ):
        """Workflow executes successfully with valid inputs"""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir_path = Path(tmpdir)

            input_json = {
                "fast_typing.query_sequences": str(query_contigs_salmonella_st1),
                "fast_typing.allele_database": str(allele_database_fasta),
                "fast_typing.profiles_table": str(profiles_table_tsv)
            }

            input_json_path = tmpdir_path / "inputs.json"
            with open(input_json_path, "w") as f:
                json.dump(input_json, f)

            result = subprocess.run(
                ["miniwdl", "run", str(fast_workflow_path),
                 "-i", str(input_json_path),
                 "-d", str(tmpdir_path)],
                capture_output=True,
                text=True,
                timeout=600
            )

            assert result.returncode == 0, \
                f"Workflow execution failed: {result.stderr}"

    def test_workflow_produces_outputs_json(
        self, fast_workflow_path, query_contigs_salmonella_st1,
        allele_database_fasta, profiles_table_tsv
    ):
        """Workflow produces outputs.json file"""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir_path = Path(tmpdir)

            input_json = {
                "fast_typing.query_sequences": str(query_contigs_salmonella_st1),
                "fast_typing.allele_database": str(allele_database_fasta),
                "fast_typing.profiles_table": str(profiles_table_tsv)
            }

            input_json_path = tmpdir_path / "inputs.json"
            with open(input_json_path, "w") as f:
                json.dump(input_json, f)

            result = subprocess.run(
                ["miniwdl", "run", str(fast_workflow_path),
                 "-i", str(input_json_path),
                 "-d", str(tmpdir_path)],
                capture_output=True,
                text=True,
                timeout=600
            )

            assert result.returncode == 0, f"Workflow execution failed: {result.stderr}"

            # Find outputs.json
            output_files = list(tmpdir_path.glob("**/outputs.json"))
            assert len(output_files) > 0, "No outputs.json found"

    def test_workflow_output_has_typing_result(
        self, fast_workflow_path, query_contigs_salmonella_st1,
        allele_database_fasta, profiles_table_tsv
    ):
        """Workflow outputs.json contains typing_result"""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir_path = Path(tmpdir)

            input_json = {
                "fast_typing.query_sequences": str(query_contigs_salmonella_st1),
                "fast_typing.allele_database": str(allele_database_fasta),
                "fast_typing.profiles_table": str(profiles_table_tsv)
            }

            input_json_path = tmpdir_path / "inputs.json"
            with open(input_json_path, "w") as f:
                json.dump(input_json, f)

            result = subprocess.run(
                ["miniwdl", "run", str(fast_workflow_path),
                 "-i", str(input_json_path),
                 "-d", str(tmpdir_path)],
                capture_output=True,
                text=True,
                timeout=600
            )

            assert result.returncode == 0, f"Workflow execution failed: {result.stderr}"

            output_files = list(tmpdir_path.glob("**/outputs.json"))
            with open(output_files[0]) as f:
                outputs = json.load(f)

            assert "fast_typing.typing_result" in outputs, \
                "Outputs do not contain typing_result"


@pytest.mark.miniwdl
class TestFastWorkflowResultFormat:
    """Test result format includes required fields and strategy metadata."""

    def test_result_is_valid_json(
        self, fast_workflow_path, query_contigs_salmonella_st1,
        allele_database_fasta, profiles_table_tsv
    ):
        """Result file is valid JSON"""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir_path = Path(tmpdir)

            input_json = {
                "fast_typing.query_sequences": str(query_contigs_salmonella_st1),
                "fast_typing.allele_database": str(allele_database_fasta),
                "fast_typing.profiles_table": str(profiles_table_tsv)
            }

            input_json_path = tmpdir_path / "inputs.json"
            with open(input_json_path, "w") as f:
                json.dump(input_json, f)

            result = subprocess.run(
                ["miniwdl", "run", str(fast_workflow_path),
                 "-i", str(input_json_path),
                 "-d", str(tmpdir_path)],
                capture_output=True,
                text=True,
                timeout=600
            )

            assert result.returncode == 0, f"Workflow execution failed: {result.stderr}"

            output_files = list(tmpdir_path.glob("**/outputs.json"))
            with open(output_files[0]) as f:
                outputs = json.load(f)

            result_path = Path(outputs["fast_typing.typing_result"])
            with open(result_path) as f:
                result_data = json.load(f)

            assert isinstance(result_data, dict), "Result is not a valid JSON object"

    def test_result_has_profile_id(
        self, fast_workflow_path, query_contigs_salmonella_st1,
        allele_database_fasta, profiles_table_tsv
    ):
        """Result includes profile_id field"""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir_path = Path(tmpdir)

            input_json = {
                "fast_typing.query_sequences": str(query_contigs_salmonella_st1),
                "fast_typing.allele_database": str(allele_database_fasta),
                "fast_typing.profiles_table": str(profiles_table_tsv)
            }

            input_json_path = tmpdir_path / "inputs.json"
            with open(input_json_path, "w") as f:
                json.dump(input_json, f)

            result = subprocess.run(
                ["miniwdl", "run", str(fast_workflow_path),
                 "-i", str(input_json_path),
                 "-d", str(tmpdir_path)],
                capture_output=True,
                text=True,
                timeout=600
            )

            assert result.returncode == 0, f"Workflow execution failed: {result.stderr}"

            output_files = list(tmpdir_path.glob("**/outputs.json"))
            with open(output_files[0]) as f:
                outputs = json.load(f)

            result_path = Path(outputs["fast_typing.typing_result"])
            with open(result_path) as f:
                result_data = json.load(f)

            has_profile = "profile_id" in result_data or "st" in result_data
            has_profile = has_profile or "sequence_type" in result_data
            assert has_profile, "Result does not include profile_id or ST"

    def test_result_has_status(
        self, fast_workflow_path, query_contigs_salmonella_st1,
        allele_database_fasta, profiles_table_tsv
    ):
        """Result includes status field"""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir_path = Path(tmpdir)

            input_json = {
                "fast_typing.query_sequences": str(query_contigs_salmonella_st1),
                "fast_typing.allele_database": str(allele_database_fasta),
                "fast_typing.profiles_table": str(profiles_table_tsv)
            }

            input_json_path = tmpdir_path / "inputs.json"
            with open(input_json_path, "w") as f:
                json.dump(input_json, f)

            result = subprocess.run(
                ["miniwdl", "run", str(fast_workflow_path),
                 "-i", str(input_json_path),
                 "-d", str(tmpdir_path)],
                capture_output=True,
                text=True,
                timeout=600
            )

            assert result.returncode == 0, f"Workflow execution failed: {result.stderr}"

            output_files = list(tmpdir_path.glob("**/outputs.json"))
            with open(output_files[0]) as f:
                outputs = json.load(f)

            result_path = Path(outputs["fast_typing.typing_result"])
            with open(result_path) as f:
                result_data = json.load(f)

            assert "status" in result_data, "Result does not include status"
            assert result_data["status"] in ["known", "novel_profile", "novel_allele"], \
                f"Invalid status: {result_data['status']}"

    def test_result_has_confidence(
        self, fast_workflow_path, query_contigs_salmonella_st1,
        allele_database_fasta, profiles_table_tsv
    ):
        """Result includes confidence field"""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir_path = Path(tmpdir)

            input_json = {
                "fast_typing.query_sequences": str(query_contigs_salmonella_st1),
                "fast_typing.allele_database": str(allele_database_fasta),
                "fast_typing.profiles_table": str(profiles_table_tsv)
            }

            input_json_path = tmpdir_path / "inputs.json"
            with open(input_json_path, "w") as f:
                json.dump(input_json, f)

            result = subprocess.run(
                ["miniwdl", "run", str(fast_workflow_path),
                 "-i", str(input_json_path),
                 "-d", str(tmpdir_path)],
                capture_output=True,
                text=True,
                timeout=600
            )

            assert result.returncode == 0, f"Workflow execution failed: {result.stderr}"

            output_files = list(tmpdir_path.glob("**/outputs.json"))
            with open(output_files[0]) as f:
                outputs = json.load(f)

            result_path = Path(outputs["fast_typing.typing_result"])
            with open(result_path) as f:
                result_data = json.load(f)

            assert "confidence" in result_data, "Result does not include confidence"
            assert isinstance(result_data["confidence"], (int, float)), \
                "Confidence is not numeric"

    def test_result_has_method_section(
        self, fast_workflow_path, query_contigs_salmonella_st1,
        allele_database_fasta, profiles_table_tsv
    ):
        """Result includes method section with strategy metadata"""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir_path = Path(tmpdir)

            input_json = {
                "fast_typing.query_sequences": str(query_contigs_salmonella_st1),
                "fast_typing.allele_database": str(allele_database_fasta),
                "fast_typing.profiles_table": str(profiles_table_tsv)
            }

            input_json_path = tmpdir_path / "inputs.json"
            with open(input_json_path, "w") as f:
                json.dump(input_json, f)

            result = subprocess.run(
                ["miniwdl", "run", str(fast_workflow_path),
                 "-i", str(input_json_path),
                 "-d", str(tmpdir_path)],
                capture_output=True,
                text=True,
                timeout=600
            )

            assert result.returncode == 0, f"Workflow execution failed: {result.stderr}"

            output_files = list(tmpdir_path.glob("**/outputs.json"))
            with open(output_files[0]) as f:
                outputs = json.load(f)

            result_path = Path(outputs["fast_typing.typing_result"])
            with open(result_path) as f:
                result_data = json.load(f)

            assert "method" in result_data, "Result does not include method section"

    def test_result_method_strategy_is_fast(
        self, fast_workflow_path, query_contigs_salmonella_st1,
        allele_database_fasta, profiles_table_tsv
    ):
        """Result method.strategy is 'fast'"""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir_path = Path(tmpdir)

            input_json = {
                "fast_typing.query_sequences": str(query_contigs_salmonella_st1),
                "fast_typing.allele_database": str(allele_database_fasta),
                "fast_typing.profiles_table": str(profiles_table_tsv)
            }

            input_json_path = tmpdir_path / "inputs.json"
            with open(input_json_path, "w") as f:
                json.dump(input_json, f)

            result = subprocess.run(
                ["miniwdl", "run", str(fast_workflow_path),
                 "-i", str(input_json_path),
                 "-d", str(tmpdir_path)],
                capture_output=True,
                text=True,
                timeout=600
            )

            assert result.returncode == 0, f"Workflow execution failed: {result.stderr}"

            output_files = list(tmpdir_path.glob("**/outputs.json"))
            with open(output_files[0]) as f:
                outputs = json.load(f)

            result_path = Path(outputs["fast_typing.typing_result"])
            with open(result_path) as f:
                result_data = json.load(f)

            assert "method" in result_data and "strategy" in result_data["method"], \
                "Result method does not include strategy"
            assert result_data["method"]["strategy"] == "fast", \
                f"Expected strategy='fast', got {result_data['method']['strategy']}"

    def test_result_method_alignment_used_is_false(
        self, fast_workflow_path, query_contigs_salmonella_st1,
        allele_database_fasta, profiles_table_tsv
    ):
        """Result method.alignment_used is false (fast strategy skips alignment)"""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir_path = Path(tmpdir)

            input_json = {
                "fast_typing.query_sequences": str(query_contigs_salmonella_st1),
                "fast_typing.allele_database": str(allele_database_fasta),
                "fast_typing.profiles_table": str(profiles_table_tsv)
            }

            input_json_path = tmpdir_path / "inputs.json"
            with open(input_json_path, "w") as f:
                json.dump(input_json, f)

            result = subprocess.run(
                ["miniwdl", "run", str(fast_workflow_path),
                 "-i", str(input_json_path),
                 "-d", str(tmpdir_path)],
                capture_output=True,
                text=True,
                timeout=600
            )

            assert result.returncode == 0, f"Workflow execution failed: {result.stderr}"

            output_files = list(tmpdir_path.glob("**/outputs.json"))
            with open(output_files[0]) as f:
                outputs = json.load(f)

            result_path = Path(outputs["fast_typing.typing_result"])
            with open(result_path) as f:
                result_data = json.load(f)

            assert "method" in result_data and "alignment_used" in result_data["method"], \
                "Result method does not include alignment_used"
            assert result_data["method"]["alignment_used"] is False, \
                "Fast strategy should have alignment_used=false"


@pytest.mark.miniwdl
class TestFastWorkflowTypingAccuracy:
    """Test workflow produces correct typing results."""

    def test_workflow_identifies_known_st(
        self, fast_workflow_path, query_contigs_salmonella_st1,
        allele_database_fasta, profiles_table_tsv
    ):
        """Workflow correctly identifies known ST=1"""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir_path = Path(tmpdir)

            input_json = {
                "fast_typing.query_sequences": str(query_contigs_salmonella_st1),
                "fast_typing.allele_database": str(allele_database_fasta),
                "fast_typing.profiles_table": str(profiles_table_tsv)
            }

            input_json_path = tmpdir_path / "inputs.json"
            with open(input_json_path, "w") as f:
                json.dump(input_json, f)

            result = subprocess.run(
                ["miniwdl", "run", str(fast_workflow_path),
                 "-i", str(input_json_path),
                 "-d", str(tmpdir_path)],
                capture_output=True,
                text=True,
                timeout=600
            )

            assert result.returncode == 0, f"Workflow execution failed: {result.stderr}"

            output_files = list(tmpdir_path.glob("**/outputs.json"))
            with open(output_files[0]) as f:
                outputs = json.load(f)

            result_path = Path(outputs["fast_typing.typing_result"])
            with open(result_path) as f:
                result_data = json.load(f)

            # Should identify ST=1
            profile_id = (result_data.get("profile_id") or
                          result_data.get("st") or
                          result_data.get("sequence_type"))
            assert profile_id == "1" or profile_id == 1, \
                f"Expected ST=1, got {profile_id}"

    def test_workflow_identifies_novel_profile(
        self, fast_workflow_path, query_contigs_novel_profile,
        allele_database_fasta, profiles_table_tsv
    ):
        """Workflow correctly identifies novel profile"""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir_path = Path(tmpdir)

            input_json = {
                "fast_typing.query_sequences": str(query_contigs_novel_profile),
                "fast_typing.allele_database": str(allele_database_fasta),
                "fast_typing.profiles_table": str(profiles_table_tsv)
            }

            input_json_path = tmpdir_path / "inputs.json"
            with open(input_json_path, "w") as f:
                json.dump(input_json, f)

            result = subprocess.run(
                ["miniwdl", "run", str(fast_workflow_path),
                 "-i", str(input_json_path),
                 "-d", str(tmpdir_path)],
                capture_output=True,
                text=True,
                timeout=600
            )

            assert result.returncode == 0, f"Workflow execution failed: {result.stderr}"

            output_files = list(tmpdir_path.glob("**/outputs.json"))
            with open(output_files[0]) as f:
                outputs = json.load(f)

            result_path = Path(outputs["fast_typing.typing_result"])
            with open(result_path) as f:
                result_data = json.load(f)

            # Should identify as novel profile
            assert result_data["status"] == "novel_profile", \
                f"Expected status=novel_profile, got {result_data['status']}"

    def test_workflow_result_has_allele_profile(
        self, fast_workflow_path, query_contigs_salmonella_st1,
        allele_database_fasta, profiles_table_tsv
    ):
        """Workflow result includes allele_profile field"""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir_path = Path(tmpdir)

            input_json = {
                "fast_typing.query_sequences": str(query_contigs_salmonella_st1),
                "fast_typing.allele_database": str(allele_database_fasta),
                "fast_typing.profiles_table": str(profiles_table_tsv)
            }

            input_json_path = tmpdir_path / "inputs.json"
            with open(input_json_path, "w") as f:
                json.dump(input_json, f)

            result = subprocess.run(
                ["miniwdl", "run", str(fast_workflow_path),
                 "-i", str(input_json_path),
                 "-d", str(tmpdir_path)],
                capture_output=True,
                text=True,
                timeout=600
            )

            assert result.returncode == 0, f"Workflow execution failed: {result.stderr}"

            output_files = list(tmpdir_path.glob("**/outputs.json"))
            with open(output_files[0]) as f:
                outputs = json.load(f)

            result_path = Path(outputs["fast_typing.typing_result"])
            with open(result_path) as f:
                result_data = json.load(f)

            assert "allele_profile" in result_data or "allele_calls" in result_data, \
                "Result does not include allele profile information"


@pytest.mark.miniwdl
class TestFastWorkflowEdgeCases:
    """Test edge cases and error handling."""

    def test_workflow_handles_empty_query(
        self, fast_workflow_path, allele_database_fasta, profiles_table_tsv
    ):
        """Workflow handles empty query sequences gracefully"""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir_path = Path(tmpdir)

            # Create empty query file
            empty_query = tmpdir_path / "empty.fasta"
            empty_query.touch()

            input_json = {
                "fast_typing.query_sequences": str(empty_query),
                "fast_typing.allele_database": str(allele_database_fasta),
                "fast_typing.profiles_table": str(profiles_table_tsv)
            }

            input_json_path = tmpdir_path / "inputs.json"
            with open(input_json_path, "w") as f:
                json.dump(input_json, f)

            result = subprocess.run(
                ["miniwdl", "run", str(fast_workflow_path),
                 "-i", str(input_json_path),
                 "-d", str(tmpdir_path)],
                capture_output=True,
                text=True,
                timeout=600
            )

            # Should either succeed with empty result or fail gracefully
            assert result.returncode in [0, 1], \
                "Unexpected return code for empty query"

    def test_workflow_handles_partial_profile(
        self, fast_workflow_path, allele_database_fasta, profiles_table_tsv
    ):
        """Workflow handles partial profile (missing loci)"""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir_path = Path(tmpdir)

            # Create query with only 2 out of 3 loci
            partial_query = tmpdir_path / "partial.fasta"
            with open(partial_query, "w") as f:
                f.write(">contig1_adk_1\n")
                f.write("ATGAATATTAACAACGCACTGGGCGACGTGCTGAAAACCCACGGCCAGATGACGAAAGAAGTGATGCAA\n")
                f.write(">contig2_fumC_1\n")
                f.write("CTGACCCAAGGTGCAACCCACGCCTTTGTGACCGCCGTGGGCGACTCGCCCGAAGAAACGCACCACGGA\n")

            input_json = {
                "fast_typing.query_sequences": str(partial_query),
                "fast_typing.allele_database": str(allele_database_fasta),
                "fast_typing.profiles_table": str(profiles_table_tsv)
            }

            input_json_path = tmpdir_path / "inputs.json"
            with open(input_json_path, "w") as f:
                json.dump(input_json, f)

            result = subprocess.run(
                ["miniwdl", "run", str(fast_workflow_path),
                 "-i", str(input_json_path),
                 "-d", str(tmpdir_path)],
                capture_output=True,
                text=True,
                timeout=600
            )

            # Should succeed and report partial profile
            assert result.returncode == 0, \
                f"Workflow should handle partial profile: {result.stderr}"
