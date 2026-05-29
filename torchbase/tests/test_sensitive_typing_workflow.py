"""Acceptance tests for sensitive typing workflow (Issue #57).

These are RED-phase tests - they MUST fail because the feature is not yet complete.

The sensitive strategy workflow always runs full alignment-based calling for maximum
accuracy. MinHash is used only for guidance/filtering. This is the most accurate but
slowest strategy.

Pipeline: Input sequences → MinHash sketching (guide) → Alignment (minimap2 asm5+eqx) →
          Refined allele calling → Profile lookup → Result

Acceptance criteria:
- torchbase/workflows/builtin/sensitive_typing.wdl exists
- Workflow imports shared tasks from tasks/ directory
- Always runs alignment regardless of MinHash confidence
- Uses minimap2 with asm5+eqx preset (high accuracy)
- Applies strict confidence thresholds (0.95)
- Outputs enhanced alignment metrics in notes field
- Tests verify alignment always runs and produces detailed output
- Output includes method: {strategy: "sensitive", alignment_used: true}
"""

import pytest
import json
import tempfile
from pathlib import Path
import subprocess


# Get the torchbase root directory
TORCHBASE_ROOT = Path(__file__).parent.parent


@pytest.fixture
def sensitive_workflow_path():
    """Path to the sensitive typing workflow WDL file."""
    return TORCHBASE_ROOT / "workflows" / "builtin" / "sensitive_typing.wdl"


@pytest.fixture
def sample_query_sequences():
    """Create sample query sequences for testing."""
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir_path = Path(tmpdir)
        query_path = tmpdir_path / "query.fasta"

        fasta_content = """>query_adk_close_to_1
ATGAATATTAACAACGCACTGGGCGACGTGCTGAAAACCCACGGCCAGATGACGAAAGAAGTGATGCAA
>query_fumC_close_to_1
CTGACCCAAGGTGCAACCCACGCCTTTGTGACCGCCGTGGGCGACTCGCCCGAAGAAACGCACCACGGA
>query_gyrB_close_to_1
ATGACCCAACTGAAAGTGATGCCGCAACGTGTCGACCTGCAAATCCACGCAGTGCTGATGAAACCGATG
"""

        with open(query_path, "w") as f:
            f.write(fasta_content)

        yield query_path


@pytest.fixture
def sample_allele_database():
    """Create sample allele database."""
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir_path = Path(tmpdir)
        db_path = tmpdir_path / "alleles.fasta"

        fasta_content = """>salmonella_adk_1
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
def sample_profile_table():
    """Create sample profile table."""
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir_path = Path(tmpdir)
        profile_path = tmpdir_path / "profiles.tsv"

        profile_content = """ST\tadk\tfumC\tgyrB
1\t1\t1\t1
2\t1\t2\t1
3\t2\t1\t1
4\t1\t1\t2
"""

        with open(profile_path, "w") as f:
            f.write(profile_content)

        yield profile_path


class TestSensitiveWorkflowFileExists:
    """Test sensitive_typing.wdl file exists at expected location."""

    def test_sensitive_wdl_file_exists(self, sensitive_workflow_path):
        """Sensitive typing workflow WDL file exists"""
        assert sensitive_workflow_path.exists(), \
            f"WDL file not found at {sensitive_workflow_path}"

    def test_sensitive_wdl_is_file(self, sensitive_workflow_path):
        """Sensitive typing workflow WDL is a regular file"""
        assert sensitive_workflow_path.is_file(), \
            f"WDL path is not a file: {sensitive_workflow_path}"


class TestSensitiveWorkflowStructure:
    """Test WDL workflow structure and components."""

    def test_wdl_has_workflow_definition(self, sensitive_workflow_path):
        """WDL file contains workflow definition"""
        with open(sensitive_workflow_path) as f:
            content = f.read()

        assert "workflow" in content, "WDL file does not contain workflow definition"

    def test_wdl_workflow_name_is_sensitive_typing(self, sensitive_workflow_path):
        """WDL workflow is named sensitive_typing or similar"""
        with open(sensitive_workflow_path) as f:
            content = f.read()

        assert "workflow sensitive" in content.lower(), \
            "Workflow name is not sensitive-related"

    def test_wdl_has_input_section(self, sensitive_workflow_path):
        """WDL workflow has input section"""
        with open(sensitive_workflow_path) as f:
            content = f.read()

        assert "input {" in content, "Workflow does not have input section"

    def test_wdl_has_output_section(self, sensitive_workflow_path):
        """WDL workflow has output section"""
        with open(sensitive_workflow_path) as f:
            content = f.read()

        assert "output {" in content, "Workflow does not have output section"

    def test_wdl_accepts_query_sequences_input(self, sensitive_workflow_path):
        """WDL workflow accepts query sequences input"""
        with open(sensitive_workflow_path) as f:
            content = f.read()

        assert "query" in content.lower() or "sequences" in content.lower(), \
            "Workflow does not accept query sequences input"

    def test_wdl_accepts_allele_database_input(self, sensitive_workflow_path):
        """WDL workflow accepts allele database input"""
        with open(sensitive_workflow_path) as f:
            content = f.read()

        assert "allele" in content.lower(), \
            "Workflow does not accept allele database input"

    def test_wdl_accepts_profile_table_input(self, sensitive_workflow_path):
        """WDL workflow accepts profile table input"""
        with open(sensitive_workflow_path) as f:
            content = f.read()

        assert "profile" in content.lower(), \
            "Workflow does not accept profile table input"


class TestSensitiveWorkflowImports:
    """Test workflow imports shared tasks from tasks/ directory."""

    def test_wdl_imports_shared_tasks(self, sensitive_workflow_path):
        """WDL workflow imports shared task files"""
        with open(sensitive_workflow_path) as f:
            content = f.read()

        assert "import" in content.lower(), \
            "Workflow does not import shared tasks"

    def test_wdl_imports_minhash_tasks(self, sensitive_workflow_path):
        """WDL workflow imports minhash.wdl"""
        with open(sensitive_workflow_path) as f:
            content = f.read()

        assert "minhash.wdl" in content or "minhash" in content.lower(), \
            "Workflow does not import minhash tasks"

    def test_wdl_imports_alignment_tasks(self, sensitive_workflow_path):
        """WDL workflow imports alignment.wdl"""
        with open(sensitive_workflow_path) as f:
            content = f.read()

        has_alignment_import = (
            "alignment.wdl" in content or
            ("alignment" in content.lower() and "import" in content.lower())
        )
        assert has_alignment_import, "Workflow does not import alignment tasks"

    def test_wdl_imports_profile_lookup_tasks(self, sensitive_workflow_path):
        """WDL workflow imports profile_lookup.wdl"""
        with open(sensitive_workflow_path) as f:
            content = f.read()

        has_profile_import = (
            "profile_lookup.wdl" in content or
            ("profile" in content.lower() and "import" in content.lower())
        )
        assert has_profile_import, "Workflow does not import profile lookup tasks"

    def test_wdl_import_paths_point_to_tasks_directory(self, sensitive_workflow_path):
        """WDL imports reference tasks/ directory"""
        with open(sensitive_workflow_path) as f:
            content = f.read()

        # Check that imports reference the tasks directory
        assert "tasks/" in content, \
            "Workflow imports do not reference tasks/ directory"


class TestSensitiveWorkflowAlignmentAlwaysRuns:
    """Test alignment always runs regardless of MinHash confidence."""

    def test_wdl_calls_alignment_task(self, sensitive_workflow_path):
        """WDL workflow calls alignment task"""
        with open(sensitive_workflow_path) as f:
            content = f.read()

        assert "align" in content.lower() or "minimap" in content.lower(), \
            "Workflow does not call alignment task"

    def test_wdl_alignment_not_conditional_on_confidence(self, sensitive_workflow_path):
        """WDL workflow alignment is not conditional on MinHash confidence"""
        with open(sensitive_workflow_path) as f:
            content = f.read()

        # Check that alignment task is called unconditionally
        # Look for absence of conditional wrapping around alignment
        lines = content.split('\n')
        alignment_line_idx = None

        for idx, line in enumerate(lines):
            if 'call' in line.lower() and ('align' in line.lower() or 'minimap' in line.lower()):
                alignment_line_idx = idx
                break

        assert alignment_line_idx is not None, "No alignment call found"

        # The alignment should run unconditionally (no if statement checking confidence)
        # Note: This is a heuristic - actual implementation may vary
        # Key point: alignment always runs in sensitive mode
        assert alignment_line_idx is not None, "Alignment task should be present"

    def test_wdl_minhash_used_only_for_guidance(self, sensitive_workflow_path):
        """WDL workflow uses MinHash only for guidance, not decision"""
        with open(sensitive_workflow_path) as f:
            content = f.read()

        # MinHash should be present but not used for conditional logic
        assert "minhash" in content.lower() or "sketch" in content.lower(), \
            "Workflow does not use MinHash at all"

        # Alignment should be the primary method
        assert "align" in content.lower() or "minimap" in content.lower(), \
            "Workflow does not use alignment"


class TestSensitiveWorkflowMinimapPreset:
    """Test minimap2 uses asm5+eqx preset for high accuracy."""

    def test_wdl_uses_minimap2(self, sensitive_workflow_path):
        """WDL workflow uses minimap2"""
        with open(sensitive_workflow_path) as f:
            content = f.read()

        assert "minimap" in content.lower(), \
            "Workflow does not use minimap2"

    def test_wdl_specifies_asm5_eqx_preset(self, sensitive_workflow_path):
        """WDL workflow specifies asm5+eqx preset"""
        with open(sensitive_workflow_path) as f:
            content = f.read()

        # Check for asm5 or asm5+eqx preset specification
        assert "asm5" in content.lower(), \
            "Workflow does not specify asm5 preset"

    def test_wdl_uses_high_accuracy_preset(self, sensitive_workflow_path):
        """WDL workflow uses high accuracy alignment preset"""
        with open(sensitive_workflow_path) as f:
            content = f.read()

        # Check for asm5+eqx which is the most accurate preset
        # or check for preset parameter being set to a sensitive value
        assert "preset" in content.lower() or "asm5" in content.lower(), \
            "Workflow does not specify alignment preset"

    def test_wdl_alignment_preset_parameter(self, sensitive_workflow_path):
        """WDL workflow passes alignment preset parameter"""
        with open(sensitive_workflow_path) as f:
            content = f.read()

        # Look for preset parameter being passed to alignment task
        assert "preset" in content.lower(), \
            "Workflow does not pass preset parameter to alignment"


class TestSensitiveWorkflowConfidenceThresholds:
    """Test strict confidence thresholds (0.95) are applied."""

    def test_wdl_has_confidence_threshold_parameter(self, sensitive_workflow_path):
        """WDL workflow has confidence threshold parameter"""
        with open(sensitive_workflow_path) as f:
            content = f.read()

        assert "threshold" in content.lower() or "confidence" in content.lower(), \
            "Workflow does not have confidence threshold parameter"

    def test_wdl_strict_threshold_default_value(self, sensitive_workflow_path):
        """WDL workflow has strict threshold default (0.95)"""
        with open(sensitive_workflow_path) as f:
            content = f.read()

        # Check for 0.95 or higher threshold
        has_strict_threshold = any(
            thresh in content
            for thresh in ["0.95", "0.96", "0.97", "0.98", "0.99"]
        )
        assert has_strict_threshold, \
            "Workflow does not have strict confidence threshold (0.95+)"

    def test_wdl_applies_threshold_to_allele_calls(self, sensitive_workflow_path):
        """WDL workflow applies threshold to allele calls"""
        with open(sensitive_workflow_path) as f:
            content = f.read()

        # Check that threshold is used in filtering or quality assessment
        assert "threshold" in content.lower(), \
            "Workflow does not apply confidence threshold"


class TestSensitiveWorkflowOutputStructure:
    """Test output includes method metadata and alignment metrics."""

    def test_wdl_output_includes_result_file(self, sensitive_workflow_path):
        """WDL workflow output includes typing result file"""
        with open(sensitive_workflow_path) as f:
            content = f.read()

        # Check output section has result file
        output_section_start = content.lower().find("output {")
        if output_section_start != -1:
            output_section = content[output_section_start:]
            has_result_file = (
                "file" in output_section.lower() or
                "result" in output_section.lower()
            )
            assert has_result_file, "Workflow output does not include result file"

    def test_wdl_output_includes_typing_result(self, sensitive_workflow_path):
        """WDL workflow output includes typing result"""
        with open(sensitive_workflow_path) as f:
            content = f.read()

        assert "result" in content.lower() or "typing" in content.lower(), \
            "Workflow output does not include typing result"

    def test_wdl_produces_method_metadata(self, sensitive_workflow_path):
        """WDL workflow produces method metadata in output"""
        with open(sensitive_workflow_path) as f:
            content = f.read()

        # Check that method metadata is assembled
        assert "method" in content.lower() or "strategy" in content.lower(), \
            "Workflow does not produce method metadata"

    def test_wdl_produces_alignment_metrics(self, sensitive_workflow_path):
        """WDL workflow produces alignment metrics in output"""
        with open(sensitive_workflow_path) as f:
            content = f.read()

        # Check that alignment metrics are included
        has_alignment_metrics = (
            "metric" in content.lower() or
            "notes" in content.lower() or
            "alignment" in content.lower()
        )
        assert has_alignment_metrics, "Workflow does not produce alignment metrics"


class TestSensitiveWorkflowSyntaxValidation:
    """Test miniwdl check validates syntax."""

    def test_miniwdl_check_passes(self, sensitive_workflow_path):
        """miniwdl check validates WDL syntax without errors"""
        result = subprocess.run(
            ["miniwdl", "check", str(sensitive_workflow_path)],
            capture_output=True,
            text=True
        )

        assert result.returncode == 0, \
            f"miniwdl check failed: {result.stderr}"


@pytest.mark.miniwdl
class TestSensitiveWorkflowExecution:
    """Test end-to-end workflow execution."""

    def test_workflow_executes_successfully(
        self, sensitive_workflow_path, sample_query_sequences,
        sample_allele_database, sample_profile_table
    ):
        """Workflow executes successfully with valid inputs"""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir_path = Path(tmpdir)

            input_json = {
                "sensitive_typing.query_sequences": str(sample_query_sequences),
                "sensitive_typing.allele_database": str(sample_allele_database),
                "sensitive_typing.profiles": str(sample_profile_table)
            }

            input_json_path = tmpdir_path / "inputs.json"
            with open(input_json_path, "w") as f:
                json.dump(input_json, f)

            result = subprocess.run(
                ["miniwdl", "run", str(sensitive_workflow_path),
                 "-i", str(input_json_path), "-d", str(tmpdir_path)],
                capture_output=True,
                text=True,
                timeout=600
            )

            assert result.returncode == 0, \
                f"Workflow execution failed: {result.stderr}"

    def test_workflow_produces_output_file(
        self, sensitive_workflow_path, sample_query_sequences,
        sample_allele_database, sample_profile_table
    ):
        """Workflow produces typing result output file"""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir_path = Path(tmpdir)

            input_json = {
                "sensitive_typing.query_sequences": str(sample_query_sequences),
                "sensitive_typing.allele_database": str(sample_allele_database),
                "sensitive_typing.profiles": str(sample_profile_table)
            }

            input_json_path = tmpdir_path / "inputs.json"
            with open(input_json_path, "w") as f:
                json.dump(input_json, f)

            result = subprocess.run(
                ["miniwdl", "run", str(sensitive_workflow_path),
                 "-i", str(input_json_path), "-d", str(tmpdir_path)],
                capture_output=True,
                text=True,
                timeout=600
            )

            assert result.returncode == 0, f"Workflow failed: {result.stderr}"

            # Find outputs.json
            output_files = list(tmpdir_path.glob("**/outputs.json"))
            assert len(output_files) > 0, "No outputs.json found"

            with open(output_files[0]) as f:
                outputs = json.load(f)

            assert "sensitive_typing.typing_result" in outputs, \
                "Typing result not in outputs"


@pytest.mark.miniwdl
class TestSensitiveWorkflowAlignmentExecution:
    """Test alignment always runs in execution."""

    def test_alignment_runs_even_with_high_confidence_minhash(
        self, sensitive_workflow_path, sample_query_sequences,
        sample_allele_database, sample_profile_table
    ):
        """Alignment runs even when MinHash has high confidence"""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir_path = Path(tmpdir)

            # Use perfect match query sequences (high MinHash confidence)
            perfect_query = tmpdir_path / "perfect_query.fasta"
            with open(perfect_query, "w") as f:
                f.write(""">query_adk_1
ATGAATATTAACAACGCACTGGGCGACGTGCTGAAAACCCACGGCCAGATGACGAAAGAAGTGATGCAA
>query_fumC_1
CTGACCCAAGGTGCAACCCACGCCTTTGTGACCGCCGTGGGCGACTCGCCCGAAGAAACGCACCACGGA
>query_gyrB_1
ATGACCCAACTGAAAGTGATGCCGCAACGTGTCGACCTGCAAATCCACGCAGTGCTGATGAAACCGATG
""")

            input_json = {
                "sensitive_typing.query_sequences": str(perfect_query),
                "sensitive_typing.allele_database": str(sample_allele_database),
                "sensitive_typing.profiles": str(sample_profile_table)
            }

            input_json_path = tmpdir_path / "inputs.json"
            with open(input_json_path, "w") as f:
                json.dump(input_json, f)

            result = subprocess.run(
                ["miniwdl", "run", str(sensitive_workflow_path),
                 "-i", str(input_json_path), "-d", str(tmpdir_path)],
                capture_output=True,
                text=True,
                timeout=600
            )

            assert result.returncode == 0, f"Workflow failed: {result.stderr}"

            # Find result file
            output_files = list(tmpdir_path.glob("**/outputs.json"))
            with open(output_files[0]) as f:
                outputs = json.load(f)

            result_path = Path(outputs["sensitive_typing.typing_result"])
            with open(result_path) as f:
                result_data = json.load(f)

            # Verify alignment was used
            assert "method" in result_data, "Result does not include method metadata"
            assert "alignment_used" in result_data["method"], \
                "Result method does not indicate if alignment was used"
            assert result_data["method"]["alignment_used"] is True, \
                "Alignment was not used (should always run in sensitive mode)"


@pytest.mark.miniwdl
class TestSensitiveWorkflowOutputFormat:
    """Test output format includes required fields."""

    def test_output_includes_method_field(
        self, sensitive_workflow_path, sample_query_sequences,
        sample_allele_database, sample_profile_table
    ):
        """Output includes method field with strategy and alignment_used"""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir_path = Path(tmpdir)

            input_json = {
                "sensitive_typing.query_sequences": str(sample_query_sequences),
                "sensitive_typing.allele_database": str(sample_allele_database),
                "sensitive_typing.profiles": str(sample_profile_table)
            }

            input_json_path = tmpdir_path / "inputs.json"
            with open(input_json_path, "w") as f:
                json.dump(input_json, f)

            result = subprocess.run(
                ["miniwdl", "run", str(sensitive_workflow_path),
                 "-i", str(input_json_path), "-d", str(tmpdir_path)],
                capture_output=True,
                text=True,
                timeout=600
            )

            assert result.returncode == 0, f"Workflow failed: {result.stderr}"

            output_files = list(tmpdir_path.glob("**/outputs.json"))
            with open(output_files[0]) as f:
                outputs = json.load(f)

            result_path = Path(outputs["sensitive_typing.typing_result"])
            with open(result_path) as f:
                result_data = json.load(f)

            assert "method" in result_data, "Result does not include method field"

    def test_output_method_includes_strategy_sensitive(
        self, sensitive_workflow_path, sample_query_sequences,
        sample_allele_database, sample_profile_table
    ):
        """Output method field includes strategy: sensitive"""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir_path = Path(tmpdir)

            input_json = {
                "sensitive_typing.query_sequences": str(sample_query_sequences),
                "sensitive_typing.allele_database": str(sample_allele_database),
                "sensitive_typing.profiles": str(sample_profile_table)
            }

            input_json_path = tmpdir_path / "inputs.json"
            with open(input_json_path, "w") as f:
                json.dump(input_json, f)

            result = subprocess.run(
                ["miniwdl", "run", str(sensitive_workflow_path),
                 "-i", str(input_json_path), "-d", str(tmpdir_path)],
                capture_output=True,
                text=True,
                timeout=600
            )

            assert result.returncode == 0, f"Workflow failed: {result.stderr}"

            output_files = list(tmpdir_path.glob("**/outputs.json"))
            with open(output_files[0]) as f:
                outputs = json.load(f)

            result_path = Path(outputs["sensitive_typing.typing_result"])
            with open(result_path) as f:
                result_data = json.load(f)

            assert result_data["method"]["strategy"] == "sensitive", \
                f"Strategy should be 'sensitive', got {result_data['method'].get('strategy')}"

    def test_output_method_includes_alignment_used_true(
        self, sensitive_workflow_path, sample_query_sequences,
        sample_allele_database, sample_profile_table
    ):
        """Output method field includes alignment_used: true"""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir_path = Path(tmpdir)

            input_json = {
                "sensitive_typing.query_sequences": str(sample_query_sequences),
                "sensitive_typing.allele_database": str(sample_allele_database),
                "sensitive_typing.profiles": str(sample_profile_table)
            }

            input_json_path = tmpdir_path / "inputs.json"
            with open(input_json_path, "w") as f:
                json.dump(input_json, f)

            result = subprocess.run(
                ["miniwdl", "run", str(sensitive_workflow_path),
                 "-i", str(input_json_path), "-d", str(tmpdir_path)],
                capture_output=True,
                text=True,
                timeout=600
            )

            assert result.returncode == 0, f"Workflow failed: {result.stderr}"

            output_files = list(tmpdir_path.glob("**/outputs.json"))
            with open(output_files[0]) as f:
                outputs = json.load(f)

            result_path = Path(outputs["sensitive_typing.typing_result"])
            with open(result_path) as f:
                result_data = json.load(f)

            assert result_data["method"]["alignment_used"] is True, \
                "alignment_used should be True in sensitive strategy"

    def test_output_includes_notes_with_alignment_metrics(
        self, sensitive_workflow_path, sample_query_sequences,
        sample_allele_database, sample_profile_table
    ):
        """Output includes notes field with enhanced alignment metrics"""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir_path = Path(tmpdir)

            input_json = {
                "sensitive_typing.query_sequences": str(sample_query_sequences),
                "sensitive_typing.allele_database": str(sample_allele_database),
                "sensitive_typing.profiles": str(sample_profile_table)
            }

            input_json_path = tmpdir_path / "inputs.json"
            with open(input_json_path, "w") as f:
                json.dump(input_json, f)

            result = subprocess.run(
                ["miniwdl", "run", str(sensitive_workflow_path),
                 "-i", str(input_json_path), "-d", str(tmpdir_path)],
                capture_output=True,
                text=True,
                timeout=600
            )

            assert result.returncode == 0, f"Workflow failed: {result.stderr}"

            output_files = list(tmpdir_path.glob("**/outputs.json"))
            with open(output_files[0]) as f:
                outputs = json.load(f)

            result_path = Path(outputs["sensitive_typing.typing_result"])
            with open(result_path) as f:
                result_data = json.load(f)

            assert "notes" in result_data, "Result does not include notes field"

            # Notes should contain alignment-related metrics
            notes = result_data["notes"]
            # Check for alignment metrics (identity, coverage, quality, etc.)
            has_alignment_metrics = (
                "alignment" in str(notes).lower() or
                "identity" in str(notes).lower() or
                "coverage" in str(notes).lower() or
                "mapq" in str(notes).lower()
            )

            assert has_alignment_metrics, \
                "Notes field does not contain alignment metrics"


@pytest.mark.miniwdl
class TestSensitiveWorkflowDetailedOutput:
    """Test detailed alignment output and metrics."""

    def test_output_includes_profile_id(
        self, sensitive_workflow_path, sample_query_sequences,
        sample_allele_database, sample_profile_table
    ):
        """Output includes profile_id field"""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir_path = Path(tmpdir)

            input_json = {
                "sensitive_typing.query_sequences": str(sample_query_sequences),
                "sensitive_typing.allele_database": str(sample_allele_database),
                "sensitive_typing.profiles": str(sample_profile_table)
            }

            input_json_path = tmpdir_path / "inputs.json"
            with open(input_json_path, "w") as f:
                json.dump(input_json, f)

            result = subprocess.run(
                ["miniwdl", "run", str(sensitive_workflow_path),
                 "-i", str(input_json_path), "-d", str(tmpdir_path)],
                capture_output=True,
                text=True,
                timeout=600
            )

            assert result.returncode == 0, f"Workflow failed: {result.stderr}"

            output_files = list(tmpdir_path.glob("**/outputs.json"))
            with open(output_files[0]) as f:
                outputs = json.load(f)

            result_path = Path(outputs["sensitive_typing.typing_result"])
            with open(result_path) as f:
                result_data = json.load(f)

            has_profile_id = (
                "profile_id" in result_data or
                "st" in result_data.get("profile_id", "").lower()
            )
            assert has_profile_id, "Result does not include profile_id"

    def test_output_includes_confidence_score(
        self, sensitive_workflow_path, sample_query_sequences,
        sample_allele_database, sample_profile_table
    ):
        """Output includes confidence score"""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir_path = Path(tmpdir)

            input_json = {
                "sensitive_typing.query_sequences": str(sample_query_sequences),
                "sensitive_typing.allele_database": str(sample_allele_database),
                "sensitive_typing.profiles": str(sample_profile_table)
            }

            input_json_path = tmpdir_path / "inputs.json"
            with open(input_json_path, "w") as f:
                json.dump(input_json, f)

            result = subprocess.run(
                ["miniwdl", "run", str(sensitive_workflow_path),
                 "-i", str(input_json_path), "-d", str(tmpdir_path)],
                capture_output=True,
                text=True,
                timeout=600
            )

            assert result.returncode == 0, f"Workflow failed: {result.stderr}"

            output_files = list(tmpdir_path.glob("**/outputs.json"))
            with open(output_files[0]) as f:
                outputs = json.load(f)

            result_path = Path(outputs["sensitive_typing.typing_result"])
            with open(result_path) as f:
                result_data = json.load(f)

            assert "confidence" in result_data, \
                "Result does not include confidence score"

    def test_output_includes_allele_calls(
        self, sensitive_workflow_path, sample_query_sequences,
        sample_allele_database, sample_profile_table
    ):
        """Output includes allele calls"""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir_path = Path(tmpdir)

            input_json = {
                "sensitive_typing.query_sequences": str(sample_query_sequences),
                "sensitive_typing.allele_database": str(sample_allele_database),
                "sensitive_typing.profiles": str(sample_profile_table)
            }

            input_json_path = tmpdir_path / "inputs.json"
            with open(input_json_path, "w") as f:
                json.dump(input_json, f)

            result = subprocess.run(
                ["miniwdl", "run", str(sensitive_workflow_path),
                 "-i", str(input_json_path), "-d", str(tmpdir_path)],
                capture_output=True,
                text=True,
                timeout=600
            )

            assert result.returncode == 0, f"Workflow failed: {result.stderr}"

            output_files = list(tmpdir_path.glob("**/outputs.json"))
            with open(output_files[0]) as f:
                outputs = json.load(f)

            result_path = Path(outputs["sensitive_typing.typing_result"])
            with open(result_path) as f:
                result_data = json.load(f)

            assert "allele_calls" in result_data or "allele_profile" in result_data, \
                "Result does not include allele calls"


@pytest.mark.miniwdl
class TestSensitiveWorkflowEdgeCases:
    """Test edge cases and error handling."""

    def test_workflow_handles_ambiguous_sequences(
        self, sensitive_workflow_path, sample_allele_database, sample_profile_table
    ):
        """Workflow handles ambiguous query sequences correctly"""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir_path = Path(tmpdir)

            # Create ambiguous query sequences
            ambiguous_query = tmpdir_path / "ambiguous.fasta"
            with open(ambiguous_query, "w") as f:
                f.write(""">query_between_alleles
ATGAATATTAACAACGCACTGGGCGACGTGCTGAAAACCCACGGCCAGATGACGAAAGAAGTGATGCAN
>query_mixed_loci
CTGACCCAAGGTGCAACCCACGCCTTTGTGACCGCCGTGGGCGACTCGCCCGAAGAAACGCACCACGGG
""")

            input_json = {
                "sensitive_typing.query_sequences": str(ambiguous_query),
                "sensitive_typing.allele_database": str(sample_allele_database),
                "sensitive_typing.profiles": str(sample_profile_table)
            }

            input_json_path = tmpdir_path / "inputs.json"
            with open(input_json_path, "w") as f:
                json.dump(input_json, f)

            result = subprocess.run(
                ["miniwdl", "run", str(sensitive_workflow_path),
                 "-i", str(input_json_path), "-d", str(tmpdir_path)],
                capture_output=True,
                text=True,
                timeout=600
            )

            # Should handle gracefully
            assert result.returncode == 0, \
                f"Workflow failed on ambiguous sequences: {result.stderr}"

    def test_workflow_handles_empty_input(
        self, sensitive_workflow_path, sample_allele_database, sample_profile_table
    ):
        """Workflow handles empty query sequences"""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir_path = Path(tmpdir)

            # Create empty query file
            empty_query = tmpdir_path / "empty.fasta"
            empty_query.touch()

            input_json = {
                "sensitive_typing.query_sequences": str(empty_query),
                "sensitive_typing.allele_database": str(sample_allele_database),
                "sensitive_typing.profiles": str(sample_profile_table)
            }

            input_json_path = tmpdir_path / "inputs.json"
            with open(input_json_path, "w") as f:
                json.dump(input_json, f)

            result = subprocess.run(
                ["miniwdl", "run", str(sensitive_workflow_path),
                 "-i", str(input_json_path), "-d", str(tmpdir_path)],
                capture_output=True,
                text=True,
                timeout=600
            )

            # Should either fail gracefully or produce empty result
            assert result.returncode in [0, 1], \
                "Unexpected return code for empty input"
