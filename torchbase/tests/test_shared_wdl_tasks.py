"""Acceptance tests for shared WDL tasks extraction (Issue #54).

These are RED-phase tests - they MUST fail because the feature is not yet complete.

Acceptance criteria:
- Tasks exist in `torchbase/workflows/builtin/tasks/` directory
- `minhash.wdl` contains sketch and compare tasks
- `alignment.wdl` contains minimap2 task with preset selection
- `profile_lookup.wdl` contains profile matching logic
- Tasks can be imported via WDL import statements
- Minimap2 task supports presets: asm20 (fast), asm5 (balanced), asm5+eqx (sensitive)
"""

import subprocess
from pathlib import Path


# Get the torchbase root directory
TORCHBASE_ROOT = Path(__file__).parent.parent


class TestTasksDirectoryStructure:
    """Test that tasks directory exists with expected structure."""

    def test_builtin_workflows_directory_exists(self):
        """builtin workflows directory exists"""
        builtin_dir = TORCHBASE_ROOT / "workflows" / "builtin"
        assert builtin_dir.exists(), f"builtin directory not found at {builtin_dir}"

    def test_tasks_directory_exists(self):
        """tasks directory exists under builtin/"""
        tasks_dir = TORCHBASE_ROOT / "workflows" / "builtin" / "tasks"
        assert tasks_dir.exists(), f"tasks directory not found at {tasks_dir}"

    def test_tasks_directory_is_directory(self):
        """tasks path is a directory"""
        tasks_dir = TORCHBASE_ROOT / "workflows" / "builtin" / "tasks"
        assert tasks_dir.is_dir(), f"tasks path is not a directory: {tasks_dir}"


class TestMinHashTaskFile:
    """Test minhash.wdl contains sketch and compare tasks."""

    def test_minhash_wdl_exists(self):
        """minhash.wdl file exists in tasks directory"""
        minhash_wdl = TORCHBASE_ROOT / "workflows" / "builtin" / "tasks" / "minhash.wdl"
        assert minhash_wdl.exists(), f"minhash.wdl not found at {minhash_wdl}"

    def test_minhash_wdl_is_file(self):
        """minhash.wdl is a regular file"""
        minhash_wdl = TORCHBASE_ROOT / "workflows" / "builtin" / "tasks" / "minhash.wdl"
        assert minhash_wdl.is_file(), f"minhash.wdl is not a file: {minhash_wdl}"

    def test_minhash_wdl_has_sketch_task(self):
        """minhash.wdl contains sketch_sequences task"""
        minhash_wdl = TORCHBASE_ROOT / "workflows" / "builtin" / "tasks" / "minhash.wdl"
        with open(minhash_wdl) as f:
            content = f.read()

        assert "task sketch_sequences" in content, "minhash.wdl does not contain sketch_sequences task"

    def test_minhash_wdl_has_compare_task(self):
        """minhash.wdl contains compare_sketches task"""
        minhash_wdl = TORCHBASE_ROOT / "workflows" / "builtin" / "tasks" / "minhash.wdl"
        with open(minhash_wdl) as f:
            content = f.read()

        assert "task compare_sketches" in content, "minhash.wdl does not contain compare_sketches task"

    def test_minhash_wdl_sketch_task_has_sequences_input(self):
        """sketch_sequences task accepts sequences input"""
        minhash_wdl = TORCHBASE_ROOT / "workflows" / "builtin" / "tasks" / "minhash.wdl"
        with open(minhash_wdl) as f:
            content = f.read()

        # Check that sketch task has sequences input
        assert "File sequences" in content, "sketch_sequences task does not have sequences input"

    def test_minhash_wdl_sketch_task_has_ksize_parameter(self):
        """sketch_sequences task accepts ksize parameter"""
        minhash_wdl = TORCHBASE_ROOT / "workflows" / "builtin" / "tasks" / "minhash.wdl"
        with open(minhash_wdl) as f:
            content = f.read()

        assert "ksize" in content, "sketch_sequences task does not have ksize parameter"

    def test_minhash_wdl_sketch_task_has_scaled_parameter(self):
        """sketch_sequences task accepts scaled parameter"""
        minhash_wdl = TORCHBASE_ROOT / "workflows" / "builtin" / "tasks" / "minhash.wdl"
        with open(minhash_wdl) as f:
            content = f.read()

        assert "scaled" in content, "sketch_sequences task does not have scaled parameter"

    def test_minhash_wdl_compare_task_has_query_sketch_input(self):
        """compare_sketches task accepts query_sketch input"""
        minhash_wdl = TORCHBASE_ROOT / "workflows" / "builtin" / "tasks" / "minhash.wdl"
        with open(minhash_wdl) as f:
            content = f.read()

        assert "query_sketch" in content, "compare_sketches task does not have query_sketch input"

    def test_minhash_wdl_compare_task_has_allele_sketch_input(self):
        """compare_sketches task accepts allele_sketch input"""
        minhash_wdl = TORCHBASE_ROOT / "workflows" / "builtin" / "tasks" / "minhash.wdl"
        with open(minhash_wdl) as f:
            content = f.read()

        assert "allele_sketch" in content, "compare_sketches task does not have allele_sketch input"

    def test_minhash_wdl_uses_sourmash(self):
        """minhash.wdl tasks use sourmash"""
        minhash_wdl = TORCHBASE_ROOT / "workflows" / "builtin" / "tasks" / "minhash.wdl"
        with open(minhash_wdl) as f:
            content = f.read()

        assert "sourmash" in content.lower(), "minhash.wdl does not use sourmash"

    def test_minhash_wdl_syntax_valid(self):
        """minhash.wdl passes miniwdl syntax validation"""
        minhash_wdl = TORCHBASE_ROOT / "workflows" / "builtin" / "tasks" / "minhash.wdl"

        result = subprocess.run(
            ["miniwdl", "check", str(minhash_wdl)],
            capture_output=True,
            text=True
        )

        assert result.returncode == 0, f"miniwdl check failed for minhash.wdl: {result.stderr}"


class TestAlignmentTaskFile:
    """Test alignment.wdl contains minimap2 task with preset selection."""

    def test_alignment_wdl_exists(self):
        """alignment.wdl file exists in tasks directory"""
        alignment_wdl = TORCHBASE_ROOT / "workflows" / "builtin" / "tasks" / "alignment.wdl"
        assert alignment_wdl.exists(), f"alignment.wdl not found at {alignment_wdl}"

    def test_alignment_wdl_is_file(self):
        """alignment.wdl is a regular file"""
        alignment_wdl = TORCHBASE_ROOT / "workflows" / "builtin" / "tasks" / "alignment.wdl"
        assert alignment_wdl.is_file(), f"alignment.wdl is not a file: {alignment_wdl}"

    def test_alignment_wdl_has_minimap2_task(self):
        """alignment.wdl contains align_sequences task"""
        alignment_wdl = TORCHBASE_ROOT / "workflows" / "builtin" / "tasks" / "alignment.wdl"
        with open(alignment_wdl) as f:
            content = f.read()

        # Accept either align_sequences or minimap2_align as task name
        assert ("task align_sequences" in content or "task minimap2_align" in content), \
            "alignment.wdl does not contain align_sequences or minimap2_align task"

    def test_alignment_wdl_has_query_input(self):
        """align task accepts query sequences input"""
        alignment_wdl = TORCHBASE_ROOT / "workflows" / "builtin" / "tasks" / "alignment.wdl"
        with open(alignment_wdl) as f:
            content = f.read()

        assert "query" in content.lower(), "align task does not have query input"

    def test_alignment_wdl_has_reference_input(self):
        """align task accepts reference sequences input"""
        alignment_wdl = TORCHBASE_ROOT / "workflows" / "builtin" / "tasks" / "alignment.wdl"
        with open(alignment_wdl) as f:
            content = f.read()

        assert "reference" in content.lower() or "allele" in content.lower(), \
            "align task does not have reference/allele input"

    def test_alignment_wdl_has_preset_parameter(self):
        """align task accepts preset parameter"""
        alignment_wdl = TORCHBASE_ROOT / "workflows" / "builtin" / "tasks" / "alignment.wdl"
        with open(alignment_wdl) as f:
            content = f.read()

        assert "preset" in content, "align task does not have preset parameter"

    def test_alignment_wdl_supports_asm20_preset(self):
        """align task supports asm20 preset (fast)"""
        alignment_wdl = TORCHBASE_ROOT / "workflows" / "builtin" / "tasks" / "alignment.wdl"
        with open(alignment_wdl) as f:
            content = f.read()

        # Check for asm20 preset documentation or usage
        assert "asm20" in content, "align task does not support asm20 preset"

    def test_alignment_wdl_supports_asm5_preset(self):
        """align task supports asm5 preset (balanced)"""
        alignment_wdl = TORCHBASE_ROOT / "workflows" / "builtin" / "tasks" / "alignment.wdl"
        with open(alignment_wdl) as f:
            content = f.read()

        # Check for asm5 preset documentation or usage
        assert "asm5" in content, "align task does not support asm5 preset"

    def test_alignment_wdl_supports_asm5_eqx_preset(self):
        """align task supports asm5+eqx preset (sensitive)"""
        alignment_wdl = TORCHBASE_ROOT / "workflows" / "builtin" / "tasks" / "alignment.wdl"
        with open(alignment_wdl) as f:
            content = f.read()

        # Check for asm5+eqx or eqx preset documentation or usage
        assert ("eqx" in content or "asm5+eqx" in content), \
            "align task does not support asm5+eqx preset"

    def test_alignment_wdl_uses_minimap2(self):
        """alignment.wdl task uses minimap2"""
        alignment_wdl = TORCHBASE_ROOT / "workflows" / "builtin" / "tasks" / "alignment.wdl"
        with open(alignment_wdl) as f:
            content = f.read()

        assert "minimap2" in content, "alignment.wdl does not use minimap2"

    def test_alignment_wdl_has_docker_container(self):
        """alignment.wdl task specifies docker container"""
        alignment_wdl = TORCHBASE_ROOT / "workflows" / "builtin" / "tasks" / "alignment.wdl"
        with open(alignment_wdl) as f:
            content = f.read()

        assert "docker:" in content, "alignment.wdl does not specify docker container"

    def test_alignment_wdl_syntax_valid(self):
        """alignment.wdl passes miniwdl syntax validation"""
        alignment_wdl = TORCHBASE_ROOT / "workflows" / "builtin" / "tasks" / "alignment.wdl"

        result = subprocess.run(
            ["miniwdl", "check", str(alignment_wdl)],
            capture_output=True,
            text=True
        )

        assert result.returncode == 0, f"miniwdl check failed for alignment.wdl: {result.stderr}"


class TestProfileLookupTaskFile:
    """Test profile_lookup.wdl contains profile matching logic."""

    def test_profile_lookup_wdl_exists(self):
        """profile_lookup.wdl file exists in tasks directory"""
        profile_lookup_wdl = TORCHBASE_ROOT / "workflows" / "builtin" / "tasks" / "profile_lookup.wdl"
        assert profile_lookup_wdl.exists(), f"profile_lookup.wdl not found at {profile_lookup_wdl}"

    def test_profile_lookup_wdl_is_file(self):
        """profile_lookup.wdl is a regular file"""
        profile_lookup_wdl = TORCHBASE_ROOT / "workflows" / "builtin" / "tasks" / "profile_lookup.wdl"
        assert profile_lookup_wdl.is_file(), f"profile_lookup.wdl is not a file: {profile_lookup_wdl}"

    def test_profile_lookup_wdl_has_lookup_task(self):
        """profile_lookup.wdl contains profile lookup task"""
        profile_lookup_wdl = TORCHBASE_ROOT / "workflows" / "builtin" / "tasks" / "profile_lookup.wdl"
        with open(profile_lookup_wdl) as f:
            content = f.read()

        # Accept various task names: lookup_profile, match_profile, profile_lookup
        assert (
            "task lookup_profile" in content or
            "task match_profile" in content or
            "task profile_lookup" in content
        ), "profile_lookup.wdl does not contain profile lookup task"

    def test_profile_lookup_wdl_has_allele_calls_input(self):
        """profile lookup task accepts allele calls input"""
        profile_lookup_wdl = TORCHBASE_ROOT / "workflows" / "builtin" / "tasks" / "profile_lookup.wdl"
        with open(profile_lookup_wdl) as f:
            content = f.read()

        assert "allele_calls" in content or "allele_profile" in content, \
            "profile lookup task does not have allele calls input"

    def test_profile_lookup_wdl_has_profiles_table_input(self):
        """profile lookup task accepts profiles table input"""
        profile_lookup_wdl = TORCHBASE_ROOT / "workflows" / "builtin" / "tasks" / "profile_lookup.wdl"
        with open(profile_lookup_wdl) as f:
            content = f.read()

        assert "profiles" in content.lower(), "profile lookup task does not have profiles table input"

    def test_profile_lookup_wdl_has_profile_matching_logic(self):
        """profile lookup task contains profile matching logic"""
        profile_lookup_wdl = TORCHBASE_ROOT / "workflows" / "builtin" / "tasks" / "profile_lookup.wdl"
        with open(profile_lookup_wdl) as f:
            content = f.read()

        # Check for matching/lookup keywords
        assert any(keyword in content.lower() for keyword in ["match", "lookup", "search", "compare"]), \
            "profile lookup task does not contain profile matching logic"

    def test_profile_lookup_wdl_has_status_output(self):
        """profile lookup task outputs status (known/novel_profile/novel_allele)"""
        profile_lookup_wdl = TORCHBASE_ROOT / "workflows" / "builtin" / "tasks" / "profile_lookup.wdl"
        with open(profile_lookup_wdl) as f:
            content = f.read()

        assert "status" in content, "profile lookup task does not output status"

    def test_profile_lookup_wdl_handles_novel_profiles(self):
        """profile lookup task handles novel profiles"""
        profile_lookup_wdl = TORCHBASE_ROOT / "workflows" / "builtin" / "tasks" / "profile_lookup.wdl"
        with open(profile_lookup_wdl) as f:
            content = f.read()

        assert "novel" in content.lower(), "profile lookup task does not handle novel profiles"

    def test_profile_lookup_wdl_syntax_valid(self):
        """profile_lookup.wdl passes miniwdl syntax validation"""
        profile_lookup_wdl = TORCHBASE_ROOT / "workflows" / "builtin" / "tasks" / "profile_lookup.wdl"

        result = subprocess.run(
            ["miniwdl", "check", str(profile_lookup_wdl)],
            capture_output=True,
            text=True
        )

        assert result.returncode == 0, f"miniwdl check failed for profile_lookup.wdl: {result.stderr}"


class TestWDLImportStatements:
    """Test that tasks can be imported via WDL import statements."""

    def test_minhash_wdl_can_be_imported(self):
        """minhash.wdl can be imported by another WDL file"""
        # Create a temporary test workflow that imports minhash.wdl
        import tempfile
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir_path = Path(tmpdir)
            test_workflow = tmpdir_path / "test_import.wdl"

            # Create a minimal workflow that imports minhash tasks
            minhash_wdl = TORCHBASE_ROOT / "workflows" / "builtin" / "tasks" / "minhash.wdl"
            workflow_content = f"""version 1.0

import "{minhash_wdl}" as minhash

workflow test_import {{
    input {{
        File sequences
    }}

    call minhash.sketch_sequences {{
        input: sequences = sequences
    }}

    output {{
        File sketch = sketch_sequences.sketch
    }}
}}
"""
            with open(test_workflow, "w") as f:
                f.write(workflow_content)

            # Validate with miniwdl
            result = subprocess.run(
                ["miniwdl", "check", str(test_workflow)],
                capture_output=True,
                text=True
            )

            assert result.returncode == 0, f"Failed to import minhash.wdl: {result.stderr}"

    def test_alignment_wdl_can_be_imported(self):
        """alignment.wdl can be imported by another WDL file"""
        import tempfile
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir_path = Path(tmpdir)
            test_workflow = tmpdir_path / "test_import.wdl"

            # Create a minimal workflow that imports alignment tasks
            alignment_wdl = TORCHBASE_ROOT / "workflows" / "builtin" / "tasks" / "alignment.wdl"

            # Read alignment.wdl to determine the actual task name
            with open(alignment_wdl) as f:
                alignment_content = f.read()

            # Determine task name
            if "task align_sequences" in alignment_content:
                task_name = "align_sequences"
            elif "task minimap2_align" in alignment_content:
                task_name = "minimap2_align"
            else:
                task_name = "align_sequences"  # default fallback

            workflow_content = f"""version 1.0

import "{alignment_wdl}" as alignment

workflow test_import {{
    input {{
        File query
        File reference
    }}

    call alignment.{task_name} {{
        input:
            query = query,
            reference = reference
    }}

    output {{
        File alignment_output = {task_name}.alignment
    }}
}}
"""
            with open(test_workflow, "w") as f:
                f.write(workflow_content)

            # Validate with miniwdl
            result = subprocess.run(
                ["miniwdl", "check", str(test_workflow)],
                capture_output=True,
                text=True
            )

            assert result.returncode == 0, f"Failed to import alignment.wdl: {result.stderr}"

    def test_profile_lookup_wdl_can_be_imported(self):
        """profile_lookup.wdl can be imported by another WDL file"""
        import tempfile
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir_path = Path(tmpdir)
            test_workflow = tmpdir_path / "test_import.wdl"

            # Create a minimal workflow that imports profile_lookup tasks
            profile_lookup_wdl = TORCHBASE_ROOT / "workflows" / "builtin" / "tasks" / "profile_lookup.wdl"

            # Read profile_lookup.wdl to determine the actual task name
            with open(profile_lookup_wdl) as f:
                lookup_content = f.read()

            # Determine task name
            if "task lookup_profile" in lookup_content:
                task_name = "lookup_profile"
            elif "task match_profile" in lookup_content:
                task_name = "match_profile"
            elif "task profile_lookup" in lookup_content:
                task_name = "profile_lookup"
            else:
                task_name = "lookup_profile"  # default fallback

            workflow_content = f"""version 1.0

import "{profile_lookup_wdl}" as profile_lookup

workflow test_import {{
    input {{
        File allele_calls
        File profiles_table
    }}

    call profile_lookup.{task_name} {{
        input:
            allele_calls = allele_calls,
            profiles_table = profiles_table
    }}

    output {{
        String status = {task_name}.status
    }}
}}
"""
            with open(test_workflow, "w") as f:
                f.write(workflow_content)

            # Validate with miniwdl
            result = subprocess.run(
                ["miniwdl", "check", str(test_workflow)],
                capture_output=True,
                text=True
            )

            assert result.returncode == 0, f"Failed to import profile_lookup.wdl: {result.stderr}"


class TestAlignmentPresetConfiguration:
    """Test that minimap2 task correctly supports different presets."""

    def test_alignment_wdl_preset_parameter_is_string(self):
        """alignment task preset parameter is String type"""
        alignment_wdl = TORCHBASE_ROOT / "workflows" / "builtin" / "tasks" / "alignment.wdl"
        with open(alignment_wdl) as f:
            content = f.read()

        # Check for String preset parameter
        assert "String preset" in content, "preset parameter is not String type"

    def test_alignment_wdl_preset_has_default_value(self):
        """alignment task preset parameter has a default value"""
        alignment_wdl = TORCHBASE_ROOT / "workflows" / "builtin" / "tasks" / "alignment.wdl"
        with open(alignment_wdl) as f:
            content = f.read()

        # Check for default value assignment
        assert "preset =" in content or "String preset" in content, \
            "preset parameter does not have default value"

    def test_alignment_wdl_uses_preset_in_command(self):
        """alignment task uses preset parameter in minimap2 command"""
        alignment_wdl = TORCHBASE_ROOT / "workflows" / "builtin" / "tasks" / "alignment.wdl"
        with open(alignment_wdl) as f:
            content = f.read()

        # Check that preset is used in command (via -x flag for minimap2)
        assert ("-x" in content and "preset" in content) or "~{preset}" in content, \
            "preset parameter is not used in minimap2 command"

    def test_alignment_wdl_documents_preset_options(self):
        """alignment task documents available preset options"""
        alignment_wdl = TORCHBASE_ROOT / "workflows" / "builtin" / "tasks" / "alignment.wdl"
        with open(alignment_wdl) as f:
            content = f.read()

        # Check for documentation of preset options (in comments or parameter description)
        preset_count = sum(1 for preset in ["asm20", "asm5", "eqx"] if preset in content)
        assert preset_count >= 2, "alignment task does not document preset options (asm20, asm5, eqx)"


class TestTasksIntegration:
    """Test that all tasks work together for end-to-end workflow."""

    def test_all_three_task_files_exist(self):
        """All three task files (minhash, alignment, profile_lookup) exist"""
        tasks_dir = TORCHBASE_ROOT / "workflows" / "builtin" / "tasks"

        minhash_wdl = tasks_dir / "minhash.wdl"
        alignment_wdl = tasks_dir / "alignment.wdl"
        profile_lookup_wdl = tasks_dir / "profile_lookup.wdl"

        assert minhash_wdl.exists(), "minhash.wdl does not exist"
        assert alignment_wdl.exists(), "alignment.wdl does not exist"
        assert profile_lookup_wdl.exists(), "profile_lookup.wdl does not exist"

    def test_all_task_files_are_valid_wdl(self):
        """All task files pass miniwdl syntax validation"""
        tasks_dir = TORCHBASE_ROOT / "workflows" / "builtin" / "tasks"

        task_files = [
            tasks_dir / "minhash.wdl",
            tasks_dir / "alignment.wdl",
            tasks_dir / "profile_lookup.wdl"
        ]

        for task_file in task_files:
            result = subprocess.run(
                ["miniwdl", "check", str(task_file)],
                capture_output=True,
                text=True
            )

            assert result.returncode == 0, f"miniwdl check failed for {task_file.name}: {result.stderr}"

    def test_tasks_can_be_imported_together(self):
        """All task files can be imported together in a single workflow"""
        import tempfile
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir_path = Path(tmpdir)
            test_workflow = tmpdir_path / "test_combined_import.wdl"

            tasks_dir = TORCHBASE_ROOT / "workflows" / "builtin" / "tasks"
            minhash_wdl = tasks_dir / "minhash.wdl"
            alignment_wdl = tasks_dir / "alignment.wdl"
            profile_lookup_wdl = tasks_dir / "profile_lookup.wdl"

            workflow_content = f"""version 1.0

import "{minhash_wdl}" as minhash
import "{alignment_wdl}" as alignment
import "{profile_lookup_wdl}" as profile_lookup

workflow test_combined_import {{
    input {{
        File sequences
        File allele_fasta
        File profiles_table
    }}

    # This workflow just tests that all imports work together
    # It doesn't need to actually call the tasks

    output {{
        String status = "imports_successful"
    }}
}}
"""
            with open(test_workflow, "w") as f:
                f.write(workflow_content)

            # Validate with miniwdl
            result = subprocess.run(
                ["miniwdl", "check", str(test_workflow)],
                capture_output=True,
                text=True
            )

            assert result.returncode == 0, f"Failed to import all task files together: {result.stderr}"
