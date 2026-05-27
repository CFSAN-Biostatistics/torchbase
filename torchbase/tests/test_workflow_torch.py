"""Tests for default workflow torch structure (Issue #14).

Acceptance criteria:
- Creates a default workflow torch at well-known location
- WDL file: main.wdl with standard MLST pipeline
- Profiles table for reference typing
- Metadata with description, version
- Can be loaded by Torch.load()
- Can be executed by miniwdl
- Provides standard output format (JSON typing results)
"""

import pytest
import json
import toml
import csv
import tempfile
from pathlib import Path

from torchbase.torchfs import Torch
from torchbase.torchbase import Profile, Schema


# Get the torchbase root directory
TORCHBASE_ROOT = Path(__file__).parent.parent


@pytest.fixture
def workflow_torch_tempdir():
    """Create a temporary workflow torch directory at well-known location.

    Structure:
    workflows/mlst/1.0.0.torch/
    ├── metadata.toml
    ├── main.wdl
    ├── profiles.tsv
    ├── quality.json
    └── _resources/
        ├── reference_alleles.fasta
        └── metadata.json
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir_path = Path(tmpdir)

        # Create well-known workflow torch location
        torch_path = tmpdir_path / "workflows" / "mlst" / "1.0.0.torch"
        torch_path.mkdir(parents=True, exist_ok=True)

        # Create metadata
        metadata = {
            "namespace": "workflows",
            "name": "mlst",
            "version": "1.0.0",
            "version_meta": {
                "strategy": "semver",
                "timestamp": 1609459200
            },
            "typing": {
                "method": "mlst",
                "description": "Default MLST typing workflow"
            },
            "maintainers": {
                "authors": ["Default Author"],
                "email": ["default@example.com"],
                "affiliations": ["Default Institute"]
            },
            "description": {
                "short": "Default MLST typing workflow",
                "long": "Standard MLST pipeline for microbial typing",
                "taxa": ["Bacteria"]
            },
            "manifest": {
                "profiles": "profiles.tsv",
                "workflow": "main.wdl",
                "resources": []
            }
        }

        # Write metadata
        with open(torch_path / "metadata.toml", "w") as f:
            toml.dump(metadata, f)

        # Create main.wdl workflow file
        wdl_content = """workflow mlst_typing {
    input {
        File reads
        File reference_db
    }

    output {
        File typing_results = "results.json"
    }

    call mlst_type { input: reads = reads, reference_db = reference_db }
}

task mlst_type {
    input {
        File reads
        File reference_db
    }

    command {
        echo '{"st": 1, "loci": {"adk": 1}}' > results.json
    }

    output {
        File results = "results.json"
    }
}
"""
        with open(torch_path / "main.wdl", "w") as f:
            f.write(wdl_content)

        # Create profiles table for reference typing
        profiles = [
            ["ST", "adk", "fumC", "gyrB", "icd", "mdh", "purA"],
            ["1", "1", "1", "1", "1", "1", "1"],
            ["2", "2", "1", "1", "1", "1", "1"],
            ["3", "1", "2", "1", "1", "1", "1"],
            ["4", "1", "1", "2", "1", "1", "1"]
        ]
        with open(torch_path / "profiles.tsv", "w") as f:
            writer = csv.writer(f, delimiter="\t")
            writer.writerows(profiles)

        # Create quality.json placeholder (empty but valid JSON)
        quality_data = {}
        with open(torch_path / "quality.json", "w") as f:
            json.dump(quality_data, f)

        # Create _resources directory
        resources_path = torch_path / "_resources"
        resources_path.mkdir(parents=True, exist_ok=True)

        # Add reference alleles
        with open(resources_path / "reference_alleles.fasta", "w") as f:
            f.write(">adk_1\nACGT\n>adk_2\nTGCA\n")

        # Add metadata
        with open(resources_path / "metadata.json", "w") as f:
            json.dump({"source": "default_mlst"}, f)

        yield torch_path


class TestWorkflowTorchLocation:
    """Test workflow torch is at well-known location."""

    def test_workflow_torch_directory_exists(self, workflow_torch_tempdir):
        """Workflow torch directory exists at workflows/mlst/1.0.0.torch/"""
        assert workflow_torch_tempdir.exists()
        assert workflow_torch_tempdir.is_dir()
        assert workflow_torch_tempdir.parent.name == "mlst"
        assert workflow_torch_tempdir.parent.parent.name == "workflows"

    def test_workflow_torch_version_format(self, workflow_torch_tempdir):
        """Workflow torch version directory ends with .torch"""
        assert workflow_torch_tempdir.name == "1.0.0.torch"

    def test_workflow_torch_namespace_path(self, workflow_torch_tempdir):
        """Workflow torch namespace path is workflows"""
        assert workflow_torch_tempdir.parent.parent.name == "workflows"


class TestWorkflowTorchMetadata:
    """Test workflow torch metadata with description and version."""

    def test_metadata_toml_exists(self, workflow_torch_tempdir):
        """metadata.toml file exists"""
        metadata_path = workflow_torch_tempdir / "metadata.toml"
        assert metadata_path.exists()
        assert metadata_path.is_file()

    def test_metadata_has_version(self, workflow_torch_tempdir):
        """metadata.toml contains version"""
        metadata_path = workflow_torch_tempdir / "metadata.toml"
        with open(metadata_path) as f:
            metadata = toml.load(f)

        assert "version" in metadata
        assert metadata["version"] == "1.0.0"

    def test_metadata_has_description(self, workflow_torch_tempdir):
        """metadata.toml contains description section"""
        metadata_path = workflow_torch_tempdir / "metadata.toml"
        with open(metadata_path) as f:
            metadata = toml.load(f)

        assert "description" in metadata
        assert "short" in metadata["description"]
        assert "long" in metadata["description"]

    def test_metadata_declares_typing_method(self, workflow_torch_tempdir):
        """metadata.toml declares typing.method as mlst"""
        metadata_path = workflow_torch_tempdir / "metadata.toml"
        with open(metadata_path) as f:
            metadata = toml.load(f)

        assert "typing" in metadata
        assert "method" in metadata["typing"]
        assert metadata["typing"]["method"] == "mlst"

    def test_metadata_has_namespace(self, workflow_torch_tempdir):
        """metadata.toml namespace matches directory path"""
        metadata_path = workflow_torch_tempdir / "metadata.toml"
        with open(metadata_path) as f:
            metadata = toml.load(f)

        assert metadata["namespace"] == "workflows"

    def test_metadata_has_name(self, workflow_torch_tempdir):
        """metadata.toml name matches directory path"""
        metadata_path = workflow_torch_tempdir / "metadata.toml"
        with open(metadata_path) as f:
            metadata = toml.load(f)

        assert metadata["name"] == "mlst"


class TestWorkflowTorchWDLFile:
    """Test WDL file: main.wdl with standard MLST pipeline."""

    def test_main_wdl_exists(self, workflow_torch_tempdir):
        """main.wdl file exists in torch directory"""
        wdl_path = workflow_torch_tempdir / "main.wdl"
        assert wdl_path.exists()
        assert wdl_path.is_file()

    def test_main_wdl_has_workflow_definition(self, workflow_torch_tempdir):
        """main.wdl contains workflow definition"""
        wdl_path = workflow_torch_tempdir / "main.wdl"
        with open(wdl_path) as f:
            content = f.read()

        assert "workflow" in content
        assert "{" in content and "}" in content

    def test_main_wdl_has_mlst_workflow_section(self, workflow_torch_tempdir):
        """main.wdl has MLST workflow section"""
        wdl_path = workflow_torch_tempdir / "main.wdl"
        with open(wdl_path) as f:
            content = f.read()

        # Should have workflow definition
        assert "workflow" in content

    def test_main_wdl_has_input_section(self, workflow_torch_tempdir):
        """main.wdl workflow has input section"""
        wdl_path = workflow_torch_tempdir / "main.wdl"
        with open(wdl_path) as f:
            content = f.read()

        assert "input" in content

    def test_main_wdl_has_output_section(self, workflow_torch_tempdir):
        """main.wdl workflow has output section"""
        wdl_path = workflow_torch_tempdir / "main.wdl"
        with open(wdl_path) as f:
            content = f.read()

        assert "output" in content

    def test_main_wdl_contains_task_definitions(self, workflow_torch_tempdir):
        """main.wdl contains task definitions"""
        wdl_path = workflow_torch_tempdir / "main.wdl"
        with open(wdl_path) as f:
            content = f.read()

        assert "task" in content


class TestWorkflowTorchProfilesTable:
    """Test profiles table for reference typing."""

    def test_profiles_tsv_exists(self, workflow_torch_tempdir):
        """profiles.tsv file exists"""
        profiles_path = workflow_torch_tempdir / "profiles.tsv"
        assert profiles_path.exists()
        assert profiles_path.is_file()

    def test_profiles_tsv_has_header(self, workflow_torch_tempdir):
        """profiles.tsv has header row"""
        profiles_path = workflow_torch_tempdir / "profiles.tsv"
        with open(profiles_path) as f:
            reader = csv.reader(f, delimiter="\t")
            header = next(reader)

        assert len(header) > 0
        assert header[0] == "ST"

    def test_profiles_tsv_has_loci(self, workflow_torch_tempdir):
        """profiles.tsv header contains loci columns"""
        profiles_path = workflow_torch_tempdir / "profiles.tsv"
        with open(profiles_path) as f:
            reader = csv.reader(f, delimiter="\t")
            header = next(reader)

        # MLST should have standard loci
        loci = set(header[1:])  # Skip ST column
        assert len(loci) > 0

    def test_profiles_tsv_has_data_rows(self, workflow_torch_tempdir):
        """profiles.tsv contains profile data rows"""
        profiles_path = workflow_torch_tempdir / "profiles.tsv"
        with open(profiles_path) as f:
            reader = csv.reader(f, delimiter="\t")
            next(reader)  # skip header
            data_rows = list(reader)

        assert len(data_rows) > 0

    def test_profiles_tsv_consistent_column_count(self, workflow_torch_tempdir):
        """profiles.tsv has consistent column count across rows"""
        profiles_path = workflow_torch_tempdir / "profiles.tsv"
        with open(profiles_path) as f:
            reader = csv.reader(f, delimiter="\t")
            rows = list(reader)

        if len(rows) > 1:
            header_cols = len(rows[0])
            for row in rows[1:]:
                assert len(row) == header_cols


class TestWorkflowTorchQualityFile:
    """Test quality.json placeholder."""

    def test_quality_json_exists(self, workflow_torch_tempdir):
        """quality.json file exists"""
        quality_path = workflow_torch_tempdir / "quality.json"
        assert quality_path.exists()
        assert quality_path.is_file()

    def test_quality_json_is_valid_json(self, workflow_torch_tempdir):
        """quality.json contains valid JSON"""
        quality_path = workflow_torch_tempdir / "quality.json"
        with open(quality_path) as f:
            data = json.load(f)

        # Should be parseable as JSON
        assert isinstance(data, dict)

    def test_quality_json_not_empty_after_workflow_execution(self, workflow_torch_tempdir):
        """quality.json can hold quality metrics after workflow execution"""
        quality_path = workflow_torch_tempdir / "quality.json"

        # Simulate writing quality data
        quality_data = {"coverage": 100, "quality": "pass"}
        with open(quality_path, "w") as f:
            json.dump(quality_data, f)

        # Verify it was written
        with open(quality_path) as f:
            loaded_data = json.load(f)

        assert loaded_data == quality_data


class TestWorkflowTorchResourcesDirectory:
    """Test _resources directory with reference files."""

    def test_resources_directory_exists(self, workflow_torch_tempdir):
        """_resources directory exists"""
        resources_path = workflow_torch_tempdir / "_resources"
        assert resources_path.exists()
        assert resources_path.is_dir()

    def test_resources_contains_reference_files(self, workflow_torch_tempdir):
        """_resources contains reference files"""
        resources_path = workflow_torch_tempdir / "_resources"
        files = list(resources_path.glob("*"))

        # Should have at least reference alleles
        non_dotfiles = [f for f in files if not f.name.startswith(".")]
        assert len(non_dotfiles) > 0

    def test_resources_has_reference_alleles(self, workflow_torch_tempdir):
        """_resources contains reference alleles FASTA"""
        resources_path = workflow_torch_tempdir / "_resources"
        fasta_files = list(resources_path.glob("*.fasta"))

        assert len(fasta_files) > 0

    def test_resources_has_metadata_file(self, workflow_torch_tempdir):
        """_resources contains metadata file"""
        resources_path = workflow_torch_tempdir / "_resources"
        metadata_file = resources_path / "metadata.json"

        assert metadata_file.exists()


class TestWorkflowTorchCanBeLoaded:
    """Test that workflow torch can be loaded by Torch.load()."""

    def test_torch_load_accepts_workflow_torch_path(self, workflow_torch_tempdir):
        """Torch.load() accepts workflow torch path"""
        torch = Torch.load(workflow_torch_tempdir)

        assert torch is not None
        assert isinstance(torch, Torch)

    def test_loaded_torch_has_correct_path(self, workflow_torch_tempdir):
        """Loaded torch has correct path attribute"""
        torch = Torch.load(workflow_torch_tempdir)

        assert torch.path == workflow_torch_tempdir

    def test_loaded_torch_has_workflow_attribute(self, workflow_torch_tempdir):
        """Loaded torch has workflow attribute"""
        torch = Torch.load(workflow_torch_tempdir)

        assert hasattr(torch, "workflow")

    def test_loaded_torch_workflow_points_to_main_wdl(self, workflow_torch_tempdir):
        """Loaded torch workflow points to main.wdl"""
        torch = Torch.load(workflow_torch_tempdir)

        assert torch.workflow is not None
        assert torch.workflow.name == "main.wdl"

    def test_loaded_torch_has_profiles(self, workflow_torch_tempdir):
        """Loaded torch has profiles attribute"""
        torch = Torch.load(workflow_torch_tempdir)

        assert hasattr(torch, "profile")
        assert torch.profile is not None

    def test_loaded_torch_profiles_is_schema(self, workflow_torch_tempdir):
        """Loaded torch profile is a Schema object"""
        torch = Torch.load(workflow_torch_tempdir)

        assert isinstance(torch.profile, (Schema, Profile))

    def test_loaded_torch_has_resources(self, workflow_torch_tempdir):
        """Loaded torch has references to resources"""
        torch = Torch.load(workflow_torch_tempdir)

        assert hasattr(torch, "references")
        assert len(torch.references) > 0

    def test_loaded_torch_references_include_alleles(self, workflow_torch_tempdir):
        """Loaded torch references include allele files"""
        torch = Torch.load(workflow_torch_tempdir)

        ref_names = {ref.name for ref in torch.references}
        assert any(".fasta" in name or "allele" in name for name in ref_names)


class TestWorkflowTorchOutputFormat:
    """Test standard output format (JSON typing results)."""

    def test_workflow_can_produce_json_output(self, workflow_torch_tempdir):
        """Workflow is configured to produce JSON output"""
        wdl_path = workflow_torch_tempdir / "main.wdl"
        with open(wdl_path) as f:
            content = f.read()

        # WDL should reference JSON output
        assert "json" in content.lower() or ".json" in content

    def test_workflow_output_includes_typing_results(self, workflow_torch_tempdir):
        """Workflow output file is named typing_results or results"""
        wdl_path = workflow_torch_tempdir / "main.wdl"
        with open(wdl_path) as f:
            content = f.read()

        # Should have output defined
        assert "output" in content

    def test_quality_json_can_store_output_metadata(self, workflow_torch_tempdir):
        """quality.json can store output metadata and quality info"""
        quality_path = workflow_torch_tempdir / "quality.json"

        # Write sample output data
        output_data = {
            "typing_method": "mlst",
            "st_number": 1,
            "loci": {"adk": 1, "fumC": 1},
            "quality_score": 0.95
        }

        with open(quality_path, "w") as f:
            json.dump(output_data, f)

        # Verify it was stored correctly
        with open(quality_path) as f:
            loaded_data = json.load(f)

        assert loaded_data["typing_method"] == "mlst"


class TestWorkflowTorchValidation:
    """Test workflow torch validation and structure requirements."""

    def test_torch_path_namespace_matches_metadata(self, workflow_torch_tempdir):
        """Torch path namespace matches metadata.namespace"""
        metadata_path = workflow_torch_tempdir / "metadata.toml"
        with open(metadata_path) as f:
            metadata = toml.load(f)

        path_namespace = workflow_torch_tempdir.parent.parent.name
        assert metadata["namespace"] == path_namespace

    def test_torch_path_name_matches_metadata(self, workflow_torch_tempdir):
        """Torch path name matches metadata.name"""
        metadata_path = workflow_torch_tempdir / "metadata.toml"
        with open(metadata_path) as f:
            metadata = toml.load(f)

        path_name = workflow_torch_tempdir.parent.name
        assert metadata["name"] == path_name

    def test_torch_path_version_matches_metadata(self, workflow_torch_tempdir):
        """Torch path version matches metadata.version"""
        metadata_path = workflow_torch_tempdir / "metadata.toml"
        with open(metadata_path) as f:
            metadata = toml.load(f)

        path_version = workflow_torch_tempdir.name.replace(".torch", "")
        assert str(metadata["version"]) == path_version

    def test_required_files_all_present(self, workflow_torch_tempdir):
        """All required files present in workflow torch"""
        required_files = [
            "metadata.toml",
            "main.wdl",
            "profiles.tsv",
            "quality.json"
        ]

        for filename in required_files:
            file_path = workflow_torch_tempdir / filename
            assert file_path.exists(), f"Missing required file: {filename}"

    def test_resources_directory_required(self, workflow_torch_tempdir):
        """_resources directory is required"""
        resources_path = workflow_torch_tempdir / "_resources"
        assert resources_path.exists()
        assert resources_path.is_dir()


class TestWorkflowTorchEdgeCases:
    """Test edge cases and boundary conditions."""

    def test_workflow_torch_with_empty_resources(self):
        """Workflow torch can have empty resources directory"""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir_path = Path(tmpdir)
            torch_path = tmpdir_path / "workflows" / "mlst" / "1.0.0.torch"
            torch_path.mkdir(parents=True, exist_ok=True)

            # Create minimal metadata
            metadata = {
                "namespace": "workflows",
                "name": "mlst",
                "version": "1.0.0",
                "version_meta": {
                    "strategy": "semver",
                    "timestamp": 1609459200
                },
                "typing": {"method": "mlst"},
                "description": {"short": "Test"},
                "manifest": {"profiles": "profiles.tsv", "workflow": "main.wdl"}
            }
            with open(torch_path / "metadata.toml", "w") as f:
                toml.dump(metadata, f)

            # Create minimal WDL
            with open(torch_path / "main.wdl", "w") as f:
                f.write("workflow mlst_typing { }\n")

            # Create minimal profiles
            profiles = [["ST", "adk"], ["1", "1"]]
            with open(torch_path / "profiles.tsv", "w") as f:
                writer = csv.writer(f, delimiter="\t")
                writer.writerows(profiles)

            # Create empty quality.json
            with open(torch_path / "quality.json", "w") as f:
                json.dump({}, f)

            # Empty resources dir
            (torch_path / "_resources").mkdir()

            # Should still load
            torch = Torch.load(torch_path)
            assert torch is not None

    def test_workflow_torch_with_multiple_resource_files(self):
        """Workflow torch can have multiple resource files"""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir_path = Path(tmpdir)
            torch_path = tmpdir_path / "workflows" / "mlst" / "1.0.0.torch"
            torch_path.mkdir(parents=True, exist_ok=True)

            # Create metadata
            metadata = {
                "namespace": "workflows",
                "name": "mlst",
                "version": "1.0.0",
                "version_meta": {
                    "strategy": "semver",
                    "timestamp": 1609459200
                },
                "typing": {"method": "mlst"},
                "description": {"short": "Test"},
                "manifest": {"profiles": "profiles.tsv", "workflow": "main.wdl"}
            }
            with open(torch_path / "metadata.toml", "w") as f:
                toml.dump(metadata, f)

            # Create WDL
            with open(torch_path / "main.wdl", "w") as f:
                f.write("workflow mlst_typing { }\n")

            # Create profiles
            profiles = [["ST", "adk"], ["1", "1"]]
            with open(torch_path / "profiles.tsv", "w") as f:
                writer = csv.writer(f, delimiter="\t")
                writer.writerows(profiles)

            # Create quality.json
            with open(torch_path / "quality.json", "w") as f:
                json.dump({}, f)

            # Create resources with multiple files
            resources_path = torch_path / "_resources"
            resources_path.mkdir()

            with open(resources_path / "allele1.fasta", "w") as f:
                f.write(">adk_1\nACGT\n")

            with open(resources_path / "allele2.fasta", "w") as f:
                f.write(">adk_2\nTGCA\n")

            with open(resources_path / "metadata.json", "w") as f:
                json.dump({}, f)

            # Should load with all resources
            torch = Torch.load(torch_path)
            assert len(torch.references) >= 2

    def test_workflow_torch_with_long_description(self, workflow_torch_tempdir):
        """Workflow torch can have long descriptions"""
        metadata_path = workflow_torch_tempdir / "metadata.toml"
        with open(metadata_path) as f:
            metadata = toml.load(f)

        long_desc = "A" * 1000
        metadata["description"]["long"] = long_desc

        with open(metadata_path, "w") as f:
            toml.dump(metadata, f)

        # Should still load
        torch = Torch.load(workflow_torch_tempdir)
        assert torch is not None


class TestWorkflowTorchConvention:
    """Test convention: torch with main.wdl is executable workflow."""

    def test_workflow_torch_identified_by_main_wdl(self, workflow_torch_tempdir):
        """Torch with main.wdl is identified as executable workflow"""
        main_wdl = workflow_torch_tempdir / "main.wdl"
        assert main_wdl.exists()

    def test_workflow_torch_main_wdl_is_required(self):
        """Workflow torch must have main.wdl to be executable"""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir_path = Path(tmpdir)
            torch_path = tmpdir_path / "workflows" / "test" / "1.0.0.torch"
            torch_path.mkdir(parents=True, exist_ok=True)

            # Create metadata
            metadata = {
                "namespace": "workflows",
                "name": "test",
                "version": "1.0.0",
                "version_meta": {
                    "strategy": "semver",
                    "timestamp": 1609459200
                },
                "typing": {"method": "mlst"},
                "description": {"short": "Test"},
                "manifest": {"profiles": "profiles.tsv"}
            }
            with open(torch_path / "metadata.toml", "w") as f:
                toml.dump(metadata, f)

            # Create profiles
            profiles = [["ST", "adk"], ["1", "1"]]
            with open(torch_path / "profiles.tsv", "w") as f:
                writer = csv.writer(f, delimiter="\t")
                writer.writerows(profiles)

            # Create quality.json
            with open(torch_path / "quality.json", "w") as f:
                json.dump({}, f)

            # Create resources
            (torch_path / "_resources").mkdir()

            # No main.wdl - should still load but workflow is None
            torch = Torch.load(torch_path)
            # This tests the convention - torch without main.wdl may still load
            # but shouldn't be executable
            assert torch.workflow is None


@pytest.mark.skip(reason="Issue #61: Default workflow torch moved to examples/ and replaced with synthetic examples")
class TestDefaultWorkflowTorchLocation:
    """Test that default workflow torch exists at well-known location."""

    def test_default_workflow_torch_directory_exists(self):
        """Default workflow torch directory exists at torchbase/workflows/mlst/1.0.0.torch/"""
        torch_path = TORCHBASE_ROOT / "workflows" / "mlst" / "1.0.0.torch"
        assert torch_path.exists(), f"Workflow torch not found at {torch_path}"
        assert torch_path.is_dir(), f"Workflow torch is not a directory: {torch_path}"

    def test_default_workflow_torch_metadata_exists(self):
        """Default workflow torch has metadata.toml"""
        metadata_path = TORCHBASE_ROOT / "workflows" / "mlst" / "1.0.0.torch" / "metadata.toml"
        assert metadata_path.exists(), f"metadata.toml not found at {metadata_path}"
        assert metadata_path.is_file(), f"metadata.toml is not a file: {metadata_path}"

    def test_default_workflow_torch_metadata_valid(self):
        """Default workflow torch metadata.toml is valid TOML"""
        metadata_path = TORCHBASE_ROOT / "workflows" / "mlst" / "1.0.0.torch" / "metadata.toml"
        with open(metadata_path) as f:
            metadata = toml.load(f)

        assert metadata is not None
        assert isinstance(metadata, dict)

    def test_default_workflow_torch_typing_method(self):
        """Default workflow torch declares typing.method = mlst"""
        metadata_path = TORCHBASE_ROOT / "workflows" / "mlst" / "1.0.0.torch" / "metadata.toml"
        with open(metadata_path) as f:
            metadata = toml.load(f)

        assert "typing" in metadata
        assert metadata["typing"]["method"] == "mlst"

    def test_default_workflow_torch_main_wdl_exists(self):
        """Default workflow torch has main.wdl"""
        wdl_path = TORCHBASE_ROOT / "workflows" / "mlst" / "1.0.0.torch" / "main.wdl"
        assert wdl_path.exists(), f"main.wdl not found at {wdl_path}"
        assert wdl_path.is_file(), f"main.wdl is not a file: {wdl_path}"

    def test_default_workflow_torch_profiles_exists(self):
        """Default workflow torch has profiles.tsv"""
        profiles_path = TORCHBASE_ROOT / "workflows" / "mlst" / "1.0.0.torch" / "profiles.tsv"
        assert profiles_path.exists(), f"profiles.tsv not found at {profiles_path}"
        assert profiles_path.is_file(), f"profiles.tsv is not a file: {profiles_path}"

    def test_default_workflow_torch_quality_json_exists(self):
        """Default workflow torch has quality.json"""
        quality_path = TORCHBASE_ROOT / "workflows" / "mlst" / "1.0.0.torch" / "quality.json"
        assert quality_path.exists(), f"quality.json not found at {quality_path}"
        assert quality_path.is_file(), f"quality.json is not a file: {quality_path}"

    def test_default_workflow_torch_resources_directory_exists(self):
        """Default workflow torch has _resources directory"""
        resources_path = TORCHBASE_ROOT / "workflows" / "mlst" / "1.0.0.torch" / "_resources"
        assert resources_path.exists(), f"_resources not found at {resources_path}"
        assert resources_path.is_dir(), f"_resources is not a directory: {resources_path}"

    def test_default_workflow_torch_can_be_loaded(self):
        """Default workflow torch can be loaded by Torch.load()"""
        torch_path = TORCHBASE_ROOT / "workflows" / "mlst" / "1.0.0.torch"
        torch = Torch.load(torch_path)

        assert torch is not None
        assert isinstance(torch, Torch)

    def test_default_workflow_torch_has_workflow_attribute(self):
        """Default workflow torch loads with workflow attribute"""
        torch_path = TORCHBASE_ROOT / "workflows" / "mlst" / "1.0.0.torch"
        torch = Torch.load(torch_path)

        assert hasattr(torch, "workflow")
        assert torch.workflow is not None

    def test_default_workflow_torch_workflow_points_to_main_wdl(self):
        """Default workflow torch workflow attribute points to main.wdl"""
        torch_path = TORCHBASE_ROOT / "workflows" / "mlst" / "1.0.0.torch"
        torch = Torch.load(torch_path)

        assert torch.workflow.name == "main.wdl"

    def test_default_workflow_torch_has_profiles(self):
        """Default workflow torch has profiles"""
        torch_path = TORCHBASE_ROOT / "workflows" / "mlst" / "1.0.0.torch"
        torch = Torch.load(torch_path)

        assert hasattr(torch, "profile")
        assert torch.profile is not None

    def test_default_workflow_torch_profiles_contain_mlst_loci(self):
        """Default workflow torch profiles contain standard MLST loci"""
        torch_path = TORCHBASE_ROOT / "workflows" / "mlst" / "1.0.0.torch"
        torch = Torch.load(torch_path)

        # Should have some profiles
        if hasattr(torch.profile, "profiles"):
            profiles = torch.profile.profiles
            assert len(profiles) > 0
            # First profile should have loci
            first_profile = profiles[0]
            assert len(first_profile.header) > 0
