"""Pytest configuration and shared fixtures."""

import pytest
from pathlib import Path
import tempfile
import toml
import csv
from unittest import mock
import os

# Mock ipyfs if it's not available
try:
    import ipyfs
except ImportError:
    ipyfs_mock = mock.MagicMock()
    ipyfs_mock.Cat = mock.MagicMock(return_value=mock.MagicMock())
    import sys
    sys.modules['ipyfs'] = ipyfs_mock


@pytest.fixture(scope="session", autouse=True)
def configure_miniwdl_local_backend():
    """Configure miniwdl to use local backend instead of Docker for tests.

    This allows workflow tests to run in environments without Docker.
    """
    # Create a temporary config file for miniwdl
    config_dir = Path(tempfile.gettempdir()) / "miniwdl_test_config"
    config_dir.mkdir(exist_ok=True)
    config_file = config_dir / "miniwdl.toml"

    config_content = """
[task_runtime]
docker = false

[runtime]
backends = ["local"]

[exe]
allow_docker_fallback = false
"""

    config_file.write_text(config_content)

    # Set environment variable to use this config
    os.environ["WDL_CFG_PATH"] = str(config_file)

    yield

    # Cleanup
    if config_file.exists():
        config_file.unlink()
    if config_dir.exists():
        try:
            config_dir.rmdir()
        except OSError:
            pass


@pytest.fixture
def multi_scheme_torch_tempdir():
    """Create a temporary multi-scheme torch directory structure.

    Structure:
    test_namespace/test_torch/1.0.0.torch/
    ├── metadata.toml
    ├── schemes/
    │   ├── ecoli/
    │   │   ├── profiles.tsv
    │   │   └── alleles/
    │   │       ├── dinB.fasta
    │   │       └── icdA.fasta
    │   └── salmonella/
    │       ├── profiles.tsv
    │       └── alleles/
    │           ├── adk.fasta
    │           └── fumC.fasta
    └── _resources/  (optional for single-scheme compat)  # noqa: E501
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir_path = Path(tmpdir)
        torch_path = (tmpdir_path / "test_namespace" / "test_torch" /
                      "1.0.0.torch")
        torch_path.mkdir(parents=True, exist_ok=True)

        # Create metadata
        metadata = {
            "namespace": "test_namespace",
            "name": "test_torch",
            "version": "1.0.0",
            "version_meta": {
                "strategy": "content-hash",
                "timestamp": 1609459200
            },
            "maintainers": {
                "authors": ["Test Author"],
                "email": ["test@example.com"],
                "affiliations": ["Test University"]
            },
            "description": {
                "short": "Test multi-scheme torch",
                "long": ("Multi-scheme test torch with E. coli and "  # noqa
                         "Salmonella"),
                "taxa": [
                    "Escherichia coli",
                    "Salmonella enterica"
                ]
            },
            "schemes": {
                "ecoli": {"organism": "Escherichia coli"},
                "salmonella": {"organism": "Salmonella enterica"}
            },
            "manifest": {
                "resources": []
            }
        }

        # Write metadata
        with open(torch_path / "metadata.toml", "w") as f:
            toml.dump(metadata, f)

        # Create E. coli scheme
        ecoli_path = torch_path / "schemes" / "ecoli"
        ecoli_path.mkdir(parents=True, exist_ok=True)

        # E. coli profiles
        ecoli_profiles = [
            ["ST", "dinB", "icdA"],
            ["1", "1", "1"],
            ["2", "2", "2"],
            ["3", "3", "1"]
        ]
        with open(ecoli_path / "profiles.tsv", "w") as f:
            writer = csv.writer(f, delimiter="\t")
            writer.writerows(ecoli_profiles)

        # E. coli alleles
        alleles_path = ecoli_path / "alleles"
        alleles_path.mkdir(parents=True, exist_ok=True)

        with open(alleles_path / "dinB.fasta", "w") as f:
            f.write(">dinB_1\nACGT\n>dinB_2\nTGCA\n")

        with open(alleles_path / "icdA.fasta", "w") as f:
            f.write(">icdA_1\nGATC\n")

        # Create Salmonella scheme
        salmonella_path = torch_path / "schemes" / "salmonella"
        salmonella_path.mkdir(parents=True, exist_ok=True)

        # Salmonella profiles
        salmonella_profiles = [
            ["ST", "adk", "fumC"],
            ["1", "1", "1"],
            ["2", "1", "2"],
            ["3", "2", "1"]
        ]
        with open(salmonella_path / "profiles.tsv", "w") as f:
            writer = csv.writer(f, delimiter="\t")
            writer.writerows(salmonella_profiles)

        # Salmonella alleles
        alleles_path = salmonella_path / "alleles"
        alleles_path.mkdir(parents=True, exist_ok=True)

        with open(alleles_path / "adk.fasta", "w") as f:
            f.write(">adk_1\nCCCC\n>adk_2\nGGGG\n")

        with open(alleles_path / "fumC.fasta", "w") as f:
            f.write(">fumC_1\nAAAA\n>fumC_2\nTTTT\n")

        yield torch_path


@pytest.fixture
def single_scheme_torch_tempdir():
    """Create a temporary single-scheme torch (backward compatibility).

    Structure (old format without schemes/ subdirectory):
    test_namespace/legacy_torch/1.0.0.torch/
    ├── metadata.toml
    ├── legacy.profiles.tsv
    └── _resources/
        ├── dinB.fasta
        └── icdA.fasta
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir_path = Path(tmpdir)
        torch_path = (tmpdir_path / "test_namespace" / "legacy_torch" /
                      "1.0.0.torch")
        torch_path.mkdir(parents=True, exist_ok=True)

        # Create metadata (legacy format without schemes section)
        metadata = {
            "namespace": "test_namespace",
            "name": "legacy_torch",
            "version": "1.0.0",
            "version_meta": {
                "strategy": "content-hash",
                "timestamp": 1609459200
            },
            "maintainers": {
                "authors": ["Test Author"],
                "email": ["test@example.com"],
                "affiliations": ["Test University"]
            },
            "description": {
                "short": "Test legacy single-scheme torch",
                "long": ("Single-scheme test torch for backward "  # noqa
                         "compatibility"),
                "taxa": ["Escherichia coli"]
            },
            "manifest": {
                "profiles": "legacy.profiles.tsv",
                "workflow": "legacy.wdl",
                "buildfile": "legacy.build.wdl",
                "resources": []
            }
        }

        # Write metadata
        with open(torch_path / "metadata.toml", "w") as f:
            toml.dump(metadata, f)

        # Create profiles
        profiles = [
            ["ST", "dinB", "icdA"],
            ["1", "1", "1"],
            ["2", "2", "2"],
            ["3", "3", "1"]
        ]
        with open(torch_path / "legacy.profiles.tsv", "w") as f:
            writer = csv.writer(f, delimiter="\t")
            writer.writerows(profiles)

        # Create resources directory with alleles
        resources_path = torch_path / "_resources"
        resources_path.mkdir(parents=True, exist_ok=True)

        with open(resources_path / "dinB.fasta", "w") as f:
            f.write(">dinB_1\nACGT\n>dinB_2\nTGCA\n")

        with open(resources_path / "icdA.fasta", "w") as f:
            f.write(">icdA_1\nGATC\n")

        # Create workflow files (placeholders)
        with open(torch_path / "legacy.wdl", "w") as f:
            f.write("workflow legacy { }\n")

        with open(torch_path / "legacy.build.wdl", "w") as f:
            f.write("workflow build { }\n")

        yield torch_path
