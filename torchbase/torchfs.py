from os import environ
from dataclasses import dataclass, field
from pathlib import Path
from typing import Tuple, Dict, Optional, Union

from torchbase.torchbase import Profile, Schema

import ipyfs

import toml
import csv

node = environ.get("TORCHBASE_IPFS_NODE", "localhost")
port = environ.get("TORCHBASE_IPFS_PORT", 5001)

TORCHBASE_REGISTRY_HASH = ""  # IPFS hash for registry file


def handle_ipfs_errors(func):
    def ipfs_error_handler(*args, **kwargs):
        try:
            return func(*args, **kwargs)
        except Exception:
            raise

    return ipfs_error_handler


# @handle_ipfs_errors
def retrieve_manifest(
    cid=TORCHBASE_REGISTRY_HASH,
    node=node,
    port=port
):
    cat = ipyfs.Cat(host=node, port=port)
    # get the registry
    # load and return
    return toml.loads(cat(cid))


@handle_ipfs_errors
def download_torch(cid, node=node, port=port):
    pass


@handle_ipfs_errors
def register_torch(path, node=node, port=port):
    with open(path / "metadata.toml") as md_file:
        metadata = toml.load(md_file)
    metadata["manifest"]["resources"] = []
    for child in (path / "_resources").iterdir():
        if child.is_file():
            if not child.name.startswith("."):
                metadata["manifest"]["resources"].append(child.name)
    with open(path / "metadata.toml", "w") as md_file:
        toml.dump(metadata, md_file)


@handle_ipfs_errors
def exists(manifest, torch_entry):
    return False  # TODO


@dataclass
class Torch:
    """Represents a torch package with single and multi-scheme support.

    Attributes:
        path: Path to the torch directory
        schemes: Dict mapping scheme names to Schema objects
        profile: Schema object for single-scheme format (backward compat)
        scheme_references: Dict mapping scheme names to reference files
        references: Tuple of reference file paths (single-scheme format)
        workflow: Path to main workflow file
        buildfile: Path to build workflow file
    """
    path: Path
    profile: Optional[Union[Profile, Schema]] = None
    workflow: Optional[Path] = None
    buildfile: Optional[Path] = None
    references: Tuple[Path, ...] = field(default_factory=tuple)

    # Multi-scheme attributes
    schemes: Dict[str, Schema] = field(default_factory=dict)
    scheme_references: Dict[str, Tuple[Path, ...]] = field(
        default_factory=dict
    )

    @staticmethod
    def load(new_path):
        """Load torch from disk with single and multi-scheme support.

        Multi-scheme format: has schemes/ subdirectory
        Single-scheme format: legacy format with manifest.profiles

        Args:
            new_path: Path to the torch directory

        Returns:
            Torch object with either single or multi-scheme data

        Raises:
            ValueError: If validation fails
            FileNotFoundError: If required files are missing
        """
        path = Path(new_path)

        # Load and validate metadata
        metadata_path = path / "metadata.toml"
        if not metadata_path.exists():
            raise FileNotFoundError(f"metadata.toml not found at {path}")

        with open(metadata_path) as metadata_file:
            metadata = toml.load(metadata_file)

        # Sanity checks on path vs metadata
        *_, namespace_from_path, name_from_path, version_from_path = (
            path.parts
        )
        if version_from_path.endswith(".torch"):
            version_from_path = version_from_path[:-6]

        if metadata.get("namespace") != namespace_from_path:
            raise ValueError(
                f"Namespace mismatch: {metadata.get('namespace')} from "
                f"metadata vs {namespace_from_path} from path"
            )
        if metadata.get("name") != name_from_path:
            raise ValueError(
                f"Name mismatch: {metadata.get('name')} from metadata vs "
                f"{name_from_path} from path"
            )
        if str(metadata.get("version")) != version_from_path:
            raise ValueError(
                f"Version mismatch: {metadata.get('version')} from metadata "
                f"vs {version_from_path} from path"
            )

        # Detect if multi-scheme or single-scheme format
        schemes_dir = path / "schemes"
        is_multi_scheme = (
            schemes_dir.exists() and schemes_dir.is_dir()
        )

        if is_multi_scheme:
            return Torch._load_multi_scheme(path, metadata)
        else:
            return Torch._load_single_scheme(path, metadata)

    @staticmethod
    def _load_multi_scheme(path: Path, metadata: dict) -> "Torch":
        """Load multi-scheme torch format.

        Args:
            path: Path to torch directory
            metadata: Loaded metadata

        Returns:
            Torch with schemes attribute populated

        Raises:
            ValueError: If validation fails
            FileNotFoundError: If required files are missing
        """
        schemes_dir = path / "schemes"
        schemes = {}
        scheme_references = {}

        # Get declared schemes from metadata
        declared_schemes = set(metadata.get("schemes", {}).keys())

        # Discover actual schemes from filesystem
        discovered_schemes = set()
        for scheme_path in sorted(schemes_dir.iterdir()):
            if not scheme_path.is_dir() or scheme_path.name.startswith("."):
                continue

            scheme_name = scheme_path.name
            discovered_schemes.add(scheme_name)

            # Load scheme profiles
            profiles_file = scheme_path / "profiles.tsv"
            if not profiles_file.exists():
                raise FileNotFoundError(
                    f"profiles.tsv not found in scheme {scheme_name} at "
                    f"{profiles_file}"
                )

            with open(profiles_file) as f:
                schema = Profile.parse(
                    scheme_name,
                    csv.reader(f, delimiter="\t")
                )
                schemes[scheme_name] = schema

            # Load scheme alleles
            alleles_dir = scheme_path / "alleles"
            if not alleles_dir.exists():
                raise FileNotFoundError(
                    f"alleles directory not found in scheme {scheme_name} "
                    f"at {alleles_dir}"
                )

            allele_files = tuple(
                f
                for f in sorted(alleles_dir.iterdir())
                if f.is_file() and not f.name.startswith(".")
            )
            scheme_references[scheme_name] = allele_files

        # Validate declared schemes match discovered schemes
        if declared_schemes and discovered_schemes != declared_schemes:
            raise ValueError(
                f"Scheme mismatch: metadata declares "
                f"{sorted(declared_schemes)} but found "
                f"{sorted(discovered_schemes)}"
            )

        return Torch(
            path=path,
            schemes=schemes,
            scheme_references=scheme_references,
        )

    @staticmethod
    def _load_single_scheme(path: Path, metadata: dict) -> "Torch":
        """Load legacy single-scheme torch format.

        Args:
            path: Path to torch directory
            metadata: Loaded metadata

        Returns:
            Torch with profile and references for backward compatibility

        Raises:
            FileNotFoundError: If required files are missing
        """
        manifest = metadata.get("manifest", {})

        # Load profiles
        profiles_file = manifest.get("profiles")
        if not profiles_file:
            raise ValueError(
                "manifest.profiles not specified in single-scheme torch"
            )

        profiles_path = path / profiles_file
        if not profiles_path.exists():
            raise FileNotFoundError(
                f"Profiles file not found at {profiles_path}"
            )

        with open(profiles_path) as profile_file:
            profile = Profile.parse(
                f"{metadata['name']}_{metadata['version']}",
                csv.reader(profile_file, delimiter="\t"),
            )

        # Load references
        resources = path / "_resources"
        references = tuple()
        if resources.exists():
            references = tuple(
                f
                for f in sorted(resources.iterdir())
                if f.is_file() and not f.name.startswith(".")
            )

        # Get workflows if specified
        workflow = None
        if "workflow" in manifest:
            workflow = path / manifest["workflow"]

        buildfile = None
        if "buildfile" in manifest:
            buildfile = path / manifest["buildfile"]

        return Torch(
            path=path,
            profile=profile,
            references=references,
            workflow=workflow,
            buildfile=buildfile,
        )
