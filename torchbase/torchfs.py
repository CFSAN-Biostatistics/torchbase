from os import environ
from dataclasses import dataclass, field
from pathlib import Path
from typing import Tuple, Dict, Optional, Union
import tempfile

from torchbase.torchbase import Profile, Schema

try:
    import ipyfs
except ImportError:
    ipyfs = None

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

        # Workflow discovery: convention-based (main.wdl) takes precedence
        workflow = None
        main_wdl = path / "main.wdl"
        if main_wdl.exists():
            workflow = main_wdl
        elif "workflow" in manifest:
            # Fallback to manifest-specified workflow
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

    def concatenate_alleles(self) -> Path:
        """Concatenate allele files from all schemes into single FASTA.

        For multi-scheme torches: concatenates all allele files with
        scheme-prefixed headers (e.g., ">ecoli_adk_1" for E. coli adk_1).

        For single-scheme torches: concatenates allele files without
        scheme prefix (backward compatible).

        Returns:
            Path to temporary concatenated FASTA file

        Raises:
            ValueError: If a scheme has no allele files
            RuntimeError: If concatenation fails
        """
        # Create temporary file for concatenated alleles
        temp_file = tempfile.NamedTemporaryFile(
            mode='w', suffix='.fasta', delete=False, dir=self.path
        )
        temp_path = Path(temp_file.name)

        try:
            if self.schemes:
                # Multi-scheme case
                self._concatenate_multi_scheme_alleles(temp_path)
            else:
                # Single-scheme case
                self._concatenate_single_scheme_alleles(temp_path)
        finally:
            temp_file.close()

        return temp_path

    def _concatenate_multi_scheme_alleles(self, output_path: Path) -> None:
        """Concatenate alleles from all schemes with scheme prefixes.

        Args:
            output_path: Path to write concatenated FASTA

        Raises:
            ValueError: If a scheme has no allele files
        """
        with open(output_path, 'w') as out_f:
            # Process schemes in sorted order for consistency
            for scheme_name in sorted(self.scheme_references.keys()):
                allele_files = self.scheme_references[scheme_name]

                if not allele_files:
                    raise ValueError(
                        f"Scheme '{scheme_name}' has no allele files"
                    )

                # Process each allele file in the scheme
                for allele_file in sorted(allele_files):
                    with open(allele_file, 'r') as in_f:
                        for line in in_f:
                            line = line.rstrip('\n')
                            if line.startswith('>'):
                                # Prefix header with scheme name
                                # e.g., ">adk_1" becomes ">ecoli_adk_1"
                                header = line[1:]  # Remove '>'
                                prefixed_header = f">{scheme_name}_{header}"
                                out_f.write(prefixed_header + '\n')
                            else:
                                # Write sequence as-is
                                out_f.write(line + '\n')

    def _concatenate_single_scheme_alleles(self, output_path: Path) -> None:
        """Concatenate alleles from single-scheme torch (backward compatible).

        Args:
            output_path: Path to write concatenated FASTA
        """
        with open(output_path, 'w') as out_f:
            for ref_file in sorted(self.references):
                with open(ref_file, 'r') as in_f:
                    out_f.write(in_f.read())

    def transform_profiles(self) -> Path:
        """Transform profile tables with scheme column and prefixed locus names.

        For multi-scheme torches: adds 'scheme' column and prefixes locus names
        (e.g., "adk" becomes "ecoli_adk").

        For single-scheme torches: returns profiles as-is (backward compatible).

        Returns:
            Path to temporary transformed profile TSV file

        Raises:
            RuntimeError: If transformation fails
        """
        # Create temporary file for transformed profiles
        temp_file = tempfile.NamedTemporaryFile(
            mode='w', suffix='.tsv', delete=False, dir=self.path
        )
        temp_path = Path(temp_file.name)
        temp_file.close()

        try:
            if self.schemes:
                # Multi-scheme case
                self._transform_multi_scheme_profiles(temp_path)
            else:
                # Single-scheme case
                self._transform_single_scheme_profiles(temp_path)
        except Exception:
            # Clean up on error
            if temp_path.exists():
                temp_path.unlink()
            raise

        return temp_path

    def _transform_multi_scheme_profiles(self, output_path: Path) -> None:
        """Transform profiles from all schemes with scheme column.

        Args:
            output_path: Path to write transformed TSV
        """
        # Collect all loci from all schemes
        all_loci = {}  # scheme_name -> set of loci
        for scheme_name, schema in sorted(self.schemes.items()):
            if schema.profiles:
                first_profile = schema.profiles[0]
                # Get loci (all headers except 'ST')
                loci = [h for h in first_profile.header if h != 'ST']
                all_loci[scheme_name] = loci

        # Build header: ST, scheme, then all prefixed loci
        header = ['ST', 'scheme']
        for scheme_name in sorted(all_loci.keys()):
            for locus in all_loci[scheme_name]:
                header.append(f"{scheme_name}_{locus}")

        # Write header
        with open(output_path, 'w') as out_f:
            out_f.write('\t'.join(header) + '\n')

            # Write profiles from each scheme
            for scheme_name in sorted(self.schemes.keys()):
                schema = self.schemes[scheme_name]
                for profile in schema.profiles:
                    row = [profile.profile, scheme_name]  # ST and scheme

                    # Add allele values for all loci
                    for other_scheme in sorted(all_loci.keys()):
                        for locus in all_loci[other_scheme]:
                            if other_scheme == scheme_name:
                                # This scheme has this locus
                                value = getattr(profile, locus, '')
                            else:
                                # Cross-scheme locus: leave empty
                                value = ''
                            row.append(value)

                    out_f.write('\t'.join(row) + '\n')

    def _transform_single_scheme_profiles(self, output_path: Path) -> None:
        """Transform profiles from single-scheme torch (backward compatible).

        Args:
            output_path: Path to write transformed TSV
        """
        # For single-scheme, copy profiles as-is
        # Profile object has attributes matching original table structure
        if not self.profile or not self.profile.profiles:
            return

        first_profile = self.profile.profiles[0]
        header = first_profile.header

        with open(output_path, 'w') as out_f:
            # Write header
            out_f.write('\t'.join(header) + '\n')

            # Write each profile
            for profile in self.profile.profiles:
                row = [str(getattr(profile, h, '')) for h in header]
                out_f.write('\t'.join(row) + '\n')

    def get_unified_files(self) -> Tuple[Path, Path]:
        """Get both unified alleles and profiles files.

        Convenience method that returns both concatenated/transformed files.

        Returns:
            Tuple of (alleles_path, profiles_path)
        """
        alleles = self.concatenate_alleles()
        profiles = self.transform_profiles()
        return alleles, profiles
