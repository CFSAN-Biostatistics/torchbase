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
        """Concatenate allele files with scheme-prefixed headers.

        For multi-scheme torches: concatenates all alleles from all schemes
        with headers prefixed by scheme name (e.g., ">ecoli_dinB_1").

        For single-scheme torches: returns unified alleles (either unprefixed
        or torch-name-prefixed for backward compatibility).

        Returns:
            Path to concatenated FASTA file

        Raises:
            ValueError: If a scheme has no allele files
            RuntimeError: If concatenation fails
        """
        output = tempfile.NamedTemporaryFile(
            mode='w',
            suffix='.fasta',
            delete=False,
            dir=self.path
        )

        try:
            if self.schemes:
                # Multi-scheme: concatenate with scheme prefixes
                return self._concatenate_multi_scheme_alleles(output)
            else:
                # Single-scheme: concatenate without scheme prefix
                return self._concatenate_single_scheme_alleles(output)
        except Exception:
            # Clean up temp file on error
            output.close()
            Path(output.name).unlink(missing_ok=True)
            raise

    def _concatenate_multi_scheme_alleles(self, output) -> Path:
        """Concatenate multi-scheme alleles with scheme prefixes."""
        for scheme_name in sorted(self.scheme_references.keys()):
            allele_files = self.scheme_references[scheme_name]

            if not allele_files:
                raise ValueError(
                    f"Scheme {scheme_name} has no allele files"
                )

            # Concatenate alleles from this scheme
            for allele_file in sorted(allele_files):
                with open(allele_file) as f:
                    for line in f:
                        line = line.rstrip('\n')
                        if line.startswith('>'):
                            # Prefix header with scheme name
                            header = line[1:]  # Remove '>'
                            prefixed_header = f">{scheme_name}_{header}"
                            output.write(prefixed_header + '\n')
                        else:
                            # Write sequence as-is
                            output.write(line + '\n')

        output.close()
        return Path(output.name)

    def _concatenate_single_scheme_alleles(self, output) -> Path:
        """Concatenate single-scheme alleles without scheme prefix."""
        for ref_file in sorted(self.references):
            with open(ref_file) as f:
                for line in f:
                    output.write(line)

        output.close()
        return Path(output.name)

    def transform_profiles(self) -> Path:
        """Transform profiles with scheme column and prefixed locus names.

        For multi-scheme torches: creates unified TSV with:
        - 'scheme' column identifying the source scheme
        - Locus columns prefixed with scheme name (e.g., "ecoli_dinB")
        - Cross-scheme loci filled with empty strings

        For single-scheme torches: returns profiles with optional scheme column.

        Returns:
            Path to transformed TSV file

        Raises:
            RuntimeError: If transformation fails
        """
        output = tempfile.NamedTemporaryFile(
            mode='w',
            suffix='.tsv',
            delete=False,
            dir=self.path
        )

        try:
            if self.schemes:
                # Multi-scheme: transform with scheme prefixes
                return self._transform_multi_scheme_profiles(output)
            else:
                # Single-scheme: transform without scheme prefix
                return self._transform_single_scheme_profiles(output)
        except Exception:
            # Clean up temp file on error
            output.close()
            Path(output.name).unlink(missing_ok=True)
            raise

    def _transform_multi_scheme_profiles(self, output) -> Path:
        """Transform multi-scheme profiles with scheme prefixes."""
        import csv

        # Collect all loci per scheme
        scheme_loci = {}
        for scheme_name, schema in self.schemes.items():
            if schema.profiles:
                # Profile.header contains all loci (not including ST)
                loci = list(schema.profiles[0].header)
                scheme_loci[scheme_name] = sorted(loci)

        # Build header: ST, all prefixed loci per scheme (alphabetically),
        # scheme column
        header = ['ST']
        for scheme_name in sorted(scheme_loci.keys()):
            for locus in scheme_loci[scheme_name]:
                header.append(f"{scheme_name}_{locus}")
        header.append('scheme')

        # Write header
        writer = csv.writer(output, delimiter='\t')
        writer.writerow(header)

        # Write profiles from each scheme
        for scheme_name, schema in sorted(self.schemes.items()):
            for profile in schema.profiles:
                # ST is the profile ID (second parameter in Profile.__init__)
                st_value = profile.profile
                row = [str(st_value)]

                # Add locus values for this scheme
                for other_scheme in sorted(scheme_loci.keys()):
                    for locus in scheme_loci[other_scheme]:
                        if other_scheme == scheme_name:
                            # Get value from this profile
                            val = profile.get(locus, '')
                            row.append(str(val) if val is not None else '')
                        else:
                            # Empty for other schemes
                            row.append('')

                # Add scheme name
                row.append(scheme_name)
                writer.writerow(row)

        output.close()
        return Path(output.name)

    def _transform_single_scheme_profiles(self, output) -> Path:
        """Transform single-scheme profiles without scheme prefix."""
        import csv

        # Get profile header and values
        if isinstance(self.profile, Schema):
            profiles = self.profile.profiles
        else:
            profiles = [self.profile]

        if not profiles:
            output.close()
            return Path(output.name)

        # Get loci from first profile (Profile.header contains all loci)
        first_profile = profiles[0]
        loci = sorted(list(first_profile.header))

        # Write header
        header = ['ST'] + loci
        writer = csv.writer(output, delimiter='\t')
        writer.writerow(header)

        # Write profiles
        for profile in profiles:
            # ST is the profile ID
            st_value = profile.profile
            row = [str(st_value)]
            for locus in loci:
                val = profile.get(locus, '')
                row.append(str(val) if val is not None else '')
            writer.writerow(row)

        output.close()
        return Path(output.name)

    def get_unified_files(self) -> Tuple[Path, Path]:
        """Get unified alleles and profiles for multi-scheme torch.

        Convenience method that returns both concatenated alleles and
        transformed profiles in a single call.

        Returns:
            Tuple of (alleles_path, profiles_path)
        """
        alleles = self.concatenate_alleles()
        profiles = self.transform_profiles()
        return alleles, profiles
