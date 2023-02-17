
from json import loads
from os import environ
from dataclasses import dataclass
from pathlib import Path
from typing import Tuple

from torchbase.torchbase import Profile

import ipyfs

import toml
import csv

node = environ.get("TORCHBASE_IPFS_NODE", "localhost")
port = environ.get("TORCHBASE_IPFS_PORT", 5001)

TORCHBASE_REGISTRY_HASH = "" # IPFS hash for registry file DO NOT CHANGE

def handle_ipfs_errors(func):
    def ipfs_error_handler(*args, **kwargs):
        try:
            return func(*args, **kwargs)
        except:
            raise

# @handle_ipfs_errors
def retrieve_manifest(cid=TORCHBASE_REGISTRY_HASH, node=node, port=port):
    cat = ipyfs.Cat(host=node, port=port)
    # get the registry
    # load and return
    return toml.loads(cat(hash))

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



@dataclass
class Torch:

    path: Path
    profile: Profile
    workflow: Path
    buildfile: Path
    references: Tuple[Path]
    
    @staticmethod
    def load(new_path):
        path = Path(new_path)
        with open(path / "metadata.toml") as metadata_file:
            metadata = toml.load(metadata_file)
            # run some sanity checks
            *_, namespace_from_path, name_from_path, version_from_path = path.parts()
            if not metadata["namespace"] == namespace_from_path:
                raise ValueError(f"Failed sanity check, namespace {metadata['namespace']} from metadata didn't match {namespace_from_path} from path")
            if not metadata["name"] == name_from_path:
                raise ValueError(f"Failed sanity check, name {metadata['name']} from metadata didn't match {name_from_path} from path")
            if not str(metadata["version"]) == version_from_path:
                raise ValueError(f"Failed sanity check, version {metadata['version']} from metadata didn't match {version_from_path} from path")
        with open(path / metadata["manifest"]["profiles"]) as profile_file:
            profile = Profile(
                schema_name = f"{metadata['name']}_{metadata['version']}",
                rows=csv.reader(profile_file, dialect='excel', delimiter='\t')
            )
        resources = path / "_resources"
        return Torch(
            path=path,
            profile=profile,
            references=tuple([file for file in resources.iterdir() if file.is_file() and not file.name.startswith(".")]),
            workflow=path / metadata['manifest']['workflow'],
            buildfile=path / metadata['manifest']['buildfile']
        )