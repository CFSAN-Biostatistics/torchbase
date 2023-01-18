from dataclasses import dataclass
from enum import Enum
from typing import Iterator, Tuple, Set, List
from pathlib import Path

from itertools import zip_longest

import tomllib
import csv

TORCHBASE_REGISTRY_HASH = "" # IPFS hash for registry file DO NOT CHANGE

class Special(Enum):
    IGNORE = "?"
    EXCLUDE = "X"

    def __repr__(self):
        return self.value

    def __str__(self):
        return repr(self)

class Schema:
    def __init__(self, name, profiles=[]):
        self.name = name
        self.profiles = tuple(profiles)

    def __getitem__(self, key):
        for profile in self.profiles:
            if profile == key:
                return profile.profile
            raise KeyError(f"{key} not found in schema.")

    def __contains__(self, item):
        for profile in self.profiles:
            if profile == item:
                return True
            return False

    def __repr__(self):
        return f"<{self.name}: {len(profiles)} profiles>"



class Profile:
    def __init__(self, schema=None, profile=None, **kwargs):
        self.schema = schema
        self.profile = profile
        self.excluded = set()
        self.header = tuple(kwargs.keys())
        self.values = tuple(kwargs.values())
        for locus, value in kwargs.items():
            if value == Special.EXCLUDE:
                self.excluded.add(locus)
            setattr(self, locus, value)

    def __contains__(self, key):
        return key in self.header

    def __getitem__(self, key):
        return getattr(self, key)

    def __iter__(self):
        return iter(self.values)

    def __len__(self):
        return len(self.values)

    def get(self, key, default=None):
        return vars(self).get(key, default)

    def zip(self, other):
        "Special function to line up values and provide default values if missing"
        h_iter = iter(self.header)
        s_iter = iter(self)
        o_iter = iter(other)
        for h, s, o in zip(h_iter, s_iter, o_iter):
            yield h, s, o
        # one of these iterators is now exhausted
        while True:
            s_f = o_f = False # flags
            try:
                s = next(s_iter)
                h = next(h_iter)
            except StopIteration:
                s = Special.IGNORE
                h = '?'
                s_f = True
            try:
                o = next(o_iter)
            except StopIteration:
                o = Special.EXCLUDE
                o_f = True
            if s_f and o_f:
                return
            yield h, s, o

    def __eq__(self, other) -> bool:
        def compare(h, s, o):
            return o==s or str(o) == str(s) or o==f"{h}_{s}" or s is Special.IGNORE
        def correct_length(o):
            return len([v for v in self.values if v is not Special.EXCLUDE]) == len(other)
        if isinstance(other, (tuple, list)):
            return all(compare(h, s, o) for h, s, o in self.zip(other)) and correct_length(other)
        if isinstance(other, (dict, Profile)):
            return all(compare(h, self[h], other.get(h, Special.EXCLUDE)) for h in self.header) and correct_length(other)
        return False

    def __hash__(self):
        return hash(self.headers) + hash(self.values)

    def __repr__(self) -> str:
        vals = ', '.join(f"{h}={v}" for h,v in zip(self.header, self.values))
        return f"""Profile({self.schema}, {self.profile}, {vals})"""

    def __str__(self) -> str:
        return f"<{self.profile} ({len(self.values)} loci)>"

    @classmethod
    def parse(cls, schema_name: str, rows: Iterator[Tuple]) -> Schema:
        rows = iter(rows)
        _, *header = next(rows)
        profiles = []
        for profile, *row in rows:
            assert len(header) == len(row)
            profiles.append(cls(schema_name, profile, **dict(zip(header, row))))
        return Schema(schema_name, profiles)
            
@dataclass
class Torch:

    path: Path
    profile: Profile
    references: Tuple[Path]
    
    @staticmethod
    def load(new_path):
        path = Path(new_path)
        with open(path / "metadata.toml") as metadata_file:
            metadata = tomllib.load(metadata_file)
            # run some sanity checks
            *_, namespace_from_path, name_from_path, version_from_path = path.parts()
            if not metadata.namespace == namespace_from_path:
                raise ValueError(f"Failed sanity check, namespace {metadata.namespace} from metadata didn't match {namespace_from_path} from path")
            if not metadata.name == name_from_path:
                raise ValueError(f"Failed sanity check, name {metadata.name} from metadata didn't match {name_from_path} from path")
            if not str(metadata.version) == version_from_path:
                raise ValueError(f"Failed sanity check, version {metadata.version} from metadata didn't match {version_from_path} from path")
        with open(path / metadata.manifest.profiles) as profile_file:
            profile = Profile(
                schema_name = f"{metadata.name}_{metadata.version}",
                rows=csv.reader(profile_file, dialect='excel', delimiter='\t')
            )
        resources = path / "_resources"
        return Torch(
            path=path,
            profile=profile,
            references=tuple([file for file in resources.iterdir() if file.is_file()])
        )
        

    def dump(self, new_path):
        pass