
from collections import defaultdict
from enum import Enum
from typing import Iterator, Tuple, Set, List
from pathlib import Path

from itertools import zip_longest
import json
import re
from datetime import datetime


import csv

class Special(Enum):
    IGNORE = "?"
    EXCLUDE = "X"

    def __repr__(self):
        return self.value

    def __str__(self):
        return repr(self)

class Schema:
    def __init__(self, name, profiles=[], version=""):
        self.name = name
        self.version = version
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

    def to_json(self):
        return json.dumps(
            dict(
                schema = self.schema.name,
                profile = self.profile,
                **{header:value for header, value in zip(self.header, self.values)}
            )
        )

    @classmethod
    def parse(cls, schema_name: str, rows: Iterator[Tuple]) -> Schema:
        rows = iter(rows)
        _, *header = next(rows)
        profiles = []
        for profile, *row in rows:
            assert len(header) == len(row)
            profiles.append(cls(schema_name, profile, **dict(zip(header, row))))
        return Schema(schema_name, profiles=profiles)


class Version:
    """Represents a versioned torch with support for multiple versioning strategies.

    Strategies:
    - snapshot: ISO date string (YYYY-MM-DD)
    - semver: Semantic versioning (X.Y.Z)
    - content-hash: Hash-based versioning (any alphanumeric string)

    All versions are ordered by timestamp (Unix epoch) for cross-strategy comparison.
    """

    def __init__(self, version_str: str, strategy: str, timestamp: int):
        """Initialize a Version.

        Args:
            version_str: The version string (format depends on strategy)
            strategy: One of 'snapshot', 'semver', or 'content-hash'
            timestamp: Unix epoch timestamp for ordering
        """
        self.version_str = version_str
        self.strategy = strategy
        self.timestamp = timestamp

    @staticmethod
    def parse(version_str: str, metadata: dict) -> "Version":
        """Parse a version string with metadata context.

        Args:
            version_str: The version string to parse
            metadata: Dictionary containing version metadata with [version] section

        Returns:
            Version object

        Raises:
            ValueError: If metadata is invalid or version format is invalid
        """
        # Validate metadata structure
        if "version" not in metadata:
            raise ValueError("Missing version section in metadata")

        version_meta = metadata["version"]

        if "strategy" not in version_meta:
            raise ValueError("Missing strategy in version metadata")

        if "timestamp" not in version_meta:
            raise ValueError("Missing timestamp in version metadata")

        strategy = version_meta["strategy"]
        timestamp = version_meta["timestamp"]

        # Validate version string based on strategy
        if strategy == "snapshot":
            Version._validate_snapshot(version_str)
        elif strategy == "semver":
            Version._validate_semver(version_str)
        elif strategy == "content-hash":
            Version._validate_content_hash(version_str)
        else:
            raise ValueError(f"Unknown strategy: {strategy}")

        return Version(version_str, strategy, timestamp)

    @staticmethod
    def _validate_snapshot(version_str: str) -> None:
        """Validate ISO date format (YYYY-MM-DD).

        Raises:
            ValueError: If format is invalid
        """
        try:
            datetime.strptime(version_str, "%Y-%m-%d")
        except ValueError:
            raise ValueError(
                f"Invalid snapshot format: '{version_str}'. Expected ISO date (YYYY-MM-DD)"
            )

    @staticmethod
    def _validate_semver(version_str: str) -> None:
        """Validate semantic versioning format (X.Y.Z).

        Raises:
            ValueError: If format is invalid
        """
        # Match semver pattern: X.Y.Z or X.Y.Z-prerelease or X.Y.Z+build
        semver_pattern = r"^\d+\.\d+\.\d+(-[a-zA-Z0-9]+)*(\+[a-zA-Z0-9]+)*$"
        if not re.match(semver_pattern, version_str):
            raise ValueError(
                f"Invalid semver format: '{version_str}'. Expected X.Y.Z format"
            )

    @staticmethod
    def _validate_content_hash(version_str: str) -> None:
        """Validate content hash format (alphanumeric string).

        Currently accepts any non-empty alphanumeric string.

        Raises:
            ValueError: If format is invalid
        """
        if not version_str or not re.match(r"^[a-zA-Z0-9]+$", version_str):
            raise ValueError(
                f"Invalid content-hash format: '{version_str}'. "
                "Expected alphanumeric string"
            )

    @staticmethod
    def compare(v1: "Version", v2: "Version") -> int:
        """Compare two versions based on timestamps.

        Args:
            v1: First version
            v2: Second version

        Returns:
            -1 if v1 < v2, 0 if v1 == v2, 1 if v1 > v2
        """
        if v1.timestamp < v2.timestamp:
            return -1
        elif v1.timestamp > v2.timestamp:
            return 1
        else:
            return 0

    def __eq__(self, other: "Version") -> bool:
        """Check equality based on timestamp."""
        if not isinstance(other, Version):
            return False
        return self.timestamp == other.timestamp

    def __lt__(self, other: "Version") -> bool:
        """Check if this version is less than other (by timestamp)."""
        if not isinstance(other, Version):
            return NotImplemented
        return self.timestamp < other.timestamp

    def __le__(self, other: "Version") -> bool:
        """Check if version <= other (by timestamp)."""
        if not isinstance(other, Version):
            return NotImplemented
        return self.timestamp <= other.timestamp

    def __gt__(self, other: "Version") -> bool:
        """Check if this version is greater than other (by timestamp)."""
        if not isinstance(other, Version):
            return NotImplemented
        return self.timestamp > other.timestamp

    def __ge__(self, other: "Version") -> bool:
        """Check if version >= other (by timestamp)."""
        if not isinstance(other, Version):
            return NotImplemented
        return self.timestamp >= other.timestamp

    def __hash__(self) -> int:
        """Hash based on strategy and version_str.

        Note: timestamp is not included as it's mutable metadata.
        """
        return hash((self.strategy, self.version_str))

    def __repr__(self) -> str:
        """String representation of Version."""
        return (
            f"Version({self.version_str!r}, strategy={self.strategy!r}, "
            f"timestamp={self.timestamp})"
        )

    def __str__(self) -> str:
        """Human-readable string representation."""
        return f"{self.version_str} ({self.strategy})"
