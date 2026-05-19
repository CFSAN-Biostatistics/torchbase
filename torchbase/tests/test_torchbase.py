#!/usr/bin/env python

"""Tests for `torchbase` package."""

import pytest

from .bigsdb_fixture import bigsdb, profile

import torchbase, torchbase.torchfs


# def test_torchbase(response):
#     """Sample pytest test function with the pytest fixture as an argument."""
#     # from bs4 import BeautifulSoup
#     # assert 'GitHub' in BeautifulSoup(response.content).title.string

@pytest.fixture
def ut():
    p = {"a":1, "b":2, "c":3, "d":torchbase.Special.IGNORE, "e":torchbase.Special.EXCLUDE}
    return torchbase.Profile("test_schema", "TEST1", **p)


class TestProfileAPI:

    def test_contains(self, ut):
        assert "a" in ut

    def test_not_contains(self, ut):
        assert "f" not in ut

    def test_getitem(self, ut):
        assert ut["a"] == 1

    def test_attr(self, ut):
        assert ut.a == 1
 
class TestEquality:

    def test_tuple(self, ut):
        assert ut == (1, 2, 3, 5)

    def test_pmlst_style_tuple(self, ut):
        assert ut == ("a_1", "b_2", "c_3", "d_5")

    def test_nonmatching_tuple(self, ut):
        assert ut != (4, 6, 8, 10)
        assert ut != ("a_4", "b_6", "c_8", "d_10")

    def test_required_wildcard_tuple(self, ut):
        assert ut != (1,2,3)

    def test_excluded_locus_tuple(self, ut):
        assert ut != (1, 2, 3, 5, 7)

    def test_excluded_locus(self, ut):
        c = torchbase.Profile(a=1, b=2, c=3, d=10, e=20)
        assert ut != c
        

    def test_dict(self, ut):
        assert ut == {"a":1, "b":2, "c":3, "d":5}

    def test_obj(self, ut):
        assert ut == torchbase.Profile(a=1, b=2, c=3, d=5)

    def test_non_matching_obj(self, ut):
        assert ut != torchbase.Profile(a=2, b=4, c=6, d=10)

    def test_required_wildcard(self, ut):
        assert ut != torchbase.Profile(a=1, b=2, c=3)

    # @pytest.mark.skip
    # def test_raises_value_error(self, ut):
    #     with pytest.raises(ValueError):
    #         assert ut == bool()



class TestProfileParser:

    def test_parse(self, bigsdb, profile):
        schema = torchbase.Profile.parse("bigsdb", bigsdb)
        assert profile in schema

    def test_determination(self, bigsdb, profile):
        schema = torchbase.Profile.parse("bigsdb", bigsdb)
        assert schema[profile] == "1"


head = """Hello and Welcome to IPFS!"""

hash = "/ipfs/QmQPeNsJPyVWPFDVHb77w8G42Fvo15z4bG2X8D2GhfbSXc/readme"


class TestTorchFS:

    def test_ipyfs(self):
        cat = torchbase.torchfs.ipyfs.Cat()
        assert cat(hash)['result'][0:len(head)] == head

    def test_handle_ipfs_errors(self):
        pass

    def test_retrieve_manifest(self):
        pass

    def test_download_torch(self):
        pass

    def test_register_torch(self):
        pass

    def test_exists(self):
        pass

class TestTorchClass:

    def test_load(self):
        pass


class TestVersionParsing:
    """Tests for Version.parse() with different strategies."""

    def test_parse_snapshot_version(self):
        """Parse snapshot (ISO date) version."""
        metadata = {
            "version": {
                "strategy": "snapshot",
                "timestamp": 1609459200  # 2021-01-01T00:00:00Z
            }
        }
        version = torchbase.Version.parse("2021-01-01", metadata)
        assert version.strategy == "snapshot"
        assert version.version_str == "2021-01-01"
        assert version.timestamp == 1609459200

    def test_parse_semantic_version(self):
        """Parse semantic (semver) version."""
        metadata = {
            "version": {
                "strategy": "semver",
                "timestamp": 1609459200  # 2021-01-01T00:00:00Z
            }
        }
        version = torchbase.Version.parse("1.2.3", metadata)
        assert version.strategy == "semver"
        assert version.version_str == "1.2.3"
        assert version.timestamp == 1609459200

    def test_parse_content_hash_version(self):
        """Parse content-addressed (hash) version."""
        metadata = {
            "version": {
                "strategy": "content-hash",
                "timestamp": 1609459200  # 2021-01-01T00:00:00Z
            }
        }
        version = torchbase.Version.parse("abc123def456", metadata)
        assert version.strategy == "content-hash"
        assert version.version_str == "abc123def456"
        assert version.timestamp == 1609459200

    def test_parse_invalid_snapshot_format_raises_error(self):
        """Invalid ISO date format raises ValueError."""
        metadata = {
            "version": {
                "strategy": "snapshot",
                "timestamp": 1609459200
            }
        }
        with pytest.raises(ValueError, match="Invalid snapshot format"):
            torchbase.Version.parse("not-a-date", metadata)

    def test_parse_invalid_semver_format_raises_error(self):
        """Invalid semver format raises ValueError."""
        metadata = {
            "version": {
                "strategy": "semver",
                "timestamp": 1609459200
            }
        }
        with pytest.raises(ValueError, match="Invalid semver format"):
            torchbase.Version.parse("1.2", metadata)

    def test_parse_missing_timestamp_raises_error(self):
        """Missing timestamp in metadata raises ValueError."""
        metadata = {
            "version": {
                "strategy": "snapshot"
            }
        }
        with pytest.raises(ValueError, match="Missing timestamp"):
            torchbase.Version.parse("2021-01-01", metadata)

    def test_parse_missing_strategy_raises_error(self):
        """Missing strategy in metadata raises ValueError."""
        metadata = {
            "version": {
                "timestamp": 1609459200
            }
        }
        with pytest.raises(ValueError, match="Missing strategy"):
            torchbase.Version.parse("1.2.3", metadata)

    def test_parse_missing_version_section_raises_error(self):
        """Missing version section in metadata raises ValueError."""
        metadata = {}
        with pytest.raises(ValueError, match="Missing version"):
            torchbase.Version.parse("1.2.3", metadata)


class TestVersionComparison:
    """Tests for Version.compare() and version ordering."""

    def test_compare_same_versions_returns_zero(self):
        """Comparing identical versions returns 0."""
        metadata = {
            "version": {
                "strategy": "snapshot",
                "timestamp": 1609459200
            }
        }
        v1 = torchbase.Version.parse("2021-01-01", metadata)
        v2 = torchbase.Version.parse("2021-01-01", metadata)
        assert torchbase.Version.compare(v1, v2) == 0

    def test_compare_earlier_timestamp_returns_negative(self):
        """Version with earlier timestamp is less than later timestamp."""
        metadata1 = {
            "version": {
                "strategy": "snapshot",
                "timestamp": 1609459200  # 2021-01-01
            }
        }
        metadata2 = {
            "version": {
                "strategy": "semver",
                "timestamp": 1640995200  # 2022-01-01
            }
        }
        v1 = torchbase.Version.parse("2021-01-01", metadata1)
        v2 = torchbase.Version.parse("1.0.0", metadata2)
        assert torchbase.Version.compare(v1, v2) == -1

    def test_compare_later_timestamp_returns_positive(self):
        """Version with later timestamp is greater than earlier timestamp."""
        metadata1 = {
            "version": {
                "strategy": "semver",
                "timestamp": 1640995200  # 2022-01-01
            }
        }
        metadata2 = {
            "version": {
                "strategy": "snapshot",
                "timestamp": 1609459200  # 2021-01-01
            }
        }
        v1 = torchbase.Version.parse("1.0.0", metadata1)
        v2 = torchbase.Version.parse("2021-01-01", metadata2)
        assert torchbase.Version.compare(v1, v2) == 1

    def test_compare_cross_strategy_same_timestamp_returns_zero(self):
        """Versions from different strategies with same timestamp are equal."""
        metadata = {
            "version": {
                "strategy": "snapshot",
                "timestamp": 1609459200
            }
        }
        v1 = torchbase.Version.parse("2021-01-01", metadata)
        metadata["version"]["strategy"] = "semver"
        v2 = torchbase.Version.parse("1.0.0", metadata)
        assert torchbase.Version.compare(v1, v2) == 0

    def test_sorting_mixed_strategy_versions(self):
        """Can sort versions with mixed strategies."""
        v1_metadata = {
            "version": {"strategy": "semver", "timestamp": 1640995200}
        }
        v2_metadata = {
            "version": {"strategy": "snapshot", "timestamp": 1609459200}
        }
        v3_metadata = {
            "version": {"strategy": "content-hash", "timestamp": 1577836800}
        }
        versions = [
            torchbase.Version.parse("1.0.0", v1_metadata),
            torchbase.Version.parse("2021-01-01", v2_metadata),
            torchbase.Version.parse("abc123", v3_metadata),
        ]
        sorted_versions = sorted(versions, key=lambda v: v.timestamp)
        assert sorted_versions[0].timestamp == 1577836800
        assert sorted_versions[1].timestamp == 1609459200
        assert sorted_versions[2].timestamp == 1640995200

    def test_same_timestamp_different_strings(self):
        """Same timestamp, different strings order correctly."""
        metadata = {
            "version": {
                "strategy": "snapshot",
                "timestamp": 1609459200
            }
        }
        v1 = torchbase.Version.parse("2021-01-01", metadata)
        metadata["version"]["strategy"] = "semver"
        v2 = torchbase.Version.parse("1.0.0", metadata)
        # They should compare as equal (same timestamp)
        assert torchbase.Version.compare(v1, v2) == 0


class TestVersionEquality:
    """Tests for Version equality and comparison operators."""

    def test_version_equality_operator(self):
        """Versions with same timestamp are equal."""
        metadata = {
            "version": {
                "strategy": "snapshot",
                "timestamp": 1609459200
            }
        }
        v1 = torchbase.Version.parse("2021-01-01", metadata)
        v2 = torchbase.Version.parse("2021-01-01", metadata)
        assert v1 == v2

    def test_version_less_than_operator(self):
        """Version with earlier timestamp is less than later."""
        v1 = torchbase.Version.parse("2021-01-01", {
            "version": {"strategy": "snapshot", "timestamp": 1609459200}
        })
        v2 = torchbase.Version.parse("1.0.0", {
            "version": {"strategy": "semver", "timestamp": 1640995200}
        })
        assert v1 < v2

    def test_version_greater_than_operator(self):
        """Version with later timestamp is greater than earlier."""
        v1 = torchbase.Version.parse("1.0.0", {
            "version": {"strategy": "semver", "timestamp": 1640995200}
        })
        v2 = torchbase.Version.parse("2021-01-01", {
            "version": {"strategy": "snapshot", "timestamp": 1609459200}
        })
        assert v1 > v2
