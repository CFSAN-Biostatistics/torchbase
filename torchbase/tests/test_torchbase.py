#!/usr/bin/env python

"""Tests for `torchbase` package."""

import pytest

from .bigsdb_fixture import bigsdb, profile

import torchbase


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
