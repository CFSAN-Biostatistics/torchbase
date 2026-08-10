#!/usr/bin/env python

"""Tests for typing-report interpretation (merge, gate, exclusion metadata).

These pin the behaviour lifted out of the `merge_allele_calls`,
`check_confidence_for_alignment` and `add_exclusion_metadata` WDL heredocs.
"""

import pytest

from torchbase.reports import (
    DEFAULT_CONFIDENCE_THRESHOLD,
    add_exclusion_metadata,
    alignment_call_preferred,
    locus_needs_alignment,
    merge_calls,
    needs_alignment,
)


@pytest.fixture
def minhash_calls():
    """MinHash calls for two loci; only abcZ has an alignment counterpart."""
    return {
        "abcZ": {"allele": "1", "similarity": 0.90, "source": "minhash"},
        "adk": {"allele": "2", "similarity": 0.99, "source": "minhash"},
    }


@pytest.fixture
def alignment_calls():
    """Alignment calls for abcZ (shared) and tonB (MinHash never saw it)."""
    return {
        "abcZ": {"allele": "7", "identity": 0.99, "source": "alignment"},
        "tonB": {"allele": "9", "identity": 0.99, "source": "alignment"},
    }


@pytest.fixture
def no_exclusions():
    """A filter_alleles exclusions document where nothing was excluded."""
    return {
        "excluded_alleles": [],
        "excluded_loci": [],
        "num_excluded_alleles": 0,
        "num_excluded_loci": 0,
    }


@pytest.fixture
def some_exclusions():
    return {
        "excluded_alleles": ["abcZ_3"],
        "excluded_loci": ["tonB"],
        "num_excluded_alleles": 1,
        "num_excluded_loci": 1,
    }


class TestDefaultThreshold:
    def test_matches_workflow_default(self):
        """Single source of truth for the balanced workflow's 0.85 default."""
        assert DEFAULT_CONFIDENCE_THRESHOLD == 0.85


class TestMergeCalls:
    def test_alignment_wins_for_shared_locus(self, minhash_calls, alignment_calls):
        merged = merge_calls(minhash_calls, alignment_calls, True, 0.85)
        assert merged["abcZ"]["source"] == "alignment"
        assert merged["abcZ"]["allele"] == "7"

    def test_minhash_kept_where_alignment_is_absent(
        self, minhash_calls, alignment_calls
    ):
        merged = merge_calls(minhash_calls, alignment_calls, True, 0.85)
        assert merged["adk"]["source"] == "minhash"

    def test_minhash_defines_the_locus_set(self, minhash_calls, alignment_calls):
        """A locus only alignment saw is dropped; a locus in neither is absent."""
        merged = merge_calls(minhash_calls, alignment_calls, True, 0.85)
        assert set(merged) == {"abcZ", "adk"}
        assert "tonB" not in merged
        assert "pgm" not in merged

    def test_weak_alignment_falls_back_to_minhash(self, minhash_calls):
        weak = {"abcZ": {"allele": "7", "identity": 0.84, "source": "alignment"}}
        merged = merge_calls(minhash_calls, weak, True, 0.85)
        assert merged["abcZ"]["source"] == "minhash"

    def test_alignment_at_threshold_wins(self, minhash_calls):
        """The identity comparison is inclusive: exactly at threshold is good."""
        exact = {"abcZ": {"allele": "7", "identity": 0.85, "source": "alignment"}}
        merged = merge_calls(minhash_calls, exact, True, 0.85)
        assert merged["abcZ"]["source"] == "alignment"

    def test_alignment_without_identity_loses(self, minhash_calls):
        merged = merge_calls(
            minhash_calls, {"abcZ": {"allele": "7", "source": "alignment"}}, True, 0.85
        )
        assert merged["abcZ"]["source"] == "minhash"

    def test_alignment_without_identity_wins_at_zero_threshold(self, minhash_calls):
        """Missing identity is 0.0, which still clears a 0.0 threshold."""
        merged = merge_calls(
            minhash_calls, {"abcZ": {"allele": "7", "source": "alignment"}}, True, 0.0
        )
        assert merged["abcZ"]["source"] == "alignment"

    def test_alignment_ignored_when_not_requested(
        self, minhash_calls, alignment_calls
    ):
        assert merge_calls(minhash_calls, alignment_calls, False, 0.85) == minhash_calls

    def test_empty_minhash_yields_empty_merge(self, alignment_calls):
        assert merge_calls({}, alignment_calls, True, 0.85) == {}
        assert merge_calls({}, alignment_calls, False, 0.85) == {}

    def test_empty_alignment_yields_minhash(self, minhash_calls):
        assert merge_calls(minhash_calls, {}, True, 0.85) == minhash_calls

    def test_both_empty(self):
        assert merge_calls({}, {}, True, 0.85) == {}

    def test_inputs_are_not_mutated(self, minhash_calls, alignment_calls):
        merged = merge_calls(minhash_calls, alignment_calls, True, 0.85)
        merged["adk"] = {"allele": "clobbered"}
        assert minhash_calls["adk"]["allele"] == "2"
        assert set(alignment_calls) == {"abcZ", "tonB"}

    def test_threshold_defaults_to_module_default(self, minhash_calls):
        borderline = {"abcZ": {"identity": 0.85, "source": "alignment"}}
        assert merge_calls(minhash_calls, borderline, True)["abcZ"]["source"] == (
            "alignment"
        )


class TestAlignmentCallPreferred:
    @pytest.mark.parametrize(
        "identity,threshold,expected",
        [
            (0.99, 0.85, True),
            (0.85, 0.85, True),
            (0.84, 0.85, False),
            (0.0, 0.0, True),
            (1.0, 1.0, True),
            (0.999, 1.0, False),
        ],
    )
    def test_inclusive_threshold(self, identity, threshold, expected):
        assert alignment_call_preferred({"identity": identity}, threshold) is expected

    def test_missing_identity_treated_as_zero(self):
        assert alignment_call_preferred({}, 0.85) is False
        assert alignment_call_preferred({}, 0.0) is True


class TestNeedsAlignment:
    def test_empty_calls_need_no_alignment(self):
        assert needs_alignment({}, 0.85) is False

    @pytest.mark.parametrize(
        "similarity,expected",
        [(0.84, True), (0.85, False), (0.86, False)],
    )
    def test_similarity_gate_at_and_around_threshold(self, similarity, expected):
        calls = {"abcZ": {"similarity": similarity}}
        assert needs_alignment(calls, 0.85) is expected

    def test_falsey_confidence_demands_alignment(self):
        assert needs_alignment({"abcZ": {"confidence": False}}, 0.85) is True
        assert needs_alignment({"abcZ": {"confidence": 0.0}}, 0.85) is True

    def test_truthy_confidence_accepts_call(self):
        assert needs_alignment({"abcZ": {"confidence": True}}, 0.85) is False
        assert needs_alignment({"abcZ": {"confidence": 0.5}}, 0.85) is False

    def test_confidence_takes_precedence_over_similarity(self):
        """A present confidence field short-circuits the similarity check."""
        calls = {"abcZ": {"confidence": True, "similarity": 0.10}}
        assert needs_alignment(calls, 0.85) is False

    def test_call_without_confidence_or_similarity_is_accepted(self):
        assert needs_alignment({"abcZ": {"allele": "1"}}, 0.85) is False

    def test_any_weak_locus_triggers_alignment(self):
        calls = {"abcZ": {"similarity": 0.99}, "adk": {"similarity": 0.10}}
        assert needs_alignment(calls, 0.85) is True

    def test_zero_threshold_never_triggers_on_similarity(self):
        assert needs_alignment({"abcZ": {"similarity": 0.0}}, 0.0) is False

    def test_threshold_defaults_to_module_default(self):
        assert needs_alignment({"abcZ": {"similarity": 0.84}}) is True
        assert needs_alignment({"abcZ": {"similarity": 0.85}}) is False


class TestLocusNeedsAlignment:
    @pytest.mark.parametrize(
        "call_data,expected",
        [
            ({"confidence": False}, True),
            ({"confidence": True}, False),
            ({"similarity": 0.10}, True),
            ({"similarity": 0.90}, False),
            ({}, False),
        ],
    )
    def test_single_locus_decision(self, call_data, expected):
        assert locus_needs_alignment(call_data, 0.85) is expected


class TestAddExclusionMetadata:
    def test_exclusions_land_under_notes(self, some_exclusions):
        annotated = add_exclusion_metadata({"st": "11"}, some_exclusions)
        assert annotated["notes"]["exclusions"] == {
            "excluded_alleles": ["abcZ_3"],
            "excluded_loci": ["tonB"],
            "num_excluded_alleles": 1,
            "num_excluded_loci": 1,
        }

    def test_result_fields_survive(self, some_exclusions):
        annotated = add_exclusion_metadata(
            {"st": "11", "profile": {"abcZ": "1"}}, some_exclusions
        )
        assert annotated["st"] == "11"
        assert annotated["profile"] == {"abcZ": "1"}

    def test_nothing_excluded_still_records_empty_metadata(self, no_exclusions):
        """The notes block is unconditional, so consumers never branch on it."""
        annotated = add_exclusion_metadata({"st": "11"}, no_exclusions)
        assert annotated["notes"]["exclusions"]["excluded_alleles"] == []
        assert annotated["notes"]["exclusions"]["num_excluded_loci"] == 0

    def test_existing_notes_are_preserved(self, some_exclusions):
        annotated = add_exclusion_metadata(
            {"notes": {"method": "balanced"}}, some_exclusions
        )
        assert annotated["notes"]["method"] == "balanced"
        assert "exclusions" in annotated["notes"]

    def test_existing_exclusions_are_replaced(self, some_exclusions):
        annotated = add_exclusion_metadata(
            {"notes": {"exclusions": {"stale": True}}}, some_exclusions
        )
        assert "stale" not in annotated["notes"]["exclusions"]

    def test_extra_exclusion_fields_are_dropped(self, some_exclusions):
        some_exclusions["excluded_profiles"] = ["ST11"]
        annotated = add_exclusion_metadata({}, some_exclusions)
        assert "excluded_profiles" not in annotated["notes"]["exclusions"]

    def test_empty_result_gets_only_notes(self, no_exclusions):
        assert list(add_exclusion_metadata({}, no_exclusions)) == ["notes"]

    def test_malformed_exclusions_raise(self):
        with pytest.raises(KeyError):
            add_exclusion_metadata({"st": "11"}, {"excluded_alleles": []})

    def test_result_is_not_mutated(self, some_exclusions):
        result = {"st": "11", "notes": {"method": "balanced"}}
        add_exclusion_metadata(result, some_exclusions)
        assert result == {"st": "11", "notes": {"method": "balanced"}}
