#!/usr/bin/env python

"""Tests for the profile -> sequence type association module.

Fixtures mirror the parity harness (.scratch/parity_profile_match.py) that
diffed this module against the verbatim `lookup_profile` heredoc from
torchbase/workflows/builtin/tasks/profile_lookup.wdl.
"""

import json

import pytest

from torchbase.profile_match import (
    STATUS_KNOWN,
    STATUS_NOVEL_PROFILE,
    UNKNOWN_PROFILE_ID,
    build_profile_record,
    build_profile_string,
    find_nearest_st,
    load_allele_calls,
    load_profiles,
    lookup_profile,
    match_profile,
    mean_confidence,
    row_profile_string,
    scheme_name,
    st_column,
)

PROFILES_TSV = "ST\tadk\tfumC\tgyrB\n1\t10\t11\t12\n2\t10\t99\t12\n3\t?\t11\t12\n"
LOCI = ["adk", "fumC", "gyrB"]


def called(**alleles):
    """Allele calls in the pipeline's shape: allele_id plus a boolean pass."""
    return {
        locus: {"allele_id": allele, "confidence": True}
        for locus, allele in alleles.items()
    }


@pytest.fixture
def table(tmp_path):
    """A three-ST profiles table inside a scheme-named directory."""
    scheme_dir = tmp_path / "ecoli_achtman_4"
    scheme_dir.mkdir()
    path = scheme_dir / "profiles.tsv"
    path.write_text(PROFILES_TSV, newline="")
    return path


@pytest.fixture
def profiles(table):
    return load_profiles(str(table))


class TestLoadProfiles:
    def test_st_column_is_not_a_locus(self, profiles):
        rows, loci_order = profiles
        assert loci_order == LOCI
        assert len(rows) == 3
        assert rows[0]["ST"] == "1"

    def test_lowercase_st_column_still_excluded(self, tmp_path):
        path = tmp_path / "profiles.tsv"
        path.write_text("st\tadk\tfumC\n7\t10\t11\n", newline="")
        rows, loci_order = load_profiles(str(path))
        assert loci_order == ["adk", "fumC"]
        assert st_column(rows[0]) == "st"

    def test_header_only_table_has_loci_but_no_rows(self, tmp_path):
        path = tmp_path / "profiles.tsv"
        path.write_text("ST\tadk\tfumC\tgyrB\n", newline="")
        rows, loci_order = load_profiles(str(path))
        assert rows == []
        assert loci_order == LOCI

    def test_empty_file_has_neither_rows_nor_loci(self, tmp_path):
        path = tmp_path / "profiles.tsv"
        path.write_text("", newline="")
        assert load_profiles(str(path)) == ([], [])

    def test_crlf_table_parses_without_stray_carriage_returns(self, tmp_path):
        path = tmp_path / "profiles.tsv"
        path.write_bytes(b"ST\tadk\tfumC\tgyrB\r\n1\t10\t11\t12\r\n")
        rows, loci_order = load_profiles(str(path))
        assert loci_order == LOCI
        assert rows[0]["gyrB"] == "12"

    def test_no_st_column_at_all(self, tmp_path):
        path = tmp_path / "profiles.tsv"
        path.write_text("adk\tfumC\tgyrB\n10\t11\t12\n", newline="")
        rows, loci_order = load_profiles(str(path))
        assert loci_order == LOCI
        assert st_column(rows[0]) is None


class TestBuildProfileString:
    def test_loci_are_emitted_in_table_order_not_call_order(self):
        calls = called(gyrB=12, adk=10, fumC=11)
        assert build_profile_string(calls, LOCI) == "10,11,12"

    def test_uncalled_locus_becomes_a_wildcard(self):
        assert build_profile_string(called(adk=10, gyrB=12), LOCI) == "10,?,12"

    def test_call_without_allele_id_becomes_a_wildcard(self):
        calls = {"adk": {"confidence": True}, "fumC": {"allele_id": 11}}
        assert build_profile_string(calls, LOCI) == "?,11,?"

    def test_no_calls_yields_an_all_wildcard_profile(self):
        assert build_profile_string({}, LOCI) == "?,?,?"

    def test_scheme_prefixed_call_keys_do_not_satisfy_bare_columns(self):
        calls = called(**{"ecoli_adk": 10, "ecoli_fumC": 11})
        assert build_profile_string(calls, LOCI) == "?,?,?"

    def test_no_loci_yields_the_empty_profile(self):
        assert build_profile_string(called(adk=10), []) == ""


class TestMatchProfile:
    def test_exact_match_returns_the_st(self, profiles):
        rows, loci_order = profiles
        assert match_profile("10,11,12", rows, loci_order) == ("1", STATUS_KNOWN)

    def test_first_matching_row_wins(self, profiles):
        rows, loci_order = profiles
        # '?,11,12' is consistent with both ST 1 and ST 3; table order decides.
        assert match_profile("?,11,12", rows, loci_order)[0] == "1"

    def test_one_locus_mismatch_is_a_novel_profile(self, profiles):
        rows, loci_order = profiles
        assert match_profile("10,11,55", rows, loci_order) == (
            None,
            STATUS_NOVEL_PROFILE,
        )

    def test_query_wildcard_matches_any_table_value(self, profiles):
        rows, loci_order = profiles
        assert match_profile("10,?,12", rows, loci_order)[0] == "1"

    def test_table_wildcard_matches_any_query_value(self, profiles):
        rows, loci_order = profiles
        # adk=77 is in no row; only ST 3, whose adk is '?', can absorb it.
        assert match_profile("77,11,12", rows, loci_order)[0] == "3"

    def test_all_wildcard_query_matches_the_first_row(self, profiles):
        rows, loci_order = profiles
        # Documents current behaviour: an entirely uncalled isolate is typed
        # as the first ST in the table rather than reported as untypeable.
        assert match_profile("?,?,?", rows, loci_order) == ("1", STATUS_KNOWN)

    def test_empty_table_is_always_novel(self):
        assert match_profile("10,11,12", [], LOCI) == (None, STATUS_NOVEL_PROFILE)

    def test_rows_of_a_different_length_are_skipped(self):
        rows = [{"ST": "1", "adk": "10", "fumC": "11"}]
        assert match_profile("10,11,12", rows, LOCI) == (None, STATUS_NOVEL_PROFILE)

    def test_table_without_an_st_column_can_only_be_novel(self):
        rows = [{"adk": "10", "fumC": "11", "gyrB": "12"}]
        assert match_profile("10,11,12", rows, LOCI) == (None, STATUS_NOVEL_PROFILE)

    def test_unknown_allele_is_reported_as_novel_profile_not_novel_allele(
        self, profiles
    ):
        rows, loci_order = profiles
        # A never-before-seen allele at one locus is indistinguishable here from
        # a new combination of known alleles: both come back novel_profile.
        _, unknown_allele = match_profile("10,11,999", rows, loci_order)
        _, new_combination = match_profile("10,99,55", rows, loci_order)
        assert unknown_allele == new_combination == STATUS_NOVEL_PROFILE

    def test_ragged_row_raises(self):
        # csv.DictReader fills a short row with None; the heredoc's join of the
        # row's values raises. Pre-existing behaviour, preserved by the lift.
        rows = [{"ST": "1", "adk": "10", "fumC": "11", "gyrB": None}]
        with pytest.raises(TypeError):
            match_profile("10,11,12", rows, LOCI)


class TestRowProfileString:
    def test_absent_locus_shortens_the_row(self):
        assert row_profile_string({"ST": "1", "adk": "10"}, LOCI) == "10"


class TestFindNearestSt:
    def test_nearest_st_is_the_smallest_hamming_distance(self, profiles):
        rows, loci_order = profiles
        assert find_nearest_st("10,99,55", rows, loci_order) == ("2", 1)

    def test_ties_go_to_the_first_row_in_table_order(self, profiles):
        rows, loci_order = profiles
        # ST 1 and ST 3 are both distance 1 from '10,11,55'.
        assert find_nearest_st("10,11,55", rows, loci_order) == ("1", 1)

    def test_wildcards_cost_nothing_on_either_side(self, profiles):
        rows, loci_order = profiles
        assert find_nearest_st("?,?,?", rows, loci_order) == ("1", 0)

    def test_missing_locus_columns_are_wildcarded_not_dropped(self):
        # Unlike match_profile, a row lacking a locus stays comparable.
        rows = [{"ST": "4", "adk": "10", "fumC": "11"}]
        assert find_nearest_st("10,11,12", rows, LOCI) == ("4", 0)

    def test_no_st_column_gives_no_neighbour(self):
        rows = [{"adk": "10", "fumC": "11", "gyrB": "12"}]
        assert find_nearest_st("10,11,12", rows, LOCI) == (None, float("inf"))


class TestMeanConfidence:
    def test_boolean_confidence_counts_as_pass_or_fail(self):
        calls = {
            "adk": {"confidence": True},
            "fumC": {"confidence": False},
        }
        assert mean_confidence(calls) == 0.5

    def test_boolean_confidence_outranks_a_numeric_similarity(self):
        # The pipeline's own calls always carry both; the boolean wins, so the
        # reported confidence is a pass fraction, not a mean similarity.
        assert mean_confidence({"adk": {"confidence": True, "similarity": 0.42}}) == 1.0

    def test_numeric_similarity_is_used_when_no_confidence_flag(self):
        calls = {"adk": {"similarity": 0.5}, "fumC": {"identity": 1.0}}
        assert mean_confidence(calls) == 0.75

    def test_similarity_wins_over_identity(self):
        assert mean_confidence({"adk": {"similarity": 0.25, "identity": 1.0}}) == 0.25

    def test_unscored_calls_do_not_dilute_the_mean(self):
        calls = {"adk": {"confidence": True}, "fumC": {"allele_id": 11}}
        assert mean_confidence(calls) == 1.0

    def test_no_calls_scores_zero(self):
        assert mean_confidence({}) == 0.0


class TestSchemeName:
    def test_explicit_scheme_wins(self):
        assert scheme_name("salmonella_mlst", "/torch/ecoli/profiles.tsv") == (
            "salmonella_mlst"
        )

    def test_scheme_falls_back_to_the_tables_directory(self):
        assert scheme_name("", "/torch/ecoli_achtman_4/profiles.tsv") == (
            "ecoli_achtman_4"
        )

    def test_bare_filename_leaves_the_scheme_empty(self):
        assert scheme_name("", "profiles.tsv") == ""


class TestBuildProfileRecord:
    def test_known_profile_record_shape(self, profiles):
        rows, loci_order = profiles
        record = build_profile_record(
            called(adk=10, fumC=11, gyrB=12), rows, loci_order, "ecoli_achtman_4"
        )
        assert record["profile_id"] == "1"
        assert record["status"] == STATUS_KNOWN
        assert record["profile_type"] == "sequence_type"
        assert record["scheme"] == "ecoli_achtman_4"
        assert record["allele_profile"] == "10,11,12"
        assert record["confidence"] == 1.0
        assert record["notes"] == {
            "num_loci": 3,
            "num_called": 3,
            "mean_confidence": 1.0,
        }
        assert "nearest_st" not in record
        assert "nearest_st_distance" not in record

    def test_novel_profile_record_carries_its_nearest_neighbour(self, profiles):
        rows, loci_order = profiles
        record = build_profile_record(
            called(adk=10, fumC=11, gyrB=55), rows, loci_order, "ecoli_achtman_4"
        )
        assert record["status"] == STATUS_NOVEL_PROFILE
        assert record["profile_id"] == UNKNOWN_PROFILE_ID
        assert record["nearest_st"] == "1"
        assert record["nearest_st_distance"] == 1

    def test_empty_table_novel_profile_has_no_nearest_neighbour(self):
        record = build_profile_record(called(adk=10), [], LOCI, "ecoli_achtman_4")
        assert record["status"] == STATUS_NOVEL_PROFILE
        assert record["profile_id"] == UNKNOWN_PROFILE_ID
        assert "nearest_st" not in record
        assert record["notes"]["num_loci"] == 3

    def test_blank_st_value_reports_known_but_unidentified(self):
        rows = [{"ST": "", "adk": "10", "fumC": "11", "gyrB": "12"}]
        record = build_profile_record(
            called(adk=10, fumC=11, gyrB=12), rows, LOCI, "ecoli_achtman_4"
        )
        assert record["status"] == STATUS_KNOWN
        assert record["profile_id"] == UNKNOWN_PROFILE_ID

    def test_confidence_is_clamped_but_notes_keep_the_raw_mean(self, profiles):
        rows, loci_order = profiles
        record = build_profile_record(
            {"adk": {"allele_id": 10, "similarity": 1.5}},
            rows,
            loci_order,
            "ecoli_achtman_4",
        )
        assert record["confidence"] == 1.0
        assert record["notes"]["mean_confidence"] == 1.5

    def test_alignment_selects_the_alignment_toolchain(self, profiles):
        rows, loci_order = profiles
        calls = called(adk=10, fumC=11, gyrB=12)
        sketch_only = build_profile_record(calls, rows, loci_order, "s", "fast", False)
        aligned = build_profile_record(calls, rows, loci_order, "s", "sensitive", True)
        assert sketch_only["method"] == {
            "strategy": "fast",
            "alignment_used": False,
            "tools": ["sourmash"],
        }
        assert aligned["method"] == {
            "strategy": "sensitive",
            "alignment_used": True,
            "tools": ["sourmash", "minimap2"],
        }

    def test_uncalled_loci_still_type_against_the_called_ones(self, profiles):
        rows, loci_order = profiles
        record = build_profile_record(
            called(adk=10, gyrB=12), rows, loci_order, "ecoli_achtman_4"
        )
        assert record["allele_profile"] == "10,?,12"
        assert record["status"] == STATUS_KNOWN
        assert record["profile_id"] == "1"
        assert record["notes"]["num_called"] == 2

    def test_record_is_json_serialisable(self, profiles):
        rows, loci_order = profiles
        record = build_profile_record(
            called(adk=10, fumC=11, gyrB=55), rows, loci_order, "ecoli_achtman_4"
        )
        assert json.loads(json.dumps(record))["nearest_st"] == "1"


class TestLookupProfile:
    def test_end_to_end_from_disk(self, tmp_path, table):
        calls_path = tmp_path / "allele_calls.json"
        calls_path.write_text(json.dumps(called(adk=10, fumC=11, gyrB=12)))
        record = lookup_profile(str(calls_path), str(table))
        assert record["profile_id"] == "1"
        assert record["status"] == STATUS_KNOWN
        # scheme inferred from the profiles table's parent directory
        assert record["scheme"] == "ecoli_achtman_4"
        assert record["allele_calls"] == load_allele_calls(str(calls_path))

    def test_explicit_scheme_overrides_the_inferred_one(self, tmp_path, table):
        calls_path = tmp_path / "allele_calls.json"
        calls_path.write_text(json.dumps({}))
        record = lookup_profile(str(calls_path), str(table), scheme="salmonella_mlst")
        assert record["scheme"] == "salmonella_mlst"


SEROTYPE_TSV = "Serotype\twzx\twzy\n1,2a\tA\tB\n3\t?\tB\n"
SEROTYPE_LOCI = ["wzx", "wzy"]


@pytest.fixture
def serotype_table(tmp_path):
    """An ECTyper/SeqSero2-shaped table: identifier column is not named ST."""
    scheme_dir = tmp_path / "shigatyper"
    scheme_dir.mkdir()
    path = scheme_dir / "profiles.tsv"
    path.write_text(SEROTYPE_TSV, newline="")
    return path


class TestExplicitIdColumn:
    """A named non-ST identifier column (docs/adr/0003), e.g. Serotype."""

    def test_load_profiles_excludes_the_named_column_from_loci(self, serotype_table):
        rows, loci_order = load_profiles(str(serotype_table), id_column="Serotype")
        assert loci_order == SEROTYPE_LOCI
        assert rows[0]["Serotype"] == "1,2a"

    def test_without_id_column_serotype_is_treated_as_a_locus(self, serotype_table):
        # Unchanged auto-detect behaviour: only a literal ST column is special.
        rows, loci_order = load_profiles(str(serotype_table))
        assert loci_order == ["Serotype"] + SEROTYPE_LOCI

    def test_st_column_honours_the_explicit_name(self):
        assert st_column({"Serotype": "1,2a", "wzx": "A"}, id_column="Serotype") == "Serotype"

    def test_st_column_explicit_name_absent_from_row_is_none(self):
        assert st_column({"wzx": "A"}, id_column="Serotype") is None

    def test_match_profile_matches_on_the_named_column(self, serotype_table):
        rows, loci_order = load_profiles(str(serotype_table), id_column="Serotype")
        profile_id, status = match_profile("A,B", rows, loci_order, id_column="Serotype")
        assert (profile_id, status) == ("1,2a", STATUS_KNOWN)

    def test_lookup_profile_end_to_end_with_named_id_column(self, tmp_path, serotype_table):
        calls_path = tmp_path / "allele_calls.json"
        calls_path.write_text(json.dumps(called(wzx="A", wzy="B")))
        record = lookup_profile(str(calls_path), str(serotype_table), id_column="Serotype")
        assert record["profile_id"] == "1,2a"
        assert record["status"] == STATUS_KNOWN
