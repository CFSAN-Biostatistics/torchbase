#!/usr/bin/env python

"""Tests for the input-preparation filters lifted from the filter_alleles and
depth_filter WDL heredocs.

Fixtures mirror `.scratch/fix_qf`, the inputs used for the byte-for-byte
parity run against the original heredocs.
"""

import json
from textwrap import dedent

import pytest

from torchbase.quality_filters import (
    build_exclusion_sets,
    depth_filter,
    empty_quality,
    extract_locus_and_allele,
    filter_alleles,
    kmer_depth,
    load_quality,
    parse_allele_fasta,
    parse_sequence_fasta,
    quality_from_data,
    write_fasta,
    write_sequence_fasta,
)

ALLELES = dedent(
    """\
    >locusA_1
    ACGTACGTAA
    >locusA_2
    ACGTACGTAC
    >locusB_1
    TTTTTTTTTT
    >locusB_2
    TTTTTTTTTA
    >locusC_1
    GGGGGGGGGG
    >locusC_2
    GGGGGGGGGA
    >locusC_3
    GGGGGGGGGT
    >locusD_1
    CCCCCCCCCC
    >locusE_1
    AAAAAAAAAA
    >noUnderscore
    ACGT
    >weird_locus_name_7
    ACGTT
    >wrapped_locus_2
    ACGT
    ACGTAC
    """
)

QUALITY = {
    "loci": {
        # whole locus flagged
        "locusA": {"suspect": True},
        # single allele flagged
        "locusB": {"alleles": {"1": {"suspect": True}, "2": {"suspect": False}}},
        # explicit threshold: only the 94.9 pair is below it
        "locusC": {
            "threshold": 95.0,
            "similarities": {"1-2": 94.9, "1-3": 95.0, "2-3": 99.0},
        },
        # default threshold (90.0); names allele 2, absent from the FASTA
        "locusE": {"similarities": {"1-2": 89.99}},
        "weird_locus_name": {"alleles": {"7": {"suspect": True}}},
        # locus absent from the FASTA entirely
        "locusZ": {"suspect": True, "alleles": {"9": {"suspect": True}}},
        "wrapped_locus": {"alleles": {"2": {"suspect": True}}},
    },
    "profiles": {
        "p1": {"suspect": True, "loci": ["locusD", "locusB"]},
        "p2": {"suspect": False, "loci": ["locusC"]},
        "p3": {"suspect": True},  # suspect but declares no loci
        "p4": {"suspect": True, "loci": ["locusMissing"]},
    },
}


@pytest.fixture
def allele_fasta(tmp_path):
    path = tmp_path / "alleles.fasta"
    path.write_text(ALLELES)
    return path


@pytest.fixture
def quality_path(tmp_path):
    path = tmp_path / "quality.json"
    path.write_text(json.dumps(QUALITY, indent=2))
    return path


@pytest.fixture
def quality(quality_path):
    return load_quality(quality_path)


def headers(entries):
    return [header for header, _seq in entries]


# ---------------------------------------------------------------- parsing


def test_parse_allele_fasta_strips_marker_and_unwraps(allele_fasta):
    entries = parse_allele_fasta(allele_fasta)
    assert len(entries) == 12
    assert entries[0] == ("locusA_1", "ACGTACGTAA")
    # wrapped record is concatenated into one sequence
    assert entries[-1] == ("wrapped_locus_2", "ACGTACGTAC")


def test_parse_allele_fasta_empty_file(tmp_path):
    path = tmp_path / "empty.fasta"
    path.write_text("")
    assert parse_allele_fasta(path) == []


@pytest.mark.parametrize(
    "header,expected",
    [
        ("locusA_1", ("locusA", "1")),
        ("weird_locus_name_7", ("weird_locus_name", "7")),
        ("noUnderscore", ("noUnderscore", "unknown")),
        ("_1", ("", "1")),
    ],
)
def test_extract_locus_and_allele(header, expected):
    assert extract_locus_and_allele(header) == expected


# ---------------------------------------------------------------- load_quality


@pytest.mark.parametrize("missing", [None, ""])
def test_load_quality_no_path_is_empty(missing):
    assert load_quality(missing) == empty_quality()


def test_load_quality_nonexistent_path_is_empty(tmp_path):
    assert load_quality(tmp_path / "nope.json") == empty_quality()


@pytest.mark.parametrize("body", ["", "   \n\t \n"])
def test_load_quality_blank_file_is_empty(tmp_path, body):
    path = tmp_path / "quality.json"
    path.write_text(body)
    assert load_quality(path) == empty_quality()


def test_load_quality_no_keys_is_empty(tmp_path):
    path = tmp_path / "quality.json"
    path.write_text("{}")
    assert load_quality(path) == empty_quality()


def test_load_quality_malformed_json_raises(tmp_path):
    path = tmp_path / "quality.json"
    path.write_text("{not json")
    with pytest.raises(ValueError):
        load_quality(path)


def test_quality_from_data_collects_all_three_sources(quality):
    assert quality["suspect_loci"] == {"locusA", "locusZ"}
    assert quality["suspect_alleles"] == {
        "locusB_1",
        "locusC_1",
        "locusC_2",
        "locusE_1",
        "locusE_2",
        "weird_locus_name_7",
        "locusZ_9",
        "wrapped_locus_2",
    }
    # only suspect profiles that declare loci contribute
    assert quality["suspect_profiles"] == {"locusD", "locusB", "locusMissing"}


@pytest.mark.parametrize(
    "similarity,threshold,suspect",
    [
        (94.9, 95.0, True),
        (95.0, 95.0, False),  # boundary is exclusive: >= threshold is fine
        (95.1, 95.0, False),
        (89.99, None, True),  # default threshold 90.0
        (90.0, None, False),
    ],
)
def test_similarity_threshold_boundary(similarity, threshold, suspect):
    locus = {"similarities": {"1-2": similarity}}
    if threshold is not None:
        locus["threshold"] = threshold
    result = quality_from_data({"loci": {"L": locus}})
    expected = {"L_1", "L_2"} if suspect else set()
    assert result["suspect_alleles"] == expected


# ---------------------------------------------------------------- flag logic


def test_no_flags_excludes_nothing(quality):
    assert build_exclusion_sets(quality) == (set(), set())


def test_exclude_alleles_flag_only_touches_alleles(quality):
    loci, alleles = build_exclusion_sets(quality, exclude_alleles=True)
    assert loci == set()
    assert alleles == quality["suspect_alleles"]


def test_exclude_loci_flag_adds_suspect_loci(quality):
    loci, alleles = build_exclusion_sets(quality, exclude_loci=True)
    assert loci == quality["suspect_loci"]
    assert alleles == quality["suspect_alleles"]


def test_exclude_profiles_flag_adds_profile_loci(quality):
    loci, alleles = build_exclusion_sets(quality, exclude_profiles=True)
    assert loci == quality["suspect_loci"] | quality["suspect_profiles"]
    assert alleles == quality["suspect_alleles"]


def test_flags_are_a_precedence_chain_not_independent(quality):
    """The most aggressive flag set wins; weaker ones add nothing."""
    profiles_only = build_exclusion_sets(quality, exclude_profiles=True)
    assert (
        build_exclusion_sets(
            quality,
            exclude_alleles=True,
            exclude_loci=True,
            exclude_profiles=True,
        )
        == profiles_only
    )
    loci_only = build_exclusion_sets(quality, exclude_loci=True)
    assert (
        build_exclusion_sets(quality, exclude_alleles=True, exclude_loci=True)
        == loci_only
    )


# ---------------------------------------------------------------- filter_alleles


def test_filter_alleles_default_includes_everything(allele_fasta, quality):
    kept, report = filter_alleles(allele_fasta, quality)
    assert len(kept) == 12
    assert report == {
        "excluded_alleles": [],
        "excluded_loci": [],
        "num_excluded_alleles": 0,
        "num_excluded_loci": 0,
        "total_input_alleles": 12,
        "total_output_alleles": 12,
    }


def test_filter_alleles_without_quality_includes_everything_under_all_flags(
    allele_fasta,
):
    kept, report = filter_alleles(
        allele_fasta,
        None,
        exclude_alleles=True,
        exclude_loci=True,
        exclude_profiles=True,
    )
    assert len(kept) == 12
    assert report["num_excluded_alleles"] == 0
    assert report["num_excluded_loci"] == 0


def test_filter_alleles_exclude_alleles(allele_fasta, quality):
    kept, report = filter_alleles(allele_fasta, quality, exclude_alleles=True)
    assert headers(kept) == [
        "locusA_1",
        "locusA_2",
        "locusB_2",
        "locusC_3",
        "locusD_1",
        "noUnderscore",
    ]
    assert report["excluded_alleles"] == [
        "locusB_1",
        "locusC_1",
        "locusC_2",
        "locusE_1",
        "weird_locus_name_7",
        "wrapped_locus_2",
    ]
    assert report["excluded_loci"] == []
    assert report["total_output_alleles"] == 6


def test_filter_alleles_exclude_loci(allele_fasta, quality):
    kept, report = filter_alleles(allele_fasta, quality, exclude_loci=True)
    assert headers(kept) == ["locusB_2", "locusC_3", "locusD_1", "noUnderscore"]
    assert sorted(report["excluded_loci"]) == ["locusA"]
    # locus-level drops are not also listed as excluded alleles; locusA held no
    # suspect allele of its own, so all 6 suspect alleles are still counted.
    assert "locusA_1" not in report["excluded_alleles"]
    assert report["num_excluded_alleles"] == 6
    assert report["total_output_alleles"] == 4


def test_filter_alleles_exclude_profiles(allele_fasta, quality):
    kept, report = filter_alleles(allele_fasta, quality, exclude_profiles=True)
    assert headers(kept) == ["locusC_3", "noUnderscore"]
    assert sorted(report["excluded_loci"]) == ["locusA", "locusB", "locusD"]
    assert report["excluded_alleles"] == [
        "locusC_1",
        "locusC_2",
        "locusE_1",
        "weird_locus_name_7",
        "wrapped_locus_2",
    ]
    assert report["num_excluded_loci"] == 3
    assert report["total_input_alleles"] == 12
    assert report["total_output_alleles"] == 2


def test_filter_alleles_all_flags_matches_profiles_only(allele_fasta, quality):
    both = filter_alleles(
        allele_fasta,
        quality,
        exclude_alleles=True,
        exclude_loci=True,
        exclude_profiles=True,
    )
    only = filter_alleles(allele_fasta, quality, exclude_profiles=True)
    assert headers(both[0]) == headers(only[0])
    assert both[1] == only[1]


def test_filter_alleles_missing_quality_file_includes_everything(
    allele_fasta, tmp_path
):
    kept, report = filter_alleles(
        allele_fasta,
        load_quality(tmp_path / "absent.json"),
        exclude_profiles=True,
    )
    assert len(kept) == 12
    assert report["num_excluded_loci"] == 0


def test_filter_alleles_empty_quality_file_includes_everything(
    allele_fasta, tmp_path
):
    path = tmp_path / "quality.json"
    path.write_text("")
    kept, _ = filter_alleles(
        allele_fasta, load_quality(path), exclude_profiles=True
    )
    assert len(kept) == 12


def test_filter_alleles_quality_naming_absent_loci_and_alleles(
    allele_fasta, tmp_path
):
    """Suspect names not present in the FASTA are silently irrelevant."""
    path = tmp_path / "quality.json"
    path.write_text(
        json.dumps(
            {
                "loci": {
                    "ghostLocus": {
                        "suspect": True,
                        "alleles": {"1": {"suspect": True}},
                    }
                },
                "profiles": {"g1": {"suspect": True, "loci": ["otherGhost"]}},
            }
        )
    )
    loaded = load_quality(path)
    assert loaded["suspect_loci"] == {"ghostLocus"}
    kept, report = filter_alleles(allele_fasta, loaded, exclude_profiles=True)
    assert len(kept) == 12
    assert report["excluded_alleles"] == []
    assert report["excluded_loci"] == []


def test_filter_alleles_empty_fasta(tmp_path, quality):
    path = tmp_path / "empty.fasta"
    path.write_text("")
    kept, report = filter_alleles(path, quality, exclude_profiles=True)
    assert kept == []
    assert report["total_input_alleles"] == 0
    assert report["total_output_alleles"] == 0
    assert report["num_excluded_alleles"] == 0


def test_filter_alleles_headerless_locus_uses_unknown_allele_id(tmp_path):
    path = tmp_path / "one.fasta"
    path.write_text(">solo\nACGT\n")
    quality = {
        "suspect_alleles": {"solo_unknown"},
        "suspect_loci": set(),
        "suspect_profiles": set(),
    }
    kept, report = filter_alleles(path, quality, exclude_alleles=True)
    assert kept == []
    assert report["excluded_alleles"] == ["solo_unknown"]


def test_write_fasta_round_trips_through_parse(tmp_path, allele_fasta, quality):
    kept, _ = filter_alleles(allele_fasta, quality, exclude_alleles=True)
    out = tmp_path / "filtered.fasta"
    write_fasta(kept, out)
    assert out.read_text().startswith(">locusA_1\nACGTACGTAA\n")
    assert parse_allele_fasta(out) == kept


# ---------------------------------------------------------------- depth_filter

# With k=3: 10xA -> 8.0, ACGTACGTAC -> 2.0, ACGT -> 1.0, AC -> 0.0 (short)
MIXED = ">high\nAAAAAAAAAA\n>boundary\nACGTACGTAC\n>low\nACGT\n>tooshort\nAC\n"


@pytest.fixture
def mixed_fasta(tmp_path):
    path = tmp_path / "mixed.fasta"
    path.write_text(MIXED)
    return path


def test_parse_sequence_fasta_keeps_marker_on_header(mixed_fasta):
    records = parse_sequence_fasta(mixed_fasta)
    assert records[0] == (">high", "AAAAAAAAAA")
    assert len(records) == 4


def test_parse_sequence_fasta_unwraps(tmp_path):
    path = tmp_path / "wrapped.fasta"
    path.write_text(">wrapped\nAAAAA\nAAAAA\n")
    assert parse_sequence_fasta(path) == [(">wrapped", "AAAAAAAAAA")]


@pytest.mark.parametrize(
    "seq,k,expected",
    [
        ("AAAAAAAAAA", 3, 8.0),  # 8 occurrences of one distinct k-mer
        ("ACGTACGTAC", 3, 2.0),
        ("ACGT", 3, 1.0),
        ("AC", 3, 0.0),  # shorter than k
        ("ACG", 3, 1.0),  # exactly k
    ],
)
def test_kmer_depth(seq, k, expected):
    assert kmer_depth(seq, k) == expected


def test_depth_filter_contigs_is_pass_through(mixed_fasta):
    kept, stats = depth_filter(mixed_fasta, input_type="contigs", min_depth=100, ksize=3)
    assert len(kept) == 4
    assert stats == {
        "input_sequences": 4,
        "kept_sequences": 4,
        "removed_sequences": 0,
        "input_type": "contigs",
        "min_depth": 100,
        "ksize": 3,
    }


def test_depth_filter_default_input_type_is_pass_through(mixed_fasta):
    kept, stats = depth_filter(mixed_fasta)
    assert stats["input_type"] == "contigs"
    assert stats["removed_sequences"] == 0
    assert len(kept) == 4


def test_depth_filter_threshold_is_inclusive(mixed_fasta):
    kept, stats = depth_filter(mixed_fasta, input_type="reads", min_depth=2, ksize=3)
    assert [h for h, _ in kept] == [">high", ">boundary"]
    assert stats["removed_sequences"] == 2

    kept, stats = depth_filter(mixed_fasta, input_type="reads", min_depth=3, ksize=3)
    assert [h for h, _ in kept] == [">high"]
    assert stats["removed_sequences"] == 3


def test_depth_filter_any_non_contig_type_filters(mixed_fasta):
    for input_type in ("reads", "fastq", ""):
        kept, stats = depth_filter(
            mixed_fasta, input_type=input_type, min_depth=8, ksize=3
        )
        assert [h for h, _ in kept] == [">high"]
        assert stats["input_type"] == input_type


def test_depth_filter_short_sequences_removed_even_at_min_depth_one(mixed_fasta):
    kept, stats = depth_filter(mixed_fasta, input_type="reads", min_depth=1, ksize=3)
    assert [h for h, _ in kept] == [">high", ">boundary", ">low"]
    assert stats["removed_sequences"] == 1


def test_depth_filter_large_ksize_removes_everything(mixed_fasta):
    kept, stats = depth_filter(mixed_fasta, input_type="reads", min_depth=3, ksize=21)
    assert kept == []
    assert stats["kept_sequences"] == 0
    assert stats["removed_sequences"] == 4


def test_depth_filter_empty_input(tmp_path):
    path = tmp_path / "empty.fasta"
    path.write_text("")
    for input_type in ("contigs", "reads"):
        kept, stats = depth_filter(path, input_type=input_type)
        assert kept == []
        assert stats["input_sequences"] == 0
        assert stats["kept_sequences"] == 0
        assert stats["removed_sequences"] == 0


def test_write_sequence_fasta_does_not_double_the_marker(tmp_path, mixed_fasta):
    kept, _ = depth_filter(mixed_fasta, input_type="reads", min_depth=2, ksize=3)
    out = tmp_path / "filtered.fasta"
    write_sequence_fasta(kept, out)
    assert out.read_text() == ">high\nAAAAAAAAAA\n>boundary\nACGTACGTAC\n"
    assert parse_sequence_fasta(out) == kept
