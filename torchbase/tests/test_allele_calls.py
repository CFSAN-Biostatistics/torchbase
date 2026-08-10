#!/usr/bin/env python

"""Tests for allele calling: sourmash similarity matrices and minimap2 SAM.

Unit level, no Docker and no minimap2 -- the compute artifacts (a comparison
CSV, a SAM file) are synthesised, because interpreting them is the only thing
under test. Fixtures here are the ones used to prove byte parity against the
WDL heredocs this module replaced.
"""

import json

import pytest

from torchbase import allele_calls


def sam_line(query, ref, mapq="60", tlen="40", nm="0", drop_nm=False):
    """One minimap2-shaped SAM record; ``nm``/``tlen`` drive identity."""
    fields = [query, "0", ref, "1", mapq, "40M", "*", "0", tlen, "ACGT", "*"]
    if not drop_nm:
        fields.append("NM:i:" + nm)
    return "\t".join(fields)


def matrix(header, *rows):
    """Rows of a sourmash comparison CSV as the reader hands them over."""
    return [list(header)] + [list(row) for row in rows]


QUERY_LENGTHS = {"q1": "A" * 40, "q2": "C" * 40, "q3": "G" * 20}


class TestExtractLocusAndAllele:
    @pytest.mark.parametrize("header,expected", [
        ("adk_1", ("adk", "1")),
        ("thrA_ec_42", ("thrA_ec", "42")),
        ("orphan", ("orphan", "unknown")),
        ("adk_", ("adk", "")),
    ])
    def test_last_underscore_field_is_the_allele(self, header, expected):
        assert allele_calls.extract_locus_and_allele(header) == expected


class TestParseFasta:
    def test_multiline_sequences_are_joined_in_file_order(self, tmp_path):
        fasta = tmp_path / "in.fasta"
        fasta.write_text(">adk_2\nACGT\nTTTT\n>adk_1\nGGGG\n")
        assert allele_calls.parse_fasta(str(fasta)) == [
            ("adk_2", "ACGTTTTT"),
            ("adk_1", "GGGG"),
        ]

    def test_empty_file_has_no_sequences(self, tmp_path):
        fasta = tmp_path / "empty.fasta"
        fasta.write_text("")
        assert allele_calls.parse_fasta(str(fasta)) == []

    def test_header_without_sequence_is_kept(self, tmp_path):
        fasta = tmp_path / "in.fasta"
        fasta.write_text(">adk_1\n")
        assert allele_calls.parse_fasta(str(fasta)) == [("adk_1", "")]

    def test_dict_form_collapses_duplicate_headers(self, tmp_path):
        fasta = tmp_path / "in.fasta"
        fasta.write_text(">adk_1\nAAAA\n>adk_1\nCCCC\n")
        assert allele_calls.parse_fasta_dict(str(fasta)) == {"adk_1": "CCCC"}


class TestConfidenceDecision:
    def test_above_excludes_the_threshold(self):
        is_confident = allele_calls.confident_above(0.9)
        assert is_confident(0.9000001) is True
        assert is_confident(0.9) is False

    def test_at_or_above_includes_the_threshold(self):
        is_confident = allele_calls.confident_at_or_above(0.9)
        assert is_confident(0.9) is True
        assert is_confident(0.8999999) is False


class TestGroupAllelesByLocus:
    def test_alleles_keep_their_matrix_column_index(self):
        grouped = allele_calls.group_alleles_by_locus(
            [("adk_1", "A"), ("gyrB_1", "C"), ("adk_2", "G")])
        assert [a["index"] for a in grouped["adk"]] == [0, 2]
        assert [a["allele_id"] for a in grouped["adk"]] == ["1", "2"]
        assert [a["index"] for a in grouped["gyrB"]] == [1]


class TestMaxSimilarityPerAllele:
    ONE_QUERY_TWO_ALLELES = matrix(
        ("q1", "adk_1", "adk_2"),
        ("1.0", "0.4", "0.9"),
        ("0.4", "1.0", "0.2"),
        ("0.9", "0.2", "1.0"),
    )

    def test_allele_columns_follow_the_query_columns(self):
        assert allele_calls.max_similarity_per_allele(
            self.ONE_QUERY_TWO_ALLELES, 1, 2) == [0.4, 0.9]

    def test_best_query_wins_per_allele(self):
        rows = matrix(
            ("q1", "q2", "adk_1", "adk_2"),
            ("1.0", "0.3", "0.4", "0.92"),
            ("0.3", "1.0", "0.99", "0.1"),
            ("0.4", "0.99", "1.0", "0.2"),
            ("0.92", "0.1", "0.2", "1.0"),
        )
        assert allele_calls.max_similarity_per_allele(rows, 2, 2) == [0.99, 0.92]

    def test_blank_cell_counts_as_no_similarity(self):
        rows = matrix(
            ("q1", "adk_1"),
            ("1.0", ""),
            ("", "1.0"),
        )
        assert allele_calls.max_similarity_per_allele(rows, 1, 1) == [0.0]

    def test_short_row_leaves_missing_alleles_at_zero(self):
        rows = matrix(
            ("q1", "adk_1", "adk_2"),
            ("1.0", "0.7"),
            ("0.7", "1.0", "0.2"),
            ("0.0", "0.2", "1.0"),
        )
        assert allele_calls.max_similarity_per_allele(rows, 1, 2) == [0.7, 0.0]

    def test_non_numeric_cell_is_an_error(self):
        rows = matrix(("q1", "adk_1"), ("1.0", "nan-ish"), ("0.0", "1.0"))
        with pytest.raises(ValueError):
            allele_calls.max_similarity_per_allele(rows, 1, 1)

    def test_matrix_that_is_not_the_expected_square_is_an_error(self):
        rows = matrix(("q1", "adk_1", "adk_2"), ("1.0", "0.9", "0.8"))
        with pytest.raises(ValueError) as excinfo:
            allele_calls.max_similarity_per_allele(rows, 1, 2)
        assert str(excinfo.value) == "Matrix size mismatch: expected 4 rows, got 2"


class TestBestSimilarityPerLocus:
    def test_tie_goes_to_the_first_allele_in_fasta_order(self):
        grouped = allele_calls.group_alleles_by_locus(
            [("gyrB_7", "A"), ("gyrB_3", "C")])
        assert allele_calls.best_similarity_per_locus(grouped, [0.95, 0.95]) == [
            ("gyrB", "7", 0.95)]

    def test_loci_are_reported_in_sorted_order(self):
        grouped = allele_calls.group_alleles_by_locus(
            [("gyrB_1", "A"), ("adk_1", "C")])
        called = allele_calls.best_similarity_per_locus(grouped, [0.5, 0.6])
        assert [locus for locus, _, _ in called] == ["adk", "gyrB"]

    def test_locus_with_no_column_in_the_matrix_is_omitted(self):
        grouped = allele_calls.group_alleles_by_locus(
            [("adk_1", "A"), ("gyrB_1", "C")])
        assert allele_calls.best_similarity_per_locus(grouped, [0.5]) == [
            ("adk", "1", 0.5)]


class TestCallsFromSimilarity:
    QUERY = [("q1", "ACGT")]
    ALLELES = [("adk_1", "AAAA"), ("adk_2", "CCCC"), ("gyrB_1", "GGGG")]
    MATRIX = matrix(
        ("q1", "adk_1", "adk_2", "gyrB_1"),
        ("1.0", "1.0", "0.5", "0.0"),
        ("1.0", "1.0", "0.4", "0.3"),
        ("0.5", "0.4", "1.0", "0.2"),
        ("0.0", "0.3", "0.2", "1.0"),
    )

    def call(self, matrix_rows=None, is_confident=None):
        return allele_calls.calls_from_similarity(
            self.MATRIX if matrix_rows is None else matrix_rows,
            self.QUERY,
            self.ALLELES,
            is_confident or allele_calls.confident_at_or_above(0.85),
        )

    def test_exact_match_is_called_and_confident(self):
        assert self.call()["adk"] == {
            "allele_id": "1", "similarity": 1.0, "confidence": True}

    def test_locus_with_no_hit_is_still_called_but_not_confident(self):
        assert self.call()["gyrB"] == {
            "allele_id": "1", "similarity": 0.0, "confidence": False}

    def test_confidence_is_decided_by_the_injected_predicate(self):
        strict = self.call(is_confident=allele_calls.confident_above(1.0))
        assert strict["adk"]["confidence"] is False

    def test_similarity_is_clamped_into_the_unit_interval(self):
        rows = matrix(
            ("q1", "adk_1", "bad_1"),
            ("1.0", "1.4", "-0.2"),
            ("1.4", "1.0", "0.0"),
            ("-0.2", "0.0", "1.0"),
        )
        calls = allele_calls.calls_from_similarity(
            rows, self.QUERY, [("adk_1", "A"), ("bad_1", "C")],
            allele_calls.confident_at_or_above(1.0))
        assert calls["adk"]["similarity"] == 1.0
        assert calls["bad"]["similarity"] == 0.0
        # Confidence is decided on the raw similarity, not the clamped one.
        assert calls["adk"]["confidence"] is True

    def test_header_without_an_allele_id_is_called_unknown(self):
        rows = matrix(("q1", "orphan"), ("1.0", "0.99"), ("0.99", "1.0"))
        calls = allele_calls.calls_from_similarity(
            rows, self.QUERY, [("orphan", "A")],
            allele_calls.confident_at_or_above(0.85))
        assert calls == {"orphan": {
            "allele_id": "unknown", "similarity": 0.99, "confidence": True}}

    @pytest.mark.parametrize("rows", [[], [["q1", "adk_1"]]])
    def test_matrix_without_data_rows_yields_no_calls(self, rows):
        assert allele_calls.calls_from_similarity(
            rows, self.QUERY, self.ALLELES,
            allele_calls.confident_at_or_above(0.85)) == {}

    def test_empty_csv_file_yields_no_calls(self, tmp_path):
        csv_path = tmp_path / "similarity.csv"
        csv_path.write_text("")
        assert allele_calls.calls_from_similarity(
            str(csv_path), self.QUERY, self.ALLELES,
            allele_calls.confident_at_or_above(0.85)) == {}

    def test_paths_and_parsed_data_give_the_same_calls(self, tmp_path):
        csv_path = tmp_path / "similarity.csv"
        csv_path.write_text(
            "\n".join(",".join(row) for row in self.MATRIX) + "\n", newline="")
        query = tmp_path / "query.fasta"
        query.write_text(">q1\nACGT\n")
        alleles = tmp_path / "alleles.fasta"
        alleles.write_text(
            "".join(">{}\n{}\n".format(h, s) for h, s in self.ALLELES))
        assert allele_calls.calls_from_similarity(
            str(csv_path), str(query), str(alleles),
            allele_calls.confident_at_or_above(0.85)) == self.call()


class TestMinhashTaskEntryPoints:
    """``call_alleles`` and ``call_alleles_minhash`` must stay distinguishable."""

    QUERY = [("q1", "ACGT")]
    ALLELES = [("a_1", "A"), ("b_1", "C"), ("c_1", "G"), ("d_1", "T")]
    # Similarities: a=0.9, b just below 0.9, c=0.85, d just below 0.85.
    MATRIX = matrix(
        ("q1", "a_1", "b_1", "c_1", "d_1"),
        ("1.0", "0.9", "0.8999999", "0.85", "0.8499999"),
        ("0.9", "1.0", "0.0", "0.0", "0.0"),
        ("0.8999999", "0.0", "1.0", "0.0", "0.0"),
        ("0.85", "0.0", "0.0", "1.0", "0.0"),
        ("0.8499999", "0.0", "0.0", "0.0", "1.0"),
    )

    def confidence(self, calls):
        return {locus: call["confidence"] for locus, call in calls.items()}

    def test_call_alleles_excludes_its_hardcoded_threshold(self):
        calls = allele_calls.call_alleles(self.MATRIX, self.QUERY, self.ALLELES)
        assert self.confidence(calls) == {
            "a": False, "b": False, "c": False, "d": False}

    def test_call_alleles_minhash_includes_its_default_threshold(self):
        calls = allele_calls.call_alleles_minhash(
            self.MATRIX, self.QUERY, self.ALLELES)
        assert self.confidence(calls) == {
            "a": True, "b": True, "c": True, "d": False}

    def test_the_two_tasks_disagree_at_similarity_exactly_zero_point_nine(self):
        strict = allele_calls.call_alleles(self.MATRIX, self.QUERY, self.ALLELES)
        inclusive = allele_calls.call_alleles_minhash(
            self.MATRIX, self.QUERY, self.ALLELES, 0.9)
        assert strict["a"]["confidence"] is False
        assert inclusive["a"]["confidence"] is True

    def test_minhash_threshold_is_caller_supplied(self):
        calls = allele_calls.call_alleles_minhash(
            self.MATRIX, self.QUERY, self.ALLELES, 1.0)
        assert self.confidence(calls) == {
            "a": False, "b": False, "c": False, "d": False}

    def test_both_tasks_report_the_same_alleles_and_similarities(self):
        strict = allele_calls.call_alleles(self.MATRIX, self.QUERY, self.ALLELES)
        inclusive = allele_calls.call_alleles_minhash(
            self.MATRIX, self.QUERY, self.ALLELES)
        drop = lambda calls: {  # noqa: E731 - confidence is the only difference
            locus: {k: v for k, v in call.items() if k != "confidence"}
            for locus, call in calls.items()
        }
        assert drop(strict) == drop(inclusive)

    def test_tie_between_two_alleles_of_a_locus_takes_the_first(self):
        alleles = [("gyrB_7", "A"), ("gyrB_3", "C")]
        rows = matrix(
            ("q1", "gyrB_7", "gyrB_3"),
            ("1.0", "0.95", "0.95"),
            ("0.95", "1.0", "0.9"),
            ("0.95", "0.9", "1.0"),
        )
        assert allele_calls.call_alleles(rows, self.QUERY, alleles) == {
            "gyrB": {"allele_id": "7", "similarity": 0.95, "confidence": True}}


class TestFormatProfile:
    def test_profile_joins_locus_and_allele_in_call_order(self):
        calls = {"adk": {"allele_id": "1"}, "gyrB": {"allele_id": "12"}}
        assert allele_calls.format_profile(calls) == "adk_1,gyrB_12"

    def test_no_calls_gives_an_empty_profile(self):
        assert allele_calls.format_profile({}) == ""

    def test_unconfident_calls_still_appear_in_the_profile(self):
        calls = allele_calls.call_alleles(
            matrix(("q1", "adk_1"), ("1.0", "0.1"), ("0.1", "1.0")),
            [("q1", "ACGT")], [("adk_1", "A")])
        assert calls["adk"]["confidence"] is False
        assert allele_calls.format_profile(calls) == "adk_1"


class TestParseSamLines:
    def test_headers_and_short_records_are_skipped(self):
        lines = [
            "@HD\tVN:1.6",
            "@SQ\tSN:adk_1\tLN:40",
            "q1\t0\tadk_1\t1\t60\t40M",
            sam_line("q2", "adk_1"),
        ]
        assert [r["query_name"] for r in allele_calls.parse_sam_lines(lines)] == ["q2"]

    def test_unmapped_records_are_skipped(self):
        lines = [
            sam_line("q1", "*"),
            sam_line("q2", "adk_1", mapq="0"),
            sam_line("q3", "adk_1"),
        ]
        assert [r["query_name"] for r in allele_calls.parse_sam_lines(lines)] == ["q3"]

    def test_nm_tag_is_read_and_defaults_to_zero(self):
        records = list(allele_calls.parse_sam_lines([
            sam_line("q1", "adk_1", nm="7"),
            sam_line("q2", "adk_1", drop_nm=True),
        ]))
        assert [r["nm"] for r in records] == [7, 0]


class TestReadSamRecords:
    def test_missing_file_yields_no_records(self, tmp_path):
        assert allele_calls.read_sam_records(str(tmp_path / "absent.sam")) == []

    def test_unparsable_line_keeps_the_records_read_so_far(self, tmp_path):
        sam = tmp_path / "alignment.sam"
        sam.write_text("\n".join([
            sam_line("q1", "adk_1"),
            sam_line("q2", "adk_2", mapq="sixty"),
            sam_line("q3", "gyrB_1"),
        ]) + "\n")
        records = allele_calls.read_sam_records(str(sam))
        assert [r["query_name"] for r in records] == ["q1"]


class TestRecordIdentity:
    def record(self, nm=0, tlen=40, query="q1"):
        return {"query_name": query, "ref_name": "adk_1", "mapq": 60,
                "nm": nm, "tlen": tlen}

    def test_identity_is_matching_bases_over_length(self):
        assert allele_calls.record_identity(
            self.record(nm=4), QUERY_LENGTHS) == pytest.approx(0.9)

    def test_exact_match_is_full_identity(self):
        assert allele_calls.record_identity(self.record(), QUERY_LENGTHS) == 1.0

    def test_query_length_falls_back_to_the_sequence_when_tlen_is_absent(self):
        # q3 is 20 bases; 2 mismatches over 20 is 0.9, not 0.95.
        assert allele_calls.record_identity(
            self.record(nm=2, tlen=0, query="q3"), QUERY_LENGTHS) == pytest.approx(0.9)

    @pytest.mark.parametrize("tlen", [0, -40])
    def test_unknown_query_without_a_usable_tlen_scores_zero(self, tlen):
        assert allele_calls.record_identity(
            self.record(tlen=tlen, query="absent"), QUERY_LENGTHS) == 0.0

    def test_more_mismatches_than_bases_floors_at_zero(self):
        assert allele_calls.record_identity(
            self.record(nm=50), QUERY_LENGTHS) == 0.0


class TestCallsFromAlignment:
    def calls(self, lines, threshold=0.90):
        return allele_calls.calls_from_alignment(
            list(allele_calls.parse_sam_lines(lines)), QUERY_LENGTHS, threshold)

    def test_no_records_yields_no_calls(self):
        assert self.calls(["@HD\tVN:1.6"]) == {}

    def test_exact_match_is_called_and_confident(self):
        assert self.calls([sam_line("q1", "adk_1")]) == {
            "adk": {"allele_id": "1", "identity": 1.0, "confidence": True}}

    def test_identity_exactly_at_the_threshold_is_confident(self):
        assert self.calls([sam_line("q1", "adk_1", nm="4")])["adk"] == {
            "allele_id": "1", "identity": pytest.approx(0.9), "confidence": True}

    def test_identity_just_below_the_threshold_is_not_confident(self):
        call = self.calls([sam_line("q1", "adk_1", nm="5")])["adk"]
        assert call["identity"] == pytest.approx(0.875)
        assert call["confidence"] is False

    def test_best_alignment_of_a_query_wins(self):
        assert self.calls([
            sam_line("q1", "adk_2", nm="4"),
            sam_line("q1", "adk_1", nm="0"),
        ])["adk"]["allele_id"] == "1"

    def test_tie_between_alignments_of_a_query_takes_the_first_record(self):
        assert self.calls([
            sam_line("q1", "adk_2", nm="4"),
            sam_line("q1", "adk_1", nm="4"),
        ])["adk"]["allele_id"] == "2"

    def test_best_query_of_a_locus_wins(self):
        calls = self.calls([
            sam_line("q1", "adk_1", nm="4"),
            sam_line("q2", "adk_2", nm="0"),
        ])
        assert calls == {
            "adk": {"allele_id": "2", "identity": 1.0, "confidence": True}}

    def test_loci_are_independent(self):
        calls = self.calls([
            sam_line("q1", "adk_1", nm="0"),
            sam_line("q2", "gyrB_3", nm="5"),
        ])
        assert calls["adk"]["allele_id"] == "1"
        assert calls["gyrB"] == {
            "allele_id": "3", "identity": pytest.approx(0.875),
            "confidence": False}

    def test_zero_identity_locus_names_the_allele_it_matched(self):
        # The WDL task seeded every locus at identity 0.0 and only replaced it
        # on a strict improvement, so a best hit of exactly 0.0 lost its
        # allele_id. Reducing per locus directly keeps the evidence: the locus
        # is reported with the allele it aligned to and a failing confidence.
        assert self.calls([sam_line("q1", "zero_9", nm="40")]) == {
            "zero": {"allele_id": "9", "identity": 0.0, "confidence": False}}

    def test_identity_is_clamped_into_the_unit_interval(self):
        call = self.calls([sam_line("q1", "adk_1", nm="-5")])["adk"]
        assert call["identity"] == 1.0
        assert call["confidence"] is True

    def test_reference_without_an_allele_id_is_called_unknown(self):
        assert self.calls([sam_line("q1", "orphan")]) == {
            "orphan": {"allele_id": "unknown", "identity": 1.0,
                       "confidence": True}}

    def test_threshold_is_caller_supplied(self):
        lines = [sam_line("q1", "adk_1", nm="5")]
        assert self.calls(lines, threshold=0.875)["adk"]["confidence"] is True
        assert self.calls(lines, threshold=0.876)["adk"]["confidence"] is False

    def test_default_threshold_matches_the_wdl_input(self):
        assert allele_calls.DEFAULT_IDENTITY_THRESHOLD == 0.90
        assert allele_calls.calls_from_alignment(
            list(allele_calls.parse_sam_lines([sam_line("q1", "adk_1", nm="4")])),
            QUERY_LENGTHS)["adk"]["confidence"] is True

    def test_sam_path_and_parsed_records_give_the_same_calls(self, tmp_path):
        lines = [
            "@HD\tVN:1.6",
            sam_line("q1", "adk_1", nm="0"),
            sam_line("q2", "gyrB_1", nm="5"),
        ]
        sam = tmp_path / "alignment.sam"
        sam.write_text("\n".join(lines) + "\n")
        assert allele_calls.calls_from_alignment(
            str(sam), QUERY_LENGTHS) == self.calls(lines)

    def test_missing_sam_file_yields_no_calls(self, tmp_path):
        assert allele_calls.calls_from_alignment(
            str(tmp_path / "absent.sam"), QUERY_LENGTHS) == {}


class TestCallsAreJsonSerialisable:
    """The caller writes these straight to JSON; nothing exotic may leak in."""

    def test_similarity_calls_round_trip(self):
        calls = allele_calls.call_alleles(
            matrix(("q1", "adk_1"), ("1.0", "0.95"), ("0.95", "1.0")),
            [("q1", "ACGT")], [("adk_1", "A")])
        assert json.loads(json.dumps(calls)) == calls

    def test_alignment_calls_round_trip(self):
        calls = allele_calls.calls_from_alignment(
            list(allele_calls.parse_sam_lines([sam_line("q1", "adk_1")])),
            QUERY_LENGTHS)
        assert json.loads(json.dumps(calls)) == calls


class TestAlignedLength:
    """Identity is measured over the aligned block, not the whole query."""

    @pytest.mark.parametrize("cigar,expected", [
        ("40M", 40),
        ("500M", 500),
        ("1900S500M300S", 500),          # a gene inside a contig
        ("300H500M1900H", 500),          # hard clipping, same answer
        ("10M2I10M", 22),                # insertions consume columns
        ("10M2D10M", 22),                # so do deletions
        ("100=5X100=", 205),             # explicit match/mismatch ops
        ("", 0),
        ("*", 0),
    ])
    def test_counts_only_alignment_columns(self, cigar, expected):
        assert allele_calls.aligned_length(cigar) == expected


class TestIdentityFromCigar:
    def test_gene_in_a_contig_scores_over_the_gene(self):
        """A 500 bp locus inside a 2.7 kb contig with 15 mismatches is 97%.

        The WDL task this replaces divided by the contig length instead, which
        reported 99.4% — high enough to accept a wrong allele as confident.
        """
        record, = allele_calls.parse_sam_lines([
            "\t".join(["contig_1", "0", "adk_2", "1", "60", "1900S500M300S",
                       "*", "0", "0", "A" * 500, "*", "NM:i:15"])
        ])
        identity = allele_calls.record_identity(record, {"contig_1": "A" * 2700})
        assert identity == pytest.approx(485 / 500)

    def test_falls_back_to_query_length_without_a_cigar(self):
        record = {"query_name": "q1", "ref_name": "adk_1", "mapq": 60,
                  "nm": 4, "tlen": 0}
        assert allele_calls.record_identity(record, QUERY_LENGTHS) == pytest.approx(0.9)

    def test_falls_back_to_tlen_before_query_length(self):
        record = {"query_name": "q1", "ref_name": "adk_1", "mapq": 60,
                  "nm": 0, "tlen": 20, "cigar": "*"}
        assert allele_calls.record_identity(record, QUERY_LENGTHS) == 1.0


class TestQueryNameKeying:
    def test_query_dict_is_keyed_like_sam_qname(self, tmp_path):
        """SAM truncates a query name at whitespace; the lookup must match.

        Keying on the full header meant no assembly with a description in its
        header could ever be called from an alignment.
        """
        fasta = tmp_path / "query.fasta"
        fasta.write_text(">contig_1 length=2700 cov=42.0\nACGTACGT\n")
        assert allele_calls.parse_fasta_dict(fasta) == {"contig_1": "ACGTACGT"}

    def test_alignment_call_survives_a_described_header(self, tmp_path):
        fasta = tmp_path / "query.fasta"
        fasta.write_text(">contig_1 assembled by spades\n" + "A" * 40 + "\n")
        calls = allele_calls.calls_from_alignment(
            list(allele_calls.parse_sam_lines([sam_line("contig_1", "adk_3")])),
            allele_calls.parse_fasta_dict(fasta),
        )
        assert calls["adk"]["allele_id"] == "3"
        assert calls["adk"]["identity"] == 1.0


class TestAssemblyQueryCallsEveryLocus:
    """One contig covering many loci yields a call per locus.

    The WDL task reduced per query sequence first, so an assembly — one contig
    aligning to every locus in the scheme — could produce at most one call.
    """

    def _contig_sam(self):
        # minimap2 emits one primary and N-1 supplementary records for a contig
        # that carries several loci.
        return [
            "\t".join(["contig_1", "0", "gyrB_1", "1", "60", "1900S500M300S",
                       "*", "0", "0", "A" * 500, "*", "NM:i:0"]),
            "\t".join(["contig_1", "2048", "adk_2", "1", "60", "300H500M1900H",
                       "*", "0", "0", "A" * 500, "*", "NM:i:0"]),
            "\t".join(["contig_1", "2048", "fumC_1", "1", "60", "1100H500M1100H",
                       "*", "0", "0", "A" * 500, "*", "NM:i:0"]),
        ]

    def test_every_locus_in_the_contig_is_called(self):
        calls = allele_calls.calls_from_alignment(
            list(allele_calls.parse_sam_lines(self._contig_sam())),
            {"contig_1": "A" * 2700},
        )
        assert {locus: call["allele_id"] for locus, call in calls.items()} == {
            "adk": "2", "fumC": "1", "gyrB": "1",
        }
        assert all(call["confidence"] for call in calls.values())

    def test_best_allele_wins_within_a_locus(self):
        sam = self._contig_sam() + [
            # A worse adk alignment must not displace adk_2.
            "\t".join(["contig_1", "2048", "adk_3", "1", "60", "300H500M1900H",
                       "*", "0", "0", "A" * 500, "*", "NM:i:30"]),
        ]
        calls = allele_calls.calls_from_alignment(
            list(allele_calls.parse_sam_lines(sam)), {"contig_1": "A" * 2700}
        )
        assert calls["adk"]["allele_id"] == "2"
        assert calls["adk"]["identity"] == 1.0
