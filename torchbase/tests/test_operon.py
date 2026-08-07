"""Unit tests for torchbase.operon: the operon-typing algorithm.

Covers §8 of docs/operon-strategy-plan.md:
- threshold selection (explicit + default fallback)
- residue-table resolution, including the fallback path
- BLAST normalization, including the reference's terminal stop codon
- frameshift stitching, locus reduction, synteny pairing, and the
  containment-based selection that decides what gets reported
- status-ladder priority

`test_operon_parity.py` covers the same algorithm end to end against
StxTyper's own golden output.
"""

import pytest
from torchbase.operon import (
    HSP,
    Candidate,
    OperonConfigError,
    assemble_candidates,
    build_residue_rules,
    call_operon,
    call_operon_status,
    combined_identity,
    normalize_hsps,
    pair_operon,
    parse_fasta,
    reduce_to_loci,
    resolve_generalized_class,
    select_operons,
    select_threshold,
    status_to_coarse,
    stitch_frameshifts,
    subtype_prefix,
    superclass_of,
    validate_operon_metadata,
)


STX_THRESHOLDS = {
    "default": 0.98,
    "1a": 0.983,
    "1c": 0.983,
    "1d": 0.983,
    "1e": 0.983,
    "2k": 0.985,
    "2l": 0.985,
}

STX2ACD_RULE = {
    "class": "2",
    "positions": [
        {"subunit": "A", "index": 312},
        {"subunit": "A", "index": 318},
        {"subunit": "B", "index": 34},
    ],
    "table": [
        {"call": "2a", "residues": [["F", "S"], ["K", "E"], ["D"]]},
        {"call": "2c", "residues": [["F"], ["K", "E"], ["N"]]},
        {"call": "2d", "residues": [["S"], ["E"], ["N"]]},
    ],
    "fallback": "2",
}


class TestSelectThreshold:
    def test_explicit_class_threshold(self):
        assert select_threshold("1a", STX_THRESHOLDS) == 0.983

    def test_default_fallback(self):
        assert select_threshold("2c", STX_THRESHOLDS) == 0.98

    def test_no_default_raises(self):
        with pytest.raises(OperonConfigError):
            select_threshold("9z", {"1a": 0.983})


class TestCombinedIdentity:
    def test_combines_across_subunits(self):
        ident = combined_identity({"A": (300, 300), "B": (80, 100)})
        assert ident == pytest.approx(380 / 400)

    def test_zero_length_is_zero(self):
        assert combined_identity({}) == 0.0


class TestResidueRuleResolution:
    @pytest.fixture
    def rule(self):
        return build_residue_rules([STX2ACD_RULE])["2"]

    @pytest.mark.parametrize(
        "residues,expected",
        [
            (["F", "K", "D"], "2a"),
            (["S", "E", "D"], "2a"),
            (["F", "K", "N"], "2c"),
            (["S", "E", "N"], "2d"),
        ],
    )
    def test_every_table_row(self, rule, residues, expected):
        assert rule.resolve(residues) == expected

    def test_unresolved_falls_back(self, rule):
        assert rule.resolve(["Y", "K", "D"]) == "2"

    def test_resolve_generalized_class_passthrough_for_non_generalized(self, rule):
        rules = {"2": rule}
        assert resolve_generalized_class("1a", ["F"], {"2": ["2a", "2c", "2d"]}, rules) == "1a"

    def test_resolve_generalized_class_dispatches(self, rule):
        rules = {"2": rule}
        result = resolve_generalized_class(
            "2", ["F", "K", "D"], {"2": ["2a", "2c", "2d"]}, rules
        )
        assert result == "2a"

    def test_missing_rule_raises(self):
        with pytest.raises(OperonConfigError):
            resolve_generalized_class("2", ["F"], {"2": ["2a"]}, {})


def _hsp(subunit, contig="NODE_1", strand="+", start=0, stop=100, ref_class="1a", **kw):
    """An HSP that covers its whole reference, unless a test says otherwise."""
    ref_len = kw.pop("ref_len", stop - start)
    defaults = dict(
        subunit=subunit,
        contig=contig,
        strand=strand,
        start=start,
        stop=stop,
        ref_subtype="stx{}1a".format(subunit),
        ref_class=ref_class,
        nident=stop - start,
        aln_len=stop - start,
        ref_len=ref_len,
        reference_accession="ACC{}".format(subunit),
        contig_len=10000,
        qstart=0,
        qend=ref_len - 1,
        stop_codon=True,
    )
    defaults.update(kw)
    return HSP(**defaults)


def _candidate(**subunits):
    return Candidate(subunits=subunits, complete=True)


class TestPairOperon:
    def test_correct_pair_accepted(self):
        a = _hsp("A", start=0, stop=900)
        b = _hsp("B", start=920, stop=1200)  # gap = 920 - 900 - 1 = 19bp
        candidates = pair_operon({"A": [a], "B": [b]}, ["A", "B"], intergenic_max=36)
        assert len(candidates) == 1
        assert candidates[0] == [a, b]

    def test_gap_too_large_rejected(self):
        a = _hsp("A", start=0, stop=900)
        b = _hsp("B", start=940, stop=1200)  # gap = 940 - 900 - 1 = 39bp > 36
        candidates = pair_operon({"A": [a], "B": [b]}, ["A", "B"], intergenic_max=36)
        assert candidates == []

    def test_wrong_strand_rejected(self):
        a = _hsp("A", start=0, stop=900, strand="+")
        b = _hsp("B", start=920, stop=1200, strand="-")
        candidates = pair_operon(
            {"A": [a], "B": [b]}, ["A", "B"], intergenic_max=36, require_same_strand=True
        )
        assert candidates == []

    def test_wrong_order_rejected(self):
        # B lies upstream of A on the + strand -> not a valid A-then-B chain.
        a = _hsp("A", start=920, stop=1200)
        b = _hsp("B", start=0, stop=900)
        candidates = pair_operon({"A": [a], "B": [b]}, ["A", "B"], intergenic_max=36)
        assert candidates == []

    def test_wrong_contig_rejected(self):
        a = _hsp("A", contig="NODE_1", start=0, stop=900)
        b = _hsp("B", contig="NODE_2", start=920, stop=1200)
        candidates = pair_operon(
            {"A": [a], "B": [b]}, ["A", "B"], intergenic_max=36, require_same_contig=True
        )
        assert candidates == []

    def test_competing_pairings_are_all_enumerated_best_first(self):
        # Two A alignments over the same locus, one better: both pairings are
        # offered — which one is reported is select_operons' call — but the
        # better-scoring alignment leads.
        a1 = _hsp("A", start=0, stop=900, nident=900, ref_class="1a")
        a2 = _hsp("A", start=0, stop=900, nident=850, ref_class="1c")
        b = _hsp("B", start=920, stop=1200)
        candidates = pair_operon(
            {"A": [a1, a2], "B": [b]}, ["A", "B"], intergenic_max=36
        )
        assert [chain[0] for chain in candidates] == [a1, a2]
        assert all(chain[1] is b for chain in candidates)

    def test_overlapping_genes_rejected_by_default(self):
        # ETEC's eltA/eltB overlap; stx's stxA/stxB never do, so the default
        # floor of 0 must reject an overlap.
        a = _hsp("A", start=1, stop=900)
        b = _hsp("B", start=897, stop=1200)
        assert pair_operon({"A": [a], "B": [b]}, ["A", "B"], intergenic_max=36) == []

    def test_overlapping_genes_accepted_when_scheme_allows(self):
        a = _hsp("A", start=1, stop=900)
        b = _hsp("B", start=897, stop=1200)
        candidates = pair_operon(
            {"A": [a], "B": [b]}, ["A", "B"], intergenic_max=36, intergenic_min=-12
        )
        assert candidates == [[a, b]]

    def test_overlap_beyond_the_floor_still_rejected(self):
        a = _hsp("A", start=1, stop=900)
        b = _hsp("B", start=880, stop=1200)  # gap = -21, below -12
        assert pair_operon(
            {"A": [a], "B": [b]}, ["A", "B"], intergenic_max=36, intergenic_min=-12
        ) == []

    def test_nested_alignment_is_not_a_downstream_partner(self):
        # B inside A is not "A then B" however small the coordinate gap looks.
        a = _hsp("A", start=1, stop=1200)
        b = _hsp("B", start=100, stop=400)
        assert pair_operon(
            {"A": [a], "B": [b]}, ["A", "B"], intergenic_max=36, intergenic_min=-2000
        ) == []

    def test_requires_at_least_two_subunits(self):
        with pytest.raises(OperonConfigError):
            pair_operon({"A": []}, ["A"], intergenic_max=36)


class TestStatusLadder:
    def test_complete(self):
        assert call_operon_status(False, False, False, False, False, False, False, False) == "COMPLETE"

    def test_frameshift_dominates_everything(self):
        assert (
            call_operon_status(
                frameshift=True,
                internal_stop=True,
                contig_end=True,
                below_threshold=True,
                class_mismatch=True,
                ambiguous=True,
                extended=True,
                partial=True,
            )
            == "FRAMESHIFT"
        )

    def test_partial_contig_end_beats_partial(self):
        assert (
            call_operon_status(False, False, True, False, False, False, False, True)
            == "PARTIAL_CONTIG_END"
        )

    def test_below_threshold_ambiguous_vs_novel(self):
        assert (
            call_operon_status(False, False, False, True, False, True, False, False)
            == "AMBIGUOUS"
        )
        assert (
            call_operon_status(False, False, False, True, False, False, False, False)
            == "COMPLETE_NOVEL"
        )


class TestStatusToCoarse:
    @pytest.mark.parametrize(
        "operon_status,expected",
        [
            ("COMPLETE", "known"),
            ("COMPLETE_NOVEL", "novel_profile"),
            ("AMBIGUOUS", "novel_profile"),
            ("PARTIAL", "incomplete"),
            ("FRAMESHIFT", "incomplete"),
        ],
    )
    def test_mapping(self, operon_status, expected):
        assert status_to_coarse(operon_status) == expected


class TestCallOperon:
    CFG = {
        "subunit_order": ["A", "B"],
        "identity_thresholds": STX_THRESHOLDS,
        "generalized_classes": {"2": ["2a", "2c", "2d"]},
        "residue_rules": [STX2ACD_RULE],
        "superclass_pattern": r"^([0-9]+)",
    }

    def test_complete_call(self):
        candidate = _candidate(
            A=_hsp("A", ref_class="1a", start=1, stop=300),
            B=_hsp("B", ref_class="1a", start=320, stop=420),
        )
        result = call_operon(candidate, self.CFG)
        assert result["operon_status"] == "COMPLETE"
        assert result["status"] == "known"
        assert result["class"] == "1a"
        assert result["type_token"] == "1a"
        assert result["combined_identity"] == 1.0

    def test_class_mismatch_is_novel_at_parent_resolution(self):
        candidate = _candidate(
            A=_hsp("A", ref_class="1a", start=1, stop=300),
            B=_hsp("B", ref_class="1c", start=320, stop=420),
        )
        result = call_operon(candidate, self.CFG)
        assert result["operon_status"] == "COMPLETE_NOVEL"
        assert result["class"] is None
        assert result["type_token"] == "1"

    def test_unrelated_classes_support_no_type_claim(self):
        candidate = _candidate(
            A=_hsp("A", ref_class="1a", start=1, stop=300),
            B=_hsp("B", ref_class="2", start=320, stop=420),
        )
        result = call_operon(candidate, self.CFG)
        assert result["operon_status"] == "COMPLETE_NOVEL"
        assert result["type_token"] is None

    def test_below_threshold_is_novel(self):
        candidate = _candidate(
            A=_hsp("A", ref_class="1a", start=1, stop=1000, nident=950,
                   aln_len=1000, ref_len=1000, qend=999),
            B=_hsp("B", ref_class="1a", start=1020, stop=2020, nident=950,
                   aln_len=1000, ref_len=1000, qend=999),
        )
        result = call_operon(candidate, self.CFG)
        assert result["combined_identity"] == 0.95
        assert result["operon_status"] == "COMPLETE_NOVEL"

    def test_generalized_class_resolved_by_residues(self):
        # stx2c: A312=F, A318=K, B34=N. The subject residue at a reference
        # offset is read out of the gapped alignment, so build one long enough
        # to contain the positions the rule asks for.
        a = _hsp(
            "A", ref_class="2", start=1, stop=960, ref_len=319,
            qseq="X" * 319, sseq="X" * 312 + "F" + "X" * 5 + "K" + "X" * 5,
        )
        b = _hsp(
            "B", ref_class="2", start=980, stop=1247, ref_len=89,
            qseq="X" * 89, sseq="X" * 34 + "N" + "X" * 54,
        )
        result = call_operon(_candidate(A=a, B=b), self.CFG)
        assert result["type_token"] == "2c"
        assert result["operon_status"] == "COMPLETE"
        assert result["residue_evidence"] == {"A312": "F", "A318": "K", "B34": "N"}

    def test_unresolved_generalized_class_is_novel(self):
        a = _hsp("A", ref_class="2", start=1, stop=960, ref_len=319,
                 qseq="X" * 319, sseq="Q" * 319)
        b = _hsp("B", ref_class="2", start=980, stop=1247, ref_len=89,
                 qseq="X" * 89, sseq="Q" * 89)
        result = call_operon(_candidate(A=a, B=b), self.CFG)
        assert result["operon_status"] == "COMPLETE_NOVEL"
        assert result["type_token"] == "2"

    def test_incomplete_reference_coverage_is_partial(self):
        candidate = _candidate(
            A=_hsp("A", start=1, stop=300, ref_len=400, qend=299),
            B=_hsp("B", start=320, stop=420),
        )
        assert call_operon(candidate, self.CFG)["operon_status"] == "PARTIAL"

    def test_missing_stop_codon_is_extended(self):
        candidate = _candidate(
            A=_hsp("A", start=1, stop=300, stop_codon=False),
            B=_hsp("B", start=320, stop=420),
        )
        assert call_operon(candidate, self.CFG)["operon_status"] == "EXTENDED"

    def test_contig_end_outranks_extension(self):
        candidate = _candidate(
            A=_hsp("A", start=1, stop=300, ref_len=400, qstart=100, qend=399,
                   stop_codon=False),
            B=_hsp("B", start=320, stop=420),
        )
        assert (
            call_operon(candidate, self.CFG)["operon_status"] == "PARTIAL_CONTIG_END"
        )


class TestNormalizeHsps:
    HEADER_FORMAT = "accession|subunit_role|reference_subtype|class"

    def _row(self, **over):
        row = {
            "qseqid": "ACC1|A|stxA1a|1a", "sseqid": "NODE_1", "pident": "100.0",
            "length": "5", "mismatch": "0", "gapopen": "0", "qstart": "1",
            "qend": "5", "sstart": "1", "send": "15", "evalue": "0",
            "bitscore": "10", "nident": "5", "qlen": "5", "slen": "1000",
            "sstrand": "plus", "qseq": "MKKT*", "sseq": "MKKT*",
        }
        row.update(over)
        return [row[col] for col in (
            "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
            "qstart", "qend", "sstart", "send", "evalue", "bitscore", "nident",
            "qlen", "slen", "sstrand", "qseq", "sseq",
        )]

    REFS = {"ACC1|A|stxA1a|1a": "MKKT*"}

    def test_terminal_stop_excluded_from_protein_and_span(self):
        hsp, = normalize_hsps([self._row()], self.REFS, self.HEADER_FORMAT)
        assert hsp.ref_len == 4          # the stop codon is not a residue
        assert hsp.nident == 4
        assert hsp.aln_len == 4
        assert (hsp.start, hsp.stop) == (1, 12)  # stop codon's 3 nt dropped
        assert hsp.stop_codon is True
        assert hsp.reference_complete
        assert not hsp.internal_stop
        assert hsp.coverage == 1.0

    def test_terminal_stop_span_trimmed_upstream_on_minus_strand(self):
        hsp, = normalize_hsps(
            [self._row(sstart="15", send="1", sstrand="minus")],
            self.REFS, self.HEADER_FORMAT,
        )
        assert (hsp.start, hsp.stop) == (4, 15)

    def test_internal_stop_detected(self):
        hsp, = normalize_hsps(
            [self._row(qseq="MKKT*", sseq="MK*T*", nident="4")],
            self.REFS, self.HEADER_FORMAT,
        )
        assert hsp.internal_stop

    def test_full_protein_without_subject_stop_is_extended(self):
        hsp, = normalize_hsps(
            [self._row(qend="4", length="4", nident="4", qseq="MKKT", sseq="MKKT",
                       send="12")],
            self.REFS, self.HEADER_FORMAT,
        )
        assert hsp.stop_codon is False
        assert hsp.extended
        assert not hsp.reference_complete

    def test_header_fields_split_by_declared_format(self):
        hsp, = normalize_hsps([self._row()], self.REFS, self.HEADER_FORMAT)
        assert hsp.reference_accession == "ACC1"
        assert hsp.subunit == "A"
        assert hsp.ref_subtype == "stxA1a"
        assert hsp.ref_class == "1a"


class TestParseFasta:
    def test_headers_truncated_at_whitespace(self):
        seqs = parse_fasta(">a|b desc here\nMKK\nTL\n>c|d\nQQ\n")
        assert seqs == {"a|b": "MKKTL", "c|d": "QQ"}


class TestStitchFrameshifts:
    def _split_pair(self, second_start=902):
        """Two HSPs of one reference, split mid-protein by a frame change."""
        a = _hsp("A", start=1, stop=900, ref_len=400, qstart=0, qend=299,
                 nident=300, aln_len=300, qseq="Q" * 300, sseq="Q" * 300,
                 stop_codon=None)
        b = _hsp("A", start=second_start, stop=second_start + 299,
                 ref_len=400, qstart=300, qend=399, nident=100, aln_len=100,
                 qseq="Q" * 100, sseq="Q" * 100)
        return a, b

    def test_frame_change_is_stitched_into_one_alignment(self):
        a, b = self._split_pair()
        assert a.frame() != b.frame()
        stitched = [h for h in stitch_frameshifts([a, b]) if h.frameshift]
        assert len(stitched) == 1
        merged = stitched[0]
        assert (merged.qstart, merged.qend) == (0, 399)
        assert merged.reference_complete
        # One column is dropped at the join: the codon spanning the shift.
        assert merged.aln_len == 399

    def test_parts_are_kept_alongside_the_stitched_alignment(self):
        a, b = self._split_pair()
        assert set(stitch_frameshifts([a, b])) >= {a, b}

    def test_same_frame_blocks_are_not_stitched(self):
        a, b = self._split_pair(second_start=901)  # in frame with a
        assert a.frame() == b.frame()
        assert not any(h.frameshift for h in stitch_frameshifts([a, b]))


class TestReduceToLoci:
    def test_best_alignment_per_locus_per_class(self):
        best = _hsp("A", start=1, stop=900, nident=900, ref_class="1a")
        worse = _hsp("A", start=1, stop=900, nident=800, ref_class="1a",
                     reference_accession="ACCZ")
        other_class = _hsp("A", start=1, stop=900, nident=500, ref_class="2")
        kept = reduce_to_loci([worse, best, other_class])
        assert set(kept) == {best, other_class}

    def test_same_coordinates_on_another_contig_is_another_locus(self):
        one = _hsp("A", contig="NODE_1", start=1, stop=900)
        two = _hsp("A", contig="NODE_2", start=1, stop=900, nident=800)
        assert set(reduce_to_loci([one, two])) == {one, two}


class TestSelectOperons:
    CFG = {"subunit_order": ["A", "B"], "identity_thresholds": STX_THRESHOLDS}

    def test_weak_candidates_are_not_reported(self):
        weak = _candidate(
            A=_hsp("A", start=1, stop=300, nident=100),
            B=_hsp("B", start=320, stop=420, nident=30),
        )
        assert select_operons([weak], self.CFG) == []

    def test_lone_subunit_inside_a_reported_operon_is_redundant(self):
        operon = _candidate(
            A=_hsp("A", start=1, stop=300),
            B=_hsp("B", start=320, stop=420),
        )
        lone = Candidate(subunits={"A": _hsp("A", start=1, stop=300)}, complete=False)
        assert select_operons([operon, lone], self.CFG) == [operon]

    def test_operons_on_different_contigs_are_both_reported(self):
        first = _candidate(
            A=_hsp("A", contig="NODE_1", start=1, stop=300),
            B=_hsp("B", contig="NODE_1", start=320, stop=420),
        )
        second = _candidate(
            A=_hsp("A", contig="NODE_2", start=1, stop=300),
            B=_hsp("B", contig="NODE_2", start=320, stop=420),
        )
        assert len(select_operons([first, second], self.CFG)) == 2


class TestAssembleCandidates:
    CFG = {
        "subunit_order": ["A", "B"],
        "identity_thresholds": STX_THRESHOLDS,
        "intergenic_max": 36,
        "intergenic_relax_factor": 2,
    }

    def test_reference_set_collapses_to_one_candidate_per_operon(self):
        # Ten references hitting the same two loci must not produce 100
        # candidate operons.
        hsps = []
        for i in range(10):
            hsps.append(_hsp("A", start=1, stop=900, nident=900 - i,
                             reference_accession="A%d" % i))
            hsps.append(_hsp("B", start=920, stop=1200, nident=280 - i,
                             reference_accession="B%d" % i))
        candidates = assemble_candidates(hsps, self.CFG)
        assert len(candidates) == 1
        assert candidates[0].subunits["A"].reference_accession == "A0"

    def test_lone_subunit_still_becomes_a_candidate(self):
        candidates = assemble_candidates([_hsp("A", start=1, stop=900)], self.CFG)
        assert len(candidates) == 1
        assert not candidates[0].complete


class TestSuperclassOf:
    def test_parent_class_from_pattern(self):
        cfg = {"superclass_pattern": r"^([0-9]+)"}
        assert superclass_of("2c", cfg) == "2"
        assert superclass_of("1a", cfg) == "1"
        assert superclass_of("2", cfg) == "2"

    def test_without_pattern_a_class_is_its_own_parent(self):
        assert superclass_of("2c", {}) == "2c"


class TestSubtypePrefix:
    def test_longest_common_prefix_of_subtype_labels(self):
        rows = [{"subtype": "stx1a"}, {"subtype": "stx2c"}, {"subtype": "stx2"}]
        assert subtype_prefix(rows) == "stx"

    def test_no_rows_means_no_prefix(self):
        assert subtype_prefix([]) == ""


class TestValidateOperonMetadata:
    def _cfg(self, **overrides):
        cfg = {
            "subunit_order": ["A", "B"],
            "identity_thresholds": dict(STX_THRESHOLDS),
            "generalized_classes": {"2": ["2a", "2c", "2d"]},
            "residue_rules": [STX2ACD_RULE],
        }
        cfg.update(overrides)
        return cfg

    def test_valid_config_passes(self):
        rows = [{"class": "1a"}, {"class": "2"}]
        validate_operon_metadata(self._cfg(), rows)  # no raise

    def test_missing_threshold_for_referenced_class_no_default_raises(self):
        cfg = self._cfg(identity_thresholds={"1a": 0.983})
        rows = [{"class": "9z"}]
        with pytest.raises(OperonConfigError):
            validate_operon_metadata(cfg, rows)

    def test_generalized_class_without_residue_rule_raises(self):
        cfg = self._cfg(residue_rules=[])
        with pytest.raises(OperonConfigError):
            validate_operon_metadata(cfg, [{"class": "1a"}])

    def test_single_subunit_order_raises(self):
        cfg = self._cfg(subunit_order=["A"])
        with pytest.raises(OperonConfigError):
            validate_operon_metadata(cfg, [{"class": "1a"}])
