"""Unit tests for torchbase.operon: the operon-typing reference algorithm.

Covers §8 of docs/operon-strategy-plan.md:
- threshold selection (explicit + default fallback)
- residue-table resolution, including the fallback path
- synteny pairing: correct pair accepted; gap too large rejected; wrong
  strand rejected; wrong order rejected
- status-ladder priority
"""

import pytest
from torchbase.operon import (
    HSP,
    OperonConfigError,
    build_residue_rules,
    call_operon,
    call_operon_status,
    combined_identity,
    pair_operon,
    resolve_generalized_class,
    select_threshold,
    status_to_coarse,
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
    defaults = dict(
        subunit=subunit,
        contig=contig,
        strand=strand,
        start=start,
        stop=stop,
        ref_subtype=f"stx{subunit}1a",
        ref_class=ref_class,
        nident=stop - start,
        aln_len=stop - start,
        ref_len=stop - start,
    )
    defaults.update(kw)
    return HSP(**defaults)


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

    def test_each_hsp_claimed_once(self):
        # Two A candidates competing for one B: only one candidate survives,
        # and only the winning pair claims the B HSP.
        a1 = _hsp("A", start=0, stop=900, nident=900)
        a2 = _hsp("A", start=0, stop=900, nident=850)
        b = _hsp("B", start=920, stop=1200)
        candidates = pair_operon(
            {"A": [a1, a2], "B": [b]}, ["A", "B"], intergenic_max=36
        )
        assert len(candidates) == 1
        assert candidates[0][0] is a1  # higher-identity seed wins

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
    def test_complete_call(self):
        candidate = {
            "A": _hsp("A", ref_class="1a", nident=300, aln_len=300, ref_len=300),
            "B": _hsp("B", ref_class="1a", nident=100, aln_len=100, ref_len=100),
        }
        cfg = {"subunit_order": ["A", "B"], "identity_thresholds": STX_THRESHOLDS}
        result = call_operon(candidate, cfg)
        assert result["operon_status"] == "COMPLETE"
        assert result["status"] == "known"
        assert result["class"] == "1a"
        assert result["combined_identity"] == 1.0

    def test_class_mismatch_is_novel(self):
        candidate = {
            "A": _hsp("A", ref_class="1a", nident=300, aln_len=300, ref_len=300),
            "B": _hsp("B", ref_class="2", nident=100, aln_len=100, ref_len=100),
        }
        cfg = {"subunit_order": ["A", "B"], "identity_thresholds": STX_THRESHOLDS}
        result = call_operon(candidate, cfg)
        assert result["operon_status"] == "COMPLETE_NOVEL"
        assert result["class"] is None

    def test_below_threshold_is_novel(self):
        candidate = {
            "A": _hsp("A", ref_class="1a", nident=950, aln_len=1000, ref_len=1000),
            "B": _hsp("B", ref_class="1a", nident=950, aln_len=1000, ref_len=1000),
        }
        cfg = {"subunit_order": ["A", "B"], "identity_thresholds": STX_THRESHOLDS}
        result = call_operon(candidate, cfg)
        assert result["combined_identity"] == 0.95
        assert result["operon_status"] == "COMPLETE_NOVEL"


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
