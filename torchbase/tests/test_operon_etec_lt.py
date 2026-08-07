"""Phase 1 of docs/operon-strategy-plan.md: a second scheme, no new algorithm.

`examples/etec_lt` encodes ETEC's heat-labile toxin operon (*eltA*/*eltB*) as
an `[operon]` config by hand. If the same code types it, the config schema
generalizes beyond stx — that is the plan's own gate for Phase 2 (§7).

The fixtures are HSPs from real `tblastn` runs of that torch against:
  h10407          the complete ETEC H10407 genome (GenBank FN649414-FN649418),
                  the reference LT1 strain — the elt operon sits on a plasmid
  lt2_synthetic   a synthetic contig: the real elt operon of LT2-producing
                  strain 2781-5 (GenBank EU113255) with random 500 bp flanks,
                  so the operon is not clipped by a contig end

What this scheme exercises that stx does not: overlapping subunit genes
(`eltA`/`eltB` overlap, hence `intergenic_min` below zero), a single class with
no super-class pattern, and a residue table keyed on five positions across both
subunits.
"""

import gzip
import json
from pathlib import Path

import pytest

from torchbase.operon import HSP, report_operons, subtype_prefix
from torchbase.torchfs import Torch

FIXTURES = Path(__file__).parent / "fixtures" / "etec_lt"
TORCH = Path(__file__).parents[2] / "examples" / "etec_lt" / "1.0.0.torch"


@pytest.fixture(scope="module")
def torch():
    return Torch.load(TORCH)


def calls_for(torch, fixture):
    with gzip.open(FIXTURES / fixture, "rt") as f:
        hsps = [HSP.from_dict(raw) for raw in json.load(f)]
    return report_operons(
        hsps,
        torch.operon_config,
        scheme="etec-lt",
        subtype_prefix=subtype_prefix(torch.operon_profiles),
    )


class TestEtecLtTorch:
    def test_torch_loads_as_an_operon_torch(self, torch):
        assert torch.typing_model == "operon"
        assert torch.operon_config["subunit_order"] == ["A", "B"]
        assert {row["subtype"] for row in torch.operon_profiles} == {"LT1", "LT2"}

    def test_overlapping_genes_are_declared(self, torch):
        # eltA and eltB overlap in the elt operon, so the pairing floor must
        # admit a negative intergenic distance.
        assert torch.operon_config["intergenic_min"] < 0

    def test_residue_positions_are_the_documented_substitutions(self, torch):
        rule, = torch.operon_config["residue_rules"]
        positions = [(p["subunit"], p["index"]) for p in rule["positions"]]
        # Mature-protein A190/A196/A213/A224 and B75 shifted by the 18- and
        # 21-residue signal peptides of the precursor references.
        assert positions == [("A", 207), ("A", 213), ("A", 230), ("A", 241),
                             ("B", 95)]


class TestEtecLtCalls:
    def test_h10407_is_lt1(self, torch):
        call, = calls_for(torch, "h10407_hsps.json.gz")
        assert call["profile_id"] == "LT1"
        assert call["operon_status"] == "COMPLETE"
        assert call["operon"]["residue_evidence"] == {
            "A207": "S", "A213": "G", "A230": "K", "A241": "S", "B95": "T",
        }
        assert call["operon"]["subunits"]["A"]["coverage"] == 1.0
        assert call["operon"]["subunits"]["B"]["coverage"] == 1.0
        # The operon is on a plasmid, not the chromosome, and the subunits
        # overlap rather than sitting apart.
        assert call["operon"]["contig"] == "FN649417.1"
        assert call["operon"]["intergenic_bp"] < 0

    def test_lt2_signature_operon_is_lt2(self, torch):
        call, = calls_for(torch, "lt2_synthetic_hsps.json.gz")
        assert call["profile_id"] == "LT2"
        assert call["operon_status"] == "COMPLETE"
        assert call["operon"]["residue_evidence"] == {
            "A207": "L", "A213": "D", "A230": "E", "A241": "T", "B95": "A",
        }
        assert call["operon"]["combined_identity"] == 1.0
