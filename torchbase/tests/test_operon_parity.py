"""Golden parity test: torchbase.operon against StxTyper's own output.

The fixtures are a real run's data, captured from `tblastn` output on
StxTyper's `test/basic.fa` against the reference set of the `ncbi/stxtyper`
1.0.45 torch (built by `torchtools convert stxtyper`):

    basic_hsps.json.gz  normalized HSPs from tasks/protein_search.wdl
    operon_config.json  the torch's [operon] config block
    profiles.tsv        the torch's subtype manifest
    basic.expected      StxTyper 1.0.45's own output for the same assembly

`basic.fa` is StxTyper's edge-case suite: a clean stx1a operon, a residue-
resolved stx2c, a novel stx2, a frameshift, an internal stop, a truncated
operon, and one truncated by a contig end. Asserting exact agreement on every
reported field is the Phase 0 parity criterion (docs/operon-strategy-plan.md
§7) as a test that needs neither BLAST nor a network.

Parity is pinned to StxTyper 1.0.45 (§9 risk 3: StxTyper is a moving target).
"""

import gzip
import json
from pathlib import Path

import pytest

from torchbase.operon import HSP, report_operons, subtype_prefix

FIXTURES = Path(__file__).parent / "fixtures" / "stxtyper"


def _pct(value):
    return "NA" if value is None else "%.2f" % (value * 100.0)


@pytest.fixture(scope="module")
def calls():
    with gzip.open(FIXTURES / "basic_hsps.json.gz", "rt") as f:
        hsps = [HSP.from_dict(raw) for raw in json.load(f)]
    config = json.loads((FIXTURES / "operon_config.json").read_text())
    rows = _read_tsv(FIXTURES / "profiles.tsv")
    return report_operons(hsps, config, scheme="stx", subtype_prefix=subtype_prefix(rows))


@pytest.fixture(scope="module")
def expected():
    rows = _read_tsv(FIXTURES / "basic.expected", strip_hash=True)
    return {row["target_contig"]: row for row in rows}


def _read_tsv(path, strip_hash=False):
    lines = path.read_text().splitlines()
    header = lines[0].lstrip("#" if strip_hash else "").split("\t")
    return [dict(zip(header, line.split("\t"))) for line in lines[1:] if line.strip()]


def _as_stxtyper_row(call):
    operon = call["operon"]
    a = operon["subunits"].get("A")
    b = operon["subunits"].get("B")
    return {
        "target_contig": operon["contig"],
        "stx_type": call["profile_id"],
        "operon": call["operon_status"],
        "identity": _pct(operon["combined_identity"]) if a and b else "NA",
        "target_start": str(operon["start"]),
        "target_stop": str(operon["stop"]),
        "target_strand": operon["strand"],
        "A_reference": a["reference"] if a else "NA",
        "A_reference_subtype": a["reference_subtype"] if a else "NA",
        "A_identity": _pct(a["identity"]) if a else "NA",
        "A_coverage": _pct(a["coverage"]) if a else "NA",
        "B_reference": b["reference"] if b else "NA",
        "B_reference_subtype": b["reference_subtype"] if b else "NA",
        "B_identity": _pct(b["identity"]) if b else "NA",
        "B_coverage": _pct(b["coverage"]) if b else "NA",
    }


CONTIGS = [
    "stx1a",                # clean COMPLETE operon
    "stx2c",                # generalized class resolved by residue table
    "stx2_novel",           # intact operon, unresolved subtype
    "stx2_fs",              # frameshift split across two HSPs
    "stxB2b_stop",          # internal stop codon
    "partial",              # truncated subunit, not at a contig end
    "partial_contig_end",   # truncated by the end of the contig
]


class TestStxTyperParity:
    def test_one_call_per_expected_contig(self, calls, expected):
        assert sorted(c["operon"]["contig"] for c in calls) == sorted(expected)

    @pytest.mark.parametrize("contig", CONTIGS)
    def test_call_matches_stxtyper(self, calls, expected, contig):
        ours = [c for c in calls if c["operon"]["contig"] == contig]
        assert len(ours) == 1, "expected exactly one call for %s" % contig
        assert _as_stxtyper_row(ours[0]) == expected[contig]
