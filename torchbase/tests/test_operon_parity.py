"""Golden parity test: torchbase.operon against StxTyper's own output.

The fixtures are one real run's data, captured from the `search_subunits` task
run on StxTyper's `test/basic.fa` against the reference set of the
`ncbi/stxtyper` 1.0.45 torch (built by `torchtools convert stxtyper`):

    basic_hits.tsv.gz   raw tblastn tabular output — exactly what the WDL layer
                        hands back to the package layer
    subunits.faa.gz     the torch's protein reference set
    operon_config.json  the torch's [operon] config block
    profiles.tsv        the torch's subtype manifest
    basic.expected      StxTyper 1.0.45's own output for the same assembly

Starting from BLAST's own table rather than pre-normalized HSPs means this
covers the entire package-layer pipeline, including the terminal-stop-codon
handling in `normalize_hsps` where an error costs every COMPLETE call.

`basic.fa` is StxTyper's edge-case suite: a clean stx1a operon, a residue-
resolved stx2c, a novel stx2, a frameshift, an internal stop, a truncated
operon, and one truncated by a contig end. Exact agreement on every reported
field is the Phase 0 parity criterion (docs/operon-strategy-plan.md §7), as a
test needing neither BLAST nor a network.

Parity is pinned to StxTyper 1.0.45 (§9 risk 3: StxTyper is a moving target).
"""

import gzip
import json
from pathlib import Path

import pytest

from torchbase.operon import type_assembly

FIXTURES = Path(__file__).parent / "fixtures" / "stxtyper"


def _pct(value):
    return "NA" if value is None else "%.2f" % (value * 100.0)


def _read_tsv(path, strip_hash=False):
    lines = path.read_text().splitlines()
    header = lines[0].lstrip("#" if strip_hash else "").split("\t")
    return [dict(zip(header, line.split("\t"))) for line in lines[1:] if line.strip()]


def _unpack(name, destination):
    with gzip.open(FIXTURES / name, "rb") as source:
        destination.write_bytes(source.read())
    return destination


@pytest.fixture(scope="module")
def calls(tmp_path_factory):
    """Calls from the raw search output, exactly as `torchbase run` derives them."""
    workdir = tmp_path_factory.mktemp("stx_parity")
    hits = _unpack("basic_hits.tsv.gz", workdir / "hits.tsv")
    reference = _unpack("subunits.faa.gz", workdir / "subunits.faa")
    config = json.loads((FIXTURES / "operon_config.json").read_text())
    profiles = _read_tsv(FIXTURES / "profiles.tsv")
    return type_assembly(hits, reference, config, profile_rows=profiles, scheme="stx")


@pytest.fixture(scope="module")
def expected():
    rows = _read_tsv(FIXTURES / "basic.expected", strip_hash=True)
    return {row["target_contig"]: row for row in rows}


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


class TestEmptySearchOutput:
    def test_no_hits_yields_no_calls(self, tmp_path):
        """An assembly with no operon produces an empty table, not an error."""
        hits = tmp_path / "hits.tsv"
        hits.write_text("")
        reference = _unpack("subunits.faa.gz", tmp_path / "subunits.faa")
        config = json.loads((FIXTURES / "operon_config.json").read_text())
        assert type_assembly(hits, reference, config) == []
