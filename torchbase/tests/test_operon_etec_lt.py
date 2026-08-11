"""Phase 1 of docs/operon-strategy-plan.md: a second scheme, no new algorithm.

The ETEC heat-labile toxin operon (*eltA*/*eltB*) encoded as an `[operon]`
config by hand. If the same code types it, the config schema generalizes beyond
stx — that is the plan's own gate for Phase 2 (§7).

The torch itself lives in the torches repository, not here: torch config and
reference data are distributable artifacts, separate from the package (see
CLAUDE.md, "Where torches live"). What this test needs is copied in as
fixtures, so it runs without that repository checked out:

    subunits.faa.gz     the torch's protein reference set (GenBank S60731,
                        EU113242-EU113255; Lasaro et al. 2008, PMID 18223074)
    operon_config.json  the torch's [operon] config block
    profiles.tsv        the torch's subtype manifest
    *_hits.tsv.gz       raw tblastn output from real runs against:
                          h10407         the complete ETEC H10407 genome
                                         (FN649414-FN649418), the reference LT1
                                         strain — the operon is on a plasmid
                          lt2_synthetic  the real elt operon of LT2-producing
                                         strain 2781-5 (EU113255) in a padded
                                         synthetic contig

What this scheme exercises that stx does not: overlapping subunit genes
(`eltA`/`eltB` overlap, hence `intergenic_min` below zero), a single class with
no super-class pattern, and a residue table keyed on five positions across both
subunits.
"""

import gzip
import json
from pathlib import Path

import pytest

from torchbase.operon import subtype_prefix, type_assembly

FIXTURES = Path(__file__).parent / "fixtures" / "etec_lt"


def _read_tsv(path):
    lines = path.read_text().splitlines()
    header = lines[0].split("\t")
    return [dict(zip(header, line.split("\t"))) for line in lines[1:] if line.strip()]


@pytest.fixture(scope="module")
def config():
    return json.loads((FIXTURES / "operon_config.json").read_text())


@pytest.fixture(scope="module")
def profiles():
    return _read_tsv(FIXTURES / "profiles.tsv")


def calls_for(fixture, config, profiles, tmp_path):
    hits = tmp_path / "hits.tsv"
    with gzip.open(FIXTURES / fixture, "rb") as source:
        hits.write_bytes(source.read())
    reference = tmp_path / "subunits.faa"
    with gzip.open(FIXTURES / "subunits.faa.gz", "rb") as source:
        reference.write_bytes(source.read())
    return type_assembly(
        hits, reference, config,
        profile_rows=profiles, scheme="etec-lt",
    )


class TestEtecLtScheme:
    def test_subtypes_are_lt1_and_lt2(self, profiles):
        assert {row["subtype"] for row in profiles} == {"LT1", "LT2"}
        assert subtype_prefix(profiles) == "LT"

    def test_overlapping_genes_are_declared(self, config):
        # eltA and eltB overlap in the elt operon, so the pairing floor must
        # admit a negative intergenic distance.
        assert config["intergenic_min"] < 0

    def test_residue_positions_are_the_documented_substitutions(self, config):
        rule, = config["residue_rules"]
        positions = [(p["subunit"], p["index"]) for p in rule["positions"]]
        # Mature-protein A190/A196/A213/A224 and B75, shifted by the 18- and
        # 21-residue signal peptides of the precursor references.
        assert positions == [("A", 207), ("A", 213), ("A", 230), ("A", 241),
                             ("B", 95)]


class TestEtecLtCalls:
    def test_h10407_is_lt1(self, config, profiles, tmp_path):
        call, = calls_for("h10407_hits.tsv.gz", config, profiles, tmp_path)
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

    def test_lt2_signature_operon_is_lt2(self, config, profiles, tmp_path):
        call, = calls_for("lt2_synthetic_hits.tsv.gz", config, profiles, tmp_path)
        assert call["profile_id"] == "LT2"
        assert call["operon_status"] == "COMPLETE"
        assert call["operon"]["residue_evidence"] == {
            "A207": "L", "A213": "D", "A230": "E", "A241": "T", "B95": "A",
        }
        assert call["operon"]["combined_identity"] == 1.0
