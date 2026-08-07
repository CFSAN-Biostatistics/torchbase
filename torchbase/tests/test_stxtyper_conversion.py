"""Tests for the StxTyper -> operon torch converter (docs/operon-strategy-plan.md).

Uses a small synthetic stx.prot fixture (real header format, short
sequences) rather than the network — network conversion is exercised
manually via `torchtools convert stxtyper --download`.
"""

import io

import pytest

from torchbase.conversions.stxtyper import convert_local
from torchbase.torchfs import Torch

# Minimal fixture covering: a non-generalized class (1a) and a generalized
# class (2, via 2a/2c partners) — enough to exercise both branches of
# profiles.tsv construction and [operon] validation.
STX_PROT_FIXTURE = """\
>AAS07500.1|stxA1a|stxA1a
MKCILFCVLTLTAIS
>AAS07501.1|stxB1a|stxB1a
MKKTLIAAAV
>AAS07596.1|stxA2c|stxA2
MKCFAAALLVSFTLA
>AAS07597.1|stxB2c|stxB2
MDVVSVKKR
"""


class TestConvertLocal:
    def test_creates_loadable_operon_torch(self, tmp_path):
        torch_path = convert_local(
            stx_prot_file=io.StringIO(STX_PROT_FIXTURE),
            output_path=tmp_path,
            namespace="ncbi",
            name="stxtyper",
            version="1.0.45",
            stxtyper_version="1.0.45",
        )

        torch = Torch.load(torch_path)
        assert torch.typing_model == "operon"
        assert torch.operon_config["subunit_order"] == ["A", "B"]
        assert torch.profile is None  # bypasses Profile matching (§2)

    def test_profiles_tsv_groups_subunits_by_type(self, tmp_path):
        torch_path = convert_local(
            stx_prot_file=io.StringIO(STX_PROT_FIXTURE),
            output_path=tmp_path,
            version="1.0.45",
        )
        torch = Torch.load(torch_path)

        by_subtype = {row["subtype"]: row for row in torch.operon_profiles}
        assert by_subtype["stx1a"]["class"] == "1a"
        assert by_subtype["stx1a"]["subunit_A"] == "stxA1a"
        assert by_subtype["stx1a"]["subunit_B"] == "stxB1a"

        # 2c collapses into generalized class "2" (identity alone cannot
        # separate 2a/2c/2d), and the reported reference subtype is
        # stx.prot's subclass field, not its famId.
        assert by_subtype["stx2c"]["class"] == "2"
        assert by_subtype["stx2c"]["subunit_A"] == "stxA2"
        assert by_subtype["stx2c"]["subunit_B"] == "stxB2"

    def test_reference_headers_use_declared_format(self, tmp_path):
        torch_path = convert_local(
            stx_prot_file=io.StringIO(STX_PROT_FIXTURE),
            output_path=tmp_path,
            version="1.0.45",
        )
        torch = Torch.load(torch_path)
        subunits_faa = torch.path / "_resources" / "subunits.faa"
        headers = [
            line[1:].strip()
            for line in subunits_faa.read_text().splitlines()
            if line.startswith(">")
        ]
        assert "AAS07500.1|A|stxA1a|1a" in headers
        assert "AAS07596.1|A|stxA2|2" in headers

    def test_identity_thresholds_transcribed(self, tmp_path):
        torch_path = convert_local(
            stx_prot_file=io.StringIO(STX_PROT_FIXTURE),
            output_path=tmp_path,
            version="1.0.45",
        )
        torch = Torch.load(torch_path)
        assert torch.operon_config["identity_thresholds"]["1a"] == 0.983
        assert torch.operon_config["identity_thresholds"]["default"] == 0.98

    def test_residue_rule_positions_match_stxtyper_offsets(self, tmp_path):
        torch_path = convert_local(
            stx_prot_file=io.StringIO(STX_PROT_FIXTURE),
            output_path=tmp_path,
            version="1.0.45",
        )
        torch = Torch.load(torch_path)
        rule = torch.operon_config["residue_rules"][0]
        positions = {(p["subunit"], p["index"]) for p in rule["positions"]}
        # StxTyper reads these residues out of qMap(), whose character i is
        # the subject residue aligned to 0-based reference offset i, so
        # A312/A318/B34 are transcribed verbatim (stxtyper.cpp:536-556).
        assert positions == {("A", 312), ("A", 318), ("B", 34)}

    def test_empty_stx_prot_raises(self, tmp_path):
        with pytest.raises(ValueError):
            convert_local(stx_prot_file=io.StringIO(""), output_path=tmp_path)
