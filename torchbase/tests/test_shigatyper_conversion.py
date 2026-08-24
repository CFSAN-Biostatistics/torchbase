"""Tests for torchbase.conversions.shigatyper.

Covers the mechanical conversion (regrouping ShigaTyper's consolidated
reference by gene into real wzx/wzy multi-allele loci, torch structure,
namespace default) and the canonical profiles.tsv table -- a mechanical
transcription of upstream's checkpoint decision cascade (docs/adr/0003).
"""

import pytest

from torchbase.conversions.shigatyper import (
    ABSENT, ACCESSORY_LOCI, LOCUS_ORDER, PRESENT, WILD,
    _build_canonical_profiles, convert_local, download_sources,
)
from torchbase.torchfs import Torch
from torchbase import profile_match

# One wzx-suffixed, one wzy-suffixed, one singleton -- enough to exercise
# download_sources's regrouping without needing all 95 real records.
CONSOLIDATED_FASTA = """\
>ipaH_c
ATGAAACCCGGG
>Sb11_wzx
TTTGGGCCCAAA
>Sb11_wzy
CATCATCATCAT
>Sb12_wzx
GGGGCCCCAAAA
"""


def _calls(**present_or_absent):
    """allele_calls-shaped calls: every value defaults to ABSENT."""
    return {
        locus: {"allele_id": present_or_absent.get(locus, ABSENT)}
        for locus in LOCUS_ORDER
    }


def _match(**present_or_absent):
    """(profile_id, status) the canonical table resolves these calls to."""
    calls = _calls(**present_or_absent)
    rows = _build_canonical_profiles()
    profile_string = profile_match.build_profile_string(calls, LOCUS_ORDER)
    return profile_match.match_profile(profile_string, rows, LOCUS_ORDER, id_column="Serotype")


class TestDownloadSources:
    def test_regroups_wzx_and_wzy_suffixed_records_by_gene(self, tmp_path, monkeypatch):
        def fake_fetch_file(url, dest, **kwargs):
            dest.write_text(CONSOLIDATED_FASTA, encoding="utf-8")
            return dest

        monkeypatch.setattr("torchbase.conversions.fetch_file", fake_fetch_file)
        result = download_sources(tmp_path)
        names = {f.name.split("\\")[-1].split("/")[-1] for f in result["sequences"]}
        assert names == {"wzx.fasta", "wzy.fasta", "ipaH_c.fasta"}

    def test_wzx_alleles_are_renamed_gene_first(self, tmp_path, monkeypatch):
        def fake_fetch_file(url, dest, **kwargs):
            dest.write_text(CONSOLIDATED_FASTA, encoding="utf-8")
            return dest

        monkeypatch.setattr("torchbase.conversions.fetch_file", fake_fetch_file)
        result = download_sources(tmp_path)
        wzx_file = next(f for f in result["sequences"] if f.name.endswith("wzx.fasta"))
        content = wzx_file.read()
        # Renamed "wzx_{variant}", the convention extract_locus_and_allele expects.
        assert ">wzx_Sb11" in content
        assert ">wzx_Sb12" in content
        assert ">Sb11_wzx" not in content

    def test_singleton_markers_get_one_file_and_an_explicit_allele_id(self, tmp_path, monkeypatch):
        def fake_fetch_file(url, dest, **kwargs):
            dest.write_text(CONSOLIDATED_FASTA, encoding="utf-8")
            return dest

        monkeypatch.setattr("torchbase.conversions.fetch_file", fake_fetch_file)
        result = download_sources(tmp_path)
        ipah = next(f for f in result["sequences"] if f.name.endswith("ipaH_c.fasta"))
        # "_1", not bare "ipaH_c": the header already contains an underscore
        # that would otherwise misparse as a locus/allele split.
        assert ipah.read().strip() == ">ipaH_c_1\nATGAAACCCGGG"


class TestConvertLocal:
    def _canonical_markers(self, tmp_path):
        """One-file-per-locus markers matching LOCUS_ORDER, minimal sequences."""
        files = []
        wzx_path = tmp_path / "wzx.fasta"
        wzx_path.write_text(">wzx_Sb11\nACGT\n>wzx_Ss\nTTTT\n")
        files.append(open(wzx_path))
        wzy_path = tmp_path / "wzy.fasta"
        wzy_path.write_text(">wzy_Sb11\nGGGG\n")
        files.append(open(wzy_path))
        for locus in ACCESSORY_LOCI:
            path = tmp_path / f"{locus}.fasta"
            path.write_text(f">{locus}_1\nACGTACGTACGT\n")
            files.append(open(path))
        return files

    def test_creates_a_loadable_torch_with_default_namespace(self, tmp_path):
        torch_path = convert_local(
            self._canonical_markers(tmp_path), output_path=str(tmp_path / "out"),
        )
        torch = Torch.load(torch_path)
        assert torch.typing_model == "allelic"

    def test_canonical_locus_set_gets_the_real_cascade_table(self, tmp_path):
        torch_path = convert_local(
            self._canonical_markers(tmp_path), output_path=str(tmp_path / "out"),
        )
        torch = Torch.load(torch_path)
        assert torch.profile is not None
        # A known row (boydii 11) must actually be present in the loaded table.
        assert any(
            row.get("Serotype") == "Shigella boydii serotype 11"
            for row in csv_rows(torch_path)
        )

    def test_non_canonical_locus_set_falls_back_to_a_stub(self, tmp_path):
        one_locus = tmp_path / "cadA.fasta"
        one_locus.write_text(">cadA\nACGT\n")
        torch_path = convert_local([open(one_locus)], output_path=str(tmp_path / "out"))
        rows = csv_rows(torch_path)
        assert rows == []  # header-only: no canonical row can apply to this locus set

    def test_supplied_profiles_override_the_canonical_table(self, tmp_path):
        profiles = tmp_path / "custom.tsv"
        profiles.write_text("Serotype\twzx\twzy\nCustom\tA\tB\n")
        torch_path = convert_local(
            self._canonical_markers(tmp_path),
            profiles_file=open(profiles), output_path=str(tmp_path / "out"),
        )
        rows = csv_rows(torch_path)
        assert rows == [{"Serotype": "Custom", "wzx": "A", "wzy": "B"}]


def csv_rows(torch_path):
    import csv
    from pathlib import Path
    with open(Path(torch_path) / "profiles.tsv", newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter="\t"))


class TestCanonicalProfilesCascade:
    """Spot-checks against upstream's actual checkpoint decision cascade."""

    def test_species_gate_ipah_absent_and_no_exception_is_no_match(self):
        profile_id, status = _match(ipaH_c=ABSENT, wzx="Sb12")
        assert profile_id is None

    def test_boydii_13_species_exception_survives_ipah_absence(self):
        profile_id, _ = _match(ipaH_c=ABSENT, wzx="Sb13")
        assert profile_id == "Shigella boydii serotype 13 (no longer classified as Shigella)"

    def test_eclacy_forces_eiec_except_for_boydii_9_and_15(self):
        assert _match(ipaH_c=PRESENT, EclacY=PRESENT, wzx="Sb12")[0] == "EIEC"
        assert _match(ipaH_c=PRESENT, EclacY=PRESENT, wzx="Sb9")[0] == "Shigella boydii serotype 9"
        assert _match(ipaH_c=PRESENT, EclacY=PRESENT, wzx="Sb15")[0] == "Shigella boydii serotype 15"

    def test_cada_present_gates_out_the_generic_numeric_fallback(self):
        # cadA+ with a wzx the cadA branch doesn't name specifically -> EIEC,
        # not "boydii serotype 12" (which the generic fallback would give).
        assert _match(ipaH_c=PRESENT, cadA=PRESENT, wzx="Sb12")[0] == "EIEC"
        assert _match(ipaH_c=PRESENT, cadA=ABSENT, wzx="Sb12")[0] == "Shigella boydii serotype 12"

    def test_sonnei_forms_distinguished_by_wzx_and_ipab(self):
        assert _match(ipaH_c=PRESENT, cadA=PRESENT, Ss_methylase=PRESENT, wzx="Ss")[0] == \
            "Shigella sonnei, form I"
        assert _match(ipaH_c=PRESENT, cadA=PRESENT, Ss_methylase=PRESENT, ipaB=ABSENT)[0] == \
            "Shigella sonnei form II"
        assert _match(ipaH_c=PRESENT, cadA=PRESENT, Ss_methylase=PRESENT, ipaB=PRESENT)[0] == \
            "Shigella sonnei (low levels of form I)"

    def test_dysenteriae_1_rfp_status_only_matters_when_cada_present(self):
        assert _match(ipaH_c=PRESENT, cadA=PRESENT, wzx="Sd1", Sd1_rfp=PRESENT)[0] == \
            "Shigella dysenteriae serotype 1"
        assert _match(ipaH_c=PRESENT, cadA=PRESENT, wzx="Sd1", Sd1_rfp=ABSENT)[0] == \
            "Shigella dysenteriae serotype 1, rfp- (phenotypically negative)"
        # cadA absent: generic fallback, no rfp distinction available or needed.
        assert _match(ipaH_c=PRESENT, cadA=ABSENT, wzx="Sd1")[0] == \
            "Shigella dysenteriae serotype 1"

    def test_boydii_6_and_10_are_never_split(self):
        """Documented gap: no raw-read depth access, so both collapse to one call."""
        assert _match(ipaH_c=PRESENT, wzx="Sb6")[0] == "Shigella boydii serotype 6 or 10"
        assert _match(ipaH_c=PRESENT, wzx="Sb6", wbaM=PRESENT)[0] == "Shigella boydii serotype 6 or 10"
        assert _match(ipaH_c=PRESENT, wzx="Sb6", wbaM=ABSENT)[0] == "Shigella boydii serotype 6 or 10"

    def test_boydii_1_vs_20_by_heparinase(self):
        assert _match(ipaH_c=PRESENT, wzx="Sb1", heparinase=PRESENT)[0] == "Shigella boydii serotype 20"
        assert _match(ipaH_c=PRESENT, wzx="Sb1", heparinase=ABSENT)[0] == "Shigella boydii serotype 1"

    def test_provisional_serotypes(self):
        assert _match(ipaH_c=PRESENT, wzx="SbProv")[0] == "Shigella boydii Provisional serotype E1621-54"
        assert _match(ipaH_c=PRESENT, wzx="SdProv")[0] == "Shigella dysenteriae Provisional serotype 96-265"
        assert _match(ipaH_c=PRESENT, wzx="SdProvE")[0] == "Shigella dysenteriae Provisional serotype E670-74"

    def test_generic_numeric_fallback_covers_every_remaining_boydii_and_dysenteriae_variant(self):
        for n in (2, 3, 4, 5, 7, 8, 9, 11, 12, 13, 14, 15, 16, 17, 18, 19):
            assert _match(ipaH_c=PRESENT, wzx=f"Sb{n}")[0] == f"Shigella boydii serotype {n}"
        for n in range(1, 16):
            assert _match(ipaH_c=PRESENT, wzx=f"Sd{n}")[0] == f"Shigella dysenteriae serotype {n}"

    def test_flexneri_6_is_its_own_wzx_variant(self):
        assert _match(ipaH_c=PRESENT, wzx="Sf6")[0] == "Shigella flexneri serotype 6"

    @pytest.mark.parametrize("present_loci,expected", [
        ((), "Shigella flexneri serotype Y"),
        (("Xv",), "Shigella flexneri Yv"),
        (("gtrI",), "Shigella flexneri serotype 1a"),
        (("gtrI", "Oac1b"), "Shigella flexneri serotype 1b"),
        (("gtrII",), "Shigella flexneri serotype 2a"),
        (("gtrII", "gtrX"), "Shigella flexneri serotype 2b"),
        (("gtrX", "Oac"), "Shigella flexneri serotype 3a"),
        (("Oac",), "Shigella flexneri serotype 3b"),
        (("Oac1b",), "Shigella flexneri serotype 3b"),
        (("gtrIV",), "Shigella flexneri serotype 4a"),
        (("gtrIV", "Oac"), "Shigella flexneri serotype 4b"),
        (("gtrIV", "Oac1b"), "Shigella flexneri serotype 4b"),
        (("gtrIV", "Oac", "Xv"), "Shigella flexneri 4bv"),
        (("gtrV", "Oac"), "Shigella flexneri serotype 5a"),
        (("gtrV",), "Shigella flexneri serotype 5a"),
        (("gtrV", "gtrX", "Oac"), "Shigella flexneri serotype 5b"),
        (("gtrV", "gtrX"), "Shigella flexneri serotype 5b"),
        (("gtrX",), "Shigella flexneri serotype X"),
        (("gtrX", "Xv"), "Shigella flexneri serotype Xv (4c)"),
        (("gtrI", "gtrIC"), "Shigella flexneri serotype 1c (7a)"),
        (("gtrI", "gtrIC", "Oac1b"), "Shigella flexneri serotype 7b"),
    ])
    def test_flexneri_accessory_gene_combinations(self, present_loci, expected):
        kwargs = {"ipaH_c": PRESENT, "cadA": ABSENT, "wzx": "Sf"}
        kwargs.update({locus: PRESENT for locus in present_loci})
        assert _match(**kwargs)[0] == expected

    def test_unrecognized_flexneri_accessory_combination_is_novel(self):
        # gtrI + gtrII together names no SfDic entry.
        profile_id, status = _match(
            ipaH_c=PRESENT, cadA=ABSENT, wzx="Sf", gtrI=PRESENT, gtrII=PRESENT,
        )
        assert profile_id is None
        assert status == "novel_profile"



class TestSingletonMarkerHeadersDontCollideWithAlleleSplitting:
    """extract_locus_and_allele splits at the LAST underscore; a singleton
    marker whose own name contains one (ipaH_c, Ss_methylase, Sd1_rfp) would
    misparse as {locus}_{spurious allele} without the "_1" suffix
    download_sources writes. Exercises the real regrouping, not a stub.
    """

    @pytest.mark.parametrize("marker", ["ipaH_c", "Ss_methylase", "Sd1_rfp"])
    def test_underscore_bearing_marker_name_round_trips_as_one_locus(self, tmp_path, monkeypatch, marker):
        from torchbase.allele_calls import extract_locus_and_allele

        fasta = f">ipaH_c\nAAAA\n>Sb11_wzx\nCCCC\n>{marker}\nGGGG\n"
        monkeypatch.setattr(
            "torchbase.conversions.fetch_file",
            lambda url, dest, **kw: (dest.write_text(fasta, encoding="utf-8"), dest)[1],
        )
        result = download_sources(tmp_path)
        marker_file = next(f for f in result["sequences"] if f.name.endswith(f"{marker}.fasta"))
        header = marker_file.read().splitlines()[0][1:]
        assert extract_locus_and_allele(header) == (marker, "1")