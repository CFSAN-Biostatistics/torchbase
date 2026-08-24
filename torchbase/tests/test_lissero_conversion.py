"""Tests for torchbase.conversions.lissero.

Covers the mechanical conversion (splitting LisSero's combined
sequences.fasta into one file per gene, torch structure) and the canonical
profiles.tsv table -- a mechanical transcription of upstream's
Serotype.report_maker() decision tree (docs/adr/0003).
"""

import pytest

from torchbase.conversions.lissero import (
    ABSENT, CANONICAL_LOCI, PRESENT, WILD,
    _build_canonical_profiles, convert_local, download_sources,
)
from torchbase.torchfs import Torch
from torchbase import profile_match

# Headers use LisSero's real "~~" separator between gene name and accession.
CONSOLIDATED_FASTA = """\
>Prs~~FR733643.1:206011-206380 Listeria monocytogenes strain, serotype 4d
GCTGAAGAGATTGCGAAAGAAGTTGGT
>lmo0737~~HG421741.1:749402-750092 Listeria monocytogenes EGD
AGGGCTTCAAGGACTTACCCTCGAAGA
"""


def _calls(**present_or_absent):
    return {
        locus: {"allele_id": present_or_absent.get(locus, ABSENT)}
        for locus in CANONICAL_LOCI
    }


def _match(**present_or_absent):
    calls = _calls(**present_or_absent)
    rows = _build_canonical_profiles()
    profile_string = profile_match.build_profile_string(calls, CANONICAL_LOCI)
    return profile_match.match_profile(profile_string, rows, CANONICAL_LOCI, id_column="Serogroup")


class TestDownloadSources:
    def test_splits_combined_fasta_into_one_file_per_gene(self, tmp_path, monkeypatch):
        def fake_fetch_file(url, dest, **kwargs):
            dest.write_text(CONSOLIDATED_FASTA, encoding="utf-8")
            return dest

        monkeypatch.setattr("torchbase.conversions.fetch_file", fake_fetch_file)
        result = download_sources(tmp_path)
        names = {f.name.split("\\")[-1].split("/")[-1] for f in result["sequences"]}
        assert names == {"Prs.fasta", "lmo0737.fasta"}

    def test_header_is_truncated_at_the_double_tilde_separator(self, tmp_path, monkeypatch):
        def fake_fetch_file(url, dest, **kwargs):
            dest.write_text(CONSOLIDATED_FASTA, encoding="utf-8")
            return dest

        monkeypatch.setattr("torchbase.conversions.fetch_file", fake_fetch_file)
        result = download_sources(tmp_path)
        prs_file = next(f for f in result["sequences"] if f.name.endswith("Prs.fasta"))
        content = prs_file.read()
        assert content.startswith(">Prs_1\n")
        assert "~~" not in content
        assert "FR733643" not in content


class TestConvertLocal:
    def _canonical_markers(self, tmp_path):
        files = []
        for locus in CANONICAL_LOCI:
            path = tmp_path / f"{locus}.fasta"
            path.write_text(f">{locus}\nACGTACGTACGT\n")
            files.append(open(path))
        return files

    def test_creates_a_loadable_torch(self, tmp_path):
        torch_path = convert_local(
            self._canonical_markers(tmp_path), output_path=str(tmp_path / "out"),
        )
        torch = Torch.load(torch_path)
        assert torch.typing_model == "allelic"
        assert {ref.stem for ref in torch.references} == set(CANONICAL_LOCI)

    def test_canonical_locus_set_gets_the_real_cascade_table(self, tmp_path):
        torch_path = convert_local(
            self._canonical_markers(tmp_path), output_path=str(tmp_path / "out"),
        )
        rows = csv_rows(torch_path)
        assert any(row.get("Serogroup") == "4b, 4d, 4e*" for row in rows)

    def test_non_canonical_locus_set_falls_back_to_a_stub(self, tmp_path):
        one_locus = tmp_path / "Prs.fasta"
        one_locus.write_text(">Prs\nACGT\n")
        torch_path = convert_local([open(one_locus)], output_path=str(tmp_path / "out"))
        assert csv_rows(torch_path) == []

    def test_supplied_profiles_override_the_canonical_table(self, tmp_path):
        profiles = tmp_path / "custom.tsv"
        profiles.write_text("Serogroup\tPrs\nCustom\tA\n")
        torch_path = convert_local(
            self._canonical_markers(tmp_path),
            profiles_file=open(profiles), output_path=str(tmp_path / "out"),
        )
        assert csv_rows(torch_path) == [{"Serogroup": "Custom", "Prs": "A"}]


def csv_rows(torch_path):
    import csv
    from pathlib import Path
    with open(Path(torch_path) / "profiles.tsv", newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter="\t"))


class TestCanonicalProfilesCascade:
    """Spot-checks against upstream's actual Serotype.report_maker()."""

    def test_prs_absent_is_nontypeable_unconditionally(self):
        assert _match(Prs=ABSENT, lmo0737=PRESENT, ORF2819=PRESENT)[0] == "Nontypeable"

    def test_lmo0737_branch_split_by_lmo1118(self):
        assert _match(Prs=PRESENT, lmo0737=PRESENT, lmo1118=PRESENT)[0] == "1/2c, 3c"
        assert _match(Prs=PRESENT, lmo0737=PRESENT, lmo1118=ABSENT)[0] == "1/2a, 3a"

    def test_orf2819_branch_split_by_orf2110(self):
        assert _match(Prs=PRESENT, ORF2819=PRESENT, ORF2110=PRESENT)[0] == "4b, 4d, 4e"
        assert _match(Prs=PRESENT, ORF2819=PRESENT, ORF2110=ABSENT)[0] == "1/2b, 3b, 7"

    def test_all_four_accessory_genes_present_split_by_lmo1118(self):
        assert _match(
            Prs=PRESENT, lmo0737=PRESENT, ORF2819=PRESENT, ORF2110=PRESENT, lmo1118=PRESENT,
        )[0] == "Nontypeable"
        assert _match(
            Prs=PRESENT, lmo0737=PRESENT, ORF2819=PRESENT, ORF2110=PRESENT, lmo1118=ABSENT,
        )[0] == "4b, 4d, 4e*"

    def test_unrecognized_combination_falls_through_as_novel(self):
        # Prs+ with only lmo1118 present matches no branch upstream defines.
        profile_id, status = _match(Prs=PRESENT, lmo1118=PRESENT)
        assert profile_id is None
        assert status == "novel_profile"
