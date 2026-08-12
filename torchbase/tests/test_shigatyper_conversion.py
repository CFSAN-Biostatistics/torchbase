"""Tests for torchbase.conversions.shigatyper.

Covers the mechanical conversion (splitting ShigaTyper's consolidated
reference into one-locus-per-marker resource files, torch structure,
namespace default) and pins the documented gap: this converter does not, and
should not, fabricate ShigaTyper's serotype decision cascade as profiles.tsv
rows.
"""

import io

import pytest

from torchbase.conversions.shigatyper import _parse_fasta_records, convert_local
from torchbase.torchfs import Torch

MULTI_MARKER_FASTA = """\
>ipaH_c
ATGAAACCCGGG
>Ss_wzx
TTTGGGCCCAAA
>Sb19_wzy
CATCATCATCAT
"""


class TestParseFastaRecords:
    def test_splits_a_multi_record_fasta(self):
        records = list(_parse_fasta_records(MULTI_MARKER_FASTA))
        assert records == [
            ("ipaH_c", "ATGAAACCCGGG"),
            ("Ss_wzx", "TTTGGGCCCAAA"),
            ("Sb19_wzy", "CATCATCATCAT"),
        ]

    def test_header_is_truncated_at_first_whitespace(self):
        records = list(_parse_fasta_records(">ipaH_c some description\nACGT\n"))
        assert records == [("ipaH_c", "ACGT")]

    def test_empty_input_yields_nothing(self):
        assert list(_parse_fasta_records("")) == []


class TestConvertLocal:
    def _markers(self, tmp_path, names):
        files = []
        for name in names:
            path = tmp_path / f"{name}.fasta"
            path.write_text(f">{name}\nACGTACGTACGT\n")
            files.append(open(path))
        return files

    def test_creates_a_loadable_torch_with_default_namespace(self, tmp_path):
        torch_path = convert_local(
            sequence_files=self._markers(tmp_path, ["ipaH_c", "Ss_wzx"]),
            output_path=tmp_path,
            version="1.0.0",
        )
        torch = Torch.load(torch_path)
        # hfp: a torch's namespace names the authority for the data, and
        # ShigaTyper is an FDA Human Foods Program product.
        assert torch.path.parent.parent.name == "hfp"
        assert torch.path.parent.name == "shigatyper"
        assert torch.typing_model == "allelic"

    def test_one_locus_file_per_marker(self, tmp_path):
        torch_path = convert_local(
            sequence_files=self._markers(tmp_path, ["ipaH_c", "Ss_wzx", "Sb19_wzy"]),
            output_path=tmp_path,
        )
        torch = Torch.load(torch_path)
        assert {ref.stem for ref in torch.references} == {"ipaH_c", "Ss_wzx", "Sb19_wzy"}

    def test_no_profiles_file_writes_a_loadable_stub(self, tmp_path):
        torch_path = convert_local(
            sequence_files=self._markers(tmp_path, ["ipaH_c"]),
            output_path=tmp_path,
        )
        torch = Torch.load(torch_path)
        assert torch.profile is not None  # header-only, but present and parseable

    def test_supplied_profiles_are_used_verbatim(self, tmp_path):
        profiles = tmp_path / "serotypes.tsv"
        profiles.write_text("Serotype\tO\tH\nA\t1\t2\n")
        torch_path = convert_local(
            sequence_files=self._markers(tmp_path, ["ipaH_c"]),
            profiles_file=open(profiles),
            output_path=tmp_path,
        )
        written = (torch_path and __import__("pathlib").Path(torch_path) / "profiles.tsv").read_text()
        assert "A\t1\t2" in written


class TestKnownGapIsDocumented:
    """ShigaTyper's serotype call is a decision cascade, not a lookup table.

    This converter must not synthesize profiles.tsv rows that pretend to
    encode that cascade -- doing so would produce a torch that types
    confidently and wrong. The gap is recorded in the module docstring and in
    the torch's own metadata so it survives being read out of context.
    """

    def test_metadata_records_that_no_typing_model_consumes_this_yet(self, tmp_path):
        (tmp_path / "ipaH_c.fasta").write_text(">ipaH_c\nACGT\n")
        torch_path = convert_local(
            sequence_files=[open(tmp_path / "ipaH_c.fasta")],
            output_path=tmp_path,
        )
        import toml
        from pathlib import Path
        metadata = toml.load(Path(torch_path) / "metadata.toml")
        long_description = metadata["description"]["long"]
        assert "decision cascade" in long_description
        assert "no typing model" in long_description.lower()
