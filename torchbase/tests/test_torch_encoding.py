"""Non-ASCII text must survive a torch write/read round-trip on any platform.

torchbase writes metadata.toml, signature.toml, profiles.tsv and
_resources/*.fasta and later reads them back. Any writer that relies on the
ambient locale encoding corrupts non-ASCII text as soon as the reader assumes
UTF-8 -- which `toml.load(path)` always does, and which every POSIX box does by
default. `ambient_encoding` reproduces that mismatch deterministically: it
supplies a chosen default to every text handle opened without an explicit
`encoding=`, separately for reads and writes, so "written on a cp1252 Windows
box, read anywhere else" is a plain unit test rather than a platform accident.
"""

import builtins
import io
from contextlib import contextmanager
from pathlib import Path

import pytest
import toml

from torchbase.conversions import pubmlst, stxtyper
from torchbase.torchfs import Torch, register_torch

EM_DASH = "\u2014"
AUTHOR = "Jos\u00e9 Mu\u00f1oz"
SPECIES = "Br\u00fccella melitensis"
LONG_DESCRIPTION = f"Brucella typing scheme {EM_DASH} snapshot of the 2024 release"


@contextmanager
def ambient_encoding(write="utf-8", read="utf-8"):
    """Force default encodings on text handles opened without an explicit one.

    Stands in for the platform locale, which is what unqualified `open()`,
    `Path.read_text()` and `Path.write_text()` would otherwise pick up.
    Reads and writes are set independently so a cp1252 writer can be paired
    with the UTF-8 reader that actually consumes the file.
    """
    real_open = builtins.open
    real_read_text = Path.read_text
    real_write_text = Path.write_text

    def fake_open(file, mode="r", *args, **kwargs):
        # args = (buffering, encoding, ...): an encoding given positionally is
        # already explicit, so leave it alone.
        if "b" not in mode and len(args) < 2 and kwargs.get("encoding") is None:
            writing = any(c in mode for c in "wax+")
            kwargs["encoding"] = write if writing else read
        return real_open(file, mode, *args, **kwargs)

    def fake_read_text(self, enc=None, errors=None):
        return real_read_text(self, enc or read, errors)

    def fake_write_text(self, data, enc=None, errors=None):
        return real_write_text(self, data, enc or write, errors)

    builtins.open = fake_open
    Path.read_text = fake_read_text
    Path.write_text = fake_write_text
    try:
        yield
    finally:
        builtins.open = real_open
        Path.read_text = real_read_text
        Path.write_text = real_write_text


def _seed_metadata(torch_dir, metadata):
    """Lay down the metadata.toml that the production writer will re-read."""
    (torch_dir / "metadata.toml").write_text(
        toml.dumps(metadata), encoding="utf-8"
    )


def test_ambient_encoding_shim_is_effective(tmp_path):
    """Guard the guard: the shim must really change unqualified handle defaults."""
    path = tmp_path / "probe.txt"
    with ambient_encoding(write="cp1252"):
        with open(path, "w") as f:
            f.write(EM_DASH)
    assert path.read_bytes() == b"\x97"
    with pytest.raises(UnicodeDecodeError):
        path.read_text(encoding="utf-8")
    with ambient_encoding(read="cp1252"):
        assert path.read_text() == EM_DASH


def test_metadata_non_ascii_survives_write_then_load(tmp_path):
    """metadata.toml written under a cp1252 locale still loads as UTF-8."""
    torch_dir = tmp_path / "brucella" / "mlst" / "1.0.0.torch"
    (torch_dir / "_resources").mkdir(parents=True)

    metadata = {
        "namespace": "brucella",
        "name": "mlst",
        "version": "1.0.0",
        "description": {"short": SPECIES, "long": LONG_DESCRIPTION},
        "authors": [AUTHOR],
        "manifest": {"profiles": "profiles.tsv", "resources": []},
    }
    _seed_metadata(torch_dir, metadata)
    (torch_dir / "profiles.tsv").write_text("ST\tadk\n1\t1\n", encoding="utf-8")
    (torch_dir / "_resources" / "adk.fasta").write_text(
        ">adk_1\nACGT\n", encoding="utf-8"
    )

    # register_torch is the production metadata.toml read-modify-write: the
    # same `toml.dump` into a text handle that every converter performs.
    # Reads are UTF-8 (the file was authored correctly); only the write picks
    # up the cp1252 locale, exactly as it would on a Windows workstation.
    with ambient_encoding(write="cp1252", read="utf-8"):
        register_torch(torch_dir)

    # ...and now read the torch on a machine whose locale is UTF-8.
    with ambient_encoding():
        torch = Torch.load(torch_dir)
        reloaded = toml.load(torch_dir / "metadata.toml")

    assert reloaded["description"]["long"] == LONG_DESCRIPTION
    assert EM_DASH in reloaded["description"]["long"]
    assert reloaded["description"]["short"] == SPECIES
    assert reloaded["authors"] == [AUTHOR]
    assert torch.path == torch_dir
    # metadata.toml on disk is UTF-8 regardless of the writing locale.
    assert LONG_DESCRIPTION.encode("utf-8") in (
        torch_dir / "metadata.toml"
    ).read_bytes()


def test_profiles_tsv_non_ascii_survives_write_then_load(tmp_path):
    """A non-ASCII profiles.tsv value written by a converter reloads intact."""
    torch_dir = tmp_path / "brucella" / "mlst" / "1.0.0.torch"
    (torch_dir / "_resources").mkdir(parents=True)

    _seed_metadata(
        torch_dir,
        {
            "namespace": "brucella",
            "name": "mlst",
            "version": "1.0.0",
            "manifest": {"profiles": "profiles.tsv", "resources": []},
        },
    )
    (torch_dir / "_resources" / "adk.fasta").write_text(
        ">adk_1\nACGT\n", encoding="utf-8"
    )

    profiles = [
        {"ST": "1", "adk": "1", "biovar": SPECIES},
        {"ST": "2", "adk": "2", "biovar": f"suis{EM_DASH}like"},
    ]
    with ambient_encoding(write="cp1252", read="utf-8"):
        pubmlst._write_profiles_tsv(torch_dir / "profiles.tsv", profiles)

    with ambient_encoding():
        torch = Torch.load(torch_dir)

    assert [p.profile for p in torch.profile.profiles] == ["1", "2"]
    assert torch.profile.profiles[0]["biovar"] == SPECIES
    assert torch.profile.profiles[1]["biovar"] == f"suis{EM_DASH}like"
    assert SPECIES.encode("utf-8") in (torch_dir / "profiles.tsv").read_bytes()


def test_stxtyper_conversion_round_trips_non_ascii_locus_names(tmp_path):
    """End-to-end converter: non-ASCII flows into profiles.tsv and subunits.faa."""
    accession = f"AAS07596.1{EM_DASH}v2"
    stx_prot = io.StringIO(
        f">{accession}|stxA2c|stxA2\nMKCILF\n"
        f">BAB{EM_DASH}1|stxB2c|stxB2\nMKKTLL\n"
    )

    with ambient_encoding(write="cp1252", read="utf-8"):
        torch_dir = Path(
            stxtyper.convert_local(
                stx_prot, output_path=str(tmp_path), version="1.0.0"
            )
        )

    with ambient_encoding():
        torch = Torch.load(torch_dir)

    subunits = (torch_dir / "_resources" / "subunits.faa").read_text(
        encoding="utf-8"
    )
    assert accession in subunits
    assert f"BAB{EM_DASH}1" in subunits
    row = torch.operon_profiles[0]
    assert row["subunit_A"] == "stxA2"
    assert (torch_dir / "profiles.tsv").read_bytes().count(b"\xe2\x80\x94") == 0
    assert accession.encode("utf-8") in subunits.encode("utf-8")


def test_signing_reads_converter_written_metadata(tmp_path):
    """The reported failure mode, end to end.

    A converter writes metadata.toml through a text handle; `signing.sign_torch`
    reads it back with `toml.load(path)`, which always decodes UTF-8. Under a
    cp1252 locale the em dash goes out as byte 0x97 and the read blows up with
    `UnicodeDecodeError: 'utf-8' codec can't decode byte 0x97`.
    """
    from torchbase import signing

    torch_dir = tmp_path / "brucella" / "mlst" / "1.0.0.torch"
    (torch_dir / "_resources").mkdir(parents=True)
    _seed_metadata(
        torch_dir,
        {
            "namespace": "brucella",
            "name": "mlst",
            "version": "1.0.0",
            "description": {"short": SPECIES, "long": LONG_DESCRIPTION},
            "authors": [AUTHOR],
            "manifest": {"profiles": "profiles.tsv", "resources": []},
        },
    )
    (torch_dir / "profiles.tsv").write_text("ST\tadk\n1\t1\n", encoding="utf-8")
    (torch_dir / "_resources" / "adk.fasta").write_text(
        ">adk_1\nACGT\n", encoding="utf-8"
    )

    priv_path, _ = signing.generate_software_keypair("brucella", tmp_path)
    signer = signing.FileKeySigner(priv_path)

    with ambient_encoding(write="cp1252", read="utf-8"):
        register_torch(torch_dir)
        sig_path = signing.sign_torch(torch_dir, signer)

    with ambient_encoding():
        torch = Torch.load(torch_dir)
        sig_data = toml.load(sig_path)
        reloaded = toml.load(torch_dir / "metadata.toml")

    assert reloaded["description"]["long"] == LONG_DESCRIPTION
    assert reloaded["authors"] == [AUTHOR]
    assert sig_data["signature"]["namespace"] == "brucella"
    assert torch.signature["signature"]["content_hash"] == (
        sig_data["signature"]["content_hash"]
    )
