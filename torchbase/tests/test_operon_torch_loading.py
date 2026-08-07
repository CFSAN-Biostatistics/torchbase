"""Tests for `typing_model = "operon"` torch loading (Torch._load_single_scheme).

Covers docs/operon-strategy-plan.md §2, §3, §6:
- typing_model defaults to "allelic" (backward compatible)
- operon torches parse the [operon] block and bypass Profile matching
- operon config is validated against profiles.tsv at load time
"""

import csv

import pytest
import toml

from torchbase.operon import OperonConfigError
from torchbase.torchfs import Torch


def _write_metadata(torch_path, extra=None):
    metadata = {
        "namespace": "test_ns",
        "name": "stx",
        "version": "1.0.0",
        "version_meta": {"strategy": "semver", "timestamp": 1609459200},
        "manifest": {"profiles": "profiles.tsv"},
        "description": {"short": "stx operon torch"},
    }
    if extra:
        metadata.update(extra)
    with open(torch_path / "metadata.toml", "w") as f:
        toml.dump(metadata, f)


def _write_profiles(torch_path, rows):
    with open(torch_path / "profiles.tsv", "w", newline="") as f:
        writer = csv.DictWriter(
            f, fieldnames=["subtype", "class", "subunit_A", "subunit_B"], delimiter="\t"
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


VALID_OPERON_CFG = {
    "subunit_order": ["A", "B"],
    "intergenic_max": 36,
    "reference": {"file": "_resources/subunits.faa", "header_format": "accession|subunit_tag|class_label"},
    "identity_thresholds": {"default": 0.98, "1a": 0.983},
    "generalized_classes": {"2": ["2a", "2c", "2d"]},
    "residue_rules": [
        {
            "class": "2",
            "positions": [
                {"subunit": "A", "index": 312},
                {"subunit": "A", "index": 318},
                {"subunit": "B", "index": 34},
            ],
            "table": [
                {"call": "2a", "residues": [["F", "S"], ["K", "E"], ["D"]]},
            ],
            "fallback": "2",
        }
    ],
}

VALID_ROWS = [
    {"subtype": "stx1a", "class": "1a", "subunit_A": "stxA1a", "subunit_B": "stxB1a"},
    {"subtype": "stx2a", "class": "2", "subunit_A": "stxA2", "subunit_B": "stxB2a"},
]


@pytest.fixture
def torch_path(tmp_path):
    path = tmp_path / "test_ns" / "stx" / "1.0.0.torch"
    path.mkdir(parents=True)
    (path / "_resources").mkdir()
    (path / "_resources" / "subunits.faa").write_text(">AAS07596.1|stxA2|stxA2\nMKC\n")
    return path


class TestTypingModelDefault:
    def test_allelic_default_when_absent(self, torch_path):
        _write_metadata(torch_path)
        _write_profiles(torch_path, [{"subtype": "1", "class": "", "subunit_A": "1", "subunit_B": "1"}])
        torch = Torch.load(torch_path)
        assert torch.typing_model == "allelic"
        assert torch.operon_config is None
        assert torch.profile is not None


class TestOperonTorchLoading:
    def test_loads_operon_config(self, torch_path):
        _write_metadata(torch_path, {"typing_model": "operon", "operon": VALID_OPERON_CFG})
        _write_profiles(torch_path, VALID_ROWS)
        torch = Torch.load(torch_path)

        assert torch.typing_model == "operon"
        assert torch.operon_config["subunit_order"] == ["A", "B"]
        assert torch.operon_profiles == VALID_ROWS

    def test_bypasses_profile_matching(self, torch_path):
        _write_metadata(torch_path, {"typing_model": "operon", "operon": VALID_OPERON_CFG})
        _write_profiles(torch_path, VALID_ROWS)
        torch = Torch.load(torch_path)

        assert torch.profile is None

    def test_missing_threshold_fails_validation(self, torch_path):
        bad_cfg = dict(VALID_OPERON_CFG)
        bad_cfg["identity_thresholds"] = {"1a": 0.983}  # no default, no "2"
        _write_metadata(torch_path, {"typing_model": "operon", "operon": bad_cfg})
        _write_profiles(torch_path, VALID_ROWS)

        with pytest.raises(OperonConfigError):
            Torch.load(torch_path)

    def test_generalized_class_without_residue_rule_fails(self, torch_path):
        bad_cfg = dict(VALID_OPERON_CFG)
        bad_cfg["residue_rules"] = []
        _write_metadata(torch_path, {"typing_model": "operon", "operon": bad_cfg})
        _write_profiles(torch_path, VALID_ROWS)

        with pytest.raises(OperonConfigError):
            Torch.load(torch_path)
