"""StxTyper converter — NCBI's stx.prot reference set to an operon torch.

StxTyper (https://github.com/ncbi/stxtyper) types the Shiga toxin (stx)
operon via protein-space (tblastn) search, syntenic pairing of the A and B
subunits, per-class combined-identity thresholds, and residue-table
tie-breaking for the stx2a/2c/2d generalized class. This converter
transcribes its bundled `stx.prot` reference set and decision rules
(source-derived from `stxtyper.cpp`, see docs/operon-strategy-plan.md and
the session handoff) into a `typing_model = "operon"` torch (§2, §3.1).

Usage:
    torchtools convert stxtyper --download
    torchtools convert stxtyper path/to/stx.prot --stxtyper-version 1.0.45
"""

from datetime import datetime, timezone
from pathlib import Path
from typing import IO, Optional

import toml

from torchbase.conversions.log import get_logger

_log = get_logger("stxtyper")

TAXA = ["Escherichia coli", "Shigella"]

_STXTYPER_RAW = "https://raw.githubusercontent.com/ncbi/stxtyper/master"

DOWNLOAD_SOURCES = {
    "repo": "https://github.com/ncbi/stxtyper",
    "stx_prot": f"{_STXTYPER_RAW}/stx.prot",
    "version": f"{_STXTYPER_RAW}/version.txt",
}

# Per-class combined-identity thresholds (stxtyper.cpp:894-911).
IDENTITY_THRESHOLDS = {
    "default": 0.98,
    "1a": 0.983,
    "1c": 0.983,
    "1d": 0.983,
    "1e": 0.983,
    "2k": 0.985,
    "2l": 0.985,
}

# stx2a/2c/2d are too close for whole-sequence identity to separate; they
# collapse to generalized class "2" and are resolved by residue table.
GENERALIZED_CLASSES = {"2": ["2a", "2c", "2d"]}

# Residue decision table (stxtyper.cpp:536-556). `index` is 0-based into the
# reference protein sequence — StxTyper's own A312/A318/B34 documentation is
# 1-based, so the source coordinates are shifted by -1 here.
RESIDUE_RULES = [
    {
        "class": "2",
        "positions": [
            {"subunit": "A", "index": 311},  # A312 (1-based)
            {"subunit": "A", "index": 317},  # A318 (1-based)
            {"subunit": "B", "index": 33},   # B34  (1-based)
        ],
        "table": [
            {"call": "2a", "residues": [["F", "S"], ["K", "E"], ["D"]]},
            {"call": "2c", "residues": [["F"], ["K", "E"], ["N"]]},
            {"call": "2d", "residues": [["S"], ["E"], ["N"]]},
        ],
        "fallback": "2",
    }
]

INTERGENIC_MAX = 36  # bp (stxtyper.cpp:147: "max intergenic region in the reference set + 2")
INTERGENIC_RELAX_FACTOR = 2

HEADER_FORMAT = "accession|subunit_role|reference_subtype|class"


def download_sources(dest_dir: Path) -> dict:
    """Download stx.prot and version.txt from ncbi/stxtyper to dest_dir.

    Returns {'stx_prot': open file handle, 'stxtyper_version': str}.
    """
    from torchbase.conversions import fetch_file

    dest_dir = Path(dest_dir)
    stx_prot_path = fetch_file(DOWNLOAD_SOURCES["stx_prot"], dest_dir / "stx.prot")
    version_path = fetch_file(DOWNLOAD_SOURCES["version"], dest_dir / "version.txt")
    stxtyper_version = version_path.read_text().strip()
    return {"stx_prot": open(stx_prot_path), "stxtyper_version": stxtyper_version}


def _parse_stx_prot(fh: IO) -> list:
    """Parse stx.prot into a list of {accession, subunit_role, type, class,
    sequence} records.

    Headers are `>accession|famId|subclass`, e.g.
    `AAS07582.1|stxA2c|stxA2` — famId is "stx" + subunit letter + 2-char
    type (e.g. "stxA2c"); subclass is famId with the type collapsed for
    classes that need residue-table resolution (e.g. "stxA2").
    """
    records = []
    header = None
    seq_lines = []

    def flush():
        if header is None:
            return
        accession, fam_id, subclass = header.split("|")
        subunit_role = fam_id[3]  # "stx" + role + type
        type_code = fam_id[4:]    # e.g. "2c", "1a"
        class_label = subclass[4:]  # e.g. "2", "1a"
        records.append({
            "accession": accession,
            "subunit_role": subunit_role,
            "reference_subtype": fam_id,
            "type_code": type_code,
            "class_label": class_label,
            "sequence": "".join(seq_lines),
        })

    for line in fh:
        line = line.rstrip("\n")
        if not line:
            continue
        if line.startswith(">"):
            flush()
            header = line[1:]
            seq_lines = []
        else:
            seq_lines.append(line)
    flush()
    return records


def _build_profiles(records: list) -> list:
    """Group per-subunit records by operon type into profiles.tsv rows:
    subtype, class, subunit_A, subunit_B."""
    by_type = {}
    for rec in records:
        by_type.setdefault(rec["type_code"], {})[rec["subunit_role"]] = rec

    rows = []
    for type_code in sorted(by_type):
        roles = by_type[type_code]
        subunit_a = roles.get("A")
        subunit_b = roles.get("B")
        class_label = (subunit_a or subunit_b)["class_label"]
        rows.append({
            "subtype": f"stx{type_code}",
            "class": class_label,
            "subunit_A": subunit_a["reference_subtype"] if subunit_a else "",
            "subunit_B": subunit_b["reference_subtype"] if subunit_b else "",
        })
    return rows


def convert_local(
    stx_prot_file: IO,
    output_path: str = ".",
    namespace: str = "ncbi",
    name: str = "stxtyper",
    version: Optional[str] = None,
    stxtyper_version: Optional[str] = None,
) -> str:
    """Convert a local stx.prot file into an operon torch.

    Args:
        stx_prot_file: Open file handle for stx.prot.
        output_path: Directory in which to create the torch.
        namespace: Torch namespace (default: "ncbi").
        name: Torch name (default: "stxtyper").
        version: Torch version string (default: derived from stxtyper_version,
            or "0.0.0" if that isn't known either).
        stxtyper_version: Upstream StxTyper version this torch claims parity
            against (§7 risk 3: "StxTyper is a moving target"). Recorded in
            metadata.toml so parity claims are pinned to a snapshot.

    Returns:
        Path to the created torch directory.
    """
    output_path = Path(output_path)
    records = _parse_stx_prot(stx_prot_file)
    if not records:
        raise ValueError("stx.prot contained no records.")

    resolved_version = version or stxtyper_version or "0.0.0"
    _log.info(
        "Starting stxtyper conversion → %s/%s %s (%d reference sequences)",
        namespace, name, resolved_version, len(records),
    )

    torch_dir = output_path / namespace / name / f"{resolved_version}.torch"
    torch_dir.mkdir(parents=True, exist_ok=True)

    resources_dir = torch_dir / "_resources"
    resources_dir.mkdir(exist_ok=True)

    subunits_path = resources_dir / "subunits.faa"
    with open(subunits_path, "w") as f:
        for rec in records:
            header = "|".join([
                rec["accession"], rec["subunit_role"],
                rec["reference_subtype"], rec["class_label"],
            ])
            f.write(f">{header}\n{rec['sequence']}\n")
    _log.debug("  wrote %s (%d sequences)", subunits_path, len(records))

    profile_rows = _build_profiles(records)
    profiles_path = torch_dir / "profiles.tsv"
    with open(profiles_path, "w") as f:
        f.write("subtype\tclass\tsubunit_A\tsubunit_B\n")
        for row in profile_rows:
            f.write(f"{row['subtype']}\t{row['class']}\t{row['subunit_A']}\t{row['subunit_B']}\n")
    _log.info("  profiles: %d subtypes", len(profile_rows))

    now = datetime.now(timezone.utc).isoformat()
    metadata = {
        "namespace": namespace,
        "name": name,
        "version": resolved_version,
        "version_info": {"strategy": "snapshot", "timestamp": now},
        "typing_model": "operon",
        "typing": {
            "method": "operon",
            "scheme": "stx (Shiga toxin)",
            "subtypes_count": len(profile_rows),
        },
        "description": {
            "short": "StxTyper-parity Shiga toxin (stx) operon torch",
            "long": (
                "Protein-space, synteny-aware typing of the stx operon, "
                "generalized from NCBI StxTyper. Transcribed from "
                f"stx.prot (upstream stxtyper {stxtyper_version or 'unknown'})."
            ),
            "taxa": TAXA,
        },
        "provenance": {
            "source": DOWNLOAD_SOURCES["repo"],
            "stxtyper_version": stxtyper_version or "unknown",
            "license": "Public Domain (NCBI)",
        },
        "manifest": {"profiles": "profiles.tsv"},
        "operon": {
            "subunit_order": ["A", "B"],
            "intergenic_max": INTERGENIC_MAX,
            "intergenic_relax_factor": INTERGENIC_RELAX_FACTOR,
            "require_same_strand": True,
            "require_same_contig": True,
            "reference": {
                "file": "_resources/subunits.faa",
                "header_format": HEADER_FORMAT,
            },
            "identity_thresholds": IDENTITY_THRESHOLDS,
            "generalized_classes": GENERALIZED_CLASSES,
            "residue_rules": RESIDUE_RULES,
        },
    }
    with open(torch_dir / "metadata.toml", "w") as f:
        toml.dump(metadata, f)
    _log.debug("  metadata.toml written")

    _log.info("Torch written: %s", torch_dir)
    return str(torch_dir)
