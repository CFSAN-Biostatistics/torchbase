"""ShigaTyper converter — local FASTA files to torch format.

ShigaTyper identifies Shigella serotypes from presence/absence of marker
genes (species markers `ipaH_c`/`ipaB`/`cadA`, plus per-serotype O-antigen
`wzx`/`wzy` and accessory loci mapped in `docs/adr/0003-profile-matching-
value-domain-agnostic.md`). If a profiles TSV is provided it is used as-is;
otherwise this converter writes the canonical table below, a mechanical
transcription of upstream's checkpoint decision cascade (see
`_build_canonical_profiles`), so a torch built with `--download` is loadable
*and* typeable, not just reference data pending a typing model.

torchtools convert shigatyper --download --output torches/

Known gap: the boydii serotype 6 vs 10 split needs raw-read pileup depth at
one bp junction (upstream's own fallback when that data is unavailable is
the same joint call this converter always makes); see
`_build_canonical_profiles`'s docstring.
"""

import csv
import json
import shutil
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, IO, Optional

import toml

from torchbase.conversions import parse_fasta_records
from torchbase.quality.kmer_analysis import analyze_locus
from torchbase.conversions.log import get_logger, TRACE

_log = get_logger("shigatyper")


TAXA = ["Shigella"]

_SHIGATYPER_RAW = "https://raw.githubusercontent.com/CFSAN-Biostatistics/shigatyper/master"

# Every accessory/species-marker locus outside wzx/wzy, in profiles.tsv
# column order. `wbaM` is cataloged (a real reference sequence exists) but
# never discriminates a row -- see the boydii 6/10 gap above.
ACCESSORY_LOCI = [
    "ipaH_c", "cadA", "ipaB", "EclacY", "Ss_methylase", "Sd1_rfp",
    "heparinase", "wbaM", "gtrI", "gtrIC", "gtrII", "gtrIV", "gtrV",
    "gtrX", "Oac", "Oac1b", "Xv",
]
LOCUS_ORDER = ["wzx", "wzy"] + ACCESSORY_LOCI

# The nine flexneri accessory loci whose presence/absence combination (never
# wzx/wzy diversity) distinguishes flexneri serotypes; see _SF_DIC.
_FLEXNERI_LOCI = ["gtrI", "gtrII", "gtrIV", "gtrV", "gtrX", "gtrIC", "Oac", "Oac1b", "Xv"]

# A single-allele locus's one FASTA record is always named "{locus}_1" (see
# download_sources), so calls_from_presence reports allele_id="1" on a hit --
# not a generic sentinel. Matching that here (rather than some other token)
# means "cadA: 1" in profiles.tsv, exactly the allele-ID convention every
# other profiles.tsv in torchbase already uses; ABSENT has no real allele to
# name and stays a true sentinel.
PRESENT, ABSENT, WILD = "1", "absent", "?"

DOWNLOAD_SOURCES = {
    "repo": "https://github.com/CFSAN-Biostatistics/shigatyper",
    # Upstream consolidated its per-locus FASTA files (wzx.fasta, wzy.fasta,
    # fliC.fasta -- none of which exist any more) into one multi-marker
    # reference, headed "{variant}_{gene}" (e.g. "Sb11_wzx"). download_sources
    # below regroups that by GENE: every wzx-suffixed record becomes one
    # allele of the single "wzx" locus (header renamed "wzx_{variant}", the
    # convention torchbase.allele_calls.extract_locus_and_allele expects),
    # likewise wzy; every other marker is already its own single-allele
    # locus and is left as-is.
    "reference": (
        "ShigellaRef5.fasta",
        f"{_SHIGATYPER_RAW}/shigatyper/resources/ShigellaRef5.fasta",
    ),
}


def download_sources(dest_dir: Path) -> dict:
    """Download ShigaTyper's consolidated reference and regroup it by gene.

    Returns {'sequences': [list of open file handles]}: one `wzx.fasta`
    (every `{variant}_wzx` record, renamed `wzx_{variant}`), one `wzy.fasta`
    (likewise), and one file per remaining singleton marker.
    """
    from torchbase.conversions import fetch_file

    dest_dir = Path(dest_dir)
    filename, url = DOWNLOAD_SOURCES["reference"]
    combined = fetch_file(url, dest_dir / filename)

    markers_dir = dest_dir / "markers"
    markers_dir.mkdir(exist_ok=True)

    genes = {"wzx": [], "wzy": []}  # type: Dict[str, list]
    singletons = []
    for header, sequence in parse_fasta_records(combined.read_text(encoding="utf-8")):
        if header.endswith("_wzx"):
            genes["wzx"].append((header[: -len("_wzx")], sequence))
        elif header.endswith("_wzy"):
            genes["wzy"].append((header[: -len("_wzy")], sequence))
        else:
            singletons.append((header, sequence))

    sequence_files = []
    for gene, variants in genes.items():
        gene_path = markers_dir / f"{gene}.fasta"
        with open(gene_path, "w", encoding="utf-8") as f:
            for variant, sequence in variants:
                f.write(f">{gene}_{variant}\n{sequence}\n")
        sequence_files.append(open(gene_path, encoding="utf-8"))

    for marker, sequence in singletons:
        marker_path = markers_dir / f"{marker}.fasta"
        # "_1" explicit allele id: several marker names (ipaH_c,
        # Ss_methylase, Sd1_rfp) already contain an underscore, which would
        # otherwise collide with extract_locus_and_allele's "split at the
        # last underscore" convention (e.g. "ipaH_c" would misparse as
        # locus "ipaH", allele "c"). Applied uniformly for one convention.
        marker_path.write_text(f">{marker}_1\n{sequence}\n", encoding="utf-8")
        sequence_files.append(open(marker_path, encoding="utf-8"))

    _log.info(
        "regrouped %s into %d loci (wzx: %d alleles, wzy: %d alleles, %d singleton markers)",
        filename, len(sequence_files), len(genes["wzx"]), len(genes["wzy"]), len(singletons),
    )
    return {"sequences": sequence_files}


# --------------------------------------------------------------------------
# Canonical serotype table: a faithful, mechanical transcription of
# upstream's checkpoint decision cascade (shigatyper.run(), "## 4. Shigella
# serotype prediction" -- https://github.com/CFSAN-Biostatistics/shigatyper)
# into profiles.tsv rows. torchbase.allele_calls.calls_from_presence's
# present/absent/ambiguous tokens plus torchbase.profile_match's
# first-matching-row-wins semantics reproduce the cascade exactly, given the
# same row priority upstream's if/elif chain uses -- with one documented
# exception: upstream's boydii 6 vs 10 split reads raw-read pileup depth at
# one bp junction (wbaM:252-253) that torchbase's calling layer, built on
# per-locus identity/coverage, does not have. Upstream's own fallback when
# that depth data is unavailable is the ambiguous joint call ("boydii
# serotype 6 or 10"); every boydii-6-shaped profile gets that joint call
# here, since depth data is *never* available to this calling layer.
# --------------------------------------------------------------------------

# Every boydii-numbered wzx variant (Sb10 does not exist -- boydii 10 shares
# Sb6's wzx, disambiguated only by the depth check above).
_ALL_SB_NUMS = [1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 12, 13, 14, 15, 16, 17, 18, 19]
_GENERIC_SB_NUMS = [n for n in _ALL_SB_NUMS if n not in (1, 6)]  # Sb1, Sb6: special-cased
_ALL_SD_NUMS = list(range(1, 16))

# Upstream's SfDic, verbatim (the wzx == "Sf_wzx" branch): a bare list is
# "this exact accessory-gene combination names this serotype" (upstream's
# `Targets == Hits`); a tuple of lists is "either combination does"
# (upstream's `Hits in Targets`, tried against each list in turn).
_SF_DIC = {
    "Shigella flexneri Yv": (["Xv"],),
    "Shigella flexneri serotype 1a": (["gtrI"],),
    "Shigella flexneri serotype 1b": (["gtrI", "Oac1b"],),
    "Shigella flexneri serotype 2a": (["gtrII"],),
    "Shigella flexneri 2av": (["gtrII", "Xv"],),
    "Shigella flexneri serotype 2b": (["gtrII", "gtrX"],),
    "Shigella flexneri serotype 3a": (["gtrX", "Oac"],),
    "Shigella flexneri serotype 3b": (["Oac"], ["Oac1b"]),
    "Shigella flexneri serotype 4a": (["gtrIV"],),
    "Shigella flexneri serotype 4av": (["gtrIV", "Xv"],),
    "Shigella flexneri serotype 4b": (["gtrIV", "Oac"], ["gtrIV", "Oac1b"]),
    "Shigella flexneri 4bv": (["gtrIV", "Oac", "Xv"],),
    "Shigella flexneri serotype 5a": (["gtrV", "Oac"], ["gtrV"]),
    "Shigella flexneri serotype 5b": (["gtrV", "gtrX", "Oac"], ["gtrV", "gtrX"]),
    "Shigella flexneri serotype X": (["gtrX"],),
    "Shigella flexneri serotype Xv (4c)": (["gtrX", "Xv"],),
    "Shigella flexneri serotype 1c (7a)": (["gtrI", "gtrIC"],),
    "Shigella flexneri serotype 7b": (["gtrI", "gtrIC", "Oac1b"],),
}


def _row(profile_id: str, wzx: str = WILD, wzy: str = WILD, **locus_values: str) -> Dict[str, str]:
    row = {"Serotype": profile_id, "wzx": wzx, "wzy": wzy}
    for locus in ACCESSORY_LOCI:
        row[locus] = locus_values.pop(locus, WILD)
    if locus_values:
        raise ValueError(f"unknown locus column(s): {sorted(locus_values)}")
    return row


def _flexneri_row(profile_id: str, present_loci) -> Dict[str, str]:
    """One SfDic row: wzx=Sf, ipaH_c present, cadA absent, exact gtr*/Oac*/Xv set."""
    values = {locus: (PRESENT if locus in present_loci else ABSENT) for locus in _FLEXNERI_LOCI}
    return _row(profile_id, wzx="Sf", ipaH_c=PRESENT, cadA=ABSENT, **values)


def _build_canonical_profiles() -> List[Dict[str, str]]:
    """Build the canonical profiles.tsv rows, in upstream's own priority order."""
    rows = []

    # Species checkpoint exception: ipaH_c- + Sb13_wzx -- upstream stops
    # here, before any of the checks below run.
    rows.append(_row(
        "Shigella boydii serotype 13 (no longer classified as Shigella)",
        wzx="Sb13", ipaH_c=ABSENT,
    ))

    # EclacY+ is an EIEC signal unless the strain is boydii 9 or 15 (upstream's
    # explicit carve-out); those two rows must out-rank the EIEC catch-all.
    rows.append(_row("Shigella boydii serotype 9", wzx="Sb9", EclacY=PRESENT, ipaH_c=PRESENT))
    rows.append(_row("Shigella boydii serotype 15", wzx="Sb15", EclacY=PRESENT, ipaH_c=PRESENT))
    rows.append(_row("EIEC", EclacY=PRESENT, ipaH_c=PRESENT))

    # cadA+ branch: Sonnei, dysenteriae 1/8, or boydii 11 by name; anything
    # else cadA+ is EIEC (this out-ranks the generic numeric fallback below,
    # so cadA+/wzx=Sb12 -- say -- is EIEC, not "boydii serotype 12").
    rows.append(_row("Shigella sonnei form II",
                      cadA=PRESENT, Ss_methylase=PRESENT, wzx=ABSENT, ipaB=ABSENT, ipaH_c=PRESENT))
    rows.append(_row("Shigella sonnei (low levels of form I)",
                      cadA=PRESENT, Ss_methylase=PRESENT, wzx=ABSENT, ipaB=PRESENT, ipaH_c=PRESENT))
    rows.append(_row("Shigella sonnei, form I",
                      cadA=PRESENT, Ss_methylase=PRESENT, wzx="Ss", ipaH_c=PRESENT))
    rows.append(_row("EIEC", cadA=PRESENT, Ss_methylase=PRESENT, ipaH_c=PRESENT))
    rows.append(_row("Shigella dysenteriae serotype 1",
                      cadA=PRESENT, wzx="Sd1", Sd1_rfp=PRESENT, ipaH_c=PRESENT))
    rows.append(_row("Shigella dysenteriae serotype 1, rfp- (phenotypically negative)",
                      cadA=PRESENT, wzx="Sd1", Sd1_rfp=ABSENT, ipaH_c=PRESENT))
    rows.append(_row("Shigella dysenteriae serotype 8", cadA=PRESENT, wzx="Sd8", ipaH_c=PRESENT))
    rows.append(_row("Shigella boydii serotype 11", cadA=PRESENT, wzx="Sb11", ipaH_c=PRESENT))
    rows.append(_row("EIEC", cadA=PRESENT, ipaH_c=PRESENT))

    # cadA- branch: boydii 6/10 (undecidable, see module docstring), then
    # heparinase-disambiguated boydii 1/20, provisional serotypes, and the
    # generic numeric fallback for every other boydii/dysenteriae variant.
    rows.append(_row("Shigella boydii serotype 6 or 10", wzx="Sb6", ipaH_c=PRESENT))
    rows.append(_row("Shigella boydii serotype 20", wzx="Sb1", heparinase=PRESENT, ipaH_c=PRESENT))
    rows.append(_row("Shigella boydii serotype 1", wzx="Sb1", heparinase=ABSENT, ipaH_c=PRESENT))
    rows.append(_row("Shigella boydii Provisional serotype E1621-54", wzx="SbProv", ipaH_c=PRESENT))
    rows.append(_row("Shigella dysenteriae Provisional serotype 96-265", wzx="SdProv", ipaH_c=PRESENT))
    rows.append(_row("Shigella dysenteriae Provisional serotype E670-74", wzx="SdProvE", ipaH_c=PRESENT))
    for n in _GENERIC_SB_NUMS:
        rows.append(_row(f"Shigella boydii serotype {n}", wzx=f"Sb{n}", ipaH_c=PRESENT))
    for n in _ALL_SD_NUMS:
        rows.append(_row(f"Shigella dysenteriae serotype {n}", wzx=f"Sd{n}", ipaH_c=PRESENT))

    # flexneri 6 (its own wzx/wzy variant); every other flexneri serotype
    # shares base wzx=Sf/wzy=Sf and is distinguished entirely by which
    # gtr*/Oac*/Xv accessory genes it carries.
    rows.append(_row("Shigella flexneri serotype 6", wzx="Sf6", ipaH_c=PRESENT))
    rows.append(_flexneri_row("Shigella flexneri serotype Y", present_loci=()))
    for name, alternatives in _SF_DIC.items():
        for combo in alternatives:
            rows.append(_flexneri_row(name, present_loci=combo))

    return rows


def convert_local(
    sequence_files: List[IO],
    profiles_file: Optional[IO] = None,
    output_path: str = ".",
    namespace: str = "hfp",
    name: str = "shigatyper",
    version: str = "1.0.0",
    kmer_size: int = 13,
    overlap_threshold: float = 0.90,
    duplicate_threshold: float = 0.95,
) -> str:
    """Convert local ShigaTyper FASTA files to torch format.

    Args:
        sequence_files: Open file handles for locus FASTA files -- `wzx.fasta`
            and `wzy.fasta` (each holding every serotype's allele) plus one
            file per singleton marker. See `download_sources`.
        profiles_file: Optional open file handle for a serotype profiles TSV,
            overriding the canonical table this converter would otherwise
            write (see module docstring). Columns: Serotype, then one column
            per locus in `LOCUS_ORDER`.
        output_path: Directory in which to create the torch.
        namespace: Torch namespace. Defaults to "hfp": a torch's namespace
            names the authority for the data, and the ShigaTyper database is
            an FDA Human Foods Program product
            (github.com/CFSAN-Biostatistics/ShigaTyper, from the center's
            former name). Compare "ncbi" for the StxTyper conversion and
            "pubmlst" for PubMLST schemes.
        name: Torch name (default: "shigatyper").
        version: Torch version string.
        kmer_size: K-mer size for quality analysis.
        overlap_threshold: Overlap similarity threshold.
        duplicate_threshold: Duplicate similarity threshold.

    Returns:
        Path to the created torch directory.
    """
    output_path = Path(output_path)
    _log.info("Starting shigatyper conversion → %s/%s %s", namespace, name, version)
    torch_dir = output_path / namespace / name / f"{version}.torch"
    torch_dir.mkdir(parents=True, exist_ok=True)

    resources_dir = torch_dir / "_resources"
    resources_dir.mkdir(exist_ok=True)

    locus_names = []
    for fasta_fh in sequence_files:
        src = Path(fasta_fh.name)
        dest = resources_dir / src.name
        shutil.copy2(src, dest)
        locus_names.append(src.stem)
        _log.debug("  copying %s → %s (%d bytes)", src.name, dest, dest.stat().st_size)

    _log.info("  %d locus file(s): %s", len(locus_names), ", ".join(locus_names))

    profiles_dest = torch_dir / "profiles.tsv"
    if profiles_file is not None:
        shutil.copy2(profiles_file.name, profiles_dest)
        profiles_file.seek(0)
        rows = list(csv.reader(profiles_file, delimiter="\t"))
        profile_count = max(0, len(rows) - 1)
    elif set(locus_names) == set(LOCUS_ORDER):
        _write_canonical_profiles(profiles_dest)
        profile_count = len(_build_canonical_profiles())
    else:
        _write_stub_profiles(profiles_dest, locus_names)
        profile_count = 0

    _log.info("  profiles: %d rows", profile_count)
    _log.debug("  profiles written to %s (%d rows)", profiles_dest, profile_count)

    _log.info("  running k-mer quality analysis on %d loci (k=%d)", len(locus_names), kmer_size)
    quality_results = _run_quality_analysis(
        resources_dir, kmer_size, overlap_threshold, duplicate_threshold
    )
    _log.info("  quality summary: %d suspect loci, %d duplicate pairs",
              quality_results.get("suspect_loci", 0), len(quality_results.get("duplicate_pairs", [])))

    now = datetime.now(timezone.utc).isoformat()
    metadata = {
        "namespace": namespace,
        "name": name,
        "version": version,
        "version_info": {"strategy": "snapshot", "timestamp": now},
        "typing": {
            "method": "serotyping",
            "scheme": "Shigella O:H antigen (presence/absence marker set)",
            "loci_count": len(locus_names),
            "profiles_count": profile_count,
            "calling_mode": "presence_absence",
            "id_column": "Serotype",
        },
        "description": {
            "short": "ShigaTyper Shigella serotyping torch",
            "long": (
                "Species and O-antigen/accessory marker sequences from "
                "CFSAN-Biostatistics/shigatyper's ShigellaRef5.fasta: wzx and "
                "wzy are real multi-allele loci (one allele per serotype's "
                "variant), the remaining markers are single-allele "
                "presence/absence loci. profiles.tsv is a mechanical "
                "transcription of ShigaTyper's checkpoint decision cascade "
                "into rows, matched by presence/absence rather than allele "
                "identity (torchbase.allele_calls.calls_from_presence); see "
                "docs/adr/0003-profile-matching-value-domain-agnostic.md. "
                "One gap is intentionally not fabricated: the boydii 6 vs 10 "
                "split needs raw-read pileup depth at one bp junction, which "
                "this calling layer does not have, so both collapse to the "
                "joint call upstream itself falls back to without that "
                "data -- see torchbase/conversions/shigatyper.py."
            ),
            "taxa": TAXA,
        },
        "data_quality": {
            "kmer_size": kmer_size,
            "overlap_threshold": overlap_threshold,
            "duplicate_threshold": duplicate_threshold,
        },
        "manifest": {
            "profiles": "profiles.tsv",
            "resources": [
                f.name for f in sorted(resources_dir.iterdir())
                if f.is_file() and not f.name.startswith(".")
            ],
        },
    }
    with open(torch_dir / "metadata.toml", "w", encoding="utf-8") as f:
        toml.dump(metadata, f)
    _log.debug("  metadata.toml written")

    quality_report = _build_quality_report(locus_names, quality_results, kmer_size, overlap_threshold)
    with open(torch_dir / "quality.json", "w", encoding="utf-8") as f:
        json.dump(quality_report, f, indent=2)
    _log.debug("  quality.json written")

    _log.info("Torch written: %s", torch_dir)
    return str(torch_dir)


def _write_stub_profiles(profiles_path: Path, locus_names: List[str]) -> None:
    with open(profiles_path, "w", newline="", encoding="utf-8") as f:
        csv.writer(f, delimiter="\t").writerow(["Serotype"] + locus_names)


def _write_canonical_profiles(profiles_path: Path) -> None:
    with open(profiles_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=["Serotype"] + LOCUS_ORDER, delimiter="\t")
        writer.writeheader()
        writer.writerows(_build_canonical_profiles())


def _build_quality_report(locus_names, quality_results, kmer_size, overlap_threshold):
    return {
        "loci": {
            locus: {
                "allele_count": quality_results["loci_results"].get(locus, {}).get("allele_count", 0),
                "kmer_size": kmer_size,
                "similarity_stats": quality_results["loci_results"].get(locus, {}).get("similarity_stats", {}),
                "threshold": quality_results["loci_results"].get(locus, {}).get("threshold", overlap_threshold),
                "threshold_method": "gap_detection",
                "suspect_pairs": quality_results["loci_results"].get(locus, {}).get("suspect_pairs", []),
            }
            for locus in locus_names
        },
        "summary": {
            "total_loci": len(locus_names),
            "suspect_loci": quality_results.get("suspect_loci", 0),
            "suspect_alleles": len(quality_results.get("duplicate_pairs", [])),
        },
    }


def _run_quality_analysis(resources_dir, kmer_size, overlap_threshold, duplicate_threshold):
    results = {
        "total_loci": 0,
        "suspect_loci": 0,
        "similar_pairs": [],
        "duplicate_pairs": [],
        "loci_results": {},
    }
    for fasta_file in sorted(resources_dir.glob("*.fasta")):
        locus_name = fasta_file.stem
        report = analyze_locus(
            fasta_file,
            k_size=kmer_size,
            overlap_threshold=overlap_threshold * 100,
            duplicate_threshold=duplicate_threshold * 100,
        )
        results["total_loci"] += 1
        entry = {
            "suspect_pairs": report.suspect_pairs,
            "threshold": report.threshold,
            "allele_count": len(report.similarities) + 1 if report.similarities else 0,
            "similarity_stats": {},
        }
        if report.similarities:
            sims = sorted(report.similarities)
            n = len(sims)
            entry["similarity_stats"] = {
                "min": sims[0],
                "median": sims[n // 2],
                "percentile_99": sims[min(int(n * 0.99), n - 1)],
            }
        results["loci_results"][locus_name] = entry
        _log.debug("    %s: %d alleles, %d suspect pairs (threshold=%.3f)",
                   locus_name, entry["allele_count"], len(report.suspect_pairs), report.threshold)
        for pair in report.suspect_pairs:
            _log.log(TRACE, "      suspect: %s ↔ %s  sim=%.4f  type=%s",
                     pair.get("allele1"), pair.get("allele2"),
                     pair.get("similarity", 0), pair.get("issue_type", "?"))
        if report.suspect_pairs:
            results["suspect_loci"] += 1
            for pair in report.suspect_pairs:
                bucket = "duplicate_pairs" if pair.get("issue_type") == "duplicate" else "similar_pairs"
                results[bucket].append({"locus": locus_name, **pair})
    return results
