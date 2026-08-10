"""Typing runs: the package layer's orchestration of a torch against a sample.

This is the layer the CLI is a thin wrapper around, and the layer that owns
*meaning*. It prepares inputs, dispatches compute to the WDL layer through
`torchbase.runner`, interprets what comes back, and returns the profile.

The division of labour, per typing model:

    operon    WDL: tblastn protein search
              here: normalize, stitch frameshifts, pair subunits, score, call

    allelic   WDL: sourmash sketch/compare, minimap2 alignment
              here: quality filtering, allele calls from similarity or from
              alignment, the alignment-fallback decision, profile association

Nothing below this layer decides anything: a task runs a tool and writes a
file. Nothing above it knows a workflow engine exists.
"""

import tempfile
from pathlib import Path
from typing import Dict, List, Optional, Sequence

from torchbase import allele_calls, operon, profile_match, quality_filters, reports, runner

BUILTIN_WORKFLOWS = Path(__file__).parent / "workflows" / "builtin"

# minimap2 presets by input type. Assemblies align as near-identical contigs;
# short reads need the sr preset.
ALIGNMENT_PRESETS = {
    "contigs": "asm5",
    "reads": "sr",
    "paired": "sr",
    "interlaced": "sr",
    "longreads": "asm20",
}

# MinHash screen resolution. `scaled` is sourmash's sampling rate: it keeps
# roughly one hash per `scaled` k-mers, so the 1000 the WDL task defaults to
# leaves a 500 bp MLST allele with *zero* hashes and a similarity matrix of
# pure zeros. Allele databases are small, so sample everything: an exhaustive
# screen of a few thousand alleles is still far cheaper than aligning.
DEFAULT_KSIZE = 31
DEFAULT_SCALED = 1


class TypingError(RuntimeError):
    """A typing run could not be completed."""


def workflow_path(name: str) -> Path:
    path = BUILTIN_WORKFLOWS / name
    if not path.exists():
        raise TypingError("built-in workflow not found: {}".format(path))
    return path


def type_operon(
    torch,
    contigs: str,
    engine: str = "miniwdl",
    workflow: Optional[Path] = None,
) -> List[dict]:
    """Type an assembly against an operon torch (docs/operon-strategy-plan.md).

    The workflow searches; every decision after that is made here.
    """
    if not torch.operon_config:
        raise TypingError("Operon torch is missing its [operon] config block.")
    reference_rel = torch.operon_config.get("reference", {}).get("file")
    if not reference_rel:
        raise TypingError("operon.reference.file not set in metadata.toml.")
    reference = torch.path / reference_rel
    if not reference.exists():
        raise TypingError("Operon reference file not found: {}".format(reference))

    outputs = runner.run_workflow(
        workflow or workflow_path("operon_typing.wdl"),
        {"contigs": contigs, "subunit_reference": str(reference)},
        engine=engine,
    )
    return operon.type_assembly(
        runner.require_file(outputs, "hits"),
        reference,
        torch.operon_config,
        profile_rows=torch.operon_profiles or (),
        scheme=torch.operon_config.get("scheme") or torch.path.parent.name,
    )


def type_allelic(
    torch,
    query: str,
    input_type: str = "contigs",
    strategy: str = "balanced",
    quality_json: Optional[str] = None,
    exclude_suspect_alleles: bool = False,
    exclude_suspect_loci: bool = False,
    exclude_suspect_profiles: bool = False,
    confidence_threshold: float = reports.DEFAULT_CONFIDENCE_THRESHOLD,
    identity_threshold: float = allele_calls.DEFAULT_IDENTITY_THRESHOLD,
    engine: str = "miniwdl",
    scratch: Optional[Path] = None,
) -> dict:
    """Type a sample against an allelic torch and return its profile.

    `strategy` is the speed/accuracy tier:

        fast       MinHash screen only
        balanced   MinHash screen, then alignment for the loci it could not
                   call confidently
        sensitive  alignment only

    `balanced`'s fallback is dispatched from here, which is what makes it a
    genuine middle tier: the decision needs the screen's calls, and calls are
    this layer's business.
    """
    if strategy not in ("fast", "balanced", "sensitive"):
        raise TypingError("unknown strategy: {}".format(strategy))

    workspace = Path(scratch or tempfile.mkdtemp(prefix="torchbase-run-"))
    allele_fasta, profiles_table = torch.get_unified_files()

    try:
        query = _prepare_query(query, input_type, workspace)
        filtered_fasta, exclusions = _prepare_alleles(
            allele_fasta,
            workspace,
            quality_json=quality_json,
            exclude_alleles=exclude_suspect_alleles,
            exclude_loci=exclude_suspect_loci,
            exclude_profiles=exclude_suspect_profiles,
        )

        calls = {}
        alignment_used = False
        if strategy in ("fast", "balanced"):
            calls = _screen(
                query, filtered_fasta, strategy, confidence_threshold, engine
            )
        if strategy == "sensitive" or (
            strategy == "balanced" and reports.needs_alignment(calls, confidence_threshold)
        ):
            aligned = _align(
                query, filtered_fasta, input_type, identity_threshold, engine
            )
            alignment_used = True
            calls = (
                aligned
                if strategy == "sensitive"
                else reports.merge_calls(calls, aligned, True, confidence_threshold)
            )

        profiles, loci_order = profile_match.load_profiles(profiles_table)
        profile = profile_match.build_profile_record(
            calls,
            profiles,
            loci_order,
            torch.path.parent.name,
            strategy=strategy,
            alignment_used=alignment_used,
        )
        return reports.add_exclusion_metadata(profile, exclusions)
    finally:
        for path in (allele_fasta, profiles_table):
            if path is not None and Path(path).exists():
                Path(path).unlink()


def _prepare_query(query: str, input_type: str, workspace: Path) -> str:
    """Drop low-depth query sequences before anything expensive runs.

    A no-op for assemblies, matching the WDL task this replaces: contig depth is
    the assembler's business, and k-mer depth of a contig means nothing.
    """
    if input_type == quality_filters.PASSTHROUGH_INPUT_TYPE:
        return query
    records, _stats = quality_filters.depth_filter(query, input_type=input_type)
    filtered = workspace / "filtered_query.fasta"
    quality_filters.write_sequence_fasta(records, filtered)
    return str(filtered)


def _prepare_alleles(
    allele_fasta,
    workspace: Path,
    quality_json=None,
    exclude_alleles=False,
    exclude_loci=False,
    exclude_profiles=False,
):
    """Apply quality filtering before anything expensive runs.

    Note the flags are a precedence chain rather than independent switches
    (loci beats alleles, profiles beats both) — preserved from the WDL task
    this replaces; see torchbase/quality_filters.py.
    """
    entries, exclusions = quality_filters.filter_alleles(
        allele_fasta,
        quality=quality_filters.load_quality(quality_json),
        exclude_alleles=exclude_alleles,
        exclude_loci=exclude_loci,
        exclude_profiles=exclude_profiles,
    )
    filtered = workspace / "filtered_alleles.fasta"
    quality_filters.write_fasta(entries, filtered)
    # The WDL task emitted excluded_loci in set order, which varied per run.
    exclusions = dict(exclusions)
    if isinstance(exclusions.get("excluded_loci"), list):
        exclusions["excluded_loci"] = sorted(exclusions["excluded_loci"])
    return filtered, exclusions


def _screen(
    query: str,
    allele_fasta: Path,
    strategy: str,
    confidence_threshold: float,
    engine: str,
    ksize: int = DEFAULT_KSIZE,
    scaled: int = DEFAULT_SCALED,
) -> Dict[str, dict]:
    outputs = runner.run_workflow(
        workflow_path("minhash_screen.wdl"),
        {
            "query_sequences": query,
            "allele_fasta": str(allele_fasta),
            "ksize": ksize,
            "scaled": scaled,
        },
        engine=engine,
    )
    matrix = runner.require_file(outputs, "similarity_matrix")
    if strategy == "fast":
        # `fast` historically used a stricter, fixed similarity bar than
        # `balanced`; kept until the two are reconciled deliberately.
        return allele_calls.call_alleles(matrix, query, allele_fasta)
    return allele_calls.call_alleles_minhash(
        matrix, query, allele_fasta, confidence_threshold=confidence_threshold
    )


def _align(
    query: str,
    allele_fasta: Path,
    input_type: str,
    identity_threshold: float,
    engine: str,
) -> Dict[str, dict]:
    outputs = runner.run_workflow(
        workflow_path("allele_alignment.wdl"),
        {
            "query_sequences": query,
            "allele_fasta": str(allele_fasta),
            "preset": ALIGNMENT_PRESETS.get(input_type, "asm5"),
        },
        engine=engine,
    )
    return allele_calls.calls_from_alignment(
        runner.require_file(outputs, "alignments"),
        allele_calls.parse_fasta_dict(query),
        identity_threshold=identity_threshold,
    )
