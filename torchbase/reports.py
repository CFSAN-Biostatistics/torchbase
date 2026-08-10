#!/usr/bin/env python

"""Typing-report interpretation: merging allele calls and annotating results.

Replaces the embedded Python heredocs of these WDL tasks:

* ``merge_allele_calls`` (``workflows/builtin/balanced_typing.wdl``)
* ``check_confidence_for_alignment`` (``workflows/builtin/balanced_typing.wdl``)
* ``add_exclusion_metadata`` (``workflows/builtin/balanced_typing.wdl``,
  ``fast_typing.wdl`` and ``sensitive_typing.wdl`` -- three copies whose
  heredoc bodies are identical apart from comments)

This is interpretation, not compute: deciding which evidence source wins for a
locus, whether the alignment fallback is warranted, and what provenance notes
belong on the returned profile are torchbase policy decisions. They need no
tool binary, so they belong in the importable package layer where they can be
unit-tested without Docker, rather than in three drifting container heredocs.
"""

from typing import Any, Dict, Mapping

#: Default MinHash/alignment confidence threshold, mirroring the
#: ``confidence_threshold`` default of the balanced typing workflow. Single
#: source of truth for the alignment-fallback gate.
DEFAULT_CONFIDENCE_THRESHOLD = 0.85

#: Exclusion fields lifted from a ``filter_alleles`` exclusions document onto
#: the typing result, in the order the WDL heredocs emitted them.
EXCLUSION_FIELDS = (
    "excluded_alleles",
    "excluded_loci",
    "num_excluded_alleles",
    "num_excluded_loci",
)


def locus_needs_alignment(call_data, threshold=DEFAULT_CONFIDENCE_THRESHOLD):
    # type: (Mapping[str, Any], float) -> bool
    """Return True if one locus call is too weak to trust without alignment.

    A ``confidence`` field takes precedence and is read for truthiness: a
    falsey value (``False``, ``0``, ``0.0``, ``None``, ``""``) demands
    alignment. Only when ``confidence`` is absent entirely does ``similarity``
    decide, and then strictly below ``threshold`` demands alignment (so a call
    exactly at the threshold is accepted). A call carrying neither field never
    demands alignment.
    """
    if "confidence" in call_data:
        return not call_data["confidence"]
    if "similarity" in call_data:
        return call_data["similarity"] < threshold
    return False


def needs_alignment(calls, threshold=DEFAULT_CONFIDENCE_THRESHOLD):
    # type: (Mapping[str, Mapping[str, Any]], float) -> bool
    """Return True if any allele call warrants the alignment fallback.

    ``calls`` maps locus name to call data as produced by the MinHash allele
    caller. Empty input means nothing is in doubt, hence False.
    """
    for call_data in calls.values():
        if locus_needs_alignment(call_data, threshold):
            return True
    return False


def alignment_call_preferred(alignment_call, threshold=DEFAULT_CONFIDENCE_THRESHOLD):
    # type: (Mapping[str, Any], float) -> bool
    """Return True if an alignment call is strong enough to beat MinHash.

    An alignment call wins when its ``identity`` is at or above ``threshold``.
    A missing ``identity`` is treated as 0.0 and therefore loses unless the
    threshold is 0.0 itself.
    """
    return alignment_call.get("identity", 0.0) >= threshold


def merge_calls(
    minhash_calls,
    alignment_calls,
    use_alignment,
    threshold=DEFAULT_CONFIDENCE_THRESHOLD,
):
    # type: (Mapping[str, Any], Mapping[str, Any], bool, float) -> Dict[str, Any]
    """Merge MinHash and alignment allele calls into the final call set.

    MinHash calls define the locus set: a locus present only in
    ``alignment_calls`` is dropped, and a locus present in neither simply does
    not appear. When ``use_alignment`` is false the MinHash calls are returned
    unchanged. Otherwise, for each MinHash locus the alignment call wins if it
    exists and passes :func:`alignment_call_preferred`; otherwise the MinHash
    call is kept.
    """
    if not use_alignment:
        return dict(minhash_calls)

    merged = {}  # type: Dict[str, Any]
    for locus, minhash_call in minhash_calls.items():
        if locus in alignment_calls and alignment_call_preferred(
            alignment_calls[locus], threshold
        ):
            merged[locus] = alignment_calls[locus]
        else:
            merged[locus] = minhash_call
    return merged


def add_exclusion_metadata(result, exclusions):
    # type: (Mapping[str, Any], Mapping[str, Any]) -> Dict[str, Any]
    """Return ``result`` with allele-exclusion provenance under ``notes``.

    The four :data:`EXCLUSION_FIELDS` are required; a missing one is a
    malformed exclusions document and raises ``KeyError``, as in the WDL
    heredocs. Neither argument is mutated. An existing ``notes`` mapping keeps
    its position and its other keys; ``notes.exclusions`` is overwritten.
    """
    annotated = dict(result)
    notes = dict(annotated.get("notes", {}))
    notes["exclusions"] = dict(
        (field, exclusions[field]) for field in EXCLUSION_FIELDS
    )
    annotated["notes"] = notes
    return annotated
