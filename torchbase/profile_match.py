"""Allelic profile -> sequence type association.

Replaces the embedded Python heredoc of the WDL task `lookup_profile`
(``torchbase/workflows/builtin/tasks/profile_lookup.wdl``).

This is interpretation, not compute: no sketching, alignment or search happens
here, only the decision of which row of a scheme's profiles table the allele
calls correspond to, and how to label the answer when no row fits. It therefore
belongs in the package layer, where it is importable and unit-testable, rather
than in a container whose only job is to run sourmash/minimap2/BLAST+.

Wildcard semantics preserved from the heredoc: an uncalled locus renders as
``?`` in the query profile, and a ``?`` on either side of a comparison is
skipped rather than counted as a mismatch, so a partially called isolate
matches any row consistent with the loci that *were* called.
"""

import csv
import json
import os
from typing import Any, Dict, List, Optional, Sequence, Tuple

#: Value used for an uncalled locus, on either side of a comparison.
WILDCARD = "?"

#: `profile_id` emitted when no row matched (or the matched row has no ST).
UNKNOWN_PROFILE_ID = "unknown"

#: Statuses the heredoc can produce. NB: the documented vocabulary in
#: CLAUDE.md/docs is ``known | novel_profile | novel_allele``; the task has only
#: ever emitted these two, and never distinguishes a novel allele.
STATUS_KNOWN = "known_profile"
STATUS_NOVEL_PROFILE = "novel_profile"

Profiles = List[Dict[str, str]]


def load_allele_calls(json_path):
    # type: (str) -> Dict[str, Any]
    """Load the per-locus allele calls JSON emitted by minhash or alignment."""
    with open(json_path) as handle:
        return json.load(handle)


def load_profiles(tsv_path, id_column=None):
    # type: (str, Optional[str]) -> Tuple[Profiles, List[str]]
    """Read a profiles TSV into rows plus the locus column order.

    `id_column` names the row identifier column verbatim (e.g. ``"Serotype"``
    for an ECTyper/SeqSero2 torch); unset auto-detects any column named ``ST``
    (any casing) as the identifier, the original behaviour. Either way the
    identifier column is excluded from the returned locus order while staying
    on each row. A header-only table yields no rows but a full locus order.
    """
    profiles = []  # type: Profiles
    loci_order = []  # type: List[str]
    with open(tsv_path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            profiles.append(row)
        if reader.fieldnames:
            if id_column is not None:
                loci_order = [col for col in reader.fieldnames if col != id_column]
            else:
                loci_order = [col for col in reader.fieldnames if col.upper() != "ST"]
    return profiles, loci_order


def st_column(profile, id_column=None):
    # type: (Dict[str, str], Optional[str]) -> Optional[str]
    """Name of the row's identifier column, else None.

    `id_column`, when given, is used verbatim (and must be present on the
    row); unset falls back to any column named ``ST`` (any casing) -- the
    original, `load_profiles`-independent auto-detection every caller used
    before an explicit identifier column existed.
    """
    if id_column is not None:
        return id_column if id_column in profile else None
    for column in profile.keys():
        if column.upper() == "ST":
            return column
    return None


def build_profile_string(allele_calls, loci_order):
    # type: (Dict[str, Any], Sequence[str]) -> str
    """Render the called alleles as a comma-joined profile in table order.

    A locus absent from the calls, or a call carrying no ``allele_id``, becomes
    ``?``. Locus names are matched exactly: a scheme-prefixed call key such as
    ``salmonella_adk`` does not satisfy a bare ``adk`` column.
    """
    parts = []
    for locus in loci_order:
        if locus in allele_calls:
            parts.append(str(allele_calls[locus].get("allele_id", WILDCARD)))
        else:
            parts.append(WILDCARD)
    return ",".join(parts)


def row_profile_string(profile, loci_order):
    # type: (Dict[str, str], Sequence[str]) -> str
    """Render a table row as a comma-joined profile in locus order.

    Loci the row does not carry at all are dropped rather than wildcarded, so a
    row from a differently-shaped table simply fails the length check in
    :func:`match_profile`. Mirroring the heredoc, the row is joined and later
    re-split on commas: a ragged TSV row (fewer cells than header columns)
    yields ``None`` values from ``csv.DictReader`` and raises ``TypeError``
    here, and an allele value containing a comma inflates the part count.
    """
    return ",".join(profile[locus] for locus in loci_order if locus in profile)


def match_profile(profile_string, profiles, loci_order, id_column=None):
    # type: (str, Profiles, Sequence[str], Optional[str]) -> Tuple[Optional[str], str]
    """First table row consistent with the query profile, and its status.

    Returns ``(st, STATUS_KNOWN)`` for the first row where every position that
    is a wildcard on neither side agrees, else ``(None, STATUS_NOVEL_PROFILE)``.
    Rows whose allele count differs from the query's are skipped, as are rows
    from a table with no identifier column -- such a table can only ever be
    novel. ``id_column`` matches :func:`load_profiles`'s.
    """
    st_col = None
    query_parts = profile_string.split(",")
    for profile in profiles:
        if st_col is None:
            st_col = st_column(profile, id_column)

        table_parts = row_profile_string(profile, loci_order).split(",")
        if len(query_parts) != len(table_parts):
            continue

        match = True
        for query_part, table_part in zip(query_parts, table_parts):
            if query_part == WILDCARD or table_part == WILDCARD:
                continue
            if query_part != table_part:
                match = False
                break

        if match and st_col:
            return profile[st_col], STATUS_KNOWN

    return None, STATUS_NOVEL_PROFILE


def find_nearest_st(profile_string, profiles, loci_order, id_column=None):
    # type: (str, Profiles, Sequence[str], Optional[str]) -> Tuple[Optional[str], float]
    """Nearest ST to a novel profile by Hamming distance; wildcards cost 0.

    Ties go to the first row in table order. Missing loci are filled with
    ``?`` rather than shortening the row, so unlike :func:`match_profile` a row
    lacking a locus column is still comparable. Returns ``(None, inf)`` when no
    row is usable (no identifier column, or every row a different length).
    ``id_column`` matches :func:`load_profiles`'s.
    """
    query_parts = profile_string.split(",")
    best_st = None  # type: Optional[str]
    best_distance = float("inf")
    st_col = None

    for profile in profiles:
        if st_col is None:
            st_col = st_column(profile, id_column)
        if st_col is None:
            continue

        table_parts = [profile.get(locus, WILDCARD) for locus in loci_order]
        if len(table_parts) != len(query_parts):
            continue

        distance = sum(
            0 if q == WILDCARD or t == WILDCARD or q == t else 1
            for q, t in zip(query_parts, table_parts)
        )
        if distance < best_distance:
            best_distance = distance
            best_st = profile[st_col]

    return best_st, best_distance


def mean_confidence(allele_calls):
    # type: (Dict[str, Any]) -> float
    """Mean per-locus confidence across the calls; 0.0 when nothing is scored.

    A boolean-ish ``confidence`` key wins and contributes 1.0/0.0, so for the
    pipeline's own calls (which always carry it) this is the fraction of loci
    that passed their threshold, not a mean similarity. Only calls with neither
    ``confidence`` nor ``similarity``/``identity`` are left out of the mean.
    """
    scores = []
    for call in allele_calls.values():
        if "confidence" in call:
            scores.append(1.0 if call["confidence"] else 0.0)
        elif "similarity" in call or "identity" in call:
            scores.append(float(call.get("similarity", call.get("identity", 0.0))))
    return sum(scores) / len(scores) if scores else 0.0


def scheme_name(scheme, profiles_table):
    # type: (str, str) -> str
    """The scheme label: the explicit one, else the table's parent directory."""
    if scheme:
        return scheme
    return os.path.basename(os.path.dirname(profiles_table))


def build_profile_record(
    allele_calls,
    profiles,
    loci_order,
    scheme,
    strategy="balanced",
    alignment_used=False,
    id_column=None,
):
    # type: (Dict[str, Any], Profiles, Sequence[str], str, str, bool, Optional[str]) -> Dict[str, Any]
    """Assemble the typing record the heredoc wrote to `profile_result.json`.

    ``nearest_st``/``nearest_st_distance`` appear only for a novel profile
    matched against a non-empty table with a usable identifier column.
    ``id_column`` matches :func:`load_profiles`'s.
    """
    profile_string = build_profile_string(allele_calls, loci_order)
    profile_id, status = match_profile(profile_string, profiles, loci_order, id_column)

    nearest_st = None
    nearest_st_distance = None
    if status == STATUS_NOVEL_PROFILE and profiles:
        nearest_st, nearest_st_distance = find_nearest_st(
            profile_string, profiles, loci_order, id_column
        )

    overall_confidence = mean_confidence(allele_calls)

    result = {
        "profile_id": profile_id if profile_id else UNKNOWN_PROFILE_ID,
        "profile_type": "sequence_type",
        "scheme": scheme,
        "status": status,
        "confidence": max(0.0, min(1.0, overall_confidence)),
        "allele_profile": profile_string,
        "allele_calls": allele_calls,
        "method": {
            "strategy": strategy,
            "alignment_used": alignment_used,
            "tools": ["sourmash", "minimap2"] if alignment_used else ["sourmash"],
        },
        "notes": {
            "num_loci": len(loci_order),
            "num_called": len(allele_calls),
            "mean_confidence": overall_confidence,
        },
    }  # type: Dict[str, Any]

    if nearest_st is not None:
        result["nearest_st"] = nearest_st
        result["nearest_st_distance"] = nearest_st_distance

    return result


def lookup_profile(
    allele_calls_path,
    profiles_table_path,
    strategy="balanced",
    alignment_used=False,
    scheme="",
    id_column=None,
):
    # type: (str, str, str, bool, str, Optional[str]) -> Dict[str, Any]
    """Read calls and a profiles table from disk and return the typing record.

    Convenience entry point for the caller that consumes workflow outputs; the
    caller decides where (and whether) to serialise the record. ``id_column``
    matches :func:`load_profiles`'s.
    """
    allele_calls = load_allele_calls(allele_calls_path)
    profiles, loci_order = load_profiles(profiles_table_path, id_column)
    return build_profile_record(
        allele_calls,
        profiles,
        loci_order,
        scheme_name(scheme, profiles_table_path),
        strategy,
        alignment_used,
        id_column,
    )
