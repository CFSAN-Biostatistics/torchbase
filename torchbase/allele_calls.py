"""Allele calling: turn compute artifacts into per-locus allele calls.

Package-layer replacement for the embedded Python heredocs of three built-in
WDL tasks:

* ``workflows/builtin/tasks/minhash.wdl`` task ``call_alleles`` -- see
  :func:`call_alleles`.
* ``workflows/builtin/tasks/minhash.wdl`` task ``call_alleles_minhash`` -- see
  :func:`call_alleles_minhash`.
* ``workflows/builtin/tasks/alignment.wdl`` task ``align_and_call`` (and its
  alias ``align_sequences``), *calling half only* -- see
  :func:`calls_from_alignment`.

This belongs at the package layer because it is interpretation, not compute:
sourmash and minimap2 do the expensive work and emit a similarity CSV or a SAM
file, and everything here is the judgement applied afterwards -- which allele of
a locus wins, how identity is derived from a mismatch count, and whether a call
clears its confidence threshold. Containers hold tools; deciding what a tool's
numbers mean is torchbase's job, and only here can it be unit-tested without
Docker.

Similarity matrix layout: ``sourmash compare --csv --ani`` over ``--singleton``
sketches emits an all-vs-all NxN matrix whose header row names every signature
(queries first, then alleles) and whose data rows carry no row labels. Query
``i`` is therefore row ``i + 1``, and allele ``j`` is column
``num_queries + j``.
"""

import csv
import os
from typing import Any, Callable, Dict, Iterable, Iterator, List, Optional, Sequence, Tuple, Union

#: ``call_alleles`` hard-codes this similarity and compares *strictly* above it.
CALL_ALLELES_CONFIDENCE_SIMILARITY = 0.9

#: Default of ``call_alleles_minhash``'s ``confidence_threshold`` WDL input.
DEFAULT_CONFIDENCE_THRESHOLD = 0.85

#: Default of ``align_and_call``'s ``identity_threshold`` WDL input.
DEFAULT_IDENTITY_THRESHOLD = 0.90

#: Presence/absence calling (docs/adr/0003-profile-matching-value-domain-agnostic.md):
#: a marker "hits" when both clear these defaults.
DEFAULT_PRESENCE_IDENTITY_THRESHOLD = 0.90
DEFAULT_PRESENCE_COVERAGE_THRESHOLD = 0.90

#: Tokens `calls_from_presence` reports as `allele_id`/`status`. Chosen so a
#: presence/absence locus's call is an ordinary string token to
#: `torchbase.profile_match` -- ABSENT/AMBIGUOUS are just tokens no
#: `profiles.tsv` row lists, so they fall through to `novel_profile` on their
#: own; a caller wanting to surface "ambiguous" as a distinct report status
#: reads `call["status"]` before handing calls to the matcher.
PRESENT = "present"
ABSENT = "absent"
AMBIGUOUS = "ambiguous"

Fasta = List[Tuple[str, str]]
FastaSource = Union[str, bytes, "os.PathLike", Fasta]
AlleleRecord = Dict[str, Any]
SamRecord = Dict[str, Any]
Call = Dict[str, Any]
Calls = Dict[str, Call]


def _is_path(value: Any) -> bool:
    return isinstance(value, (str, bytes, os.PathLike))


# --------------------------------------------------------------------------
# Shared name parsing and the confidence decision
# --------------------------------------------------------------------------


def parse_fasta(fasta_path: Union[str, bytes, "os.PathLike"]) -> Fasta:
    """Read a FASTA file into ``(header, sequence)`` pairs, in file order.

    Multi-line sequences are joined; a header with no sequence yields an empty
    string. Anything before the first ``>`` is ignored.
    """
    sequences = []  # type: Fasta
    with open(fasta_path) as f:
        current_header = None  # type: Optional[str]
        current_seq = []  # type: List[str]
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_header is not None:
                    sequences.append((current_header, "".join(current_seq)))
                current_header = line[1:]
                current_seq = []
            else:
                current_seq.append(line)
        if current_header is not None:
            sequences.append((current_header, "".join(current_seq)))
    return sequences


def parse_fasta_dict(fasta_path: Union[str, bytes, "os.PathLike"]) -> Dict[str, str]:
    """Read a FASTA file into ``{sequence name: sequence}``.

    The key is the header up to its first whitespace, which is the name aligners
    put in SAM's QNAME field. Keying on the full header — as the WDL task this
    replaces did — means no lookup ever matches for a query whose header carries
    a description, i.e. for essentially every real assembly, and every
    alignment-based call silently comes back unconfident.

    Duplicate names collapse, last one winning.
    """
    return dict(
        (header.split()[0] if header.split() else header, sequence)
        for header, sequence in parse_fasta(fasta_path)
    )


def extract_locus_and_allele(header: str) -> Tuple[str, str]:
    """Split a ``<locus>_<allele_id>`` FASTA header on its *last* underscore.

    ``"thrA_ec_42"`` is ``("thrA_ec", "42")``. A header with no underscore has
    no allele id and is reported as ``(header, "unknown")``.
    """
    parts = header.split("_")
    if len(parts) >= 2:
        allele_id = parts[-1]
        locus = "_".join(parts[:-1])
        return locus, allele_id
    return header, "unknown"


def confident_above(threshold: float) -> Callable[[float], bool]:
    """Confidence decision that *excludes* the threshold (``score > t``)."""

    def is_confident(score: float) -> bool:
        return score > threshold

    return is_confident


def confident_at_or_above(threshold: float) -> Callable[[float], bool]:
    """Confidence decision that *includes* the threshold (``score >= t``)."""

    def is_confident(score: float) -> bool:
        return score >= threshold

    return is_confident


def _call(allele_id: Optional[str], score: float, score_key: str,
          is_confident: Callable[[float], bool]) -> Call:
    """One call record: score clamped to [0, 1], confidence from the raw score."""
    return {
        "allele_id": allele_id,
        score_key: max(0.0, min(1.0, score)),
        "confidence": is_confident(score),
    }


# --------------------------------------------------------------------------
# MinHash path: sourmash comparison matrix -> calls
# --------------------------------------------------------------------------


def read_similarity_matrix(matrix_path: Union[str, bytes, "os.PathLike"]) -> List[List[str]]:
    """Read a sourmash comparison CSV into raw rows, header row included."""
    with open(matrix_path, newline="") as f:
        return list(csv.reader(f))


def _as_matrix_rows(matrix_or_path: Any) -> List[List[str]]:
    if _is_path(matrix_or_path):
        return read_similarity_matrix(matrix_or_path)
    return [list(row) for row in matrix_or_path]


def _as_fasta(source: FastaSource) -> Fasta:
    if _is_path(source):
        return parse_fasta(source)
    return [(header, seq) for header, seq in source]


def _as_fasta_dict(source: Union[FastaSource, Dict[str, str]]) -> Dict[str, str]:
    if _is_path(source):
        return parse_fasta_dict(source)
    if isinstance(source, dict):
        return source
    return dict(source)


def is_empty_matrix(rows: Sequence[Sequence[str]]) -> bool:
    """True when the CSV holds no data rows: nothing was compared."""
    return len(rows) <= 1


def group_alleles_by_locus(allele_seqs: Fasta) -> Dict[str, List[AlleleRecord]]:
    """Group allele records by locus, keeping each allele's matrix column index.

    The index is the allele's position in the FASTA file, which is what fixes
    its column in the sourmash comparison matrix.
    """
    alleles_by_locus = {}  # type: Dict[str, List[AlleleRecord]]
    for idx, (header, _seq) in enumerate(allele_seqs):
        locus, allele_id = extract_locus_and_allele(header)
        alleles_by_locus.setdefault(locus, []).append(
            {"allele_id": allele_id, "header": header, "index": idx}
        )
    return alleles_by_locus


def max_similarity_per_allele(
    rows: Sequence[Sequence[str]], num_queries: int, num_alleles: int
) -> List[float]:
    """Best similarity any query reached, per allele column.

    Blank cells and columns absent from a short row count as 0.0; a
    non-numeric cell raises ``ValueError``. Raises ``ValueError`` when the
    matrix is not the ``num_queries + num_alleles`` square the sketches imply.
    """
    expected_size = num_queries + num_alleles
    if len(rows) != expected_size + 1:  # +1 for header
        raise ValueError(
            "Matrix size mismatch: expected {} rows, got {}".format(
                expected_size + 1, len(rows)
            )
        )

    max_similarities = [0.0] * num_alleles
    for query_idx in range(num_queries):
        data_row_idx = query_idx + 1  # Skip header row
        if data_row_idx < len(rows):
            row = rows[data_row_idx]
            for allele_idx in range(num_alleles):
                col_idx = num_queries + allele_idx  # Allele columns follow queries
                if col_idx < len(row):
                    sim = float(row[col_idx]) if row[col_idx] else 0.0
                    max_similarities[allele_idx] = max(
                        max_similarities[allele_idx], sim
                    )
    return max_similarities


def best_similarity_per_locus(
    alleles_by_locus: Dict[str, List[AlleleRecord]],
    max_similarities: Sequence[float],
) -> List[Tuple[str, str, float]]:
    """Highest-similarity allele of each locus, loci in sorted order.

    Ties break on FASTA order: the first allele of the locus reaching the
    winning similarity wins. A locus whose every allele lies past the end of
    the matrix is omitted entirely.
    """
    best = []  # type: List[Tuple[str, str, float]]
    for locus, alleles in sorted(alleles_by_locus.items()):
        best_match = None  # type: Optional[str]
        best_similarity = -1.0
        for allele in alleles:
            idx = allele["index"]
            if idx < len(max_similarities):
                sim = max_similarities[idx]
                if sim > best_similarity:
                    best_similarity = sim
                    best_match = allele["allele_id"]
        if best_match is not None:
            best.append((locus, best_match, best_similarity))
    return best


def calls_from_similarity(
    matrix_or_path: Any,
    query_sequences: FastaSource,
    allele_fasta: FastaSource,
    is_confident: Callable[[float], bool],
) -> Calls:
    """Call one allele per locus from a sourmash comparison matrix.

    ``matrix_or_path`` is a CSV path or already-parsed rows; the two FASTA
    arguments are paths or already-parsed ``(header, sequence)`` pairs. Only
    the *count* of query sequences matters (it fixes the allele columns);
    allele headers supply the locus and allele ids.

    An empty matrix (header row only, or no rows at all) yields no calls.
    ``is_confident`` decides confidence from the raw similarity -- the two
    minhash tasks disagree about it, so it is injected rather than assumed.
    """
    query_seqs = _as_fasta(query_sequences)
    allele_seqs = _as_fasta(allele_fasta)
    rows = _as_matrix_rows(matrix_or_path)

    if is_empty_matrix(rows):
        return {}

    max_similarities = max_similarity_per_allele(
        rows, len(query_seqs), len(allele_seqs)
    )
    return {
        locus: _call(allele_id, similarity, "similarity", is_confident)
        for locus, allele_id, similarity in best_similarity_per_locus(
            group_alleles_by_locus(allele_seqs), max_similarities
        )
    }


def call_alleles(
    matrix_or_path: Any,
    query_sequences: FastaSource,
    allele_fasta: FastaSource,
) -> Calls:
    """WDL task ``call_alleles``: fixed 0.9 threshold, *exclusive*.

    Kept separate from :func:`call_alleles_minhash` -- a different workflow
    calls it, and it differs observably: the threshold is hard-coded and the
    comparison excludes it, so a locus at exactly 0.9 is *not* confident. The
    task also emits a profile string; see :func:`format_profile`.
    """
    return calls_from_similarity(
        matrix_or_path,
        query_sequences,
        allele_fasta,
        confident_above(CALL_ALLELES_CONFIDENCE_SIMILARITY),
    )


def call_alleles_minhash(
    matrix_or_path: Any,
    query_sequences: FastaSource,
    allele_fasta: FastaSource,
    confidence_threshold: float = DEFAULT_CONFIDENCE_THRESHOLD,
) -> Calls:
    """WDL task ``call_alleles_minhash``: caller's threshold, *inclusive*.

    Differs from :func:`call_alleles` in exactly two ways: the threshold is an
    input defaulting to 0.85 rather than a hard-coded 0.9, the comparison
    includes it (a locus at exactly the threshold *is* confident), and no
    profile string is produced.
    """
    return calls_from_similarity(
        matrix_or_path,
        query_sequences,
        allele_fasta,
        confident_at_or_above(confidence_threshold),
    )


def format_profile(calls: Calls) -> str:
    """Render calls as the comma-separated ``<locus>_<allele_id>`` profile.

    Only ``call_alleles`` produces this; the string is empty when nothing was
    called. Confidence is *not* considered -- every call contributes.
    """
    return ",".join(
        "{}_{}".format(locus, call["allele_id"]) for locus, call in calls.items()
    )


# --------------------------------------------------------------------------
# Alignment path: minimap2 SAM records -> calls
# --------------------------------------------------------------------------


def parse_sam_lines(lines: Iterable[str]) -> Iterator[SamRecord]:
    """Yield mapped alignment records from SAM text lines.

    Skipped: ``@`` header lines, records with fewer than 11 fields, records
    whose reference is ``*``, and records with ``MAPQ == 0`` (minimap2's
    multi-mapping / unmapped marker). A record with no ``NM:i:`` tag is treated
    as having zero mismatches. Each record is a dict with ``query_name``,
    ``ref_name``, ``mapq``, ``nm`` and ``tlen``.
    """
    for line in lines:
        if line.startswith("@"):
            continue
        fields = line.strip().split("\t")
        if len(fields) < 11:
            continue

        query_name = fields[0]
        ref_name = fields[2]
        mapq = int(fields[4])

        if mapq == 0 or ref_name == "*":  # Skip unmapped
            continue

        nm = 0
        for tag in fields[11:]:
            if tag.startswith("NM:i:"):
                nm = int(tag.split(":")[-1])
                break

        yield {
            "query_name": query_name,
            "ref_name": ref_name,
            "mapq": mapq,
            "nm": nm,
            "tlen": int(fields[8]),
            "cigar": fields[5],
        }


def read_sam_records(sam_path: Union[str, bytes, "os.PathLike"]) -> List[SamRecord]:
    """Read alignment records from a SAM file, tolerating a broken one.

    A missing file yields no records, and a line that cannot be parsed (a
    non-integer MAPQ, say) stops parsing and returns the records collected so
    far. That tolerance is the WDL task's behaviour and is deliberately
    preserved: minimap2 failing must not fail the whole job.
    """
    records = []  # type: List[SamRecord]
    try:
        with open(sam_path) as f:
            for record in parse_sam_lines(f):
                records.append(record)
    except Exception:
        pass
    return records


def _as_sam_records(records_or_path: Any) -> List[SamRecord]:
    if _is_path(records_or_path):
        return read_sam_records(records_or_path)
    return list(records_or_path)


def aligned_length(cigar: str) -> int:
    """Bases the alignment actually covers, from a CIGAR string.

    Counts the operations that consume alignment columns (``M``, ``=``, ``X``,
    ``I``, ``D``) and ignores clipping (``S``, ``H``) and padding, so a gene
    found inside a long contig measures the gene, not the contig.
    """
    length = 0
    digits = ""
    for character in cigar:
        if character.isdigit():
            digits += character
            continue
        if character in "M=XID" and digits:
            length += int(digits)
        digits = ""
    return length


def record_identity(record: SamRecord, query_sequences: Dict[str, str]) -> float:
    """Fractional identity of an alignment: ``(aligned - mismatches) / aligned``.

    Length comes from the CIGAR's aligned block. The WDL task this replaces used
    the whole query sequence's length, which is only right when a query holds a
    single locus: for an assembly it divides a gene's mismatches by megabases
    and reports near-perfect identity for anything that aligns at all.

    Falls back to TLEN and then to the query's length when there is no CIGAR, so
    records from an aligner that omits it still score. An unknown query or a
    zero length scores 0.0, and more mismatches than bases floors at 0.0.
    """
    length = aligned_length(record.get("cigar", "") or "")
    if length <= 0:
        tlen = record["tlen"]
        length = tlen if tlen > 0 else len(
            query_sequences.get(record["query_name"], "")
        )
    if length <= 0:
        return 0.0
    return max(0.0, (length - record["nm"]) / length)


def best_allele_per_locus(
    records: Iterable[SamRecord], query_sequences: Dict[str, str]
) -> Dict[str, Dict[str, Any]]:
    """Reduce alignment records to ``{locus: {allele_id, identity}}``.

    Best-wins per locus over every record, ties going to the earlier record.

    The WDL task this replaces reduced per *query* first — the best allele for
    each query sequence, then the best query per locus — which silently assumes
    one locus per query. Given an assembly, where a single contig aligns to
    every locus in the scheme, that assumption throws away all but one locus, so
    `torchbase run -c assembly.fasta` could return at most one allele call. The
    locus a reference belongs to is a property of the reference, not of the
    query it matched.
    """
    locus_results = {}  # type: Dict[str, Dict[str, Any]]
    for record in records:
        identity = record_identity(record, query_sequences)
        locus, allele_id = extract_locus_and_allele(record["ref_name"])
        incumbent = locus_results.get(locus)
        if incumbent is None or identity > incumbent["identity"]:
            locus_results[locus] = {
                "allele_id": allele_id,
                "identity": identity,
            }
    return locus_results


def calls_from_alignment(
    records_or_path: Any,
    query_sequences: Dict[str, str],
    identity_threshold: float = DEFAULT_IDENTITY_THRESHOLD,
) -> Dict[str, Dict[str, Any]]:
    """Per-locus allele calls from alignment records.

    The calling half of WDL task ``align_and_call`` (and of its alias
    ``align_sequences``, whose calling half is identical -- the two differ only
    in how the minimap2 preset is chosen, which stays in the container).
    ``records_or_path`` is a SAM path or already-parsed records;
    ``query_sequences`` is ``{sequence name: sequence}``, consulted only for the
    identity fallback when a record carries neither CIGAR nor TLEN.

    ``confidence`` is ``identity >= identity_threshold``, *inclusive*, decided
    on the raw identity; the reported ``identity`` is clamped to [0, 1].
    """
    is_confident = confident_at_or_above(identity_threshold)
    return {
        locus: _call(result["allele_id"], result["identity"], "identity", is_confident)
        for locus, result in best_allele_per_locus(
            _as_sam_records(records_or_path), query_sequences
        ).items()
    }



# --------------------------------------------------------------------------
# Presence/absence path: alignment records -> present/absent/ambiguous calls
# --------------------------------------------------------------------------
#
# See docs/adr/0003-profile-matching-value-domain-agnostic.md. A presence/
# absence marker panel (LisSero's five serogroup-determinant genes;
# ShigaTyper's O-antigen and accessory markers) is called the same way an
# allelic locus is -- best candidate per locus above a threshold -- with two
# differences an identity-only reduction does not need: no qualifying
# candidate is an explicit ABSENT, not a low-confidence best guess, and more
# than one qualifying candidate at once is AMBIGUOUS (contamination, or two
# members of a mutually-exclusive marker family both present), not silently
# resolved by picking the higher-scoring one. Once called, the token feeds
# `torchbase.profile_match` exactly like any other locus's `allele_id`.


def record_coverage(record: SamRecord, reference_sequences: Dict[str, str]) -> float:
    """Fractional reference coverage of an alignment: aligned length / reference length.

    `reference_sequences` is `{reference name: sequence}` -- the marker/allele
    FASTA aligned against, keyed the same way as `record["ref_name"]`. Coverage
    is clamped to 1.0 (an alignment can run past its reference's ends) and is
    0.0 for a reference this dict does not know.
    """
    length = aligned_length(record.get("cigar", "") or "")
    ref_len = len(reference_sequences.get(record["ref_name"], ""))
    if ref_len <= 0:
        return 0.0
    return min(1.0, length / ref_len)


def calls_from_presence(
    records_or_path: Any,
    query_sequences: Union[FastaSource, Dict[str, str]],
    reference_sequences: Union[FastaSource, Dict[str, str]],
    identity_threshold: float = DEFAULT_PRESENCE_IDENTITY_THRESHOLD,
    coverage_threshold: float = DEFAULT_PRESENCE_COVERAGE_THRESHOLD,
    ambiguity_gap: Optional[float] = None,
) -> Calls:
    """Presence/absence allele calls: one PRESENT/ABSENT/AMBIGUOUS token per locus.

    A candidate qualifies when its identity and coverage both clear their
    thresholds (`record_identity`, `record_coverage`). Per locus:

    +- zero qualifying candidates -> `allele_id`/`status` = ABSENT
    +- exactly one qualifying candidate -> its allele_id, status PRESENT
    +- more than one, with `ambiguity_gap` unset (the default) -- LisSero's and
      ShigaTyper's own rule, "more than one hit is contamination," with no
      tolerance -- -> AMBIGUOUS, and the call carries `candidates` (every
      qualifying allele_id)

    `ambiguity_gap`, when set, tolerates ties: two qualifying candidates are
    only AMBIGUOUS if the best does not beat the second by more than that
    much identity. Every locus present in `reference_sequences` is reported,
    even one with zero alignment records at all, so a marker minimap2 never
    even attempted to place is reported ABSENT rather than silently omitted.
    """
    reference_sequences = _as_fasta_dict(reference_sequences)
    query_sequences = _as_fasta_dict(query_sequences)
    records = _as_sam_records(records_or_path)

    qualifying_by_locus = {}  # type: Dict[str, List[Tuple[str, float]]]
    for record in records:
        identity = record_identity(record, query_sequences)
        coverage = record_coverage(record, reference_sequences)
        if identity < identity_threshold or coverage < coverage_threshold:
            continue
        locus, allele_id = extract_locus_and_allele(record["ref_name"])
        qualifying_by_locus.setdefault(locus, []).append((allele_id, identity))

    all_loci = {
        extract_locus_and_allele(ref_name)[0] for ref_name in reference_sequences
    }

    calls = {}  # type: Calls
    for locus in all_loci:
        qualifying = sorted(
            qualifying_by_locus.get(locus, []), key=lambda c: -c[1]
        )
        if not qualifying:
            calls[locus] = {
                "allele_id": ABSENT, "status": ABSENT, "confidence": False,
            }
            continue

        tied = (
            len(qualifying) > 1 if ambiguity_gap is None
            else len(qualifying) > 1 and qualifying[0][1] - qualifying[1][1] <= ambiguity_gap
        )
        if tied:
            if ambiguity_gap is None:
                candidates = [allele_id for allele_id, _ in qualifying]
            else:
                candidates = [
                    allele_id for allele_id, identity in qualifying
                    if qualifying[0][1] - identity <= ambiguity_gap
                ]
            calls[locus] = {
                "allele_id": AMBIGUOUS, "status": AMBIGUOUS, "confidence": False,
                "candidates": candidates,
            }
        else:
            allele_id, identity = qualifying[0]
            calls[locus] = {
                "allele_id": allele_id, "status": PRESENT,
                "identity": max(0.0, min(1.0, identity)), "confidence": True,
            }
    return calls