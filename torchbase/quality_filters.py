#!/usr/bin/env python

"""Input-preparation filters lifted out of two built-in WDL task heredocs.

This module replaces the embedded ``python3`` heredocs of:

* ``torchbase/workflows/builtin/tasks/filter_alleles.wdl``, task
  ``filter_alleles`` -- quality.json-driven suspect filtering under the
  ``exclude_suspect_alleles`` / ``exclude_suspect_loci`` /
  ``exclude_suspect_profiles`` flags.
* ``torchbase/workflows/builtin/tasks/depth_filter.wdl``, task
  ``depth_filter`` -- k-mer-depth filtering of read inputs (a pass-through
  for contigs).

Why the package layer and not a container: neither task is compute. Both are
*interpretation of torch metadata* -- they read a quality report the torch
already carries and decide which alleles the caller is willing to trust, and
they decide which input sequences are worth aligning. Sourmash, minimap2 and
BLAST+ do the compute; these two only shape their inputs. Concretely they run
at the same preparation stage as ``Torch.get_unified_files()``, which already
materializes the unified allele FASTA and profiles table inside the package,
so keeping the filters here means the whole "assemble the reference inputs"
step is importable and unit-testable without Docker.

Pure library code: paths and parsed data in, Python objects out. Callers
serialize the reports (``exclusions.json``, ``depth_filter_stats.json``) and
write the FASTAs; ``write_fasta`` is provided so the on-disk byte layout the
downstream compute tasks expect stays in one place.

Behaviour is preserved verbatim from the heredocs, including two warts noted
inline: the two tasks disagree on whether a parsed FASTA header keeps its
``>``, and ``excluded_loci`` is serialized from a set.
"""

import json
import os
from typing import Dict, Iterable, List, Optional, Set, Tuple

FastaEntry = Tuple[str, str]

#: Similarity threshold assumed when a locus in quality.json declares none.
DEFAULT_SIMILARITY_THRESHOLD = 90.0

#: Input type that skips depth filtering entirely.
PASSTHROUGH_INPUT_TYPE = "contigs"

DEFAULT_MIN_DEPTH = 3
DEFAULT_KSIZE = 21


# --------------------------------------------------------------------------
# filter_alleles
# --------------------------------------------------------------------------


def parse_allele_fasta(fasta_path) -> List[FastaEntry]:
    """Parse an allele FASTA into ``(header, sequence)`` entries.

    Headers have their leading ``>`` stripped and lines are whitespace-
    stripped, matching the ``filter_alleles`` heredoc. Sequence lines are
    concatenated, so wrapped records collapse to a single line on output.
    """
    entries = []  # type: List[FastaEntry]
    with open(fasta_path) as f:
        current_header = None  # type: Optional[str]
        current_seq = []  # type: List[str]
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_header is not None:
                    entries.append((current_header, ''.join(current_seq)))
                current_header = line[1:]
                current_seq = []
            else:
                current_seq.append(line)
        if current_header is not None:
            entries.append((current_header, ''.join(current_seq)))
    return entries


def extract_locus_and_allele(header: str) -> Tuple[str, str]:
    """Split a ``locus_name_allele_id`` header into ``(locus, allele_id)``.

    A header with no underscore has no allele id; it is returned unchanged
    as the locus with the literal allele id ``"unknown"``.
    """
    parts = header.split('_')
    if len(parts) >= 2:
        allele_id = parts[-1]
        locus = '_'.join(parts[:-1])
        return locus, allele_id
    return header, "unknown"


def empty_quality() -> Dict[str, Set[str]]:
    """Suspect sets meaning "no quality report available"."""
    return {
        'suspect_alleles': set(),
        'suspect_loci': set(),
        'suspect_profiles': set(),
    }


def quality_from_data(data: dict) -> Dict[str, Set[str]]:
    """Derive suspect alleles, loci and profile loci from quality report data.

    Three independent sources of suspicion, all merged into one payload:

    * ``loci[locus]['suspect']`` marks the whole locus;
    * ``loci[locus]['alleles'][allele]['suspect']`` marks ``locus_allele``;
    * ``loci[locus]['similarities']`` marks *both* alleles of any
      ``"a-b"`` pair scoring strictly below the locus ``threshold``
      (default :data:`DEFAULT_SIMILARITY_THRESHOLD`).

    ``profiles[pid]['suspect']`` contributes that profile's ``loci`` list,
    kept separate because only the ``exclude_profiles`` flag honours it.
    """
    suspect_alleles = set()  # type: Set[str]
    suspect_loci = set()  # type: Set[str]

    if 'loci' in data:
        for locus, locus_data in data['loci'].items():
            if locus_data.get('suspect', False):
                suspect_loci.add(locus)

            if 'alleles' in locus_data:
                for allele, allele_data in locus_data['alleles'].items():
                    if allele_data.get('suspect', False):
                        suspect_alleles.add("{}_{}".format(locus, allele))

            if 'similarities' in locus_data:
                threshold = locus_data.get(
                    'threshold', DEFAULT_SIMILARITY_THRESHOLD
                )
                for pair, similarity in locus_data['similarities'].items():
                    if similarity < threshold:
                        allele1, allele2 = pair.split('-')
                        suspect_alleles.add("{}_{}".format(locus, allele1))
                        suspect_alleles.add("{}_{}".format(locus, allele2))

    suspect_profiles = set()  # type: Set[str]
    if 'profiles' in data:
        for _profile_id, profile_data in data['profiles'].items():
            if profile_data.get('suspect', False):
                if 'loci' in profile_data:
                    suspect_profiles.update(profile_data['loci'])

    return {
        'suspect_alleles': suspect_alleles,
        'suspect_loci': suspect_loci,
        'suspect_profiles': suspect_profiles,
    }


def load_quality(quality_path) -> Dict[str, Set[str]]:
    """Load a quality report, tolerating "there isn't one".

    ``None``, the empty string, a path that does not exist, and a file whose
    contents are empty or whitespace-only all yield :func:`empty_quality`
    -- i.e. nothing is suspect and every flag becomes a no-op.
    """
    if not quality_path or not os.path.exists(str(quality_path)):
        return empty_quality()

    with open(quality_path) as f:
        raw = f.read().strip()
    if not raw:
        return empty_quality()
    return quality_from_data(json.loads(raw))


def build_exclusion_sets(
    quality: Dict[str, Set[str]],
    exclude_alleles: bool = False,
    exclude_loci: bool = False,
    exclude_profiles: bool = False,
) -> Tuple[Set[str], Set[str]]:
    """Turn suspect data plus the three flags into ``(loci, alleles)`` to drop.

    The flags are a strict precedence chain, not independent switches -- the
    most aggressive one set wins and the others are ignored:

    * ``exclude_profiles``: suspect profile loci + suspect loci + suspect
      alleles;
    * else ``exclude_loci``: suspect loci + suspect alleles;
    * else ``exclude_alleles``: suspect alleles only;
    * else: nothing is excluded (the default include-everything path).
    """
    to_exclude_loci = set()  # type: Set[str]
    to_exclude_alleles = set()  # type: Set[str]

    if exclude_profiles:
        to_exclude_loci.update(quality['suspect_profiles'])
        to_exclude_loci.update(quality['suspect_loci'])
        to_exclude_alleles.update(quality['suspect_alleles'])
    elif exclude_loci:
        to_exclude_loci.update(quality['suspect_loci'])
        to_exclude_alleles.update(quality['suspect_alleles'])
    elif exclude_alleles:
        to_exclude_alleles.update(quality['suspect_alleles'])

    return to_exclude_loci, to_exclude_alleles


def filter_allele_entries(
    entries: Iterable[FastaEntry],
    to_exclude_loci: Set[str],
    to_exclude_alleles: Set[str],
) -> Tuple[List[FastaEntry], List[str], Set[str]]:
    """Split entries into ``(kept, excluded allele names, excluded loci)``.

    A locus-level exclusion short-circuits the allele-level check, so an
    allele dropped because its locus is suspect is *not* listed in the
    excluded-allele names. Excluded allele names are the normalized
    ``locus_alleleid`` form, which for a header without an underscore is
    ``header_unknown``.
    """
    excluded_alleles = []  # type: List[str]
    excluded_loci = set()  # type: Set[str]
    kept = []  # type: List[FastaEntry]

    for header, sequence in entries:
        locus, allele_id = extract_locus_and_allele(header)
        full_allele_name = "{}_{}".format(locus, allele_id)

        exclude = False
        if locus in to_exclude_loci:
            exclude = True
            excluded_loci.add(locus)
        elif full_allele_name in to_exclude_alleles:
            exclude = True
            excluded_alleles.append(full_allele_name)

        if not exclude:
            kept.append((header, sequence))

    return kept, excluded_alleles, excluded_loci


def exclusion_report(
    total_input: int,
    kept: List[FastaEntry],
    excluded_alleles: List[str],
    excluded_loci: Set[str],
) -> dict:
    """Build the payload downstream tasks read as ``exclusions.json``.

    ``excluded_loci`` is serialized from a set, so its order is not stable
    across interpreters; callers that care must sort it.
    """
    return {
        'excluded_alleles': list(excluded_alleles),
        'excluded_loci': list(excluded_loci),
        'num_excluded_alleles': len(excluded_alleles),
        'num_excluded_loci': len(excluded_loci),
        'total_input_alleles': total_input,
        'total_output_alleles': len(kept),
    }


def filter_alleles(
    allele_fasta,
    quality: Optional[Dict[str, Set[str]]] = None,
    exclude_alleles: bool = False,
    exclude_loci: bool = False,
    exclude_profiles: bool = False,
) -> Tuple[List[FastaEntry], dict]:
    """Filter an allele FASTA against a quality report.

    ``quality`` is the payload from :func:`load_quality`; ``None`` means no
    quality report, which excludes nothing regardless of the flags.

    Returns ``(kept entries, exclusion report)``. The entries carry no ``>``;
    hand them to :func:`write_fasta` to reproduce the FASTA the compute tasks
    consume.
    """
    if quality is None:
        quality = empty_quality()

    entries = parse_allele_fasta(allele_fasta)
    to_exclude_loci, to_exclude_alleles = build_exclusion_sets(
        quality,
        exclude_alleles=exclude_alleles,
        exclude_loci=exclude_loci,
        exclude_profiles=exclude_profiles,
    )
    kept, excluded_alleles, excluded_loci = filter_allele_entries(
        entries, to_exclude_loci, to_exclude_alleles
    )
    report = exclusion_report(
        len(entries), kept, excluded_alleles, excluded_loci
    )
    return kept, report


def write_fasta(entries: Iterable[FastaEntry], out_path) -> None:
    """Write ``(header, sequence)`` pairs as unwrapped FASTA with ``>``."""
    with open(out_path, 'w') as f:
        for header, sequence in entries:
            f.write('>{}\n'.format(header))
            f.write('{}\n'.format(sequence))


# --------------------------------------------------------------------------
# depth_filter
# --------------------------------------------------------------------------


def parse_sequence_fasta(path) -> List[FastaEntry]:
    """Parse a query FASTA into ``(header, sequence)`` records.

    Note the deliberate asymmetry with :func:`parse_allele_fasta`: the
    ``depth_filter`` heredoc keeps the leading ``>`` on the header and only
    right-strips lines, so query headers round-trip byte for byte including
    any leading whitespace. Preserved as-is.
    """
    records = []  # type: List[FastaEntry]
    header = None  # type: Optional[str]
    seq_lines = []  # type: List[str]
    with open(path) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(seq_lines)))
                header, seq_lines = line, []
            else:
                seq_lines.append(line)
    if header is not None:
        records.append((header, "".join(seq_lines)))
    return records


def kmer_depth(seq: str, k: int) -> float:
    """Mean k-mer depth: k-mer occurrences over distinct k-mers.

    A sequence shorter than ``k`` has depth ``0.0`` (and so is always
    removed when filtering is active).
    """
    if len(seq) < k:
        return 0.0
    counts = {}  # type: Dict[str, int]
    for i in range(len(seq) - k + 1):
        km = seq[i:i + k]
        counts[km] = counts.get(km, 0) + 1
    return sum(counts.values()) / len(counts)


def filter_by_depth(
    records: Iterable[FastaEntry],
    input_type: str = PASSTHROUGH_INPUT_TYPE,
    min_depth: int = DEFAULT_MIN_DEPTH,
    ksize: int = DEFAULT_KSIZE,
) -> Tuple[List[FastaEntry], int]:
    """Return ``(kept records, removed count)``.

    ``input_type == "contigs"`` passes everything through untouched and
    removes nothing. Any other value (reads, fastq, ...) keeps a record when
    its mean k-mer depth is ``>= min_depth`` -- the boundary is inclusive.
    """
    records = list(records)
    if input_type == PASSTHROUGH_INPUT_TYPE:
        return records, 0

    kept = []  # type: List[FastaEntry]
    removed = 0
    for header, seq in records:
        if kmer_depth(seq, ksize) >= min_depth:
            kept.append((header, seq))
        else:
            removed += 1
    return kept, removed


def depth_filter_report(
    total_input: int,
    kept: List[FastaEntry],
    removed: int,
    input_type: str,
    min_depth: int,
    ksize: int,
) -> dict:
    """Build the payload downstream tasks read as ``depth_filter_stats.json``."""
    return {
        "input_sequences": total_input,
        "kept_sequences": len(kept),
        "removed_sequences": removed,
        "input_type": input_type,
        "min_depth": min_depth,
        "ksize": ksize,
    }


def depth_filter(
    sequences,
    input_type: str = PASSTHROUGH_INPUT_TYPE,
    min_depth: int = DEFAULT_MIN_DEPTH,
    ksize: int = DEFAULT_KSIZE,
) -> Tuple[List[FastaEntry], dict]:
    """Drop low-coverage query sequences before allele calling.

    Returns ``(kept records, stats report)``. Records keep their ``>``
    headers; use :func:`write_sequence_fasta` to write them.
    """
    records = parse_sequence_fasta(sequences)
    kept, removed = filter_by_depth(
        records, input_type=input_type, min_depth=min_depth, ksize=ksize
    )
    report = depth_filter_report(
        len(records), kept, removed, input_type, min_depth, ksize
    )
    return kept, report


def write_sequence_fasta(records: Iterable[FastaEntry], out_path) -> None:
    """Write records whose headers already carry their ``>``."""
    with open(out_path, 'w') as out:
        for header, seq in records:
            out.write("{}\n{}\n".format(header, seq))
