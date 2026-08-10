"""Core operon-typing algorithm: HSP normalization, frameshift stitching,
locus reduction, synteny pairing, and status-ladder scoring.

Reference implementation for the `operon` typing model — see
docs/operon-strategy-plan.md §2-§5 for the design rationale (protein-space,
synteny-aware calling of multi-subunit operons, generalized from NCBI
StxTyper).

This module is the *only* implementation: `operon_typing.wdl`'s tasks receive
it as a `File` input and import it inside their containers, so the algorithm
that runs in a workflow is the algorithm covered by `tests/test_operon.py`.
It therefore has to stay dependency-free (standard library only) and 3.8+
compatible.

Everything here is driven by the torch's `[operon]` config block; there is no
stx-specific logic. Where StxTyper hard-codes stx behaviour, the
corresponding knob is named in the docstring:

    subunit_order            transcription order, e.g. ["A", "B"]
    intergenic_max           max intergenic distance, bp (stx: 36)
    intergenic_min           min intergenic distance, bp; negative allows
                             overlapping genes (ETEC LT: eltA/eltB overlap 4)
    intergenic_relax_factor  multiplier for the recovery passes (stx: 2)
    identity_thresholds      per-class combined-identity floor + `default`
    generalized_classes      classes identity alone cannot separate
    residue_rules            residue decision table per generalized class
    superclass_pattern       regex whose group 1 is a class's parent class,
                             used to report a coarser type for operons that
                             are not COMPLETE (stx: "2c" -> "2")
    min_operon_identity      floor below which a candidate is not reported
                             (StxTyper's identity_min, 0.8)
    overlap_slack            bp slack for "this operon is inside that one"
                             containment tests (StxTyper's slack, 30)
"""

import re
from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Tuple

# Disruption/completeness status ladder, most to least complete. This is also
# the priority order used to resolve overlapping candidate operons, and it
# matches StxTyper's own ordering (stxtyper.cpp changelog 1.0.37/1.0.40:
# ... < PARTIAL_CONTIG_END < EXTENDED < PARTIAL < AMBIGUOUS < COMPLETE_NOVEL
# < COMPLETE, read most-preferred-first here).
STATUS_PRIORITY = (
    "COMPLETE",
    "COMPLETE_NOVEL",
    "AMBIGUOUS",
    "PARTIAL",
    "EXTENDED",
    "PARTIAL_CONTIG_END",
    "INTERNAL_STOP",
    "FRAMESHIFT",
)

STATUS_RANK = {status: i for i, status in enumerate(STATUS_PRIORITY)}

# Nucleotides per aligned residue in a translated search.
CODON = 3

# Defaults mirroring StxTyper's compile-time parameters.
DEFAULT_MIN_OPERON_IDENTITY = 0.8   # stxtyper.cpp identity_min
DEFAULT_OVERLAP_SLACK = 30          # stxtyper.cpp slack
DEFAULT_INTERGENIC_MAX = 36         # stxtyper.cpp intergenic_max

# Residue window BLAST may report twice (or skip) where a frameshift splits
# one alignment into two HSPs (StxTyper's Hsp::Merge window).
MERGE_WINDOW = 20

# A lone subunit is called PARTIAL_CONTIG_END when its missing partner would
# not have fit on the contig: intergenic_max + a minimal domain length
# (StxTyper's BlastAlignment::otherTruncated, 20 residues).
MIN_PARTNER_LENGTH = 20 * CODON


class OperonConfigError(ValueError):
    """Raised when an `[operon]` metadata.toml block fails validation."""


def status_to_coarse(operon_status: str) -> str:
    """Map the eight-value operon status ladder to the coarse,
    cross-typing-model `status` vocabulary (§5)."""
    if operon_status == "COMPLETE":
        return "known"
    if operon_status in ("COMPLETE_NOVEL", "AMBIGUOUS"):
        return "novel_profile"
    return "incomplete"


def select_threshold(class_label: str, thresholds: Dict[str, float]) -> float:
    """Per-class identity threshold, falling back to `default` (§3.1)."""
    if class_label in thresholds:
        return thresholds[class_label]
    if "default" in thresholds:
        return thresholds["default"]
    raise OperonConfigError(
        f"No identity threshold for class {class_label!r} and no "
        f"[operon.identity_thresholds] default set."
    )


def combined_identity(subunit_idents: Dict[str, Tuple[int, int]]) -> float:
    """Combined identity = sum(nident) / sum(aln_len) across subunits (§3.1).

    `subunit_idents` maps subunit tag -> (nident, aln_len). Denominator is the
    alignment length, matching StxTyper's `Operon::relIdentity()`; reference
    coverage is reported separately.
    """
    nident = sum(n for n, _ in subunit_idents.values())
    length = sum(n for _, n in subunit_idents.values())
    if length == 0:
        return 0.0
    return nident / length


def superclass_of(class_label: Optional[str], operon_cfg: dict) -> Optional[str]:
    """Parent class of `class_label` under `[operon] superclass_pattern`.

    The pattern's group 1 is the parent (stx: `^([0-9]+)` maps "2c" -> "2",
    "1a" -> "1"). Without a pattern a class is its own parent, which disables
    parent-class reporting rather than guessing at scheme semantics.
    """
    if class_label is None:
        return None
    pattern = operon_cfg.get("superclass_pattern")
    if not pattern:
        return class_label
    match = re.match(pattern, class_label)
    if not match or not match.group(1):
        return class_label
    return match.group(1)


@dataclass
class ResidueRule:
    class_label: str
    positions: List[dict]  # [{"subunit": "A", "index": 311}, ...]
    table: List[dict]      # [{"call": "2a", "residues": [["F","S"], ...]}, ...]
    fallback: str

    def resolve(self, residues: Sequence[Optional[str]]) -> str:
        """Resolve observed `residues` (same order as `positions`) against
        the decision table; unresolved rows fall through to `fallback`."""
        for row in self.table:
            row_residues = row["residues"]
            if len(row_residues) != len(residues):
                continue
            if all(
                observed is not None and observed in alternatives
                for observed, alternatives in zip(residues, row_residues)
            ):
                return row["call"]
        return self.fallback


def build_residue_rules(raw_rules: List[dict]) -> Dict[str, ResidueRule]:
    """Build `class -> ResidueRule` from a parsed `[[operon.residue_rules]]` list."""
    rules = {}
    for raw in raw_rules:
        class_label = raw["class"]
        rules[class_label] = ResidueRule(
            class_label=class_label,
            positions=raw["positions"],
            table=raw["table"],
            fallback=raw.get("fallback", class_label),
        )
    return rules


def resolve_generalized_class(
    class_label: str,
    residues: Sequence[Optional[str]],
    generalized_classes: Dict[str, List[str]],
    residue_rules: Dict[str, ResidueRule],
) -> str:
    """Resolve a generalized class (e.g. "2") to a specific subtype (e.g.
    "2a") using its residue rule. Classes that aren't generalized pass
    through unchanged."""
    if class_label not in generalized_classes:
        return class_label
    rule = residue_rules.get(class_label)
    if rule is None:
        raise OperonConfigError(
            f"Generalized class {class_label!r} has no matching residue rule."
        )
    return rule.resolve(residues)


# Identity comparison, not field comparison: two alignments of different
# references can be field-identical yet are distinct hits, and the pipeline
# tracks them in sets.
@dataclass(eq=False)
class HSP:
    """A single translated-search hit, normalized at the protein_search.wdl
    task boundary (§4) — never leak BLAST tabular format past this shape.

    Coordinates: `start`/`stop` are 1-based inclusive nucleotide coordinates
    on the contig with `start <= stop` regardless of strand. `qstart`/`qend`
    are 0-based inclusive offsets into the reference protein, with the
    reference's trailing stop codon excluded from `ref_len` (StxTyper trims
    it in `Hsp::finishHsp`).
    """

    subunit: str
    contig: str
    strand: str
    start: int
    stop: int
    ref_subtype: str
    ref_class: str
    nident: int
    aln_len: int
    ref_len: int
    reference_accession: str = ""
    contig_len: int = 0
    qstart: int = 0
    qend: int = 0
    qseq: str = ""
    sseq: str = ""
    frameshift: bool = False
    internal_stop: bool = False
    # True/False when the alignment reaches the reference's stop codon and the
    # subject does/doesn't have one; None when the alignment stops earlier
    # (StxTyper's tri-state `c_complete`).
    stop_codon: Optional[bool] = None
    ambiguous: int = 0

    @property
    def identity(self) -> float:
        return self.nident / self.aln_len if self.aln_len else 0.0

    @property
    def ref_covered(self) -> int:
        return self.qend - self.qstart + 1

    @property
    def coverage(self) -> float:
        return self.ref_covered / self.ref_len if self.ref_len else 0.0

    @property
    def reference_complete(self) -> bool:
        """Whole reference protein aligned, stop codon included where the
        reference has one (StxTyper's `qComplete()`)."""
        return self.ref_covered == self.ref_len and self.stop_codon is not False

    @property
    def extended(self) -> bool:
        """Reference aligned from its first residue but the subject runs past
        the reference's stop codon (StxTyper's `c_extended()`)."""
        return self.qstart == 0 and self.stop_codon is False

    @property
    def contig_end(self) -> bool:
        """The alignment is clipped by the end of the contig on the side where
        the reference is incomplete (StxTyper's `sTruncated()`)."""
        head_room = self.start - 1
        tail_room = self.contig_len - self.stop if self.contig_len else CODON
        ref_head_missing = self.qstart > 0
        ref_tail_missing = self.qend + 1 < self.ref_len
        if self.strand == "-":
            ref_head_missing, ref_tail_missing = ref_tail_missing, ref_head_missing
        return (
            (head_room < CODON and ref_head_missing)
            or (tail_room < CODON and ref_tail_missing)
        )

    def frame(self) -> int:
        """Reading frame implied by this alignment, relative to the reference.

        Two HSPs of the same reference on the same contig/strand with
        different frames are the signature of a frameshift.
        """
        if self.strand == "+":
            return (self.start - CODON * self.qstart) % CODON
        return (self.stop + CODON * self.qstart) % CODON

    def residue_at(self, index: int) -> Optional[str]:
        """Subject residue aligned to 0-based reference offset `index`, or None
        when that reference position falls in a gap or outside the alignment
        (§3.1, §7 risk 1-2: the projection is where off-by-ones hide)."""
        qpos = self.qstart
        for qchar, schar in zip(self.qseq, self.sseq):
            if qchar == "-":
                continue
            if qpos == index:
                return None if schar == "-" else schar
            qpos += 1
        return None

    def to_dict(self) -> dict:
        return {
            "subunit": self.subunit,
            "contig": self.contig,
            "strand": self.strand,
            "start": self.start,
            "stop": self.stop,
            "ref_subtype": self.ref_subtype,
            "ref_class": self.ref_class,
            "nident": self.nident,
            "aln_len": self.aln_len,
            "ref_len": self.ref_len,
            "reference_accession": self.reference_accession,
            "contig_len": self.contig_len,
            "qstart": self.qstart,
            "qend": self.qend,
            "qseq": self.qseq,
            "sseq": self.sseq,
            "frameshift": self.frameshift,
            "internal_stop": self.internal_stop,
            "stop_codon": self.stop_codon,
            "ambiguous": self.ambiguous,
        }

    @classmethod
    def from_dict(cls, raw: dict) -> "HSP":
        fields = {
            "subunit", "contig", "strand", "start", "stop", "ref_subtype",
            "ref_class", "nident", "aln_len", "ref_len", "reference_accession",
            "contig_len", "qstart", "qend", "qseq", "sseq", "frameshift",
            "internal_stop", "stop_codon", "ambiguous",
        }
        return cls(**{k: v for k, v in raw.items() if k in fields})


# BLAST tabular columns `protein_search.wdl` requests, in order. This is the
# only place BLAST's format is understood; everything downstream sees HSPs.
BLAST_COLUMNS = (
    "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
    "qstart", "qend", "sstart", "send", "evalue", "bitscore", "nident",
    "qlen", "slen", "sstrand", "qseq", "sseq",
)


def parse_fasta(text: str) -> Dict[str, str]:
    """Minimal FASTA reader: header up to the first whitespace -> sequence."""
    sequences: Dict[str, str] = {}
    name = None
    chunks: List[str] = []
    for line in text.splitlines():
        if line.startswith(">"):
            if name is not None:
                sequences[name] = "".join(chunks)
            name = line[1:].split()[0] if line[1:].strip() else ""
            chunks = []
        elif name is not None:
            chunks.append(line.strip())
    if name is not None:
        sequences[name] = "".join(chunks)
    return sequences


def parse_ref_header(qseqid: str, header_format: str) -> Dict[str, str]:
    """Split a reference header into its declared fields (§3.1).

    Default format `accession|subunit_role|reference_subtype|class`:
    `subunit_role` is the structural position ("A"/"B", matched against
    `subunit_order`), `class` is the already-collapsed operon-level class
    token used for threshold lookup and cross-subunit class agreement.
    """
    return dict(zip(header_format.split("|"), qseqid.split("|")))


def normalize_hsps(
    rows: Sequence[Sequence[str]],
    references: Dict[str, str],
    header_format: str,
) -> List[HSP]:
    """Turn BLAST tabular rows into HSPs (§4's task boundary).

    Reference proteins conventionally carry their terminal stop codon (stx's
    `stx.prot` does). That column is not part of the protein, so it is trimmed
    here and recorded as `stop_codon` — otherwise every full-length hit reads
    as an internal stop and no operon is ever COMPLETE.
    """
    hsps: List[HSP] = []
    for row in rows:
        rec = dict(zip(BLAST_COLUMNS, row))
        qseqid = rec["qseqid"]
        ref = parse_ref_header(qseqid, header_format)
        ref_seq = references.get(qseqid, "")
        has_terminal_stop = ref_seq.endswith("*")
        ref_len = int(rec["qlen"]) - (1 if has_terminal_stop else 0)

        qstart = int(rec["qstart"]) - 1
        qend = int(rec["qend"]) - 1
        qseq, sseq = rec["qseq"], rec["sseq"]
        nident, aln_len = int(rec["nident"]), int(rec["length"])
        sstart, send = int(rec["sstart"]), int(rec["send"])
        strand = "+" if rec["sstrand"] == "plus" else "-"
        start, stop = min(sstart, send), max(sstart, send)

        stop_codon = None
        if has_terminal_stop and qend == ref_len:
            # The alignment reaches the reference's stop codon: whether the
            # subject has one there decides COMPLETE vs EXTENDED.
            qseq, sseq, trimmed = _trim_terminal_stop(qseq, sseq)
            stop_codon = trimmed == "*"
            qend -= 1
            aln_len -= 1
            if trimmed == "*":
                nident -= 1
            if trimmed != "-":
                # The stop codon's three nucleotides are not part of the
                # reported operon span either.
                if strand == "+":
                    stop -= CODON
                else:
                    start += CODON
        elif has_terminal_stop and qend == ref_len - 1:
            # Whole protein aligned but the stop codon did not: if the contig
            # continues, the subject's coding sequence runs past the
            # reference's terminator (StxTyper's `c_complete = efalse`).
            tail = (int(rec["slen"]) - stop) if strand == "+" else (start - 1)
            if tail >= CODON:
                stop_codon = False

        hsps.append(HSP(
            subunit=ref.get("subunit_role", qseqid),
            contig=rec["sseqid"],
            strand=strand,
            start=start,
            stop=stop,
            ref_subtype=ref.get("reference_subtype", qseqid),
            ref_class=ref.get("class", ref.get("reference_subtype", qseqid)),
            nident=nident,
            aln_len=aln_len,
            ref_len=ref_len,
            reference_accession=ref.get("accession", qseqid),
            contig_len=int(rec["slen"]),
            qstart=qstart,
            qend=qend,
            qseq=qseq,
            sseq=sseq,
            internal_stop="*" in sseq,
            stop_codon=stop_codon,
            ambiguous=sseq.count("X"),
        ))
    return hsps


def _trim_terminal_stop(qseq: str, sseq: str) -> Tuple[str, str, str]:
    """Drop the alignment column holding the reference's terminal stop codon.

    Returns the trimmed alignment plus the subject character that was aligned
    to the stop ("*" when the subject has its own stop codon, "-" when the
    alignment gapped over it).
    """
    for i in range(len(qseq) - 1, -1, -1):
        if qseq[i] != "-":
            return qseq[:i] + qseq[i + 1:], sseq[:i] + sseq[i + 1:], sseq[i]
    return qseq, sseq, "-"


def _upstream_gap(a: HSP, b: HSP) -> Optional[int]:
    """Intergenic distance from `a`'s downstream end to `b`'s upstream end.

    Strand-aware, and signed: a negative value means the two genes overlap,
    which is real — ETEC's `eltA`/`eltB` overlap by 4 bp — so callers decide
    how much overlap a scheme tolerates via `intergenic_min`. None when `b`
    does not follow `a` at all (wrong strand, or upstream of it).
    """
    if a.strand != b.strand:
        return None
    if a.strand == "+":
        if b.start <= a.start or b.stop <= a.stop:
            return None
        return b.start - a.stop - 1
    if b.stop >= a.stop or b.start >= a.start:
        return None
    return a.start - b.stop - 1


def stitch_frameshifts(hsps: Sequence[HSP], max_gap: int = DEFAULT_INTERGENIC_MAX) -> List[HSP]:
    """Add stitched alignments for HSP runs split by a frame change.

    A frameshift in the subject splits one protein alignment into two HSPs of
    the same reference on the same contig and strand: consecutive in reference
    coordinates, consecutive in contig coordinates, and in *different* reading
    frames. Stitching them recovers one alignment across the whole reference
    and flags it `frameshift=True` — the FRAMESHIFT arm of the status ladder
    (StxTyper does this in `processDisruptions`).

    The stitched alignment is *offered alongside* its parts rather than
    replacing them, so locus reduction decides on identity and coverage
    whether the frameshift reading is the better explanation. That is what
    keeps a short terminal frame change reported as EXTENDED, matching
    StxTyper, instead of overriding a full-length ungapped alignment.

    HSPs in the same frame are left alone: those are ordinary gapped
    alignments, and locus reduction picks the better one.
    """
    by_reference: Dict[Tuple[str, str, str, str], List[HSP]] = {}
    for hsp in hsps:
        key = (hsp.contig, hsp.strand, hsp.subunit, hsp.reference_accession)
        by_reference.setdefault(key, []).append(hsp)

    stitched: List[HSP] = []
    for group in by_reference.values():
        if len(group) < 2:
            continue
        # Reference order, which is contig order once strand is accounted for.
        group = sorted(group, key=lambda h: h.qstart)
        current = None
        for nxt in group:
            if current is not None and _stitchable(current, nxt, max_gap):
                current = _stitch(current, nxt)
                stitched.append(current)
            else:
                current = nxt
    return list(hsps) + stitched


def _stitchable(a: HSP, b: HSP, max_gap: int) -> bool:
    """Whether `b` continues `a`'s alignment across a frame change.

    BLAST usually reports the two blocks with a residue or two of overlap in
    both reference and contig coordinates, so exact adjacency is the wrong
    test; StxTyper allows a 20-unit window when merging.
    """
    if b.qstart <= a.qstart or b.qend <= a.qend:
        return False
    if b.qstart - a.qend - 1 < -MERGE_WINDOW:
        return False
    if a.strand == "+":
        offset = b.stop - a.stop
        gap = b.start - a.stop - 1
    else:
        offset = a.start - b.start
        gap = a.start - b.stop - 1
    if offset <= 0 or gap < -CODON * MERGE_WINDOW or gap > max_gap:
        return False
    return a.frame() != b.frame()


def _stitch(a: HSP, b: HSP) -> HSP:
    """Concatenate two co-linear HSPs into one frameshifted alignment.

    Three adjustments make the join honest:
      * `b`'s leading columns are dropped where they re-cover reference
        positions `a` already aligned;
      * reference positions neither block covers are carried as query-only
        columns, so `residue_at` keeps projecting reference coordinates
        correctly across the join;
      * one column is dropped at the join, because the codon straddling the
        frame change is not interpretable in either frame — it is taken from
        the downstream block, whose first codon is the shifted one. StxTyper's
        merge does the same, which is why its reported alignment length is one
        shorter than the reference span it covers.
    """
    qseq, sseq = _drop_query_prefix(b.qseq, b.sseq, b.qstart, a.qend)
    skipped = max(0, b.qstart - a.qend - 1)
    qseq = a.qseq + "X" * skipped + qseq[1:]
    sseq = a.sseq + "-" * skipped + sseq[1:]
    return HSP(
        subunit=a.subunit,
        contig=a.contig,
        strand=a.strand,
        start=min(a.start, b.start),
        stop=max(a.stop, b.stop),
        ref_subtype=a.ref_subtype,
        ref_class=a.ref_class,
        nident=sum(1 for q, s in zip(qseq, sseq) if q == s and q != "-"),
        aln_len=len(qseq),
        ref_len=a.ref_len,
        reference_accession=a.reference_accession,
        contig_len=a.contig_len,
        qstart=a.qstart,
        qend=b.qend,
        qseq=qseq,
        sseq=sseq,
        frameshift=True,
        internal_stop=a.internal_stop or b.internal_stop,
        stop_codon=b.stop_codon,
        ambiguous=a.ambiguous + b.ambiguous,
    )


def _drop_query_prefix(
    qseq: str, sseq: str, qstart: int, keep_after: int
) -> Tuple[str, str]:
    """Trim alignment columns covering reference positions <= `keep_after`."""
    qpos = qstart
    for i, qchar in enumerate(qseq):
        if qpos > keep_after:
            return qseq[i:], sseq[i:]
        if qchar != "-":
            qpos += 1
    return "", ""


def _hsp_rank(hsp: HSP) -> tuple:
    """Best-first ordering within a locus (StxTyper's `BlastAlignment::less`):
    identity, then relative reference coverage, then how many reference
    residues were actually aligned — with equally identical, equally covered
    references of different lengths, the longer alignment explains more of the
    contig — then accession for determinism."""
    return (
        -hsp.identity, -hsp.coverage, -hsp.ref_covered, hsp.reference_accession
    )


def reduce_to_loci(hsps: Sequence[HSP]) -> List[HSP]:
    """Collapse a reference set's worth of HSPs to one per locus per class.

    A reference set holds many accessions per subtype (stx: 160 proteins for
    19 subtypes), so every subunit locus produces dozens of overlapping HSPs.
    Keeping the best HSP per (contig, strand, subunit, class) locus is what
    makes downstream pairing produce one candidate per real operon instead of
    the cross product; keeping it *per class* rather than globally preserves
    the class-agreeing candidates the first pairing pass needs (StxTyper's
    `paretoBest` over the same grouping key).
    """
    groups: Dict[Tuple[str, str, str, str], List[HSP]] = {}
    for hsp in hsps:
        key = (hsp.contig, hsp.strand, hsp.subunit, hsp.ref_class)
        groups.setdefault(key, []).append(hsp)

    kept: List[HSP] = []
    for group in groups.values():
        for locus in _cluster_by_overlap(group):
            kept.append(min(locus, key=_hsp_rank))
    return kept


def _cluster_by_overlap(hsps: Sequence[HSP]) -> List[List[HSP]]:
    """Group HSPs into loci of overlapping intervals on the same contig strand.

    Coordinates only mean anything within one contig and strand, so those are
    part of the grouping key — clustering on coordinates alone would merge
    unrelated loci that happen to sit at the same offset on different contigs.
    """
    by_placement: Dict[Tuple[str, str], List[HSP]] = {}
    for hsp in hsps:
        by_placement.setdefault((hsp.contig, hsp.strand), []).append(hsp)

    clusters: List[List[HSP]] = []
    for placed in by_placement.values():
        span = None
        for hsp in sorted(placed, key=lambda h: (h.start, h.stop)):
            if span is not None and hsp.start <= span:
                clusters[-1].append(hsp)
                span = max(span, hsp.stop)
            else:
                clusters.append([hsp])
                span = hsp.stop
    return clusters


@dataclass
class Candidate:
    """A candidate operon: one HSP per subunit tag it could place."""

    subunits: Dict[str, HSP]
    complete: bool  # every subunit in subunit_order is present
    pass_name: str = "strict"

    @property
    def hsps(self) -> List[HSP]:
        return list(self.subunits.values())

    @property
    def contig(self) -> str:
        return self.hsps[0].contig

    @property
    def strand(self) -> str:
        return self.hsps[0].strand

    @property
    def start(self) -> int:
        return min(h.start for h in self.hsps)

    @property
    def stop(self) -> int:
        return max(h.stop for h in self.hsps)

    @property
    def identity(self) -> float:
        return combined_identity(
            {tag: (h.nident, h.aln_len) for tag, h in self.subunits.items()}
        )

    @property
    def coverage(self) -> float:
        covered = sum(h.ref_covered for h in self.hsps)
        ref_len = sum(h.ref_len for h in self.hsps)
        return covered / ref_len if ref_len else 0.0

    @property
    def perfect(self) -> bool:
        return all(
            h.reference_complete and not h.frameshift for h in self.hsps
        )

    def inside(self, other: "Candidate", slack: int) -> bool:
        return (
            self.contig == other.contig
            and self.strand == other.strand
            and self.start + slack >= other.start
            and self.stop <= other.stop + slack
        )

    def to_dict(self) -> dict:
        return {
            "subunits": {tag: h.to_dict() for tag, h in self.subunits.items()},
            "complete": self.complete,
            "pass": self.pass_name,
        }

    @classmethod
    def from_dict(cls, raw: dict) -> "Candidate":
        return cls(
            subunits={
                tag: HSP.from_dict(h) for tag, h in raw["subunits"].items()
            },
            complete=raw["complete"],
            pass_name=raw.get("pass", "strict"),
        )


def pair_operon(
    hsps_by_subunit: Dict[str, List[HSP]],
    subunit_order: Sequence[str],
    intergenic_max: int,
    require_same_strand: bool = True,
    require_same_contig: bool = True,
    intergenic_min: int = 0,
) -> List[List[HSP]]:
    """Every syntenic chain across the subunits in `subunit_order` (§4).

    Returns candidate operons — lists of HSPs, one per subunit, ordered as
    `subunit_order` — each satisfying same contig/strand (if required),
    correct transcription order, and an intergenic distance within
    [`intergenic_min`, `intergenic_max`] between every pair of consecutive
    subunits. `intergenic_min` is negative for schemes whose genes overlap.

    All combinations are enumerated rather than greedily claimed: which
    reference pairing best explains a locus is an identity question, and it is
    `select_operons` that answers it. Locus reduction has already collapsed
    the reference set to one HSP per locus per class, so the enumeration is
    bounded by the number of classes, not the size of the reference set.
    Chains are returned best-first by the leading subunit's rank.
    """
    if len(subunit_order) < 2:
        raise OperonConfigError("subunit_order must have at least 2 subunits.")

    chains: List[List[HSP]] = [
        [seed] for seed in sorted(hsps_by_subunit.get(subunit_order[0], []), key=_hsp_rank)
    ]
    for tag in subunit_order[1:]:
        extended: List[List[HSP]] = []
        for chain in chains:
            prev = chain[-1]
            for cand in sorted(hsps_by_subunit.get(tag, []), key=_hsp_rank):
                if require_same_contig and cand.contig != prev.contig:
                    continue
                if require_same_strand and cand.strand != prev.strand:
                    continue
                gap = _upstream_gap(prev, cand)
                if gap is None or gap > intergenic_max or gap < intergenic_min:
                    continue
                extended.append(chain + [cand])
        chains = extended
    return chains


# Pairing passes, in StxTyper's order (`goodBlasts2operons` calls in
# stxtyper.cpp): class-agreeing operons over the identity threshold first,
# then threshold-passing operons of mixed class, then perfect alignments, and
# finally anything syntenic — the last pass with a relaxed intergenic
# distance, which is what recovers partial operons.
#   name, same_class, gate, relax
PAIRING_PASSES = (
    ("same_class", True, "threshold", False),
    ("strong", False, "threshold", False),
    ("perfect", False, "perfect", False),
    ("any", False, "none", True),
)


def assemble_candidates(hsps: Sequence[HSP], operon_cfg: dict) -> List[Candidate]:
    """Frameshift-stitch, reduce to loci, then pair into candidate operons.

    Each pass only sees HSPs no earlier pass used, so a locus is explained at
    the highest stringency that accepts it and is not re-reported by a looser
    pass. Nothing below `min_operon_identity` is ever accepted: claiming an
    HSP for a candidate that will not be reported would silently lose the
    lone-subunit call it should have produced instead.
    """
    subunit_order = list(operon_cfg["subunit_order"])
    if len(subunit_order) < 2:
        raise OperonConfigError("operon.subunit_order must list at least 2 subunits.")

    intergenic_max = operon_cfg.get("intergenic_max", DEFAULT_INTERGENIC_MAX)
    # Negative for schemes whose genes overlap (ETEC's eltA/eltB overlap 4 bp).
    intergenic_min = operon_cfg.get("intergenic_min", 0)
    relax_factor = operon_cfg.get("intergenic_relax_factor", 2)
    thresholds = operon_cfg.get("identity_thresholds", {})
    min_identity = operon_cfg.get("min_operon_identity", DEFAULT_MIN_OPERON_IDENTITY)
    require_same_strand = operon_cfg.get("require_same_strand", True)
    require_same_contig = operon_cfg.get("require_same_contig", True)

    available = reduce_to_loci(
        stitch_frameshifts(hsps, max_gap=intergenic_max * relax_factor)
    )

    candidates: List[Candidate] = []
    for pass_name, same_class, gate, relax in PAIRING_PASSES:
        gap_max = intergenic_max * (relax_factor if relax else 1)
        accepted: List[Candidate] = []
        for chain in _pairs_for_pass(
            available, subunit_order, gap_max, same_class,
            require_same_strand, require_same_contig, intergenic_min,
        ):
            candidate = Candidate(
                subunits=dict(zip(subunit_order, chain)),
                complete=True,
                pass_name=pass_name,
            )
            if candidate.identity < min_identity:
                continue
            if not _passes_gate(candidate, gate, thresholds):
                continue
            accepted.append(candidate)
        candidates.extend(accepted)
        available = [
            h for h in available
            if not any(_spans(c, h) for c in accepted)
        ]

    # Leftover lone subunits: still reported, as single-subunit candidates.
    for tag in subunit_order:
        for hsp in _best_per_locus([h for h in available if h.subunit == tag]):
            candidates.append(
                Candidate(subunits={tag: hsp}, complete=False, pass_name="single")
            )
    return candidates


def _spans(candidate: Candidate, hsp: HSP) -> bool:
    return (
        hsp.contig == candidate.contig
        and hsp.strand == candidate.strand
        and hsp.start >= candidate.start
        and hsp.stop <= candidate.stop
    )


def _pairs_for_pass(
    available: Sequence[HSP],
    subunit_order: Sequence[str],
    gap_max: int,
    same_class: bool,
    require_same_strand: bool,
    require_same_contig: bool,
    gap_min: int = 0,
) -> List[List[HSP]]:
    by_subunit: Dict[str, List[HSP]] = {tag: [] for tag in subunit_order}
    for hsp in available:
        if hsp.subunit in by_subunit:
            by_subunit[hsp.subunit].append(hsp)

    if not same_class:
        return pair_operon(
            by_subunit, subunit_order, gap_max,
            require_same_strand, require_same_contig, gap_min,
        )

    chains: List[List[HSP]] = []
    classes = {h.ref_class for h in available}
    for class_label in sorted(classes):
        per_class = {
            tag: [h for h in hsps if h.ref_class == class_label]
            for tag, hsps in by_subunit.items()
        }
        chains.extend(
            pair_operon(
                per_class, subunit_order, gap_max,
                require_same_strand, require_same_contig, gap_min,
            )
        )
    return chains


def _passes_gate(candidate: Candidate, gate: str, thresholds: Dict[str, float]) -> bool:
    if gate == "none":
        return True
    if gate == "perfect":
        return candidate.perfect
    identity = candidate.identity
    return all(
        identity >= select_threshold(h.ref_class, thresholds)
        for h in candidate.hsps
    )


def _best_per_locus(hsps: Sequence[HSP]) -> List[HSP]:
    return [
        min(locus, key=_hsp_rank) for locus in _cluster_by_overlap(hsps)
    ]


def select_operons(candidates: Sequence[Candidate], operon_cfg: dict) -> List[Candidate]:
    """Drop weak and redundant candidates (StxTyper's `goodOperons` pass).

    Candidates below `min_operon_identity` are not reported at all, and a
    candidate contained within an already-kept, at-least-as-identical operon
    is redundant — that is what stops one locus from being reported once per
    reference class and once per lone subunit.
    """
    min_identity = operon_cfg.get("min_operon_identity", DEFAULT_MIN_OPERON_IDENTITY)
    slack = operon_cfg.get("overlap_slack", DEFAULT_OVERLAP_SLACK)

    ranked = sorted(
        (c for c in candidates if c.identity >= min_identity),
        key=lambda c: (not c.complete, -c.identity, -c.coverage, c.contig, c.start),
    )

    kept: List[Candidate] = []
    for candidate in ranked:
        if any(
            better.inside(candidate, 3 * slack) or candidate.inside(better, 3 * slack)
            for better in kept
        ):
            continue
        kept.append(candidate)
    return kept


def call_operon_status(
    frameshift: bool,
    internal_stop: bool,
    contig_end: bool,
    below_threshold: bool,
    class_mismatch: bool,
    ambiguous: bool,
    extended: bool,
    partial: bool,
) -> str:
    """Priority-ordered disruption/completeness status (§5).

    Order mirrors StxTyper's `Operon::saveTsvOut` quality ladder: a
    frameshift or internal stop outranks a truncation, truncation at a contig
    end outranks a C-terminal extension, and any incompleteness outranks the
    identity-based COMPLETE/COMPLETE_NOVEL/AMBIGUOUS distinction.
    """
    if frameshift:
        return "FRAMESHIFT"
    if internal_stop:
        return "INTERNAL_STOP"
    if contig_end:
        return "PARTIAL_CONTIG_END"
    if extended:
        return "EXTENDED"
    if partial:
        return "PARTIAL"
    if class_mismatch or below_threshold:
        return "AMBIGUOUS" if ambiguous else "COMPLETE_NOVEL"
    return "COMPLETE"


def call_operon(candidate: Candidate, operon_cfg: dict) -> dict:
    """Score one candidate operon into the §5 output shape.

    Type resolution follows StxTyper's `Operon::getStxType`: classes must
    agree, or the call falls back to their shared parent class; a generalized
    class is resolved through its residue decision table; and anything short
    of COMPLETE is reported at parent-class resolution, because a disrupted
    operon does not support a subtype-level claim.
    """
    subunit_order = list(operon_cfg["subunit_order"])
    thresholds = operon_cfg.get("identity_thresholds", {})
    generalized = operon_cfg.get("generalized_classes", {})
    residue_rules = build_residue_rules(operon_cfg.get("residue_rules", []))

    subunits = candidate.subunits
    hsps = candidate.hsps
    complete_operon = candidate.complete and set(subunits) == set(subunit_order)

    classes = {h.ref_class for h in hsps}
    class_mismatch = len(classes) != 1
    class_label = next(iter(classes)) if not class_mismatch else None
    superclasses = {superclass_of(c, operon_cfg) for c in classes}
    superclass = next(iter(superclasses)) if len(superclasses) == 1 else None

    identity = candidate.identity
    residue_evidence: Dict[str, Optional[str]] = {}
    resolved = class_label
    unresolved_generalized = False

    if class_label in generalized:
        rule = residue_rules.get(class_label)
        if rule is None:
            raise OperonConfigError(
                f"Generalized class {class_label!r} has no matching residue rule."
            )
        residues = []
        for position in rule.positions:
            hsp = subunits.get(position["subunit"])
            index = position["index"]
            residue = hsp.residue_at(index) if hsp else None
            residues.append(residue)
            residue_evidence["{}{}".format(position["subunit"], index)] = residue
        resolved = rule.resolve(residues)
        unresolved_generalized = resolved == rule.fallback and resolved in generalized

    threshold_applied = None
    below_threshold = False
    if class_label is not None:
        threshold_applied = select_threshold(class_label, thresholds)
        below_threshold = identity < threshold_applied

    status = call_operon_status(
        frameshift=any(h.frameshift for h in hsps),
        internal_stop=any(h.internal_stop for h in hsps),
        contig_end=any(h.contig_end for h in hsps)
        or _partner_would_not_fit(candidate, operon_cfg),
        # An unresolved generalized class is a novel call in its own right:
        # the operon is intact but its subtype is not in the scheme
        # (StxTyper's `stxType.size() <= 1` novelty test).
        below_threshold=below_threshold or unresolved_generalized,
        class_mismatch=class_mismatch,
        ambiguous=any(h.ambiguous for h in hsps),
        extended=any(h.extended for h in hsps),
        partial=not complete_operon
        or not all(h.reference_complete for h in hsps),
    )

    # Only a COMPLETE operon supports a subtype-level claim; everything else
    # is reported at parent-class resolution.
    if status == "COMPLETE":
        type_token = resolved
    else:
        type_token = superclass if superclass is not None else resolved
    if class_mismatch and superclass is None:
        type_token = None

    intergenic_bp = None
    if complete_operon and len(subunit_order) >= 2:
        ordered = [subunits[tag] for tag in subunit_order]
        gaps = [
            _upstream_gap(a, b) for a, b in zip(ordered, ordered[1:])
        ]
        if all(gap is not None for gap in gaps):
            intergenic_bp = sum(gaps)

    return {
        "type_token": type_token,
        "class": class_label,
        "operon_status": status,
        "status": status_to_coarse(status),
        "combined_identity": identity,
        "threshold_applied": threshold_applied,
        "residue_evidence": residue_evidence,
        "unresolved_generalized": unresolved_generalized,
        "contig": candidate.contig,
        "start": candidate.start,
        "stop": candidate.stop,
        "strand": candidate.strand,
        "intergenic_bp": intergenic_bp,
        "subunits": {
            tag: {
                "reference": hsp.reference_accession,
                "reference_subtype": hsp.ref_subtype,
                "identity": hsp.identity,
                "coverage": hsp.coverage,
            }
            for tag, hsp in subunits.items()
        },
    }


def _partner_would_not_fit(candidate: Candidate, operon_cfg: dict) -> bool:
    """For a lone subunit, whether its missing partners had no room on the
    contig — StxTyper's `otherTruncated()`, which is what separates
    PARTIAL_CONTIG_END from PARTIAL for single-subunit calls."""
    subunit_order = list(operon_cfg["subunit_order"])
    if candidate.complete and set(candidate.subunits) == set(subunit_order):
        return False
    intergenic_max = operon_cfg.get("intergenic_max", DEFAULT_INTERGENIC_MAX)
    missed_max = intergenic_max + MIN_PARTNER_LENGTH
    for tag, hsp in candidate.subunits.items():
        position = subunit_order.index(tag)
        upstream_missing = position > 0
        downstream_missing = position < len(subunit_order) - 1
        head_room = hsp.start - 1
        tail_room = hsp.contig_len - hsp.stop if hsp.contig_len else missed_max + 1
        if hsp.strand == "-":
            head_room, tail_room = tail_room, head_room
        if upstream_missing and head_room <= missed_max:
            return True
        if downstream_missing and tail_room <= missed_max:
            return True
    return False


def report_operons(
    hsps: Sequence[HSP],
    operon_cfg: dict,
    scheme: str = "",
    subtype_prefix: str = "",
) -> List[dict]:
    """Full pipeline: HSPs -> reported operon calls in the §5 output shape."""
    return format_calls(
        assemble_candidates(hsps, operon_cfg), operon_cfg, scheme, subtype_prefix
    )


def format_calls(
    candidates: Sequence[Candidate],
    operon_cfg: dict,
    scheme: str = "",
    subtype_prefix: str = "",
) -> List[dict]:
    """Select the reportable candidates and render them in the §5 shape.

    `subtype_prefix` is the common prefix of the scheme's subtype labels (stx:
    "stx"), so a resolved class token becomes the profile id a curator would
    recognise ("2c" -> "stx2c") without the algorithm knowing the scheme.
    """
    candidates = select_operons(candidates, operon_cfg)
    results = []
    for candidate in candidates:
        call = call_operon(candidate, operon_cfg)
        token = call["type_token"]
        results.append({
            # Subunits from unrelated classes support no type claim beyond the
            # scheme itself, but the operon is still reported.
            "profile_id": subtype_prefix + (token or ""),
            "profile_type": "operon_subtype",
            "scheme": scheme,
            "status": call["status"],
            "operon_status": call["operon_status"],
            "confidence": max(0.0, min(1.0, call["combined_identity"])),
            "operon": {
                "contig": call["contig"],
                "start": call["start"],
                "stop": call["stop"],
                "strand": call["strand"],
                "intergenic_bp": call["intergenic_bp"],
                "combined_identity": call["combined_identity"],
                "threshold_applied": call["threshold_applied"],
                "subunits": call["subunits"],
                "residue_evidence": call["residue_evidence"],
            },
            "method": {"typing_model": "operon", "tools": ["tblastn"]},
        })
    results.sort(
        key=lambda r: (
            STATUS_RANK.get(r["operon_status"], len(STATUS_PRIORITY)),
            r["operon"]["contig"],
            r["operon"]["start"],
        )
    )
    return results


def subtype_prefix(profile_rows: Sequence[dict]) -> str:
    """Longest common prefix of a scheme's subtype labels (stx: "stx")."""
    labels = [row["subtype"] for row in profile_rows if row.get("subtype")]
    if not labels:
        return ""
    prefix = labels[0]
    for label in labels[1:]:
        while prefix and not label.startswith(prefix):
            prefix = prefix[:-1]
    return prefix


def validate_operon_metadata(operon_cfg: dict, profile_rows: List[dict]) -> None:
    """Validate an `[operon]` metadata.toml block against a torch's
    profiles.tsv rows (§3.2):

    - `subunit_order` lists at least 2 subunits
    - every `class` referenced by a profiles.tsv row has an identity
      threshold, explicit or via `default`
    - every class in `[operon.generalized_classes]` has a matching
      `[[operon.residue_rules]]` entry
    """
    subunit_order = operon_cfg.get("subunit_order", [])
    if len(subunit_order) < 2:
        raise OperonConfigError("operon.subunit_order must list at least 2 subunits.")

    thresholds = operon_cfg.get("identity_thresholds", {})
    generalized = operon_cfg.get("generalized_classes", {})
    residue_rule_classes = {r["class"] for r in operon_cfg.get("residue_rules", [])}

    classes_seen = {row["class"] for row in profile_rows if row.get("class")}
    for class_label in classes_seen:
        select_threshold(class_label, thresholds)  # raises if unresolvable

    for class_label in generalized:
        if class_label not in residue_rule_classes:
            raise OperonConfigError(
                f"Generalized class {class_label!r} has no matching "
                f"[[operon.residue_rules]] entry."
            )
