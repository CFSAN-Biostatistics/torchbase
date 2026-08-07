"""Core operon-typing algorithm: threshold selection, residue-table
resolution, HSP synteny pairing, and status-ladder scoring.

Reference implementation for the `operon` typing model — see
docs/operon-strategy-plan.md §2-§5 for the design rationale (protein-space,
synteny-aware calling of multi-subunit operons, generalized from NCBI
StxTyper).

`torchfs.Torch._load_single_scheme` uses `validate_operon_metadata` to check
a torch's `[operon]` metadata.toml block against its profiles.tsv at load
time. `torchbase/workflows/builtin/operon_typing.wdl` tasks reimplement this
same algorithm inline in Python heredocs, because WDL tasks run in isolated
containers with no torchbase install (the same convention already used by
`tasks/minhash.wdl`, `tasks/alignment.wdl`, etc.) — keep the two in lockstep
when either changes.
"""

from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Tuple

# Disruption/completeness status ladder, most to least complete. This is
# also the priority order used to resolve overlapping candidate operons.
STATUS_PRIORITY = (
    "COMPLETE",
    "COMPLETE_NOVEL",
    "AMBIGUOUS",
    "PARTIAL",
    "PARTIAL_CONTIG_END",
    "EXTENDED",
    "INTERNAL_STOP",
    "FRAMESHIFT",
)

STATUS_RANK = {status: i for i, status in enumerate(STATUS_PRIORITY)}


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

    `subunit_idents` maps subunit tag -> (nident, aln_len).
    """
    nident = sum(n for n, _ in subunit_idents.values())
    length = sum(n for _, n in subunit_idents.values())
    if length == 0:
        return 0.0
    return nident / length


@dataclass
class ResidueRule:
    class_label: str
    positions: List[dict]  # [{"subunit": "A", "index": 312}, ...]
    table: List[dict]      # [{"call": "2a", "residues": [["F","S"], ...]}, ...]
    fallback: str

    def resolve(self, residues: Sequence[str]) -> str:
        """Resolve observed `residues` (same order as `positions`) against
        the decision table; unresolved rows fall through to `fallback`."""
        for row in self.table:
            row_residues = row["residues"]
            if len(row_residues) != len(residues):
                continue
            if all(
                observed in alternatives
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
    residues: Sequence[str],
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


@dataclass
class HSP:
    """A single translated-search hit, normalized at the protein_search.wdl
    task boundary (§4) — never leak BLAST tabular format past this shape."""

    subunit: str
    contig: str
    strand: str
    start: int  # nucleotide coordinates, min(start, stop)
    stop: int
    ref_subtype: str
    ref_class: str
    nident: int
    aln_len: int
    ref_len: int
    frameshift: bool = False
    internal_stop: bool = False
    contig_end: bool = False


def _upstream_gap(a: HSP, b: HSP) -> Optional[int]:
    """Intergenic distance from `a`'s downstream end to `b`'s upstream end,
    respecting strand direction. None if `b` does not follow `a`."""
    if a.strand != b.strand:
        return None
    if a.strand == "+":
        if b.start <= a.stop:
            return None
        return b.start - a.stop - 1
    else:
        if a.start <= b.stop:
            return None
        return a.start - b.stop - 1


def pair_operon(
    hsps_by_subunit: Dict[str, List[HSP]],
    subunit_order: Sequence[str],
    intergenic_max: int,
    require_same_strand: bool = True,
    require_same_contig: bool = True,
) -> List[List[HSP]]:
    """Greedy synteny pairing across subunits in `subunit_order` (§4).

    Returns candidate operons — lists of HSPs, one per subunit, ordered as
    `subunit_order` — each satisfying same contig/strand (if required),
    correct transcription order, and intergenic distance <= `intergenic_max`
    between every pair of consecutive subunits. Each HSP is claimed by at
    most one candidate; first-subunit seeds are tried highest-identity
    first, mirroring StxTyper's `reported`-flag greedy claim.
    """
    if len(subunit_order) < 2:
        raise OperonConfigError("subunit_order must have at least 2 subunits.")

    claimed = set()
    candidates: List[List[HSP]] = []

    first_tag = subunit_order[0]
    seeds = sorted(
        hsps_by_subunit.get(first_tag, []),
        key=lambda h: (h.nident / h.aln_len) if h.aln_len else 0.0,
        reverse=True,
    )

    for seed in seeds:
        if id(seed) in claimed:
            continue
        chain = [seed]
        ok = True
        for tag in subunit_order[1:]:
            prev = chain[-1]
            best, best_gap = None, None
            for cand in hsps_by_subunit.get(tag, []):
                if id(cand) in claimed:
                    continue
                if require_same_contig and cand.contig != prev.contig:
                    continue
                if require_same_strand and cand.strand != prev.strand:
                    continue
                gap = _upstream_gap(prev, cand)
                if gap is None or gap > intergenic_max:
                    continue
                if best_gap is None or gap < best_gap:
                    best, best_gap = cand, gap
            if best is None:
                ok = False
                break
            chain.append(best)
        if ok:
            for hsp in chain:
                claimed.add(id(hsp))
            candidates.append(chain)

    return candidates


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
    """Priority-ordered disruption/completeness status (§5)."""
    if frameshift:
        return "FRAMESHIFT"
    if internal_stop:
        return "INTERNAL_STOP"
    if extended:
        return "EXTENDED"
    if partial and contig_end:
        return "PARTIAL_CONTIG_END"
    if partial:
        return "PARTIAL"
    if class_mismatch or below_threshold:
        return "AMBIGUOUS" if ambiguous else "COMPLETE_NOVEL"
    return "COMPLETE"


def call_operon(candidate: Dict[str, HSP], operon_cfg: dict) -> dict:
    """Score one paired candidate operon (one HSP per subunit tag) into the
    §5 output shape. Residue-table resolution for generalized classes is
    intentionally out of scope here — it requires the reference-coordinate
    alignment projection (§3.1, §7 risk 1), which only exists once an
    alignment is available; use `resolve_generalized_class` once residues
    have been read."""
    thresholds = operon_cfg.get("identity_thresholds", {})
    subunit_order = operon_cfg["subunit_order"]

    ref_classes = {tag: candidate[tag].ref_class for tag in subunit_order}
    distinct_classes = set(ref_classes.values())
    class_mismatch = len(distinct_classes) != 1
    class_label = next(iter(distinct_classes)) if not class_mismatch else None

    ident = combined_identity(
        {tag: (h.nident, h.aln_len) for tag, h in candidate.items()}
    )

    frameshift = any(h.frameshift for h in candidate.values())
    internal_stop = any(h.internal_stop for h in candidate.values())
    contig_end = any(h.contig_end for h in candidate.values())
    extended = any(h.aln_len > h.ref_len for h in candidate.values())
    partial = any(h.aln_len < h.ref_len for h in candidate.values())

    below_threshold = False
    threshold_applied = None
    if not class_mismatch:
        threshold_applied = select_threshold(class_label, thresholds)
        below_threshold = ident < threshold_applied

    status = call_operon_status(
        frameshift=frameshift,
        internal_stop=internal_stop,
        contig_end=contig_end,
        below_threshold=below_threshold,
        class_mismatch=class_mismatch,
        ambiguous=False,
        extended=extended,
        partial=partial,
    )

    return {
        "class": class_label,
        "operon_status": status,
        "status": status_to_coarse(status),
        "combined_identity": ident,
        "threshold_applied": threshold_applied,
    }


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
        select_threshold(class_label, thresholds)  # raises OperonConfigError if unresolvable

    for class_label in generalized:
        if class_label not in residue_rule_classes:
            raise OperonConfigError(
                f"Generalized class {class_label!r} has no matching "
                f"[[operon.residue_rules]] entry."
            )
