version 1.0

# Scores each candidate operon from operon_assembly.wdl into the §5 output
# shape: combined identity vs per-class threshold, residue-table
# resolution for generalized classes, and the eight-value disruption status
# ladder. Mirrors torchbase/operon.py (kept inline per the tasks/*.wdl
# container-isolation convention — see operon_assembly.wdl header).
task call_operons {
    input {
        File candidates
        File operon_config_json
        File profiles_table
        String scheme = ""
    }

    command <<<
        python3 <<'PYTHON_SCRIPT'
import json
import os

STATUS_PRIORITY = (
    "COMPLETE", "COMPLETE_NOVEL", "AMBIGUOUS", "PARTIAL",
    "PARTIAL_CONTIG_END", "EXTENDED", "INTERNAL_STOP", "FRAMESHIFT",
)


def status_to_coarse(operon_status):
    if operon_status == "COMPLETE":
        return "known"
    if operon_status in ("COMPLETE_NOVEL", "AMBIGUOUS"):
        return "novel_profile"
    return "incomplete"


def select_threshold(class_label, thresholds):
    if class_label in thresholds:
        return thresholds[class_label]
    return thresholds.get("default")


def combined_identity(subunits):
    nident = sum(h["nident"] for h in subunits.values())
    length = sum(h["aln_len"] for h in subunits.values())
    return nident / length if length else 0.0


def project_residue(hsp, index):
    """Reference-anchored residue projection (§3.1, §7 risk 1-2): walk the
    gapped query/subject alignment together, tracking the reference (query)
    offset, and return the subject residue observed at `index` — or None if
    that reference position isn't covered by this HSP (gap, or outside the
    aligned block)."""
    qpos = hsp.get("qstart")
    qseq, sseq = hsp.get("qseq"), hsp.get("sseq")
    if qpos is None or qseq is None or sseq is None:
        return None
    for qc, sc in zip(qseq, sseq):
        if qc != "-":
            if qpos == index:
                return None if sc == "-" else sc
            qpos += 1
    return None


def resolve_generalized_class(class_label, subunits, generalized_classes, residue_rules):
    """Resolve a generalized class (e.g. "2") to a specific subtype (e.g.
    "2a") via its residue decision table. Returns (resolved_label, evidence,
    resolved) — `resolved` is False when the table fell through to the
    fallback value, which the caller uses to flag COMPLETE_NOVEL."""
    if class_label not in generalized_classes:
        return class_label, {}, True

    rule = residue_rules.get(class_label)
    if rule is None:
        return class_label, {}, False

    residues = []
    evidence = {}
    for pos in rule["positions"]:
        subunit_tag, index = pos["subunit"], pos["index"]
        hsp = subunits.get(subunit_tag)
        residue = project_residue(hsp, index) if hsp else None
        residues.append(residue)
        evidence[f"{subunit_tag}{index}"] = residue

    for row in rule["table"]:
        row_residues = row["residues"]
        if len(row_residues) != len(residues):
            continue
        if all(
            r is not None and r in alts
            for r, alts in zip(residues, row_residues)
        ):
            return row["call"], evidence, True

    return rule.get("fallback", class_label), evidence, False


def call_operon_status(frameshift, internal_stop, contig_end, below_threshold,
                        class_mismatch, unresolved_generalized, extended, partial):
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
    if class_mismatch or below_threshold or unresolved_generalized:
        return "COMPLETE_NOVEL"
    return "COMPLETE"


with open("~{candidates}") as f:
    candidates = json.load(f)

with open("~{operon_config_json}") as f:
    cfg = json.load(f)

subunit_order = cfg["subunit_order"]
thresholds = cfg.get("identity_thresholds", {})
generalized_classes = cfg.get("generalized_classes", {})
residue_rules = {r["class"]: r for r in cfg.get("residue_rules", [])}

provided_scheme = "~{scheme}"
scheme = provided_scheme if provided_scheme else os.path.basename(
    os.path.dirname("~{profiles_table}")
)

results = []

for entry in candidates:
    subunit_list = entry["subunits"]
    subunits = {h["subunit"]: h for h in subunit_list}
    complete = entry["complete"] and set(subunits) == set(subunit_order)

    frameshift = any(h.get("frameshift") for h in subunit_list)
    internal_stop = any(h.get("internal_stop") for h in subunit_list)
    contig_end = any(h.get("contig_end") for h in subunit_list)
    extended = any(h["aln_len"] > h["ref_len"] for h in subunit_list)
    partial = (not complete) or any(h["aln_len"] < h["ref_len"] for h in subunit_list)

    distinct_classes = {h["ref_class"] for h in subunit_list}
    class_mismatch = complete and len(distinct_classes) != 1
    class_label = next(iter(distinct_classes)) if len(distinct_classes) == 1 else None

    ident = combined_identity(subunits) if subunit_list else 0.0
    threshold_applied = None
    below_threshold = False
    if complete and not class_mismatch:
        threshold_applied = select_threshold(class_label, thresholds)
        below_threshold = threshold_applied is None or ident < threshold_applied

    resolved_subtype, residue_evidence, resolved = (class_label, {}, True)
    unresolved_generalized = False
    if complete and not class_mismatch and class_label is not None:
        resolved_subtype, residue_evidence, resolved = resolve_generalized_class(
            class_label, subunits, generalized_classes, residue_rules
        )
        unresolved_generalized = not resolved

    status = call_operon_status(
        frameshift, internal_stop, contig_end, below_threshold,
        class_mismatch, unresolved_generalized, extended, partial,
    )

    # Intergenic distance across the paired chain (sum of consecutive gaps).
    intergenic_bp = None
    if complete and len(subunit_list) >= 2:
        gaps = []
        ordered = [subunits[tag] for tag in subunit_order if tag in subunits]
        for a, b in zip(ordered, ordered[1:]):
            if a["strand"] == "+":
                gaps.append(b["start"] - a["stop"] - 1)
            else:
                gaps.append(a["start"] - b["stop"] - 1)
        intergenic_bp = sum(gaps)

    contigs_involved = {h["contig"] for h in subunit_list}
    strands_involved = {h["strand"] for h in subunit_list}

    subunit_out = {
        tag: {
            "reference": h.get("reference_accession"),
            "reference_subtype": h.get("ref_subtype"),
            "identity": (h["nident"] / h["aln_len"]) if h["aln_len"] else 0.0,
            "coverage": (h["aln_len"] / h["ref_len"]) if h["ref_len"] else 0.0,
        }
        for tag, h in subunits.items()
    }

    result = {
        "profile_id": resolved_subtype if complete else "unknown",
        "profile_type": "operon_subtype",
        "scheme": scheme,
        "status": status_to_coarse(status),
        "operon_status": status,
        "confidence": max(0.0, min(1.0, ident)),
        "operon": {
            "contig": next(iter(contigs_involved)) if len(contigs_involved) == 1 else None,
            "start": min(h["start"] for h in subunit_list),
            "stop": max(h["stop"] for h in subunit_list),
            "strand": next(iter(strands_involved)) if len(strands_involved) == 1 else None,
            "intergenic_bp": intergenic_bp,
            "combined_identity": ident,
            "threshold_applied": threshold_applied,
            "subunits": subunit_out,
            "residue_evidence": residue_evidence,
        },
        "method": {"typing_model": "operon", "tools": ["tblastn"]},
    }
    results.append(result)

# Priority order: COMPLETE first, FRAMESHIFT last (§5); stable within a tier.
rank = {s: i for i, s in enumerate(STATUS_PRIORITY)}
results.sort(key=lambda r: rank.get(r["operon_status"], len(STATUS_PRIORITY)))

with open("operon_calls.json", "w") as f:
    json.dump(results, f, indent=2)
PYTHON_SCRIPT
    >>>

    output {
        File result = "operon_calls.json"
    }

    runtime {
        docker: "python:3.12-slim"
        cpu: 1
        memory: "2 GB"
    }
}
