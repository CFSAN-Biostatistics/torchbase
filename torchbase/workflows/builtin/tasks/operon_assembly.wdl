version 1.0

# Pairs subunit HSPs into candidate operons under the syntenic constraints
# declared in [operon] metadata (docs/operon-strategy-plan.md §4): same
# contig, same strand, subunits in `subunit_order`, intergenic gap <=
# `intergenic_max`. Greedy claim so no HSP is double-counted (mirrors
# StxTyper's `reported` flag). Three passes at relaxing stringency:
#   1. strict intergenic_max — full operons
#   2. intergenic_max * intergenic_relax_factor — partial-operon recovery
#   3. unpaired leftover HSPs, reported as single-subunit PARTIAL candidates
#
# This mirrors torchbase/operon.py's pair_operon(); kept inline (rather than
# importing torchbase) because this task runs in an isolated container with
# no torchbase install, matching the existing tasks/*.wdl convention.
task assemble_operons {
    input {
        File hsps
        File operon_config_json
    }

    command <<<
        python3 <<'PYTHON_SCRIPT'
import json


def upstream_gap(a, b):
    if a["strand"] != b["strand"]:
        return None
    if a["strand"] == "+":
        if b["start"] <= a["stop"]:
            return None
        return b["start"] - a["stop"] - 1
    if a["start"] <= b["stop"]:
        return None
    return a["start"] - b["stop"] - 1


def pair_operon(hsps_by_subunit, subunit_order, intergenic_max,
                 require_same_strand=True, require_same_contig=True):
    claimed = set()
    candidates = []

    first_tag = subunit_order[0]
    seeds = sorted(
        hsps_by_subunit.get(first_tag, []),
        key=lambda h: (h["nident"] / h["aln_len"]) if h["aln_len"] else 0.0,
        reverse=True,
    )

    for seed in seeds:
        seed_id = id(seed)
        if seed_id in claimed:
            continue
        chain = [seed]
        ok = True
        for tag in subunit_order[1:]:
            prev = chain[-1]
            best, best_gap = None, None
            for cand in hsps_by_subunit.get(tag, []):
                if id(cand) in claimed:
                    continue
                if require_same_contig and cand["contig"] != prev["contig"]:
                    continue
                if require_same_strand and cand["strand"] != prev["strand"]:
                    continue
                gap = upstream_gap(prev, cand)
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

    return candidates, claimed


with open("~{hsps}") as f:
    all_hsps = json.load(f)

with open("~{operon_config_json}") as f:
    operon_cfg = json.load(f)

subunit_order = operon_cfg["subunit_order"]
intergenic_max = operon_cfg.get("intergenic_max", 36)
relax_factor = operon_cfg.get("intergenic_relax_factor", 2)
require_same_strand = operon_cfg.get("require_same_strand", True)
require_same_contig = operon_cfg.get("require_same_contig", True)

hsps_by_subunit = {}
for hsp in all_hsps:
    hsps_by_subunit.setdefault(hsp["subunit"], []).append(hsp)

results = []

# Pass 1: strict threshold, full operons.
strict_candidates, claimed = pair_operon(
    hsps_by_subunit, subunit_order, intergenic_max,
    require_same_strand, require_same_contig,
)
for chain in strict_candidates:
    results.append({"subunits": chain, "pass": "strict", "complete": True})

# Pass 2: relaxed intergenic distance, partial-operon recovery, over
# whatever HSPs pass 1 left unclaimed.
remaining_by_subunit = {
    tag: [h for h in hsps if id(h) not in claimed]
    for tag, hsps in hsps_by_subunit.items()
}
relaxed_candidates, relaxed_claimed = pair_operon(
    remaining_by_subunit, subunit_order, intergenic_max * relax_factor,
    require_same_strand, require_same_contig,
)
for chain in relaxed_candidates:
    results.append({"subunits": chain, "pass": "relaxed", "complete": True})
claimed |= relaxed_claimed

# Pass 3: leftover unpaired HSPs — single-subunit PARTIAL candidates, one
# per remaining HSP, so a lone subunit hit is still reported rather than
# silently dropped.
for tag in subunit_order:
    for hsp in hsps_by_subunit.get(tag, []):
        if id(hsp) not in claimed:
            results.append({"subunits": [hsp], "pass": "any", "complete": False})
            claimed.add(id(hsp))

with open("candidates.json", "w") as f:
    json.dump(results, f, indent=2)
PYTHON_SCRIPT
    >>>

    output {
        File candidates = "candidates.json"
    }

    runtime {
        docker: "python:3.12-slim"
        cpu: 1
        memory: "2 GB"
    }
}
