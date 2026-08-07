version 1.0

# Pairs subunit HSPs into candidate operons under the syntenic constraints
# declared in [operon] metadata (docs/operon-strategy-plan.md §4): same
# contig, same strand, subunits in `subunit_order`, intergenic gap <=
# `intergenic_max`. Co-linear HSPs split by a frame change are stitched first,
# then a reference set's many accessions per subtype are reduced to one HSP
# per locus per class, then four pairing passes at relaxing stringency claim
# HSPs greedily (mirroring StxTyper's `reported` flag).
#
# The algorithm is torchbase/operon.py, localized into the container as
# `operon_module`: one implementation, unit-tested in tests/test_operon.py.
task assemble_operons {
    input {
        File hsps
        File operon_config_json
        File operon_module
    }

    command <<<
        python3 <<'PYTHON_SCRIPT'
import json
import os
import sys

sys.path.insert(0, os.path.dirname("~{operon_module}"))
from operon import HSP, assemble_candidates

with open("~{hsps}") as f:
    hsps = [HSP.from_dict(raw) for raw in json.load(f)]

with open("~{operon_config_json}") as f:
    cfg = json.load(f)

candidates = assemble_candidates(hsps, cfg)

with open("candidates.json", "w") as f:
    json.dump([c.to_dict() for c in candidates], f, indent=2)
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
