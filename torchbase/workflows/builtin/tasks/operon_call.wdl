version 1.0

# Scores each candidate operon from operon_assembly.wdl into the §5 output
# shape: combined identity vs the per-class threshold, residue-table
# resolution for generalized classes, the disruption status ladder, and
# suppression of candidates that are redundant with a better overlapping one.
#
# The algorithm is torchbase/operon.py, localized into the container as
# `operon_module`: one implementation, unit-tested in tests/test_operon.py.
task call_operons {
    input {
        File candidates
        File operon_config_json
        File profiles_table
        File operon_module
        String scheme = ""
    }

    command <<<
        python3 <<'PYTHON_SCRIPT'
import csv
import json
import os
import sys

sys.path.insert(0, os.path.dirname("~{operon_module}"))
from operon import Candidate, format_calls, subtype_prefix

with open("~{candidates}") as f:
    candidates = [Candidate.from_dict(raw) for raw in json.load(f)]

with open("~{operon_config_json}") as f:
    cfg = json.load(f)

with open("~{profiles_table}", newline="") as f:
    profile_rows = list(csv.DictReader(f, delimiter="\t"))

scheme = "~{scheme}" or os.path.basename(os.path.dirname("~{profiles_table}"))

calls = format_calls(candidates, cfg, scheme, subtype_prefix(profile_rows))

with open("operon_calls.json", "w") as f:
    json.dump(calls, f, indent=2)
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
