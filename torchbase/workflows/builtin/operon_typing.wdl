version 1.0

# Built-in workflow for the "operon" typing model (docs/operon-strategy-plan.md).
# This is the *only* operon implementation for v1 — `--strategy` does not
# apply here (torchbase/cli.py rejects it for operon torches) because the
# typing model and the speed/accuracy axis are orthogonal (§2). Output is a
# list, not a single profile: multiple operons per assembly are normal
# (e.g. an isolate carrying both stx1 and stx2).

import "tasks/protein_search.wdl" as search
import "tasks/operon_assembly.wdl" as assembly
import "tasks/operon_call.wdl" as scoring

workflow operon_typing {
    input {
        File contigs
        File subunit_reference
        File profiles_table
        File operon_config_json
        # torchbase/operon.py — the typing algorithm itself, localized into
        # each task's container so the workflow and the unit tests run the
        # same code (see tasks/*.wdl headers).
        File operon_module
        String scheme = ""
    }

    call search.search_subunits {
        input:
            contigs = contigs,
            subunit_reference = subunit_reference,
            operon_config_json = operon_config_json,
            operon_module = operon_module,
            gapextend = 2
    }

    call assembly.assemble_operons {
        input:
            hsps = search_subunits.hsps,
            operon_config_json = operon_config_json,
            operon_module = operon_module
    }

    call scoring.call_operons {
        input:
            candidates = assemble_operons.candidates,
            operon_config_json = operon_config_json,
            profiles_table = profiles_table,
            operon_module = operon_module,
            scheme = scheme
    }

    output {
        File result = call_operons.result
    }
}
