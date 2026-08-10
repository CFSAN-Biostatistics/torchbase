"""Workflow dispatch: the seam between the package layer and the WDL layer.

Torchbase's package layer dispatches compute to a workflow engine and then
interprets what comes back. This module owns exactly that boundary: build the
engine invocation, run it, and hand back the workflow's declared outputs as
data. It knows nothing about typing models, and nothing above it should know
how miniwdl is invoked.

Keeping the boundary here is what lets the WDL layer stay compute-only. A task
runs a tool in a container and writes a file; deciding what that file *means*
is the package layer's job, in `operon.py`, `allele_calls.py`,
`profile_match.py` and friends.
"""

import json
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Union

Inputs = Dict[str, Union[str, int, float, bool, Path]]


class WorkflowError(RuntimeError):
    """Raised when the workflow engine fails or returns something unreadable."""


def format_input(name: str, value) -> str:
    """One `name=value` argument in the engine's command line.

    Booleans need WDL spelling, and paths must be strings — a `Path` or an open
    file object stringifies to something no engine can resolve.
    """
    if isinstance(value, bool):
        rendered = "true" if value else "false"
    elif isinstance(value, Path):
        rendered = str(value)
    else:
        rendered = str(value)
    return "{}={}".format(name, rendered)


def build_command(
    workflow: Union[str, Path],
    inputs: Inputs,
    engine: str = "miniwdl",
    extra_args: Sequence[str] = (),
) -> List[str]:
    command = [engine, "run", str(workflow)]
    command.extend(extra_args)
    command.extend(format_input(name, value) for name, value in inputs.items())
    return command


def run_workflow(
    workflow: Union[str, Path],
    inputs: Inputs,
    engine: str = "miniwdl",
    extra_args: Sequence[str] = (),
    cwd: Optional[Union[str, Path]] = None,
) -> Dict[str, object]:
    """Run `workflow` and return its outputs as `{output_name: value}`.

    The engine's progress log goes to stderr and is left alone so the user
    still sees it; its result JSON comes back on stdout and is parsed here.
    File outputs come back as paths into the run directory.

    Raises WorkflowError on a non-zero exit or unparseable result.
    """
    command = build_command(workflow, inputs, engine=engine, extra_args=extra_args)
    if shutil.which(engine) is None:
        raise WorkflowError(
            "workflow engine {!r} not found on PATH".format(engine)
        )
    completed = subprocess.run(
        command,
        cwd=str(cwd) if cwd else None,
        stdout=subprocess.PIPE,
        stderr=None,  # stream the engine's log straight through
        universal_newlines=True,
    )
    if completed.returncode != 0:
        raise WorkflowError(
            "workflow execution failed with code {}".format(completed.returncode)
        )
    return parse_outputs(completed.stdout)


def parse_outputs(stdout: str) -> Dict[str, object]:
    """Pull the outputs mapping out of an engine's result JSON.

    miniwdl prints `{"dir": ..., "outputs": {"wf.name": value, ...}}`. Output
    names are returned unqualified — callers ask for `result`, not
    `operon_typing.result`.
    """
    payload = stdout.strip()
    if not payload:
        raise WorkflowError("workflow produced no result JSON")
    start = payload.find("{")
    if start > 0:  # tolerate anything the engine printed before the JSON
        payload = payload[start:]
    try:
        parsed = json.loads(payload)
    except ValueError as error:
        raise WorkflowError("could not parse workflow result: {}".format(error))

    outputs = parsed.get("outputs", parsed)
    if not isinstance(outputs, dict):
        raise WorkflowError("workflow result had no outputs mapping")
    return {name.split(".")[-1]: value for name, value in outputs.items()}


def require_file(outputs: Dict[str, object], name: str) -> Path:
    """The path of a declared `File` output, or a clear error naming what ran."""
    value = outputs.get(name)
    if not value:
        raise WorkflowError(
            "workflow did not produce the expected output {!r} (got: {})".format(
                name, ", ".join(sorted(outputs)) or "nothing"
            )
        )
    path = Path(value if isinstance(value, str) else str(value))
    if not path.exists():
        raise WorkflowError("workflow output {!r} is missing at {}".format(name, path))
    return path


def emit(payload, destination: Optional[Union[str, Path]] = None) -> None:
    """Write a result document to `destination`, or to stdout when unset."""
    rendered = json.dumps(payload, indent=2)
    if destination:
        Path(destination).write_text(rendered + "\n")
    else:
        sys.stdout.write(rendered + "\n")
