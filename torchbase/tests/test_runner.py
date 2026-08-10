"""Tests for torchbase.runner — the package/WDL layer boundary.

The runner's contract is narrow and worth pinning precisely: it builds an
engine invocation from data inputs, and it turns the engine's result into
outputs the package layer can interpret. Everything about *what* the outputs
mean lives elsewhere.
"""

import json
from pathlib import Path

import pytest

from torchbase import runner


class TestFormatInput:
    def test_booleans_use_wdl_spelling(self):
        assert runner.format_input("flag", True) == "flag=true"
        assert runner.format_input("flag", False) == "flag=false"

    def test_paths_become_strings(self):
        rendered = runner.format_input("contigs", Path("a") / "b.fasta")
        assert rendered.startswith("contigs=")
        assert "b.fasta" in rendered
        assert "PosixPath" not in rendered and "WindowsPath" not in rendered

    def test_numbers_pass_through(self):
        assert runner.format_input("threads", 4) == "threads=4"
        assert runner.format_input("identity", 0.98) == "identity=0.98"


class TestBuildCommand:
    def test_engine_run_workflow_then_inputs(self):
        command = runner.build_command("wf.wdl", {"contigs": "a.fa", "threads": 2})
        assert command[:3] == ["miniwdl", "run", "wf.wdl"]
        assert command[3:] == ["contigs=a.fa", "threads=2"]

    def test_extra_args_precede_inputs(self):
        command = runner.build_command(
            "wf.wdl", {"contigs": "a.fa"}, extra_args=["--verbose"]
        )
        assert command == ["miniwdl", "run", "wf.wdl", "--verbose", "contigs=a.fa"]


class TestParseOutputs:
    def test_output_names_are_unqualified(self):
        stdout = json.dumps({"dir": "/runs/1", "outputs": {"operon_typing.hits": "/runs/1/hits.tsv"}})
        assert runner.parse_outputs(stdout) == {"hits": "/runs/1/hits.tsv"}

    def test_engine_chatter_before_the_json_is_tolerated(self):
        stdout = 'pulling image...\n{"outputs": {"wf.result": "r.json"}}\n'
        assert runner.parse_outputs(stdout) == {"result": "r.json"}

    def test_empty_result_is_an_error(self):
        with pytest.raises(runner.WorkflowError):
            runner.parse_outputs("   \n")

    def test_unparseable_result_is_an_error(self):
        with pytest.raises(runner.WorkflowError):
            runner.parse_outputs("not json at all")

    def test_result_without_outputs_mapping_is_an_error(self):
        with pytest.raises(runner.WorkflowError):
            runner.parse_outputs(json.dumps({"outputs": ["a", "b"]}))


class TestRequireFile:
    def test_returns_the_path_when_present(self, tmp_path):
        hits = tmp_path / "hits.tsv"
        hits.write_text("")
        assert runner.require_file({"hits": str(hits)}, "hits") == hits

    def test_names_the_missing_output(self):
        with pytest.raises(runner.WorkflowError) as raised:
            runner.require_file({"other": "x"}, "hits")
        assert "hits" in str(raised.value)
        assert "other" in str(raised.value)

    def test_declared_but_absent_file_is_an_error(self, tmp_path):
        with pytest.raises(runner.WorkflowError):
            runner.require_file({"hits": str(tmp_path / "gone.tsv")}, "hits")


class TestEmit:
    def test_writes_json_to_a_destination(self, tmp_path):
        destination = tmp_path / "calls.json"
        runner.emit([{"profile_id": "stx2c"}], destination)
        assert json.loads(destination.read_text()) == [{"profile_id": "stx2c"}]

    def test_writes_json_to_stdout_without_a_destination(self, capsys):
        runner.emit({"profile_id": "stx1a"})
        assert json.loads(capsys.readouterr().out) == {"profile_id": "stx1a"}


class TestRunWorkflow:
    def test_missing_engine_is_reported_clearly(self):
        with pytest.raises(runner.WorkflowError) as raised:
            runner.run_workflow("wf.wdl", {}, engine="no-such-engine-xyz")
        assert "no-such-engine-xyz" in str(raised.value)
