# Torchbase Codebase Audit

**Date**: 2026-06-06
**Method**: Code inspection + test execution. Spec document: `docs/PRD-MLST-Implementation.md`.
All conclusions derived from source only; documentation, comments, and commit messages ignored.
**Test run**: `pytest -v` (after `pip install toml`): 697 passed, 12 failed, 13 skipped, 82 deselected in 12.42 seconds.

---

## Executive Summary

The core data structures (`Profile`, `Schema`, `Version`, `Torch`, multi-scheme loader, registry manager, versioning) are implemented and have genuine test coverage. The quality analysis modules (`kmer_analysis.py`, `report.py`) and the BIGSdb client are implemented, though the k-mer module uses pure-Python Jaccard similarity instead of Jellyfish (contradicting the PRD spec). The most critical gap is that **no actual sourmash, jellyfish, or minimap2 calls occur anywhere in the codebase**: the MinHash WDL tasks are explicit mock stubs, the allele-fetching step in the PubMLST converter writes synthetic placeholder FASTA, and IPFS/IPNS integration returns a hardcoded `/tmp/ipfs/<CID>` path rather than fetching from a real node. Tests labeled "RED-phase—they MUST fail" are in fact passing (220 tests passing across 7 such files), because those tests inspect WDL file content as text rather than executing the workflows. The system cannot produce a real typing result end-to-end.

---

## 1. What Is Genuinely Implemented

### 1.1 Core Data Model (`torchbase/torchbase.py`)

`Profile`, `Schema`, and their parsing logic are genuinely implemented. The `Profile.__eq__` wildcard and PubMLST-format matching logic works. `Version` supports all three strategies (snapshot, semver, content-hash) with cross-strategy timestamp comparison. Tests in `test_torchbase.py` exercise real behavior with non-trivial data.

### 1.2 Multi-Scheme Torch Loader (`torchbase/torchfs.py`)

`Torch.load()` handles both single- and multi-scheme formats, validates namespace/name/version against path, discovers `schemes/*/` directories, loads profiles TSV, and scans allele FASTA files. `concatenate_alleles()` and `transform_profiles()` apply scheme-prefixed naming. The 19 passing tests in `test_multi_scheme_torch.py` and 21 in `test_multi_scheme_concatenation.py` exercise real file-based behavior.

### 1.3 K-mer Quality Analysis (`torchbase/quality/kmer_analysis.py`)

`analyze_locus()` is genuinely implemented using pure-Python k-mer extraction (not Jellyfish). Gap detection (`>2%` jump in sorted similarities) falls back to 99th percentile threshold. Overlap detection (one-way containment ≥95%) and duplicate detection (symmetric Jaccard ≥98%) are functional. 24 tests in `test_kmer_analysis.py` use synthetic FASTA data with known properties and make non-trivial assertions.

### 1.4 Quality Report Generator (`torchbase/quality/report.py`)

`generate_report()` produces all three output formats (text, JSON, both). ASCII histograms are rendered. Hierarchical suspect propagation collapses suspect loci to suspect profiles (though this mapping is simplified: `suspect_profiles = sorted(suspect_loci)`). 31 tests in `test_quality_report.py` are genuine.

### 1.5 BIGSdb REST Client (`torchbase/conversions/bigsdb_client.py`)

The client is implemented with typed response dataclasses, pagination, temporal filtering, and error handling. All tests in `test_bigsdb_client.py` (25 tests) mock HTTP responses and test real parsing logic. This module is complete relative to the PRD's API contract section.

### 1.6 Registry Manager (`torchbase/registry.py`)

`RegistryManager.resolve()`, `fetch_torch()`, and `pin_torch()` are implemented. Registry fallback ordering, pin precedence, workflow dependency pinning, and atomic config writes work. Tests in `test_registry.py` (25 tests) mock network calls and verify real resolution logic. **Critical caveat**: `_cid_to_local_path()` at line 169 is explicitly a mock: it returns `Path("/tmp/ipfs") / cid` with the comment "Currently returns a mock path. Real IPFS integration would be added here."

### 1.7 Registry Config (`torchbase/config.py`)

`RegistryConfig.load()` correctly implements hierarchical configuration: project `.torchbase.toml` overrides user `~/.torchbase/config.toml`. 9 tests in `test_registry_config.py` are genuine.

### 1.8 Versioning Support (`torchbase/torchbase.py`)

`Version.parse()`, `Version.compare()`, and all comparison operators are implemented for all three strategies. Genuine tests in `test_torchbase.py` cover parsing and comparison.

### 1.9 Workflow Discovery (`torchbase/torchfs.py` + `torchbase/cli.py`)

Convention-based discovery (`main.wdl` takes precedence over manifest `workflow` field) is implemented in `_load_single_scheme()` at lines 276–281. The `--strategy` vs. embedded-workflow conflict check in `_run()` (lines 544–548) is implemented. 31 tests in `test_workflow_discovery.py` are genuine (they build temp torch directories and check `Torch.workflow`).

### 1.10 CLI Strategy Routing

The `--strategy` flag with choices `fast/balanced/sensitive/auto` is implemented in `cli.py`. Strategy-to-workflow mapping at lines 594–598 resolves built-in workflow files. The auto-strategy analysis (`_analyze_sequences()`) detects FASTA/FASTQ format and uses mean sequence length to select between fast (contigs >1000 bp) and balanced (reads <500 bp). The `--quality-json`, `--include-suspect-alleles`, `--exclude-suspect-alleles`, `--exclude-suspect-loci`, `--exclude-suspect-profiles` flags exist in the `run` command definition (lines 516–520).

### 1.11 WDL Task Files (Structural)

Four shared WDL tasks exist and pass `miniwdl check` syntax validation:
- `tasks/minhash.wdl` — sketch, compare, call_alleles, call_alleles_minhash
- `tasks/alignment.wdl` — align_and_call, align_sequences
- `tasks/profile_lookup.wdl` — lookup_profile
- `tasks/filter_alleles.wdl` — filter_alleles

Three built-in strategy workflows exist: `fast_typing.wdl`, `balanced_typing.wdl`, `sensitive_typing.wdl`. The MLST workflow torch at `workflows/mlst/1.0.0.torch/main.wdl` exists.

---

## 2. Critical Gaps

### 2.1 MinHash (sourmash) Is a Mock Stub

**Claim** (PRD §Implementation Decisions): "MinHash sketching (sourmash) against all scheme alleles."

**Reality**: `torchbase/workflows/builtin/tasks/minhash.wdl`, lines 15–27:
```
# Create a mock sketch file (for testing without sourmash)
# In a real environment, this would call sourmash
sketch_file = "sequences.sig"
with open(sketch_file, 'w') as f:
    f.write("")  # Create an empty file to represent the sketch
```
The `compare_sketches` task at lines 52–99 creates a synthetic identity matrix that simulates similarity, not a real sourmash comparison. The `call_alleles_minhash` task operates on this synthetic matrix. There is no `sourmash` binary call anywhere in the WDL tasks.

**Consequence**: Every test asserting "WDL uses sourmash" (e.g., `test_wdl_command_uses_sourmash`) passes by checking that the string `"sourmash"` appears in a code comment inside the WDL file, not by executing sourmash. No real MinHash comparison can occur at runtime.

### 2.2 PubMLST Converter Does Not Fetch Real Allele Sequences

**Claim** (PRD §Implementation Decisions): "Orchestrates conversion workflow: API fetch → directory structure → k-mer analysis → metadata generation."

**Reality**: `torchbase/conversions/pubmlst.py`, lines 86–93:
```python
# For now, we'll create stub FASTA files
# In a real implementation, these would be fetched from the API
fasta_path = alleles_dir / f"{locus.locus_id}.fasta"
# Create minimal FASTA content (stub for now)
_write_stub_fasta(fasta_path, locus.locus_id)
```
`_write_stub_fasta()` (lines 169–181) writes two identical repetitive sequences for every locus regardless of what the API returns. The BIGSdb client does not implement a `GET /db/{database}/loci/{locus}/alleles_fasta` endpoint at all. `alleles_fasta` endpoint is listed in the PRD (§API Contracts) but absent from `bigsdb_client.py`.

**Consequence**: Any torch created via `torchtools convert-pubmlst` contains synthetic placeholder alleles, making the resulting database useless for real typing.

### 2.3 IPFS Integration Is Hardcoded Mock

**Claim** (PRD §Solution): "IPNS-based torch registry with version pinning."

**Reality**: `torchbase/registry.py`, lines 151–170:
```python
def _cid_to_local_path(self, cid: str) -> Path:
    """
    Note:
        Currently returns a mock path. Real IPFS integration would
        be added here.
    """
    ipfs_cache_dir = Path("/tmp/ipfs")
    return ipfs_cache_dir / cid
```
`torchbase/torchfs.py`, line 20: `TORCHBASE_REGISTRY_HASH = ""`. The `download_torch()` function at line 46–47 is `pass`. `exists()` at lines 63–65 returns `False  # TODO`.

**Consequence**: `torchbase pull` cannot actually retrieve a torch from IPFS. `manager.fetch_torch()` returns a nonexistent path. The `torchbase list` command filters on `exists()` which always returns False, so it will always show nothing.

### 2.4 Jellyfish Is Not Used

**Claim** (PRD §Implementation Decisions): "Encapsulates all Jellyfish operations for pairwise k-mer similarity within loci."

**Reality**: `torchbase/quality/kmer_analysis.py` has zero Jellyfish references. It uses `_get_kmers()` (pure Python set comprehension, lines 64–69) and `_jaccard_similarity()` (pure Python set intersection, lines 73–82). The `pyproject.toml` dependency `"jellyfish~=1.0.0"` refers to the Python string-similarity library, not the k-mer counting tool.

**Consequence**: The k-mer analysis is O(n²·L·k) in pure Python. On real MLST databases with 100–3000 alleles per locus, this will be prohibitively slow.

### 2.5 The `torchbase` CLI Binary Is Not Installed

**Claim** (PRD §Implementation): "CLI integration via `torchtools convert pubmlst <database_url> <scheme_id>`."

**Reality**: `which torchbase` returns nothing. The package has not been installed in editable mode (`pip install -e .` not run). Tests that invoke `subprocess.run(["torchbase", "run", "--help"], ...)` all fail with `FileNotFoundError: [Errno 2] No such file or directory: 'torchbase'`. This is the root cause of 10 of the 12 failing tests.

### 2.6 Auto Strategy Never Selects `sensitive`

**Claim** (PRD §Implementation): "auto: Automatically selects based on input characteristics." The `_analyze_sequences()` docstring lists `selected_strategy: 'fast', 'balanced', or 'sensitive'`.

**Reality**: `cli.py` lines 455–466:
```python
if mean_length > 1000:
    selected_strategy = 'fast'
elif mean_length < 500:
    selected_strategy = 'balanced'
else:
    selected_strategy = 'balanced'
```
`sensitive` is never selected by `auto`. Long reads (which would benefit from sensitive) are not detected separately from short reads.

### 2.7 Output JSON Missing `scheme` Field in Built-in Workflows

**Claim** (CLAUDE.md / PRD output format): The standardized output JSON must include `"scheme": "salmonella_mlst"`.

**Reality**: `torchbase/workflows/builtin/tasks/profile_lookup.wdl`, lines 112–128 — the `result` dict has no `scheme` key. The MLST workflow torch at `workflows/mlst/1.0.0.torch/main.wdl` does set `result["scheme"]` (line 532) via inferred scheme logic, but the three built-in strategy workflows (`fast_typing.wdl`, `balanced_typing.wdl`, `sensitive_typing.wdl`) all route through `profile_lookup.wdl` which does not emit `scheme`.

### 2.8 `nearest_st` Not Implemented in Built-in Workflows

**Claim** (PRD §User Story 13, §Implementation): "Novel profile/allele reporting with nearest ST"; output includes `"nearest_st"`.

**Reality**: `torchbase/workflows/builtin/tasks/profile_lookup.wdl` produces no `nearest_st` field. Only the MLST workflow torch (`workflows/mlst/1.0.0.torch/main.wdl`, lines 471–556) implements `calculate_nearest_st()`. The three built-in strategy workflows do not call this logic.

### 2.9 `torchbase info` and `torchtools version/build` Are Stubs

`cli.py` lines 124–127:
```python
@cli.command("info")
@torch
def _info(torch):
    "Display info for the selected torch."
    pass
```
`torchtools version` (line 719) and `torchtools build` (line 726) are identical `pass` stubs.

### 2.10 `torchbase/workflows/mlst/` Directory Should Have Been Removed

Test `test_synthetic_example_torches.py::TestOldMLSTTorchRemoved::test_old_mlst_workflow_directory_removed` asserts that `torchbase/workflows/mlst/` should not exist. It exists and causes this test to fail. The test comment says it was supposed to be replaced by examples in `examples/`.

---

## 3. Missing Features

| Feature | PRD Section | Status | Notes |
|---------|-------------|--------|-------|
| Real sourmash calls in WDL tasks | §Implementation / User Story 15 | Not implemented | WDL tasks write empty .sig files and synthetic similarity matrices |
| Allele FASTA fetch from BIGSdb API | §API Contracts (`alleles_fasta` endpoint) | Not implemented | `pubmlst.py` writes stub FASTAs; `bigsdb_client.py` has no `alleles_fasta` method |
| Real IPFS download (`download_torch`) | §Solution, §Registry Manager | Stub (`pass`) | `torchfs.py:47` |
| IPNS manifest fetch (real network) | §Registry Manager | Stub (HTTP GET that returns TOML; no actual IPNS resolution) | `registry.py:_fetch_manifest` makes a raw HTTP GET, not IPNS |
| `torchbase info` command | §Workflow Discovery | Stub | `cli.py:127` |
| `torchtools version` command | §Authoring | Stub | `cli.py:720` |
| `torchtools build` command | §Authoring | Stub | `cli.py:726` |
| `torchtools convert pubcgmlst` | §Conversion | Stub (`pass`) | `cli.py:767` |
| `torchtools convert chewie-ns` | §Conversion | Stub (`pass`) | `cli.py:774` |
| `torchtools convert shigatyper` | §Conversion | Stub (`pass`); test class empty | `cli.py:781`, `test_shigatyper_conversion.py` |
| `torchtools convert pubmlst` (legacy) | §Conversion | Stub (`pass`) | `cli.py:762` |
| Jellyfish-based k-mer counting | §K-mer Analysis Module | Not implemented | Pure Python used instead |
| Auto strategy selecting `sensitive` | §Strategy Selection | Not implemented | `_analyze_sequences` never returns `sensitive` |
| `scheme` field in built-in workflow output | §Output Format | Missing in fast/balanced/sensitive workflows | Only present in MLST torch `main.wdl` |
| `nearest_st` in built-in workflow output | §User Story 13, §Output Format | Missing in built-in workflows | Only in MLST `main.wdl` |
| `torchbase list` showing real installed torches | §User Story 2 | Broken (always empty) | `exists()` always returns False |
| Metadata `[version]` section with `strategy`/`timestamp` | §Versioning Support / §Data Structures | Not generated by `pubmlst.py` | Converter generates no `[version]` section |
| Metadata `[provenance]` `last_updated` from BIGSdb | §Data Structures | Partial | `last_updated` from scheme metadata captured; locus `last_updated` not written to metadata |
| Quality JSON `allele_count`/`kmer_size`/`similarity_stats` per PRD schema | §Data Structures | Schema mismatch | `report.py` produces different keys (`similarities`, `statistics`) vs PRD spec (`allele_count`, `similarity_stats`) |
| Depth filtering for reads in built-in workflows | §Default MLST Workflow Stage 2 | Not in fast/balanced/sensitive | Only in MLST torch `main.wdl` |
| Ambiguity triggers (top-2 within 3%, best <92%, >2 schemes, coverage <80%) | §Implementation Decisions | Not in built-in workflows | Only threshold `confidence_threshold` used; specific triggers not enumerated |

---

## 4. Test Suite Integrity Summary

**Totals**: 804 collected; 82 deselected (marked `not miniwdl` in `pyproject.toml`); 722 run; 697 passed; 12 failed; 13 skipped.

**Wall time**: 12.42 seconds for 697 passing tests. This is consistent with tests that read files and do string matching rather than executing biological tools.

### What Passing Tests Actually Prove

The 697 passing tests prove:
1. **Data model correctness**: `Profile`, `Schema`, `Version` parsing and equality logic works on synthetic data.
2. **File structure**: WDL files exist at expected paths and contain expected string patterns (task names, input/output variable names, keywords like `"sourmash"` in comments).
3. **WDL syntax**: `miniwdl check` passes for all WDL files (structure, not semantics).
4. **Registry logic**: Version resolution, pin precedence, multi-registry fallback, and atomic writes work when network calls are mocked.
5. **K-mer analysis**: Pure-Python Jaccard similarity and threshold detection work on synthetic FASTA files.
6. **Quality reports**: All three output formats render correctly.
7. **Torch loading**: `Torch.load()` correctly handles single- and multi-scheme formats from temp directories.
8. **CLI routing**: Strategy flag and workflow-conflict detection work at the Python level (not end-to-end with subprocess).

The passing tests **do not** prove:
- That sourmash can be called and produces correct MinHash results.
- That minimap2 alignment produces correct allele calls.
- That Jellyfish is used for k-mer counting.
- That IPFS fetches or uploads work.
- That the end-to-end pipeline (`torchbase run`) produces a correct typing result for any input.
- That converted torches contain real allele sequences.

### Deselected Tests (82 tests, `not miniwdl` marker)

The 82 deselected tests are marked `@pytest.mark.miniwdl` and require Docker and miniwdl execution. These are the only tests that would exercise the WDL tasks end-to-end. They are excluded from CI by default (`addopts = "-m 'not miniwdl'"` in `pyproject.toml`). Their presence indicates the developers know the WDL tasks need real execution testing; their exclusion means this testing never runs automatically.

---

## Appendix A: Bad Tests

### A.1 Unconditionally Ignored Tests

**`test_workflow_torch.py::TestDefaultWorkflowTorchLocation`** (13 tests, all skipped):
```python
@pytest.mark.skip(reason="Issue #61: Default workflow torch moved to examples/ and replaced with synthetic examples")
class TestDefaultWorkflowTorchLocation:
```
These tests check that `torchbase/workflows/mlst/1.0.0.torch/` exists and is loadable. The skip was applied when the torch was supposed to move, but the directory was never actually removed (confirmed by `test_old_mlst_workflow_directory_removed` failing). These tests cover real behavior but are silenced.

### A.2 Vacuous Tests

**`test_pubmlst_conversion.py::TestTorchFSIntegration::test_resulting_torch_is_loadable`** (line 519):
```python
# This is a placeholder - actual test would use Torch.load()
assert True
```
The test creates a temp directory, writes a metadata file with missing required fields, then asserts `True` unconditionally. It does not call `Torch.load()`.

**`test_pubmlst_conversion.py::TestTorchFSIntegration::test_torch_load_validates_metadata`** (line 549):
```python
# Torch.load() should reject invalid metadata
# This is a placeholder test
assert True
```
Same pattern: the test sets up invalid metadata but then asserts unconditionally without calling the function being tested.

**`test_pubmlst_conversion.py::TestTorchFSIntegration::test_torch_load_scans_resources_directory`** (line 565):
```python
# Torch should discover these files
assert True
```
Same pattern.

**`test_suspect_data_workflow_flags.py`** — multiple test methods that are `pass` (lines 295, 316, 591+):
```python
def test_exclude_suspect_loci_implies_exclude_suspect_alleles(self):
    """Excluding suspect loci implicitly excludes all their alleles"""
    pass
```
These test functions have docstrings asserting behavior but no assertions. They pass vacuously.

**`test_pubmlst_conversion.py::TestKmerAnalysisIntegration::test_respects_kmer_size_parameter`** (line 499):
```python
def test_respects_kmer_size_parameter(self):
    assert True  # Parameter acceptance verified by CLI test
```
The assertion is unconditional; no actual verification occurs.

### A.3 Tests Named for Feature X That Test Feature Y

**All 220 tests in the 8 "RED-phase" files** (`test_alignment_fallback.py`, `test_fast_typing_workflow.py`, `test_balanced_typing_workflow.py`, `test_sensitive_typing_workflow.py`, `test_minhash_allele_calling.py`, `test_mlst_workflow_orchestration.py`, `test_shared_wdl_tasks.py`, `test_suspect_data_workflow_flags.py`):

The test class names and test names assert behavioral properties of WDL workflow execution (e.g., `test_wdl_command_uses_sourmash`, `test_wdl_detects_low_similarity_threshold`, `test_wdl_has_depth_filtering_task`). However, the test bodies open WDL files as text and search for string patterns. For example, `test_wdl_command_uses_sourmash` in `test_minhash_allele_calling.py` (line 297):

```python
def test_wdl_command_uses_sourmash(self):
    wdl_path = TORCHBASE_ROOT / "workflows" / "minhash_allele_calling.wdl"
    with open(wdl_path) as f:
        content = f.read()
    assert "sourmash" in content, "Task command does not use sourmash"
```

The word "sourmash" appears in a comment ("# Create a mock sketch file (for testing without sourmash)"). The test passes because the comment exists, not because sourmash is invoked. These tests are named to suggest the WDL pipeline correctly uses sourmash, minimap2, or other tools, but they only assert that certain strings appear somewhere in the WDL file text.

The 8 files' docstrings all say "These are RED-phase tests - they MUST fail because the feature is not yet complete." They do not fail; they pass, obscuring the implementation gap.

### A.4 Stale "Currently Fails" Comments

**All 8 RED-phase test files**: The module docstring in each file states the tests "MUST fail because the feature is not yet complete." All 220 tests in these files pass. The comment is stale. The tests do not verify execution behavior; they verify WDL file text content, so they pass even though the underlying functionality is not complete.

---

## Appendix B: Non-Implemented Features

### B.1 Core Execution Pipeline

| Feature | Spec Location | Status | Notes |
|---------|---------------|--------|-------|
| Sourmash MinHash sketching (real) | §Implementation, PRD line 112 | Not implemented | WDL stub creates empty `.sig` file |
| Sourmash sketch comparison (real) | §Implementation, PRD line 112 | Not implemented | WDL stub generates synthetic identity matrix |
| Minimap2 alignment execution via WDL | §Alignment Fallback | Code present | Will call minimap2 if available but never exercised by CI tests |
| Depth filtering for reads (WDL task) | §Default MLST Workflow Stage 2 | Present in MLST torch, absent in builtin fast/balanced/sensitive | |
| Scheme auto-detection from MinHash across 100+ schemes | §User Story 6, 15 | Not implemented | Only single profile table lookup in profile_lookup.wdl |

### B.2 Distribution Infrastructure

| Feature | Spec Location | Status | Notes |
|---------|---------------|--------|-------|
| Real IPFS download | §Solution | Stub (`pass` at `torchfs.py:47`) | |
| IPNS manifest resolution | §Registry Manager | HTTP GET only; no IPNS protocol | `/ipns/` URL is not resolved to IPFS CID |
| `torchbase pull` end-to-end | §User Story 2 | Broken | `fetch_torch` returns `/tmp/ipfs/<cid>` which does not exist |
| `torchbase list` (installed torches) | §User Story | Broken | `exists()` always returns False |

### B.3 Conversion Pipeline

| Feature | Spec Location | Status | Notes |
|---------|---------------|--------|-------|
| Allele FASTA download from BIGSdb | §API Contracts | Not implemented | Stub FASTAs written; `alleles_fasta` API endpoint not implemented in client |
| `[version]` section in generated metadata | §Data Structures | Not generated | `pubmlst.py:_generate_metadata()` does not write `[version]` |
| `torchtools convert pubcgmlst` | §Conversion | Stub | `cli.py:767` |
| `torchtools convert chewie-ns` | §Conversion | Stub | `cli.py:774` |
| `torchtools convert shigatyper` (CLI) | §Conversion | Stub | `cli.py:781` |

### B.4 Output Format

| Feature | Spec Location | Status | Notes |
|---------|---------------|--------|-------|
| `scheme` field in built-in workflow results | §Output Format | Missing | Only in MLST torch main.wdl |
| `nearest_st` in novel profile results | §User Story 13 | Missing in built-in workflows | Only in MLST torch main.wdl |
| Quality JSON schema matching PRD spec | §Data Structures | Schema mismatch | PRD: `allele_count`, `similarity_stats.min/median/percentile_99`, `threshold_method`. Actual: `similarities`, `statistics.mean/std_dev/min/max/percentile_99`, `threshold_type` |

### B.5 CLI Commands

| Feature | Spec Location | Status | Notes |
|---------|---------------|--------|-------|
| `torchbase info` | §CLI Layer | Stub `pass` | |
| `torchtools version` | §Authoring | Stub `pass` | |
| `torchtools build` | §Authoring | Stub `pass` | |
| `torchbase` binary installed | §Entry Points | Not installed | `which torchbase` fails |
| Auto strategy selecting `sensitive` | §Strategy Selection | Not implemented | Dead code path |
