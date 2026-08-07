# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Torchbase is a Python framework for generalized allelic typing from sequencing data. It works for any allelic profile-based typing system: MLST, serotyping, antimicrobial resistance prediction, and more. It addresses challenges around database distribution, versioning, reproducibility, and maintainability.

**Core Concept**: "Torches" are versioned, distributed databases containing:
- Allele reference sequences
- Allelic profile tables (schema definitions)
- Optional: WDL workflows for custom typing logic
- Build files and metadata

Torches are distributed via IPFS to enable versioned, reproducible typing across different users and institutions.

**Typing Model**: What kind of scheme this is — a torch property (`typing_model` in `metadata.toml`):
- `allelic` (default, backward compatible): independent loci matched via `Profile`.
- `operon`: multi-subunit operons (e.g. Shiga toxin *stxAB*) typed by protein-space
  synteny-aware search — combined identity, adjacency, disruption status. See
  docs/operon-strategy-plan.md.

**Typing Strategy**: Orthogonal speed/accuracy tradeoff, an invocation flag
(`--strategy`), applicable only to the `allelic` typing model — `operon`
torches have exactly one built-in workflow and reject `--strategy`:
- `fast`: MinHash-based calling only (fastest)
- `balanced`: MinHash with alignment fallback (default)
- `sensitive`: Full alignment-based calling (most accurate)
- `auto`: Automatically selects based on input characteristics

## Architecture

### Three-Layer System

1. **Torch Definition Layer** (`torchbase/torchbase.py`, `torchbase/operon.py`)
   - `Schema`: Container for typing profiles with version info
   - `Profile`: Represents allelic profiles with special handling for wildcards (`IGNORE = "?"`) and exclusions (`EXCLUDE = "X"`)
   - Profile equality supports multiple formats: tuples, dicts, PubMLST-style strings (e.g., "locus_allele")
   - `operon.py`: the algorithm for the `operon` typing model — BLAST
     normalization, frameshift stitching, locus reduction, synteny pairing,
     threshold and residue-table resolution, status-ladder scoring. Operon
     torches bypass `Profile` entirely (no adjacency/threshold concept there).
     `operon_typing.wdl`'s tasks take this module as a `File` input and import
     it inside their containers, so there is one implementation, not a copy per
     task — do not reintroduce inline copies.

2. **Filesystem/Distribution Layer** (`torchbase/torchfs.py`)
   - `Torch` dataclass: Loads and validates torch packages from disk
   - `Torch.typing_model` ("allelic" default, or "operon") and
     `Torch.operon_config`/`Torch.operon_profiles` for operon torches,
     validated at load time via `torchbase.operon.validate_operon_metadata`
   - IPFS integration for distributed torch retrieval (via `ipyfs`)
   - Manifest system for tracking available torches
   - Environment-based IPFS configuration (`TORCHBASE_IPFS_NODE`, `TORCHBASE_IPFS_PORT`)

3. **CLI Layer** (`torchbase/cli.py`)
   - Two command groups:
     - `torchbase`: User-facing commands (list, pull, info, run, workflow)
     - `torchtools`: Authoring commands (build, version, convert)
   - Dispatch on the pair `(typing_model, strategy)`:
     - `typing_model = "allelic"`: built-in workflows in `torchbase/workflows/builtin/` (fast/balanced/sensitive), or a custom `main.wdl`
     - `typing_model = "operon"`: routes to `torchbase/workflows/builtin/operon_typing.wdl`; `--strategy` is rejected (one implementation for v1)
     - Custom workflows via torch-embedded `main.wdl` (error if `--strategy` also given)
   - Execution via `miniwdl` for WDL workflows
   - Automatic file decompression/compression to zstandard format


### Torch Package Structure

**Single-Scheme Format** (simple, most common):
```
<namespace>/<torchname>/<version>.torch/
├── metadata.toml           # Package metadata, citations, maintainers
├── profiles.tsv            # Tab-separated allelic profile table
├── main.wdl                # Optional: custom workflow (overrides built-in)
└── _resources/             # Reference FASTA files for alleles
    ├── locus1.fasta
    └── locus2.fasta
```

**Multi-Scheme Format** (advanced, multiple organisms):
```
<namespace>/<torchname>/<version>.torch/
├── metadata.toml
└── schemes/
    ├── organism1/
    │   ├── profiles.tsv
    │   └── alleles/
    │       ├── locus1.fasta
    │       └── locus2.fasta
    └── organism2/
        ├── profiles.tsv
        └── alleles/
            ├── locus1.fasta
            └── locus2.fasta
```

**Operon Format** (`typing_model = "operon"`, e.g. StxTyper-parity stx torches):
```
<namespace>/<torchname>/<version>.torch/
├── metadata.toml           # typing_model = "operon" + [operon] config block
├── profiles.tsv            # declarative manifest of valid subtypes (not the matcher)
└── _resources/
    └── subunits.faa        # protein reference set, accession|subunit_role|reference_subtype|class headers
```
No `main.wdl` — routes to the built-in `operon_typing.wdl`. Worked examples:
`torchtools convert stxtyper` (stx) and `examples/etec_lt` (ETEC LT, hand-encoded).

`[operon]` keys, all scheme-agnostic (see `torchbase/operon.py`'s docstring):
`subunit_order`, `intergenic_max`/`intergenic_min` (negative where subunit
genes overlap, as ETEC's *eltA*/*eltB* do), `intergenic_relax_factor`,
`require_same_strand`/`require_same_contig`, `identity_thresholds` (per class,
plus `default`), `generalized_classes` + `residue_rules` (for classes identity
alone cannot separate — residue indices are 0-based offsets into the reference
protein), `superclass_pattern` (parent class, reported when an operon is not
COMPLETE), `min_operon_identity`, `overlap_slack`.

**Workflow Discovery**:
- If torch has `main.wdl` → use it (user cannot specify `--strategy`)
- If `typing_model == "operon"` → use built-in `operon_typing.wdl` (user cannot specify `--strategy`)
- Else → use built-in workflow with selected strategy
- CLI concatenates multi-scheme torches with scheme-prefixed locus names (e.g., `salmonella_adk_1`)

### Conversion System (`torchbase/conversions/`)

Converts external typing schemes to torch format:
- `pubmlst.py`: PubMLST MLST schemes
- `pubcgmlst.py`: PubMLST cgMLST schemes
- `shigatyper.py`: ShigaTyper database
- `stxtyper.py`: NCBI StxTyper's stx.prot → operon torch (typing_model = "operon")
- `chewie-ns`: Chewie-NS wgMLST (planned)

All conversions use cookiecutter templates in `torchbase/templates/`.

## Development Commands

### Setup
```bash
make install-dev          # Install in editable mode with dev dependencies
# or
pip install -e '.[dev]'
```

### Testing
```bash
make test                 # Run pytest
make test-all             # Run tests on all Python versions via tox
make coverage             # Generate coverage report and open in browser
pytest                    # Direct pytest invocation
```

### Linting
```bash
make lint                 # Run flake8
```

### Building/Distribution
```bash
make dist                 # Build source and wheel distributions
make release              # Upload to PyPI (requires twine)
```

### Cleanup
```bash
make clean                # Remove all build/test/Python artifacts
make clean-pyc            # Remove Python file artifacts only
make clean-test           # Remove test/coverage artifacts only
```

## Key Implementation Details

### Profile Comparison Logic

The `Profile.__eq__` method handles flexible matching:
- Compares against tuples, lists, dicts, or other Profile objects
- Special values:
  - `Special.IGNORE ("?")`: wildcard, matches any value at that locus
  - `Special.EXCLUDE ("X")`: locus should not be present in query
- Supports PubMLST naming convention where alleles are prefixed with locus name (e.g., "dinB_1")
- Length validation ensures excluded loci aren't counted

### Torch Loading and Validation

`Torch.load()` performs sanity checks:
- Validates namespace/name/version consistency between metadata.toml and directory path
- `typing_model` (default `"allelic"`): allelic torches parse the profile table via
  `Profile.parse()`; operon torches (`typing_model = "operon"`) parse the `[operon]`
  block and profiles.tsv as a raw manifest instead — see `torchbase.operon.validate_operon_metadata`
- Scans `_resources/` for reference files (ignores dotfiles)
- Returns fully-loaded Torch dataclass

### CLI File Handling

`ReadsFile` custom Click parameter type:
- Auto-detects compression (gzip, bzip2, zip, zstd) via magic bytes
- Transparently converts all input to zstandard format
- Used for `-c/--contigs`, `-r/--reads`, `-pe1/-pe2/--paired`, `-i/--interlaced`, `-l/--longreads`

## Testing Notes

- Test fixtures in `torchbase/tests/bigsdb_fixture.py` provide realistic BigsDB/MLST schema data
- Profile parsing and equality tests cover edge cases: wildcards, exclusions, PubMLST formats
- IPFS tests currently stubbed (hash = `/ipfs/QmQPeNsJPyVWPFDVHb77w8G42Fvo15z4bG2X8D2GhfbSXc/readme`)

## Workflow System

### Built-in Strategies

Located in `torchbase/workflows/builtin/`:
- `fast_typing.wdl`: MinHash-only pipeline, fastest
- `balanced_typing.wdl`: MinHash + conditional alignment (default)
- `sensitive_typing.wdl`: Always runs alignment, most accurate

All three import shared tasks from `torchbase/workflows/builtin/tasks/`:
- `minhash.wdl`: Sourmash sketching and comparison
- `alignment.wdl`: Minimap2 with preset selection (asm20/asm5/asm5+eqx)
- `profile_lookup.wdl`: Profile matching and scheme inference

`operon_typing.wdl` is the single built-in workflow for `typing_model = "operon"`
torches (`--strategy` does not apply). Tasks in `tasks/`:
- `protein_search.wdl`: tblastn of the subunit reference set against contigs
  (StxTyper's own search parameters), normalized to internal HSP JSON at the
  task boundary — BLAST format never leaks downstream
- `operon_assembly.wdl`: stitches frameshift-split HSPs, reduces a reference
  set's many accessions to one alignment per locus per class, then pairs them
  into candidate operons under synteny constraints (contig/strand/order/
  intergenic distance) across four passes of relaxing stringency
- `operon_call.wdl`: combined identity vs threshold, residue-table resolution
  for generalized classes, disruption status ladder, and suppression of
  candidates redundant with a better overlapping operon

Each task takes `torchbase/operon.py` as its `operon_module` File input and
imports it; the WDL is plumbing, the algorithm lives in one place.

### Strategy Selection

CLI routes to appropriate workflow:
```python
strategy_to_workflow = {
    "fast": "torchbase/workflows/builtin/fast_typing.wdl",
    "balanced": "torchbase/workflows/builtin/balanced_typing.wdl",
    "sensitive": "torchbase/workflows/builtin/sensitive_typing.wdl",
}
```

For `auto` strategy: CLI pre-analyzes inputs and picks fast/balanced/sensitive once.

`typing_model = "operon"` torches skip this entirely and route to `operon_typing.wdl`.

### Output Format

All workflows produce standardized JSON:
```json
{
  "profile_id": "ST1",
  "profile_type": "sequence_type",
  "scheme": "salmonella_mlst",
  "status": "known|novel_profile|novel_allele",
  "confidence": 0.98,
  "allele_profile": "salmonella_adk_1,salmonella_fumC_2",
  "allele_calls": {...},
  "method": {
    "strategy": "balanced",
    "alignment_used": false,
    "tools": ["sourmash", "minimap2"]
  },
  "notes": {
    // Strategy-specific metadata (alignment metrics, decision rationale, etc.)
  }
}
```

Operon output (`typing_model = "operon"`) is a **list** (multiple operons per
assembly are normal, e.g. an isolate with both stx1 and stx2), with an eight-
value `operon_status` ladder (most to least intact: `COMPLETE > COMPLETE_NOVEL
> AMBIGUOUS > PARTIAL > EXTENDED > PARTIAL_CONTIG_END > INTERNAL_STOP >
FRAMESHIFT`) alongside the coarse `status`. When assigning a status the
disruptions are checked in the reverse order — a frameshift outranks an
internal stop, which outranks truncation at a contig end. See
docs/operon-strategy-plan.md §5 for the full schema.

## Known Incomplete Features

- Full end-to-end validation/benchmarking (in-scope tests are unit/simple integration only)
- IPFS functionality partially implemented (error handling incomplete)
- Some conversion modules need completion (cgMLST, Chewie-NS)
- Operon typing model (docs/operon-strategy-plan.md): Phase 0 parity with
  StxTyper 1.0.45 is done (182/182 contigs of its golden suite agree on every
  field; `tests/test_operon_parity.py` pins it offline), and Phase 1 is done —
  `examples/etec_lt` is a second, hand-encoded scheme (ETEC *eltA*/*eltB*) that
  types real data through the same workflow, needing only the `intergenic_min`
  config key for its overlapping genes (`tests/test_operon_etec_lt.py`). Still
  open: broader concordance against MicroBIGG-E, and Phase 2
  (`torchtools derive-operon`).

## Entry Points

Defined in pyproject.toml:
- `torchbase` → `torchbase.cli:cli` (main user commands)
- `torchtools` → `torchbase.cli:tools` (authoring tools)

## Dependencies of Note

- `miniwdl`: Executes WDL workflows
- BLAST+ (`ncbi/blast` container): tblastn search for the operon typing model
- `zstandard`: File compression
- `ipyfs`: IPFS Python client
- `toml`: Metadata parsing
- `click`: CLI framework
- `cookiecutter`: Template system for torch generation
