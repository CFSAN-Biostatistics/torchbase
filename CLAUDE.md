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

**Typing Strategies**: Users choose speed/accuracy tradeoff:
- `fast`: MinHash-based calling only (fastest)
- `balanced`: MinHash with alignment fallback (default)
- `sensitive`: Full alignment-based calling (most accurate)
- `auto`: Automatically selects based on input characteristics

## Architecture

### Three-Layer System

1. **Torch Definition Layer** (`torchbase/torchbase.py`)
   - `Schema`: Container for typing profiles with version info
   - `Profile`: Represents allelic profiles with special handling for wildcards (`IGNORE = "?"`) and exclusions (`EXCLUDE = "X"`)
   - Profile equality supports multiple formats: tuples, dicts, PubMLST-style strings (e.g., "locus_allele")

2. **Filesystem/Distribution Layer** (`torchbase/torchfs.py`)
   - `Torch` dataclass: Loads and validates torch packages from disk
   - IPFS integration for distributed torch retrieval (via `ipyfs`)
   - Manifest system for tracking available torches
   - Environment-based IPFS configuration (`TORCHBASE_IPFS_NODE`, `TORCHBASE_IPFS_PORT`)

3. **CLI Layer** (`torchbase/cli.py`)
   - Two command groups:
     - `torchbase`: User-facing commands (list, pull, info, run, workflow)
     - `torchtools`: Authoring commands (build, version, convert)
   - Strategy-based workflow routing:
     - Built-in workflows in `torchbase/workflows/builtin/` (fast/balanced/sensitive)
     - Custom workflows via torch-embedded `main.wdl`
     - `--strategy` flag selects typing approach (error if used with embedded workflows)
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

**Workflow Discovery**:
- If torch has `main.wdl` → use it (user cannot specify `--strategy`)
- If torch has no workflow → use built-in workflow with selected strategy
- CLI concatenates multi-scheme torches with scheme-prefixed locus names (e.g., `salmonella_adk_1`)

### Conversion System (`torchbase/conversions/`)

Converts external typing schemes to torch format:
- `pubmlst.py`: PubMLST MLST schemes
- `pubcgmlst.py`: PubMLST cgMLST schemes
- `shigatyper.py`: ShigaTyper database
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
- Parses profile table using `Profile.parse()`
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

## Known Incomplete Features

- Full end-to-end validation/benchmarking (in-scope tests are unit/simple integration only)
- IPFS functionality partially implemented (error handling incomplete)
- Some conversion modules need completion (cgMLST, Chewie-NS)

## Entry Points

Defined in pyproject.toml:
- `torchbase` → `torchbase.cli:cli` (main user commands)
- `torchtools` → `torchbase.cli:tools` (authoring tools)

## Dependencies of Note

- `miniwdl`: Executes WDL workflows
- `zstandard`: File compression
- `ipyfs`: IPFS Python client
- `toml`: Metadata parsing
- `click`: CLI framework
- `cookiecutter`: Template system for torch generation
