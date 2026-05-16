# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Torchbase is a Python framework for microbial serotyping prediction from sequencing data. It addresses the problem that most microbial typing pipelines are written by domain experts without software engineering backgrounds, resulting in poor UX, reliability, maintainability, and database distribution/upkeep.

**Core Concept**: "Torches" are versioned, distributed databases containing:
- Allele reference sequences
- Allelic profile tables (schema definitions)
- WDL workflows for execution
- Build files and metadata

Torches are distributed via IPFS to enable versioned, reproducible typing across different users and institutions.

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
     - `torchbase`: User-facing commands (list, pull, info, run)
     - `torchtools`: Authoring commands (build, version, convert)
   - Execution via `miniwdl` for WDL workflows
   - Automatic file decompression/compression to zstandard format

### Torch Package Structure

```
<namespace>/<torchname>/<version>.torch/
├── metadata.toml           # Package metadata, citations, maintainers
├── <buildname>.profiles.tsv  # Tab-separated allelic profile table
├── <buildname>.wdl         # Main WDL workflow
├── <torchname>.build.wdl   # Build workflow for database
└── _resources/             # Reference FASTA files for alleles
```

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

## Known Incomplete Features

- Conversion modules in `torchbase/conversions/` are empty stubs
- IPFS functionality partially implemented (error handling incomplete)
- WDL templates need completion (see `torchbase/templates/torch/`)
- Docker image build target defined but not tested

## Entry Points

Defined in setup.cfg:
- `torchbase` → `cli:cli` (main user commands)
- `torchtools` → `cli:tools` (authoring tools)

## Dependencies of Note

- `miniwdl`: Executes WDL workflows
- `zstandard`: File compression
- `ipyfs`: IPFS Python client
- `toml`: Metadata parsing
- `click`: CLI framework
- `cookiecutter`: Template system for torch generation
