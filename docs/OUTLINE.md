# Torchbase Documentation Outline

**Target:** ReadTheDocs site with comprehensive user guide and authoring documentation.

---

## I. Introduction

### 1.1 Welcome
- What is Torchbase?
- Problem statement: reproducibility, versioning, distribution challenges in typing systems
- Core concept: "torch" as versioned, distributed typing database

### 1.2 Key Features
- Reproducible typing (version pinning)
- IPFS-based distribution
- Multi-scheme support
- Flexible typing strategies (fast/balanced/sensitive/auto)
- Quality analysis and filtering
- Cryptographic signing and verification

### 1.3 Quick Start
- Installation
- Download a torch
- Run a typing workflow
- Interpret results

### 1.4 Use Cases
- MLST typing
- Serotyping (Salmonella, E. coli, Shigella, Listeria)
- cgMLST/wgMLST
- Custom allelic profile systems

---

## II. Concepts and Theory

### 2.1 What is a Torch?
- Definition and structure
- Components: alleles, profiles, workflows, metadata
- Single-scheme vs multi-scheme torches
- Versioning semantics (semver)

### 2.2 Allelic Profile Typing
- Theory: allele calling → profile lookup
- Schema and profile definitions
- Special values: wildcards (`?`), exclusions (`X`)
- Novel alleles vs novel profiles

### 2.3 Distribution Model
- IPFS content addressing
- Content-addressable storage (CID)
- Manifest system and registries
- IPNS publishing for mutable references

### 2.4 Typing Strategies
- **Fast:** MinHash-only (sketch-based similarity)
- **Balanced:** MinHash + conditional alignment fallback
- **Sensitive:** Full alignment-based calling
- **Auto:** Input-driven strategy selection
- Tradeoffs: speed vs accuracy
- When to use each strategy

### 2.5 Quality and Provenance
- K-mer similarity analysis
- Suspect allele detection (duplicates, overlapping)
- Filtering at typing time
- Cryptographic signing (Ed25519, P-256/YubiKey PIV)
- IPLD commit chains for tamper-evident history

---

## III. User Guide

### 3.1 Installation
- PyPI install (`pip install torchbase`)
- Development install
- Dependencies and requirements
- Platform notes (Linux/Mac/WSL)

### 3.2 Configuration
- Config hierarchy: user (`~/.torchbase/config.toml`) vs project (`.torchbase.toml`)
- Registries configuration
- IPFS node configuration (`TORCHBASE_IPFS_NODE`, `TORCHBASE_IPFS_PORT`)
- Version pinning
- Trusted keys and key registries

### 3.3 Discovering and Fetching Torches
- `torchbase list` — browse available torches
- `torchbase info` — inspect torch metadata
- `torchbase pull` — download a torch
- `torchbase pin` — lock version for reproducibility
- Working with multiple registries

### 3.4 Running Typing Workflows
- Input formats: contigs, paired reads, interleaved, long reads
- Automatic compression detection and conversion (zstd)
- Choosing a strategy (`--strategy fast|balanced|sensitive|auto`)
- Interpreting result JSON
  - `profile_id`, `profile_type`, `scheme`, `status`, `confidence`
  - Allele calls and method metadata
- Multi-scheme torches: scheme selection and output

### 3.5 Quality Filtering at Runtime
- Using `quality.json` for suspect data filtering
- `--exclude-suspect-alleles` vs `--exclude-suspect-loci` vs `--exclude-suspect-profiles`
- Filtering hierarchy and implications
- Reviewing exclusion metadata in results
- When and why to filter

### 3.6 Workflows and WDL
- Built-in workflows: `fast_typing.wdl`, `balanced_typing.wdl`, `sensitive_typing.wdl`
- Torch-embedded custom workflows (`main.wdl`)
- `torchbase workflow inspect` — visualize workflow logic
- Understanding WDL execution (miniwdl backend)
- Workflow outputs and provenance

### 3.7 Verification and Trust
- `torchbase verify` — check torch signatures
- Public key resolution (flag, config, key registry, embedded)
- Requiring signatures (`--require-signature`)
- Trust model: namespace keys, chain-of-custody via IPLD

---

## IV. Authoring Torches

### 4.1 Overview
- When to create a new torch
- Namespace conventions and governance
- torch vs torchtools CLI split

### 4.2 Torch Package Structure
- Directory layout: `<namespace>/<name>/<version>.torch/`
- `metadata.toml` schema
- `profiles.tsv` format
- `_resources/` or `schemes/*/alleles/` FASTA files
- Optional `main.wdl` workflow
- `quality.json` report

### 4.3 Converting Existing Schemes
#### 4.3.1 PubMLST (MLST)
- `torchtools convert-pubmlst` — fetch from API
- `torchtools convert pubmlst` — from local files
- Scheme ID resolution, locus mapping

#### 4.3.2 PubMLST (cgMLST)
- `torchtools convert pubcgmlst` — local cgMLST files
- Handling hundreds of loci

#### 4.3.3 Salmonella Serotyping
- `torchtools convert seqsero2` — UGA/Deng lab database
- `torchtools convert seqsero2s` — LSTUGA curated database (antigen + 7-gene MLST)
- Antigen notation (White-Kauffmann-Le Minor)

#### 4.3.4 E. coli / Shigella Serotyping
- `torchtools convert ectyper` — phac-nml database
- Combined vs split FASTA input

#### 4.3.5 Shigella Serotyping
- `torchtools convert shigatyper` — CFSAN database
- O/H antigen loci

#### 4.3.6 Listeria Serogroups
- `torchtools convert lissero` — MDU-PHL database
- Presence/absence profile logic

#### 4.3.7 Chewie-NS (wgMLST)
- `torchtools convert chewie-ns` — placeholder, future feature

### 4.4 Building from Scratch
- Preparing allele FASTA files
- Writing profiles TSV
- Writing metadata.toml
- Cookiecutter templates for torch scaffolding
- `torchtools build` — validate torch structure

### 4.5 Quality Analysis
- K-mer analysis options: `--kmer-size`, `--overlap-threshold`, `--duplicate-threshold`
- Interpreting `quality.json` output
  - `loci[].suspect`, `loci[].similarities`, `loci[].alleles[].suspect`
  - `profiles[].suspect`
- Tuning thresholds for conserved gene families
- Using quality reports for torch curation

### 4.6 Versioning
- Semantic versioning (major.minor.patch)
- `torchtools version` — bump version in metadata
- When to increment major/minor/patch
- IPFS CID checkpointing

### 4.7 Custom Workflows
- When to embed a `main.wdl`
- WDL task definitions (`tasks/` directory)
- Input/output contract for `main.wdl`
- Strategy flag behavior (disabled when custom workflow present)
- Testing custom workflows locally

---

## V. Publishing and Distribution

### 5.1 Cryptographic Signing
- Key generation: `torchtools keygen`
- Software Ed25519 keys (PEM)
- YubiKey PIV-backed keys (P-256)
- Signing a torch: `torchtools sign`
- Signature format (`signature.toml`)

### 5.2 IPFS Upload
- Manual IPFS add + pinning
- `torchtools publish` — sign + upload + sign CID
- IPFS node configuration
- Pinning strategies

### 5.3 Namespace Management and IPLD Commit Chains
- Namespace registration: `torchtools namespace register`
- Genesis block creation
- IPNS key import and publishing
- Chain structure: genesis → update blocks
- Tamper-evident history (previous-CID links)

### 5.4 Publishing Updates
- `torchtools manifest add` — complete workflow
- Steps: sign torch → upload → sign CID → create update block → publish to IPNS
- YubiKey signing in CI/CD
- Version conflict resolution

### 5.5 Registries and Discovery
- Flat-manifest registries (HTTP/IPFS)
- Chain-based namespace registries (IPLD + IPNS)
- Resolution order: chain-based → flat-manifest
- Configuring multiple registries
- `torchtools manifest show` / `torchtools namespace show`

### 5.6 CT-Style Log Operators
- What log operators do (permissionless recording of genesis CIDs)
- Trust model: no approval, just transparency
- Configuring log submission (`log_operators` in config)
- `--submit-to` flag override
- Future: log auditing and gossip

---

## VI. Advanced Topics

### 6.1 Multi-Scheme Torches
- Use cases: single organism with multiple typing systems
- `schemes/` directory structure
- CLI concatenation of scheme-prefixed locus names
- Profile resolution across schemes

### 6.2 Profile Comparison Logic
- Flexible equality: tuples, dicts, PubMLST-style strings
- Special value semantics (IGNORE, EXCLUDE)
- Length validation and locus exclusion
- Profile matching algorithm

### 6.3 Registry Management
- Hierarchical resolution
- Fallback registries
- Cache TTL and key cache (`key_cache.toml`)
- Adding custom registries

### 6.4 IPFS Backend Details
- Kubo/go-ipfs node setup
- Environment variable overrides
- Multipart upload for large torches
- Pinning services (Pinata, Infura, etc.)
- IPFS gateway fallback

### 6.5 Extending Torchbase
- Adding new converters (plugin architecture, future)
- Custom workflow patterns
- Integrating torchbase in pipelines (NextFlow, Snakemake, etc.)

### 6.6 Performance Tuning
- Strategy selection for throughput
- Parallel batch execution
- MinHash k-mer size and sketch parameters
- Alignment presets (asm20/asm5/asm5+eqx)

---

## VII. Development and Contributing

### 7.1 Development Setup
- `make install-dev`
- Running tests (`make test`, `make test-all`)
- Code coverage (`make coverage`)
- Linting (`make lint`)

### 7.2 Project Architecture
- Three-layer system: definition, filesystem, CLI
- Module structure: `torchbase.py`, `torchfs.py`, `cli.py`
- Conversions module (`torchbase/conversions/`)
- Workflows (`torchbase/workflows/builtin/`)

### 7.3 Testing
- Test fixtures (`bigsdb_fixture.py`)
- Profile parsing and equality tests
- IPFS stubbing
- Workflow execution tests

### 7.4 Contributing Guidelines
- Fork and branch workflow
- Writing tests for new features
- Documentation requirements
- PR review process

### 7.5 Roadmap and Known Gaps
- Chewie-NS converter implementation
- Full end-to-end benchmarking
- IPFS error handling completion
- Plugin architecture for converters
- Log operator gossip protocol

---

## VIII. Reference

### 8.1 CLI Reference
- `torchbase` commands
  - `list`, `info`, `pull`, `pin`, `run`, `workflow`, `verify`
- `torchtools` commands
  - `convert`, `convert-pubmlst`, `build`, `version`, `keygen`, `sign`, `publish`, `namespace`, `manifest`

### 8.2 File Formats
- `metadata.toml` schema
- `profiles.tsv` format
- `quality.json` schema
- `signature.toml` format
- IPLD block schemas (genesis, update)

### 8.3 Configuration Reference
- `~/.torchbase/config.toml` options
- `.torchbase.toml` project overrides
- Environment variables

### 8.4 Workflow Output Schema
- JSON result structure
- Status codes: `known`, `novel_profile`, `novel_allele`
- Method metadata fields
- Notes and exclusions

### 8.5 API Reference (Python)
- `torchbase.torchbase` module
  - `Schema`, `Profile`, `Special`
- `torchbase.torchfs` module
  - `Torch`, `TorchLoader`
- `torchbase.registry` module
  - `RegistryManager`, `RegistryConfig`
- `torchbase.quality.kmer_analysis`
  - `analyze_locus`, quality report structure

---

## IX. Appendices

### A. Glossary
- Allele, locus, profile, scheme, torch, CID, IPFS, IPNS, IPLD, WDL, etc.

### B. FAQ
- Common issues and troubleshooting
- IPFS connectivity problems
- Version pinning conflicts
- Signature verification failures

### C. Citation and License
- How to cite Torchbase
- License terms
- Third-party database licenses (PubMLST, SeqSero2, etc.)

### D. Example Workflows
- End-to-end MLST typing
- Batch processing with shell scripting
- Integrating with bioinformatics pipelines
- CI/CD for torch publishing

---

## Notes

- Documentation generated from commit: [TBD]
- Last updated: 2026-06-10
- Target ReadTheDocs theme: sphinx_rtd_theme
- Syntax: reStructuredText (`.rst`) or MyST Markdown (`.md`)
