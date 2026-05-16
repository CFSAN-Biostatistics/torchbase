# PRD: MLST Multi-Scheme Typing Implementation

## Problem Statement

Microbial typing tools like Tseemann's MLST have no mechanism to distribute database updates independently from the tool itself. PubMLST maintains typing schemes for 100+ organisms via BIGSdb REST API, but accessing this requires users to either run queries manually or install tools that bundle stale database snapshots. Additionally, existing tools don't provide quality assessment of typing databases, making it difficult to identify ambiguous or suspect alleles. Users need a way to create versioned, distributable MLST databases that can be updated independently and include comprehensive quality reporting.

## Solution

Implement the core MLST typing infrastructure for torchbase, including:
- PubMLST database conversion via BIGSdb REST API to hierarchical multi-scheme torch packages
- K-mer-based quality analysis with auto-tuned thresholds detecting overlapping and duplicate alleles
- IPNS-based torch registry with version pinning for reproducible deployments
- Default MLST workflow using MinHash for rapid scheme detection with alignment fallback for ambiguous cases
- Quality reports (JSON + human-readable with ASCII histograms) distributed with torch packages

This provides the foundation for distributing typing databases as versioned, immutable IPFS artifacts while maintaining quality provenance.

## User Stories

1. As a bioinformatician, I want to convert PubMLST schemes to torches, so that I can distribute versioned MLST databases to my team
2. As a tool developer, I want to pull MLST databases via IPFS, so that I don't need to bundle databases with my tool
3. As a researcher, I want quality reports for typing databases, so that I can assess data reliability before use
4. As a pipeline author, I want to pin torch versions in my workflow, so that results are reproducible across runs
5. As a database curator, I want to see k-mer similarity distributions per locus, so that I can identify problematic alleles
6. As a lab technician, I want automatic scheme detection from contigs, so that I don't need to know which organism I'm typing
7. As a system administrator, I want torch registries configurable per-project, so that different teams can use different database sources
8. As a quality analyst, I want hierarchical suspect flagging (allele→locus→profile), so that I can trace quality issues to their source
9. As a workflow developer, I want snapshot-based versioning with timestamps, so that I can track when database states were captured
10. As a bioinformatics researcher, I want to query BIGSdb with temporal filters, so that I can build torches from specific time points
11. As a typing pipeline user, I want depth-filtered k-mer matching on raw reads, so that I can skip assembly
12. As a data scientist, I want quality.json files with torch packages, so that I can programmatically filter suspect data
13. As a microbiologist, I want novel profile detection with nearest ST reporting, so that I can understand novel variants epidemiologically
14. As a core facility director, I want multi-scheme torches (all PubMLST in one package), so that users can type any organism with one download
15. As a developer, I want MinHash-based scheme auto-detection, so that typing is fast across 100+ schemes
16. As a pipeline maintainer, I want alignment fallback when MinHash is ambiguous, so that edge cases are handled precisely
17. As an MLST user, I want torches to follow the `schemes/<organism>/` structure, so that I can navigate databases easily
18. As a reproducibility advocate, I want --pin to lock all torch versions on first use, so that deployments stay consistent
19. As a fork maintainer, I want to publish my own IPNS registry, so that I can distribute modified torches independently
20. As a tool author, I want the default MLST workflow distributed as a torch, so that workflow updates don't require tool updates
21. As a database publisher, I want gap detection for k-mer similarity thresholds, so that outliers are flagged automatically
22. As a bioinformatician, I want CLI-overridable k-mer parameters, so that I can tune quality analysis for specific schemes
23. As a typing service operator, I want per-torch workflow configuration in metadata.toml, so that schemes can customize thresholds
24. As a quality reviewer, I want ASCII histogram visualization in terminal, so that I can quickly assess similarity distributions
25. As a researcher, I want both snapshot and semver versioning strategies, so that I can choose between date-based and change-based versioning
26. As a registry operator, I want IPNS-based registry updates, so that the registry can evolve without code changes
27. As a CI/CD engineer, I want directory-scoped torch pins (.torchbase.toml), so that projects have isolated version locks
28. As a typing analyst, I want novel allele detection below confidence threshold, so that I can identify potential new variants
29. As a workflow developer, I want convention-based workflow discovery (main.wdl), so that custom workflows are automatically detected
30. As a database validator, I want exact substring matching for overlap detection, so that perfect containment is flagged
31. As a sequence analyst, I want metadata recording of BIGSdb last_updated timestamps, so that provenance is traceable
32. As a power user, I want multiple registry support, so that I can combine official and custom torch sources
33. As a typing pipeline user, I want positive flag names (--include-suspect-alleles), so that I don't have to reason about double negatives
34. As a bioinformatics educator, I want self-contained torches with embedded quality reports, so that students can explore data quality

## Implementation Decisions

### Module Architecture

**K-mer Analysis Module** (`torchbase/quality/kmer_analysis.py`)
- Encapsulates all Jellyfish operations for pairwise k-mer similarity within loci
- Auto-tuning algorithm: gap detection (>2% jumps) in sorted similarity distribution, falls back to 99th percentile
- Returns structured `SimilarityReport` with statistics (min/median/99th percentile, suspect pairs with reasons)
- Interface: `analyze_locus(fasta_path: Path, k_size: int = 21) → SimilarityReport`
- Thresholds: overlap detection at 95% containment (one-way), duplicate detection at 98% symmetric similarity

**Quality Report Generator** (`torchbase/quality/report.py`)
- Pure data transformation: `List[SimilarityReport]` → formatted outputs
- Three output formats: human-readable text with ASCII histograms, structured JSON, both
- Hierarchical flagging: generates suspect_alleles, suspect_loci, suspect_profiles sections for metadata.toml
- Interface: `generate_report(locus_reports: List[SimilarityReport], format: str = 'text') → Report`
- ASCII histograms use simple bar charts showing similarity distribution bins

**BIGSdb REST Client** (`torchbase/conversions/bigsdb_client.py`)
- Wraps BIGSdb REST API with typed response models
- Handles pagination, temporal filtering (added_after, updated_after), scheme/locus/profile retrieval
- Returns dataclasses: `SchemeMetadata`, `LocusData`, `ProfileTable`
- Interface: `BIGSdbClient(base_url: str).fetch_scheme(database: str, scheme_id: str) → SchemeData`
- Captures `last_updated` timestamps from API for provenance

**PubMLST Converter** (`torchbase/conversions/pubmlst.py`)
- Orchestrates conversion workflow: API fetch → directory structure → k-mer analysis → metadata generation
- Implements hierarchical torch structure: `schemes/<organism>/{profiles.tsv, alleles/*.fasta}`
- Writes metadata.toml with provenance section (source URL, database ID, scheme ID, timestamps, fetch date)
- Generates quality.json and embeds summary stats in metadata.toml [data_quality] section
- CLI integration via `torchtools convert pubmlst <database_url> <scheme_id>`
- Supports CLI overrides: --kmer-size, --overlap-threshold, --duplicate-threshold

**Multi-Scheme Torch Loader** (update `torchbase/torchfs.py`)
- Changes `Torch` dataclass: `profile: Profile` → `schemes: Dict[str, Schema]`
- Discovers schemes by scanning `schemes/*/` subdirectories
- Each scheme loads independently: profiles.tsv + alleles/*.fasta
- Validates metadata.toml declares all discovered schemes
- Interface unchanged: `Torch.load(path: Path) → Torch`

**Registry Manager** (`torchbase/registry.py`)
- Manages multiple IPNS registries with hierarchical config (directory .torchbase.toml overrides user ~/.torchbase/config.toml)
- Resolves torch references: `namespace/name` → IPFS CID, respecting version constraints and pins
- Implements --pin behavior: on first use with --pin, write current versions to config; subsequent calls use pinned versions
- Pin dependencies: when pinning data torch, also pin referenced workflow torch
- Interface: `RegistryManager.fetch_torch(name: str, version: Optional[str], pin: bool) → Path`
- Default registry: `/ipns/registry.torchbase.org` (configurable)

**Versioning Support** (update `torchbase/torchbase.py`)
- Three strategy implementations: snapshot (ISO date strings), content-addressed (hashes), semantic (semver)
- metadata.toml declares strategy in `[version]` section with `strategy` and `timestamp` fields
- Cross-strategy comparison via required `timestamp` field (Unix epoch)
- Interface: `Version.parse(version_str: str, metadata: dict) → Version`, `Version.compare(v1: Version, v2: Version) → int`

**Default MLST Workflow** (new WDL, distributed as torch)
- Workflow torch at `torchbase/workflows/mlst/1.0.0.torch/` containing main.wdl
- Pipeline stages:
  1. Input validation (reads or contigs)
  2. Depth filtering for reads (histogram-based, fallback ≥3x coverage)
  3. MinHash sketching (sourmash) against all scheme alleles
  4. Allele calling from best matches
  5. Scheme inference from called alleles
  6. Profile lookup in scheme table
  7. Alignment fallback (minimap2) if ambiguous (triggers: top 2 within 3%, best <92%, >2 schemes, coverage <80%)
  8. Novel profile/allele reporting with nearest ST
- Parameterized via metadata.toml thresholds
- Containerized tasks (Docker/Singularity) for sourmash, jellyfish, minimap2

**Workflow Discovery** (update `torchbase/cli.py`)
- Convention: torch with `main.wdl` uses that file; otherwise use default workflow
- Default workflow resolution: fetch `torchbase/default-workflow` (latest) from registry
- Validation: torch with WDL files must be named `main.wdl` or torch creation fails
- Users can override: `torch run --workflow torchbase/workflows/mlst/1.0.0`

### Data Structures

**Torch Metadata Additions** (metadata.toml):
```toml
[version]
string = "2024-03-15"
strategy = "snapshot"
timestamp = 1710489600

[provenance]
source = "pubmlst.org/bigsdb"
database = "pubmlst_senterica_seqdef"
scheme_id = "1"
last_updated = "2024-03-15T10:23:45Z"
fetched_at = "2024-03-20T14:30:00Z"

[data_quality]
report = "quality.json"
suspect_alleles_count = 8
suspect_loci_count = 2
suspect_profiles_count = 5

[typing]
method = "mlst"

[schemes]
"senterica" = "schemes/senterica"
"ecoli-1" = "schemes/ecoli-1"
```

**Quality Report JSON Schema**:
```json
{
  "loci": {
    "aroC": {
      "allele_count": 127,
      "kmer_size": 21,
      "similarity_stats": {
        "min": 0.852,
        "median": 0.943,
        "percentile_99": 0.981
      },
      "threshold": 0.985,
      "threshold_method": "gap_detection",
      "suspect_pairs": [
        {
          "allele_1": "aroC_45",
          "allele_2": "aroC_102",
          "similarity": 0.998,
          "type": "duplicate"
        }
      ]
    }
  },
  "summary": {
    "total_loci": 7,
    "total_alleles": 823,
    "suspect_alleles": 8,
    "suspect_loci": 2
  }
}
```

**Registry Config** (~/.torchbase/config.toml):
```toml
[registries]
default = "/ipns/registry.torchbase.org"
additional = ["/ipns/myorg.torches.org"]

[pins]
"pubmlst/mlst-database" = "2024-03-15"
"torchbase/default-workflow" = "1.0.0"
```

### Technical Clarifications

- K-mer size default (k=21) chosen for typical MLST gene lengths (400-600bp) with low divergence (few SNPs between alleles)
- Gap detection threshold: >2% jump in sorted similarities
- Jellyfish selected over KMC for simpler API, good Bioconda support
- Sourmash selected over Mash for Python integration, active maintenance
- MinHash confidence thresholds configurable per-torch but defaults are: ambiguity margin 3%, confidence floor 92%, coverage 80%
- IPNS registry updates managed by torchbase maintainers, not automated
- Workflow torch follows same structure as data torches (metadata.toml, quality.json placeholder)
- Pin scope: directory .torchbase.toml for project-specific, ~/.torchbase/config.toml for user-global

### Architectural Decisions

- Multi-scheme torches preferred over per-scheme torches (reduces IPFS overhead, mirrors PubMLST organization)
- Workflow-as-torch reuses entire distribution infrastructure rather than special-casing workflow updates
- Quality analysis runs at conversion time, not at typing time (static analysis where possible, runtime discoveries deferred to future scope)
- Convention over configuration for workflow discovery reduces required metadata fields
- Positive flag semantics (--include-suspect-alleles) reduces cognitive load
- Hierarchical pins (directory overrides user) matches package manager patterns (package-lock.json vs global config)

### Schema Changes

- `Torch` dataclass: add `schemes: Dict[str, Schema]`, remove `profile: Profile`
- `Schema` class: add `locus_paths: Dict[str, Path]` for allele FASTA locations
- `Profile` class: already supports wildcards (Special.IGNORE), no changes needed
- metadata.toml: add `[version]`, `[provenance]`, `[data_quality]`, `[typing]`, `[schemes]` sections

### API Contracts

BIGSdb REST API endpoints used:
- GET `/db/{database}/schemes/{scheme_id}` - scheme metadata
- GET `/db/{database}/schemes/{scheme_id}/loci` - locus list
- GET `/db/{database}/loci/{locus}/alleles_fasta` - allele sequences
- GET `/db/{database}/schemes/{scheme_id}/profiles_csv` - profile table
- Temporal filters: `added_after`, `updated_after`, `last_updated`

Registry API (IPNS resolution):
- Fetch manifest from `/ipns/<registry_name>`
- Manifest format: TOML mapping `namespace/name/version` → IPFS CID
- Registry pinning via torchbase version: each release declares default registry and default workflow version

WDL workflow interface:
- Inputs: `File? reads`, `File? contigs`, `File profiles`, `Array[File] references`, `String? scheme_hint`
- Outputs: `File results_json`, `String sequence_type`, `String status`, `String? nearest_st`
- Runtime requirements: Docker/Singularity images for sourmash, minimap2

## Testing Decisions

### What Makes a Good Test

- **Test external behavior, not implementation**: Test that `analyze_locus()` returns correct suspect pairs, not that it calls Jellyfish with specific flags
- **Use synthetic data for quality modules**: Generate FASTA files with known overlap/duplicate patterns
- **Mock external services**: BIGSdb client tests use recorded API responses, not live HTTP calls
- **Test version comparison logic thoroughly**: All three strategies, edge cases (same timestamp different strings, invalid formats)
- **Validate torch structure**: Loader tests verify scheme discovery, metadata validation, error handling for malformed torches

### Modules with Comprehensive Test Coverage

**Registry Manager** (`torchbase/registry.py`)
- Pin behavior: first use writes config, subsequent uses read config, directory overrides global
- Version resolution: latest without constraint, pinned versions, explicit versions
- Multi-registry fallback: try default, then additional registries in order
- Dependency pinning: pinning data torch also pins workflow torch
- Prior art: Similar to package manager tests (pip, npm lockfiles)

**Versioning Support** (`torchbase/torchbase.py`)
- Parsing all three strategies from version strings + metadata
- Cross-strategy comparison via timestamps
- Invalid format handling (malformed semver, non-ISO dates, missing timestamp)
- Sorting: multiple versions of same torch, mixed strategies
- Prior art: Version comparison in setup.py, existing Profile equality tests

**Multi-Scheme Torch Loader** (`torchbase/torchfs.py`)
- Hierarchical scheme discovery from directory structure
- Metadata validation: declared schemes match discovered schemes
- Error handling: missing profiles.tsv, missing alleles/, malformed FASTA
- Backward compatibility: single-scheme torches (no schemes/ directory) still load
- Prior art: Existing `Torch.load()` tests, Profile.parse() tests

**Quality Report Generator** (`torchbase/quality/report.py`)
- All three output formats (text, JSON, both)
- Hierarchical flagging: suspect pairs → suspect loci → suspect profiles propagation
- ASCII histogram rendering: various distribution shapes (normal, bimodal, uniform)
- Edge cases: no suspect data, all alleles suspect, empty loci
- Prior art: tabulate usage in CLI (similar formatting patterns)

### Testing Approach

- Unit tests for pure functions (quality report generation, version comparison)
- Integration tests for modules with external dependencies (BIGSdb client with mocked responses, Jellyfish via temp files)
- Fixture-based testing: use existing bigsdb_fixture.py patterns, expand with multi-scheme examples
- End-to-end test: convert small synthetic PubMLST scheme → torch → load → validate structure
- WDL workflow testing via miniwdl check (syntax validation), defer runtime testing to manual validation

## Out of Scope

- Runtime quality annotations (users reporting newly discovered ambiguities) - future scope
- Validation datasets distributed as torches - future scope
- Migration wizard for converting arbitrary tools - separate PRD
- SeqSero2 or other non-MLST typing systems - separate implementations
- Web UI or TUI for torch management - CLI only initially
- Automated registry submission (PR creation) - manual process
- Update tracking for source databases (re-running conversion on PubMLST changes) - manual re-conversion
- Full assembly workflow integration - reads support via depth filtering only
- Read QC (fastp, trimmomatic) - rely on depth filtering for quality
- Custom alignment parameters per-locus - use workflow defaults
- Workflow composition (importing WDL tasks from other torches) - single workflow per torch
- Data torch dependencies (one torch referencing another torch's data) - workflow dependencies only

## Further Notes

- This implementation provides the foundation for torchbase's vision: decoupling typing databases from typing tools
- The k-mer quality analysis approach is novel in this domain and may reveal database quality issues not previously documented
- Default workflow design (MinHash + alignment fallback) balances speed and accuracy, but individual torches can override with custom workflows if needed
- IPNS-based registry is chosen over git-based or HTTP-based registries for decentralization and IPFS ecosystem integration
- Multi-scheme torch structure matches PubMLST's organization, making it intuitive for users already familiar with that resource
- Snapshot versioning (dates) is recommended for database torches, but semver is appropriate for workflow torches
- The pin mechanism balances "use latest" (default) with "lock for reproducibility" (opt-in via --pin), following modern package manager patterns
- Quality report ASCII histograms make terminal-based exploratory analysis viable without requiring plotting libraries
- Convention-based workflow discovery (main.wdl) reduces configuration burden while allowing full customization when needed
- This design intentionally avoids assembly to keep the stack simple and fast, relying on depth filtering to handle read errors
