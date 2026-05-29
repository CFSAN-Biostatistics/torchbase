# Torchbase

**A Python framework for reproducible microbial serotyping from sequencing data**

Torchbase addresses challenges in microbial typing pipelines around database distribution, versioning, and reproducibility. It provides a standardized, versioned, and distributed approach to microbial typing databases and workflows.

## What is a Torch?

A "torch" is a versioned, self-contained typing database package distributed via IPFS that contains:

- **Allele reference sequences** - FASTA files in `_resources/`
- **Allelic profile tables** - Tab-separated schema definitions
- **WDL workflows** - Execution workflows for typing
- **Metadata** - Citations, maintainers, build information

Torches are distributed via IPFS to enable versioned, reproducible typing across different users and institutions.

## Features

- **Reproducible typing** - Pin specific torch versions for deterministic results
- **Distributed databases** - IPFS-based distribution eliminates single points of failure
- **Multiple scheme support** - Convert from PubMLST, cgMLST, ShigaTyper, and more
- **Flexible typing strategies** - Choose speed/accuracy tradeoff (fast/balanced/sensitive/auto)
- **Generalized allelic typing** - Works for MLST, serotyping, and other allelic profile systems
- **Multi-scheme torches** - Single torch can contain multiple typing schemes
- **Quality analysis** - K-mer analysis for detecting similar/duplicate alleles
- **Registry system** - Hierarchical configuration with fallback registries
- **WDL workflows** - Standardized execution via miniwdl

## Installation

```bash
pip install torchbase
```

For development:

```bash
git clone https://github.com/CFSAN-Biostatistics/torchbase
cd torchbase
make install-dev
```

## Quick Start

### List available torches

```bash
torchbase list
```

### Pull a torch

```bash
torchbase pull pubmlst/ecoli
```

### Get torch information

```bash
torchbase info pubmlst/ecoli
```

### Run typing on sequencing data

```bash
# Run with default balanced strategy
torchbase run pubmlst/ecoli --contigs contigs.fasta

# Choose typing strategy (fast/balanced/sensitive/auto)
torchbase run pubmlst/ecoli --reads reads.fastq --strategy fast
torchbase run pubmlst/ecoli --contigs contigs.fasta --strategy sensitive

# Auto strategy analyzes input and picks optimal approach
torchbase run pubmlst/ecoli --contigs contigs.fasta --strategy auto
```

**Typing Strategies:**
- `fast` - MinHash-based calling only, fastest (best for high-quality assemblies)
- `balanced` - MinHash with alignment fallback if needed (default, good for most cases)
- `sensitive` - Full alignment-based calling, most accurate (best for challenging samples)
- `auto` - Automatically selects strategy based on input characteristics

### Pin a torch version

```bash
torchbase pin pubmlst/ecoli 1.2.0
```

## Configuration

Torchbase uses hierarchical configuration with two levels:

1. **User config**: `~/.torchbase/config.toml` - Global settings
2. **Project config**: `.torchbase.toml` - Project-specific overrides

### Example configuration

```toml
[registries]
default = "https://registry.torchbase.org/manifest.toml"
additional = [
    "https://alt-registry.example.com/manifest.toml"
]

[pins]
"pubmlst/ecoli" = "1.2.0"
"pubmlst/salmonella" = "2.1.5"
```

Pins lock torch versions for reproducibility. Project pins override user pins.

## Creating Torches

### Convert existing schemes

```bash
# Convert PubMLST MLST scheme
torchtools convert pubmlst --scheme ecoli

# Convert PubMLST cgMLST scheme
torchtools convert pubcgmlst --scheme listeria

# Convert ShigaTyper database
torchtools convert shigatyper
```

### Build a torch

```bash
torchtools build <namespace>/<torchname>/<version>.torch
```

### Version a torch

```bash
torchtools version <namespace>/<torchname> --increment minor
```

## Architecture

Torchbase consists of three layers:

### 1. Torch Definition Layer (`torchbase/torchbase.py`)

- `Schema`: Container for typing profiles with version info
- `Profile`: Represents allelic profiles with wildcard (`?`) and exclusion (`X`) support
- Profile equality supports tuples, dicts, and PubMLST-style strings

### 2. Filesystem/Distribution Layer (`torchbase/torchfs.py`)

- `Torch` dataclass: Loads and validates torch packages
- IPFS integration via `ipyfs`
- Manifest system for tracking available torches
- Environment-based IPFS configuration

### 3. CLI Layer (`torchbase/cli.py`)

- `torchbase`: User-facing commands (list, pull, info, run, workflow)
- `torchtools`: Authoring commands (build, version, convert)
- Strategy-based workflow routing (fast/balanced/sensitive/auto)
- Automatic file compression to zstandard format

## Workflow Strategies

Torchbase provides three built-in typing strategies that balance speed and accuracy:

### Fast Strategy
- **Method**: MinHash-based similarity only
- **Speed**: Fastest (~1-2 min for typical MLST)
- **Best for**: High-quality assemblies, large batches, screening
- **Pipeline**: Sketch → Compare → Call alleles → Lookup profile

### Balanced Strategy (Default)
- **Method**: MinHash with conditional alignment fallback
- **Speed**: Moderate (~2-5 min)
- **Best for**: Most use cases, mixed quality data
- **Pipeline**: Sketch → Compare → Call alleles → If confidence <85% → Align → Refine

### Sensitive Strategy
- **Method**: Full alignment-based calling
- **Speed**: Slower (~5-15 min)
- **Best for**: Novel alleles, difficult samples, maximum accuracy
- **Pipeline**: Sketch (guide) → Align → Refine calls → Lookup profile

### Auto Strategy
- **Method**: Analyzes input and picks appropriate strategy
- **Logic**: Contigs → fast, Reads → balanced, Edge cases → balanced
- **Best for**: Unknown data quality, automated pipelines

### Inspecting Workflows

Visualize any workflow's pipeline:

```bash
# View built-in workflow
torchbase workflow inspect balanced

# View torch-embedded workflow
torchbase workflow inspect path/to/torch/
```

## Torch Package Structure

### Single-Scheme Torch (Simple)

```
<namespace>/<torchname>/<version>.torch/
├── metadata.toml              # Package metadata, citations
├── profiles.tsv               # Allelic profile table
├── main.wdl                   # Optional: custom workflow
└── _resources/                # Reference FASTA files
    ├── locus1.fasta
    └── locus2.fasta
```

### Multi-Scheme Torch (Advanced)

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

### Custom Workflows

Torches can include a `main.wdl` workflow for custom typing logic. If present, the torch's workflow is used instead of built-in strategies. Note: `--strategy` flag cannot be used with torch-embedded workflows.

## Quality Control and Suspect Data Filtering

Torchbase provides built-in quality analysis and filtering capabilities to handle suspect or low-quality alleles in typing databases.

### K-mer Analysis

Torchbase includes k-mer analysis for quality control:

```python
from torchbase.quality.kmer_analysis import analyze_locus

report = analyze_locus(fasta_path, k=21)
print(f"Found {len(report.suspect_pairs)} suspect pairs")
```

### Filtering Suspect Data

When running typing workflows, you can provide a `quality.json` file to filter suspect alleles, loci, or profiles. This is useful for databases with known quality issues or to exclude problematic sequences.

#### Using Quality Filtering

```bash
# Exclude suspect alleles only (default behavior)
torchbase run pubmlst/ecoli --contigs sample.fasta \
    --quality-json quality.json \
    --exclude-suspect-alleles

# Exclude all alleles from suspect loci
torchbase run pubmlst/ecoli --contigs sample.fasta \
    --quality-json quality.json \
    --exclude-suspect-loci

# Exclude all loci from suspect profiles (most aggressive)
torchbase run pubmlst/ecoli --contigs sample.fasta \
    --quality-json quality.json \
    --exclude-suspect-profiles
```

**Filtering Levels (hierarchical):**
1. `--exclude-suspect-alleles` - Excludes only specific flagged alleles
2. `--exclude-suspect-loci` - Excludes all alleles from flagged loci (implies level 1)
3. `--exclude-suspect-profiles` - Excludes all loci from flagged profiles (implies levels 1 & 2)

Note: If no `--quality-json` is provided, exclusion flags are silently ignored and no filtering occurs.

#### Quality.json Schema

The `quality.json` file contains quality annotations for alleles, loci, and profiles:

```json
{
  "loci": {
    "locus_name": {
      "suspect": false,
      "threshold": 90.0,
      "similarities": {
        "allele_1-allele_2": 45.5,
        "allele_1-allele_3": 98.5
      },
      "alleles": {
        "1": {
          "suspect": false,
          "length": 450,
          "gc_content": 52.3
        },
        "2": {
          "suspect": true,
          "reason": "low similarity to other alleles"
        }
      },
      "statistics": {
        "mean": 72.0,
        "std_dev": 31.2,
        "min": 45.5,
        "max": 98.5,
        "percentile_99": 97.0,
        "threshold_type": "percentile"
      }
    }
  },
  "profiles": {
    "ST1": {
      "suspect": false,
      "loci": ["locus1", "locus2", "locus3"]
    },
    "ST42": {
      "suspect": true,
      "loci": ["locus1", "locus2"],
      "reason": "incomplete profile"
    }
  }
}
```

**Key fields:**
- `loci[].suspect`: Boolean flag marking entire locus as suspect
- `loci[].similarities`: Pairwise allele similarity scores (pairs below threshold are suspect)
- `loci[].threshold`: Similarity cutoff for suspect pairs (default: 90.0)
- `loci[].alleles[].suspect`: Per-allele quality flags
- `profiles[].suspect`: Boolean flag marking entire profile as suspect
- `profiles[].loci`: List of loci in this profile

Allele pairs in the `similarities` object with scores below the threshold are automatically marked as suspect. This enables detection of duplicate or highly similar alleles that may represent sequencing errors or database quality issues.

#### Result Metadata

When filtering is applied, the typing result includes exclusion metadata:

```json
{
  "profile_id": "ST1",
  "status": "known",
  "confidence": 0.95,
  "notes": {
    "exclusions": {
      "excluded_alleles": ["locus1_42", "locus2_7"],
      "excluded_loci": ["locus3"],
      "num_excluded_alleles": 2,
      "num_excluded_loci": 1
    }
  }
}
```

This allows you to track which alleles/loci were filtered during typing for provenance and quality control purposes.

## Development

### Running tests

```bash
make test                # Run pytest
make test-all           # Test on all Python versions
make coverage           # Generate coverage report
```

### Code quality

```bash
make lint               # Run flake8
```

### Building distributions

```bash
make dist               # Build source and wheel
make release            # Upload to PyPI
```

## Registry System

The `RegistryManager` resolves torch references to IPFS CIDs:

```python
from torchbase.registry import RegistryManager
from torchbase.config import RegistryConfig

config = RegistryConfig.load()
manager = RegistryManager(config)

# Resolve torch to CID
cid = manager.resolve("pubmlst/ecoli", version="1.2.0")

# Fetch torch to local path
path = manager.fetch_torch("pubmlst/ecoli")
```

## Contributing

Contributions welcome! Please:

1. Fork the repository
2. Create a feature branch
3. Add tests for new functionality
4. Ensure `make test` and `make lint` pass
5. Submit a pull request

## License

See LICENSE file for details.

## Citation

If you use Torchbase in your research, please cite:

```
[Citation information to be added]
```

## Support

- GitHub Issues: https://github.com/CFSAN-Biostatistics/torchbase/issues
- Documentation: [Coming soon]

---

Created from Binfie-cookiecutter, https://github.com/crashfrog/binfie-cookiecutter