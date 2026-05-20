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
torchbase run pubmlst/ecoli --contigs contigs.fasta
torchbase run pubmlst/ecoli --reads reads_R1.fastq.gz reads_R2.fastq.gz
```

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

- `torchbase`: User-facing commands (list, pull, info, run)
- `torchtools`: Authoring commands (build, version, convert)
- Automatic file compression to zstandard format

## Torch Package Structure

```
<namespace>/<torchname>/<version>.torch/
├── metadata.toml              # Package metadata, citations
├── <buildname>.profiles.tsv   # Allelic profile table
├── <buildname>.wdl            # Main WDL workflow
├── <torchname>.build.wdl      # Build workflow
└── _resources/                # Reference FASTA files
    ├── locus1.fasta
    └── locus2.fasta
```

## Quality Analysis

Torchbase includes k-mer analysis for quality control:

```python
from torchbase.quality.kmer_analysis import analyze_locus

report = analyze_locus(fasta_path, k=21)
print(f"Found {len(report.suspect_pairs)} suspect pairs")
```

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