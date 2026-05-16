# ADR 0001: MLST Multi-Scheme Architecture

**Status:** Accepted  
**Date:** 2026-05-16  
**Context:** PubMLST MLST conversion for torchbase

## Decision

MLST torches will use hierarchical structure (S3) to support multiple schemes in a single torch package:

```
pubmlst/mlst-database/2024-03-15.torch/
├── metadata.toml
├── schemes/
│   ├── senterica/
│   │   ├── profiles.tsv
│   │   └── alleles/
│   ├── ecoli-1/
│   └── ecoli-2/
└── mlst.wdl
```

Auto-detection workflow:
1. MinHash query against all alleles across all schemes
2. Call alleles based on best matches
3. Infer scheme from allele calls
4. Look up profile (ST) in scheme's profile table
5. Fall back to alignment if ambiguous

Fallback triggers (any of):
- Top 2 allele matches within 3% similarity
- Best match < 92% similarity for any locus
- Loci match >2 different schemes
- Coverage < 80% of expected loci

Novel outcomes:
- Novel profile: known alleles, combination not in profile table → report nearest ST
- Novel allele: sequence doesn't match any known allele → report match quality
- Both trigger quality feedback system (future scope)

For reads input: depth-filter k-mers before MinHash sketching (histogram-based with ≥3x fallback).

Alignment fallback: BLAST for assemblies, minimap2 for reads.

## Rationale

**Why hierarchical structure:** Natural organization preventing allele file collisions across organisms. Preserves legibility. Allows single-artifact distribution of entire PubMLST while supporting per-scheme execution.

**Why MinHash first:** Fast screening across 100+ schemes. Sufficient sensitivity for typical MLST typing (low allele divergence in housekeeping genes).

**Why alignment fallback:** Handles edge cases (novel alleles, ambiguous matches) with precision at cost of speed. Only invoked when needed.

**Why depth filtering for reads:** Provides error correction similar to assembly without assembly complexity. Standard approach in k-mer-based tools.

## Consequences

- `Torch` dataclass changes from `profile: Profile` to `schemes: Dict[str, Schema]`
- WDL workflow must support scheme auto-detection logic
- Conversion tool must build MinHash sketches during torch creation
- Need to integrate Jellyfish (k-mer counting) and BLAST/minimap2 (alignment)
- Thresholds configurable per-torch in metadata.toml
