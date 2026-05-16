# Torchbase Domain Language

## Core Concepts

**Database**
A collection of FASTA files containing allele sequences, organized by locus (e.g., `aroC.fasta`, `dnaN.fasta`). Each file contains all known sequence variants for that genetic locus.

**Profile table** (also: schema)
Tab-separated file mapping named combinations of alleles to sequence types (STs). Format: first column is ST identifier, subsequent columns are locus names containing allele identifiers. The combination of alleles defines the sequence type.

**Torch**
Versioned, immutable package containing:
- Database (allele sequence FASTA files)
- Profile table
- Metadata (provenance, authorship, citations)
- WDL workflow for execution
- Build instructions

Distributed via IPFS. Structure: `<namespace>/<torchname>/<version>.torch/`

**Registry**
IPNS-based index mapping torch identifiers to IPFS CIDs. The official registry lives at `/ipns/registry.torchbase.org`. Users can configure additional registries to support independent forks and alternative distributions.

**Typing**
The process of:
1. Matching query sequences against the database to identify which alleles are present (allele calling)
2. Looking up the allele combination in the profile table
3. Returning the sequence type (ST) if the combination matches

**Novel profile**
When allele calling succeeds (all or most loci match known alleles) but the specific combination doesn't exist in the profile table. Reported with nearest known ST for epidemiological context.

**Novel allele**
When a locus sequence doesn't match any known allele in the database above confidence threshold. May indicate new variant or sequencing/assembly error.

**Scheme**
A typing system for a specific organism or lineage, defining which loci are examined and which allele combinations constitute recognized sequence types. Multi-scheme torches contain multiple schemes in hierarchical structure.

## Versioning Strategies

Torches support three versioning approaches, declared in `metadata.toml`:

**Snapshot versioning**
Version string is a date (ISO 8601 format). Represents the state of the source database at that point in time. Example: `2024-03-15`

**Content-addressed versioning**
Version string is a hash of the database + profile contents. Maximally precise but opaque to humans. Example: `a3f5b891`

**Semantic versioning**
Version string follows semver (major.minor.patch). Changes classified by impact: new alleles = patch, new loci = minor, breaking changes = major. Example: `1.2.3`

All strategies must provide a `timestamp` field in metadata for cross-strategy comparison.

## Registry Architecture

**Single official registry** at IPNS name, updatable by maintainers
**Multi-registry support** allowing users to add alternative IPNS endpoints
**Namespace ownership** enforced per-registry (no cross-registry collision resolution)

Users configure registries in `~/.torchbase/registries.toml`

## Quality Assessment

**K-mer analysis**
Conversion tools use k-mer frequency tables (via Jellyfish or similar) to detect suspect alleles within each locus. For each locus, pairwise k-mer similarities are computed to identify:
- **Overlaps**: One allele's k-mer set is largely contained within another (possible subsequence)
- **Duplicates**: Two alleles share nearly identical k-mer sets (possible redundant entries)

Auto-tuning detects outliers using gap detection in the similarity distribution, falling back to 99th percentile thresholds when no clear gap exists.

**Quality report**
Generated during conversion as both human-readable text (with ASCII histograms) and structured JSON. The JSON report (`quality.json`) distributes with the torch package to provide full provenance. Contains per-locus statistics, similarity distributions, and flagged allele pairs.

**Hierarchical flagging** 
Suspect data marked at three levels in torch metadata:
- Suspect alleles (specific sequence pairs)
- Suspect loci (genes containing problematic alleles)  
- Suspect profiles (STs affected by suspect alleles)

Default behavior includes suspect alleles in typing. Users control via positive flags like `--include-suspect-alleles` (not double-negative flags).

## Future Scope

**Runtime annotations**: Mechanism for users to report quality issues discovered during typing (e.g., newly found ambiguities). Would use IPFS/IPNS for distributed annotation log separate from immutable torch packages.

**Validation datasets**: Distributed as torches containing:
- Curated reference isolates with laboratory-confirmed types
- Synthetic controls for edge cases
- Expected results for each sample
- Both raw reads and assemblies
Validation torches used to verify typing pipeline correctness (biological accuracy) rather than tool replication.
