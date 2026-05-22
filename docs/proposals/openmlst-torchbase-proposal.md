# OpenMLST: A Proposal for Torchbase Integration
## Toward a Community Fork of PubMLST with a Decentralized Future

---

## Background and Motivation

PubMLST (hosted at the University of Oxford, powered by the BIGSdb platform) has been the
canonical registry for MLST schemes since 1998 — covering over 130 bacterial species and
serving as the authoritative nomenclature for alleles and sequence types used across 25 years
of published epidemiological literature.

In January 2025, PubMLST implemented mandatory authentication for all data added after
31 December 2024, and both PubMLST and its Institut Pasteur mirror adopted an explicit
non-redistribution policy for new data. This has broken unauthenticated database update
scripts in downstream tools (most visibly Torsten Seemann's widely-used `mlst`), created
legal ambiguity around bundling schemes in open-source tools distributed via Bioconda, and
introduced friction into automated public health surveillance pipelines that depend on
programmatic access.

The pre-2025 data is still accessible and was accumulated under a culture of open access.
The allele sequences themselves — fragments of naturally occurring housekeeping genes — are
almost certainly not copyrightable as individual works in any jurisdiction. Oxford's
protectable interest lies primarily in the database right over the *curation work product*:
the quality-controlled, error-corrected, continuously updated numbered allele registry.

This creates a well-defined opportunity and a clear design challenge:

- **Opportunity**: Snapshot the pre-2025 PubMLST data into a community-stewarded resource
  with an open governance model, freeing downstream tools from registry dependency.
- **Challenge**: Any future alleles need to be identified in a way that doesn't recreate
  the centralized registry bottleneck — while remaining interoperable with the legacy
  integer-based nomenclature the literature depends on.

The solution is a **dual-namespace allele identification scheme**: legacy PubMLST integer
IDs for pre-2025 alleles, and content-based (hash) identifiers for all new alleles
discovered thereafter. Torchbase is the right vehicle for this.

---

## The Hash-Based Identification Principle

Traditional MLST allele IDs are positional artifacts: allele 672 of *abcZ* means "the 672nd
unique sequence observed, in the order of registry submission." The number carries no
information about the sequence; it is purely a registry artifact created and owned by whoever
controls the central database. This is what ties the entire community to that registry.

Content-based addressing inverts this. A cryptographic hash of the allele sequence *is* the
identifier. Two independent labs sequencing the same allele in isolation will compute the
same hash without ever communicating. The identifier namespace is self-organizing,
decentralized by construction, and requires no central authority to mint new IDs.

This is the approach advocated by Lee Katz (CDC), independently described for *C. difficile*
cgMLST by Eyre et al. (2019), formalized in the chewieSnake workflow (Deneke et al., 2021),
and implemented as a general engine in PHAC-NML's locidex. It is an idea whose time has
clearly come, accelerated directly by the PubMLST access changes.

The tradeoff is real but manageable: hash identity is binary (exact match or no match).
There is no "nearest allele" relationship encoded. A novel allele one SNP away from a known
one is an entirely new entity. This is handled at a higher layer — the torch's WDL workflow
can implement fuzzy matching and report both the exact hash (if matched) and the nearest
legacy allele (if known), preserving the biological interpretability that the integer scheme
provides.

---

## Proposed Torchbase Features for OpenMLST Support

### 1. Dual-Namespace Allele Identification

**What**: Each allele entry in a torch's `_resources/` FASTA and profile table carries
both a legacy integer ID (where one exists) and a canonical hash ID. The hash is always
authoritative for identity; the integer is a human-readable annotation for literature
cross-referencing.

**Implementation**:

The FASTA record format in `_resources/` should be extended from:
```
>allele_number
ATGCGT...
```
to:
```
>sha256:<hash> pubmlst:<integer_id> locus:<locus_name>
ATGCGT...
```

For pre-2025 alleles imported from PubMLST, both fields are populated. For new alleles
submitted by users, only the `sha256:` field is populated; `pubmlst:` is absent or
`pubmlst:novel`.

The `Profile` dataclass should gain a `allele_id_t` union type:
```python
@dataclass
class AlleleId:
    hash: str                  # sha256:<hex>, always present
    legacy_int: int | None     # PubMLST integer, present for pre-2025 alleles
    locus: str
```

Profile comparison (`Profile.__eq__`) should continue to work on hash identity for
new alleles, and on integer identity for legacy alleles where the hash has been
back-computed from the known sequence (which it can always be, since the pre-2025
sequences are accessible).

**Why**: This allows a torch built from the pre-2025 PubMLST snapshot to report results
in the legacy `ST-11 complex` / `abcZ(672)` notation that published literature uses,
while new alleles discovered and contributed by users get stable, registry-independent
identities immediately upon discovery.

---

### 2. The `pubmlst` Conversion Module (Complete Implementation)

The `torchbase/conversions/pubmlst.py` stub should be fully implemented to perform the
initial OpenMLST torch build from the pre-2025 PubMLST snapshot.

**What it needs to do**:

1. Accept a PubMLST scheme name (e.g. `neisseria`, `saureus`) and a local snapshot
   directory (or fetch via the still-accessible unauthenticated API for pre-2025 data).
2. Download all allele FASTA files for all loci in the scheme.
3. Back-compute SHA-256 hashes for each allele sequence and store them.
4. Download the profiles table (allele integer combos → ST integers).
5. Emit a torch package with the dual-namespace FASTA and profile table, plus populated
   `metadata.toml` with PubMLST citation, scheme version/date, and a `data_cutoff` field
   set to `2024-12-31`.

**New `torchtools` command**:
```
torchtools convert pubmlst <scheme_name> [--snapshot-dir DIR] [--cutoff-date DATE]
```

The `data_cutoff` in `metadata.toml` is important — it marks the provenance boundary
between legacy-integer-indexed data and hash-only data.

---

### 3. Novel Allele Reporting and Back-Propagation

This is the "backwards" flow described in the CLAUDE.md overview — gathering new alleles
from end users. For OpenMLST this becomes a first-class workflow.

**What**: When a user runs a typing workflow and encounters a sequence that doesn't match
any known allele exactly (currently reported as `novel` or `~nearest`), the torch workflow
should:

1. Extract the novel sequence.
2. Compute its SHA-256 hash.
3. Emit a standardized `novel_alleles.tsv` record:

```
locus     hash                        nearest_legacy_id  nearest_distance  sequence
abcZ      sha256:3f2a...              672                1 SNP             ATGCGT...
```

4. Optionally push this to a community allele repository (see below).

The novel allele immediately has a stable, portable identity (its hash) even before any
community validation. If it is later validated and added to the torch's allele database,
its hash ID remains its canonical identity forever. The integer it might eventually be
assigned in a legacy registry (if anyone still runs one) becomes merely a label, not
the identity.

**New `torchbase` command**:
```
torchbase contribute <novel_alleles.tsv> --torch <torch_ref> [--repo <git_or_ipfs_ref>]
```

---

### 4. Torch Versioning and the Allele Accumulation Model

**What**: A torch version represents a snapshot of the allele database at a point in time.
The pre-2025 PubMLST import becomes version `1.0.0` of, e.g.,
`openmlst/neisseria/1.0.0.torch`. As community-contributed alleles are validated and
incorporated, new minor versions are cut: `1.1.0`, `1.2.0`, etc.

The version changelog should track which alleles were added, as hash IDs, with optional
cross-references to the submitting institution and sequencing run accession. This is the
audit trail that replaces PubMLST's curatorial record.

**Existing torch versioning** (`torchtools version`) should be extended to support an
`--add-alleles` subcommand that ingests a validated `novel_alleles.tsv`, bumps the minor
version, and rebuilds the BLAST database in `_resources/`.

**Governance implication**: Because the version history is a git repo distributed via
IPFS, the community can inspect the full provenance of every allele, fork at any point,
and run competing "registries" that are simply forks of the same git history. There is
no single point of failure or control.

---

### 5. IPFS Distribution and the Community Manifest

IPFS is the right distribution layer for OpenMLST torches because it is content-addressed
at the file level — the IPFS CID of a torch package is a function of its contents,
which means a torch at a given CID is immutable and verifiable. This aligns perfectly
with the hash-based allele identity philosophy.

**What needs to be implemented** (currently stubbed):

- `torchbase pull openmlst/neisseria` should resolve the latest version CID from the
  community manifest and fetch it via IPFS.
- The manifest itself should be a signed, versioned document distributed via a stable
  IPFS name (IPNS) or a DNS TXT record pointing to the current CID.
- `torchtools publish` should pin the new torch to IPFS, compute the CID, and submit
  a pull request (or signed update) to the community manifest.

The key point: once a torch is published to IPFS, it requires no ongoing server
infrastructure to remain accessible. Any node that has fetched it can serve it. This
directly addresses the "sustainability is potentially limited by the funding available"
problem identified in the hash-cgMLST literature.

---

### 6. Legacy Compatibility Layer

For users who depend on PubMLST integer notation in their existing pipelines and reports,
the torch workflow output should include a `--legacy-notation` flag that renders results
in traditional format (`abcZ(672)`, `ST-11`, `CC-11`) wherever the allele was imported
from the pre-2025 snapshot and thus has a known integer ID.

For novel alleles (hash-only), the output should render as `abcZ(sha256:3f2a...)` or an
abbreviated form like `abcZ(~3f2a)` — visually distinct from integer allele IDs so users
know at a glance which namespace applies.

The `Profile` comparison logic (already flexible in handling tuples, dicts, PubMLST-style
strings) should be extended to handle hash-prefixed strings as a first-class format.

---

## Data Provenance and the `metadata.toml` Schema Extension

The existing `metadata.toml` should gain an `[openmlst]` section:

```toml
[openmlst]
source = "pubmlst"                     # or "community", "enterobase", etc.
source_scheme = "neisseria"
data_cutoff = "2024-12-31"            # date after which no legacy integers exist
legacy_namespace = "pubmlst_integer"  # how to interpret integer allele IDs
hash_algorithm = "sha256"
bootstrap_citation = """
This torch incorporates data from the PubMLST website (https://pubmlst.org/)
developed by Keith Jolley (Wellcome Open Res. 2018 Sep 24:3:124) sited at the
University of Oxford. Pre-2025 data incorporated under open access terms.
"""
```

This makes the provenance of the dual-namespace explicit and machine-readable, so
downstream tools can reason about which alleles have legacy cross-references and which
do not.

---

## What This Does Not Solve (Honest Limitations)

**New sequence types for legacy schemes**: The ST integer (the combination of allele
integers into a sequence type number) is also a registry artifact. For pre-2025 STs,
the integer is preserved. For novel allele combinations, the "ST" concept either has
to be abandoned in favor of a hash of the allele-profile tuple, or the community needs
to agree on a distributed ST assignment protocol. This is an open governance question.
One pragmatic option: define the ST as the sorted, concatenated SHA-256 of all locus
hashes in the profile — deterministic, self-organizing, no registry needed.

**Curation quality**: PubMLST's curators catch errors like contaminated alleles and
PCR artefacts. An open git-based model needs a comparable review process. The torch
versioning system supports this (alleles can be yanked in a patch version with a
documented reason) but the social infrastructure for review is out of scope here.

**Retroactive compatibility**: Labs with post-2025 PubMLST data that they obtained
via authenticated API access cannot legally contribute those integer IDs to OpenMLST.
They can contribute the sequences (which hash to stable IDs), but the mapping from
those sequences to PubMLST's post-2025 integer assignments is controlled by Oxford.
This is acceptable — the hash IDs for those sequences are interoperable without the
integer mapping.

---

## Summary of Required Torchbase Changes

| Feature | Component | Status |
|---|---|---|
| Dual-namespace `AlleleId` type | `torchbase.py` | New |
| Hash back-computation for legacy alleles | `torchbase.py` | New |
| Hash-format Profile comparison | `torchbase.py` | Extend existing |
| `metadata.toml` `[openmlst]` section | `torchfs.py` | Extend existing |
| `pubmlst` conversion (full) | `conversions/pubmlst.py` | Implement stub |
| `pubcgmlst` conversion (full) | `conversions/pubcgmlst.py` | Implement stub |
| Novel allele reporting in WDL output | `templates/torch/` | New |
| `torchbase contribute` command | `cli.py` | New |
| `torchtools version --add-alleles` | `cli.py` | Extend existing |
| IPFS pull/publish (full) | `torchfs.py` | Implement stub |
| Community manifest resolution | `torchfs.py` | New |
| `--legacy-notation` output flag | `cli.py` | New |
| Hash-format rendering in output | `cli.py` | New |
