# Building Torches

Torches are built with `torchtools`, the authoring companion to `torchbase`.
All build commands live under `torchtools convert`; two utility commands
(`build` and `version`) handle validation and versioning.

```
torchtools --help
torchtools convert --help
```

---

## Common options

Every `torchtools convert` subcommand accepts these quality-analysis options:

| Option | Default | Description |
|---|---|---|
| `--kmer-size` | `13` | K-mer size used in pairwise allele similarity analysis |
| `--overlap-threshold` | `0.90` | Similarity above this value flags allele pairs as overlapping |
| `--duplicate-threshold` | `0.95` | Similarity above this value flags allele pairs as duplicates |
| `--output` | `.` | Directory under which the torch is written (`<output>/<namespace>/<name>/<version>.torch/`) |

---

## PubMLST schemes

### `torchtools convert-pubmlst` — fetch from API

Downloads a complete MLST scheme from a BIGSdb-compatible REST API, fetches
all allele FASTA sequences, runs k-mer quality analysis, and writes a torch.

```
torchtools convert-pubmlst \
    --url https://pubmlst.org/api \
    --scheme-id 1 \
    --output ./torches
```

| Option | Required | Description |
|---|---|---|
| `--url` | yes | Base URL of the BIGSdb database API |
| `--scheme-id` | yes | Numeric scheme ID |
| `--output` | yes | Output directory |

### `torchtools convert pubmlst` — from local files (legacy)

Accepts a locally-exported scheme file and pre-downloaded allele FASTA files.

```
torchtools convert pubmlst scheme.json locus1.fasta locus2.fasta ...
```

### `torchtools convert pubcgmlst` — cgMLST from local files

Same structure as the MLST local converter but sets namespace to `pubcgmlst`
and is suited to cgMLST schemes (hundreds of loci).

```
torchtools convert pubcgmlst \
    --output ./torches \
    --namespace pubcgmlst \
    --name ecoli_cgmlst \
    cgst_profiles.tsv \
    locus1.fasta locus2.fasta ...
```

| Argument | Description |
|---|---|
| `SCHEME` | Profiles TSV (cgST/ID column + per-locus allele number columns) |
| `SEQUENCES` | Per-locus FASTA files (one file per locus) |

Extra options: `--namespace`, `--name` (defaults to scheme filename stem).

---

## Salmonella

### `torchtools convert seqsero2` — SeqSero2 (UGA Center for Food Safety / Deng lab)

Serotypes *Salmonella* using wzx/wzy (O-antigen), fliC (H1), and fljB (H2)
antigen loci. Corresponds to the [UGA Center for Food Safety / Deng lab](https://github.com/denglab/SeqSero2)
database.

```
torchtools convert seqsero2 \
    --profiles serotypes.tsv \
    --output ./torches \
    wzx_wzy.fasta fliC.fasta fljB.fasta
```

**Profiles TSV columns:** `Serotype`, `O`, `H1`, `H2`

Antigen values follow White-Kauffmann-Le Minor notation; use `-` for absent
and `?` for ambiguous antigens.

If `--profiles` is omitted a stub header-only table is written; the torch is
loadable but cannot resolve serotypes until a table is supplied.

### `torchtools convert seqsero2s` — SeqSero2S (LSTUGA)

A curated Salmonella serotyping database from
[LSTUGA/SeqSero2S](https://github.com/LSTUGA/SeqSero2S) that combines a
corrected antigen allele database with 7-gene MLST sequence typing
(aroC, dnaN, hemD, hisD, purE, sucA, thrA).

```
torchtools convert seqsero2s \
    --antigen-db H_and_O_and_specific_genes.fasta \
    --mlst-profiles salmonella_profile.txt \
    --serotype-profiles serotypes.tsv \
    --output ./torches \
    aroC.tfa dnaN.tfa hemD.tfa hisD.tfa purE.tfa sucA.tfa thrA.tfa
```

| Option/Argument | Description |
|---|---|
| `[MLST_FASTAS]...` | Per-locus FASTA files for the 7 MLST loci (any extension; renamed to `.fasta`) |
| `--antigen-db` | Combined antigen FASTA (`H_and_O_and_specific_genes.fasta` from the repo) |
| `--mlst-profiles` | MLST allele-to-ST table (`salmonella_profile.txt`) |
| `--serotype-profiles` | Serotype definitions TSV — columns: `Serotype`, `O`, `H1`, `H2`, `ST` |

All three `--*-profiles` flags are optional; the torch is loadable without
them but serotype and ST resolution will be incomplete.

---

## Escherichia coli / Shigella

### `torchtools convert ectyper` — ECTyper (phac-nml)

Serotypes *E. coli* and *Shigella* by O-antigen and H-antigen loci.
Corresponds to the [phac-nml/ectyper](https://github.com/phac-nml/ectyper)
database.

**Option A — combined database FASTA** (recommended; matches ECTyper's
`ECTyperDB.fasta` directly). Headers beginning with `O` are written to
`O_antigen.fasta`; headers beginning with `H` to `H_antigen.fasta`.

```
torchtools convert ectyper \
    --db ECTyperDB.fasta \
    --profiles ectyper.tsv \
    --output ./torches
```

**Option B — pre-split FASTA files**

```
torchtools convert ectyper \
    --profiles ectyper.tsv \
    --output ./torches \
    O_antigen.fasta H_antigen.fasta
```

**Profiles TSV columns:** `Serotype`, `O`, `H`

---

## Shigella

### `torchtools convert shigatyper` — ShigaTyper

Serotypes *Shigella* using wzx/wzy (O-antigen) and fliC (H-antigen) loci.
Corresponds to the [CFSAN/ShigaTyper](https://github.com/CFSAN-Biostatistics/ShigaTyper)
database.

```
torchtools convert shigatyper \
    --profiles serotypes.tsv \
    --output ./torches \
    wzx.fasta wzy.fasta fliC.fasta
```

**Profiles TSV columns:** `Serotype`, `O`, `H`

---

## Listeria monocytogenes

### `torchtools convert lissero` — LisSero (MDU-PHL)

Determines *Listeria monocytogenes* serogroup by presence/absence of eight
PCR-target loci: `prs`, `LMOSA`, `LMOSB`, `ORF2110`, `ORF2819`, `ldh`,
`lin0764`, `lin1118`. Corresponds to the
[MDU-PHL/LisSero](https://github.com/MDU-PHL/LisSero) database.

```
torchtools convert lissero \
    --profiles serogroups.tsv \
    --output ./torches \
    prs.fasta LMOSA.fasta LMOSB.fasta ORF2110.fasta ORF2819.fasta \
    ldh.fasta lin0764.fasta lin1118.fasta
```

**Profiles TSV columns:** `Serogroup`, `prs`, `LMOSA`, `LMOSB`, `ORF2110`,
`ORF2819`, `ldh`, `lin0764`, `lin1118`

Locus values in the profiles table should be `1` (present) or `0` (absent).

---

## Chewie-NS (wgMLST)

### `torchtools convert chewie-ns`

Placeholder for [Chewie-NS](https://chewbbaca.online/) wgMLST scheme
conversion. Not yet implemented.

---

## Torch utilities

### `torchtools build` — validate a torch

Loads a torch from disk and reports its contents. Useful for verifying a
newly-built torch before distributing it.

```
torchtools build ./torches/pubcgmlst/ecoli_cgmlst/1.0.0.torch
```

### `torchtools version` — bump the version

Increments the semantic version recorded in a torch's `metadata.toml`.

```
torchtools version ./torches/pubcgmlst/ecoli_cgmlst/1.0.0.torch --bump minor
# writes 1.1.0 into metadata.toml
```

| Option | Default | Description |
|---|---|---|
| `--bump` | `patch` | Version component to increment: `patch`, `minor`, or `major` |

The optional `CHECKPOINT` argument pins the torch to a specific IPFS CID after
a successful version bump.

---

## Cryptographic signing

Torches you build and distribute can be cryptographically signed so consumers
can verify they came from you. Signatures cover the full file content of the
torch and, separately, the IPFS CID recorded in the manifest. Two key backends
are supported: a software Ed25519 key file, or a YubiKey PIV slot (key never
leaves hardware).

### Key generation

**Software key (default):**

```
torchtools keygen --namespace us-fda-hfp
```

Writes `~/.torchbase/keys/us-fda-hfp.key` (mode 0600) and `~/.torchbase/keys/us-fda-hfp.pub`.
The public key is printed at the end — add it to your key registry TOML.

**YubiKey PIV:**

```
torchtools keygen --namespace us-fda-hfp --yubikey [--slot 9c] [--pin <pin>]
```

Generates a key on the PIV slot (default `9c`) and exports the public key to
`~/.torchbase/keys/us-fda-hfp.pub`. The private key never leaves the device.

### Signing a torch

```
torchtools sign ./torches/us-fda-hfp/seqsero2/1.0.0.torch
torchtools sign ./torches/us-fda-hfp/seqsero2/1.0.0.torch --yubikey [--slot 9c] [--pin <pin>]
```

Writes `signature.toml` inside the torch directory. Signing is idempotent —
re-running produces a new signature over the same content.

### Publishing (sign + IPFS add + sign CID)

```
torchtools publish ./torches/us-fda-hfp/seqsero2/1.0.0.torch
torchtools publish ./torches/us-fda-hfp/seqsero2/1.0.0.torch --yubikey
```

Chains three steps: signs the torch content, uploads to IPFS, then signs the
resulting CID. Prints a manifest snippet ready to paste into your registry:

```toml
["us-fda-hfp/seqsero2"]
"1.0.0" = "QmABC..."
latest   = "QmABC..."

["us-fda-hfp/seqsero2".signatures]
"1.0.0" = "base64url-encoded-signature"
```

### Verifying a torch

```
torchbase verify ./torches/us-fda-hfp/seqsero2/1.0.0.torch
torchbase verify ./torches/us-fda-hfp/seqsero2/1.0.0.torch --public-key <base64url-key>
torchbase verify ./torches/us-fda-hfp/seqsero2/1.0.0.torch --require-signature
```

Recomputes the content hash, checks it matches `signature.toml`, and verifies
the Ed25519 (or P-256) signature. Key lookup order:

1. `--public-key` flag (explicit override)
2. Key registries listed in `~/.torchbase/config.toml`
3. `trusted_keys` in config
4. Public key embedded in `signature.toml`

Without `--require-signature`, a missing or unknown-namespace signature prints
a warning and exits 0. With `--require-signature`, missing signatures are
treated as failures.

### Key registry

Publish a TOML file (via IPFS/IPNS or HTTP) mapping namespaces to public keys:

```toml
[keys]
us-fda-hfp = "base64url-pubkey..."
pubmlst    = "base64url-pubkey..."
```

Reference it in `~/.torchbase/config.toml`:

```toml
key_registries = ["https://example.com/torchbase-keys.toml"]
# or
key_registries = ["/ipns/k51q...yourname"]

key_cache_ttl_hours = 24   # optional, default 24

[trusted_keys]
us-fda-hfp = "base64url-pubkey..."  # inline fallback, no network needed
```

The registry is fetched once and cached at `~/.torchbase/key_cache.toml`;
re-fetched after `key_cache_ttl_hours` hours.

---

## Publishing torches with IPLD commit chains

`torchtools manifest add` is the complete publishing workflow — it signs, uploads,
and appends a tamper-evident update block to your namespace's IPLD commit chain.
Every namespace has its own append-only history; altering any past block breaks
every forward CID link.

### One-time: generate a key and register the namespace

```
# 1. Generate an Ed25519 key for your namespace
torchtools keygen --namespace us-fda-hfp

# 2. Register the namespace (creates genesis block, publishes to IPNS)
torchtools namespace register --namespace us-fda-hfp
```

`namespace register` does the following:

1. Imports your Ed25519 PEM key into Kubo's keystore as the key named `us-fda-hfp`
2. Creates and signs a genesis block (namespace claim)
3. Uploads the genesis block to IPFS and pins it
4. Publishes the genesis CID to IPNS under your key
5. Submits the genesis CID to any configured log operators
6. Writes the IPNS address to `~/.torchbase/config.toml`:

```toml
[namespaces]
"us-fda-hfp" = "/ipns/k51q..."
```

Share that config snippet with collaborators so they can resolve your torches.

### Publishing a torch

```
torchtools manifest add ./torches/us-fda-hfp/seqsero2/2.0.0.torch
```

Steps executed:

1. Sign the torch (writes/overwrites `signature.toml`)
2. Upload the torch directory to IPFS via multipart upload; pin the CID
3. Sign the CID: `"{namespace}:{version}:{cid}"` with your namespace key
4. Resolve the current IPNS head (the previous block CID)
5. Build and upload an update block:

```toml
type      = "update"
namespace = "us-fda-hfp"
previous  = "<previous block CID>"
timestamp = "2026-06-07T12:00:00+00:00"
signature = "<base64url>"

["us-fda-hfp/seqsero2"]
"2.0.0" = "<torch CID>"
latest   = "<torch CID>"

["us-fda-hfp/seqsero2".signatures]
"2.0.0" = "<base64url CID signature>"
```

6. Publish the update block CID to IPNS (chain head advances)
7. Optionally submit the genesis CID to configured log operators

Output includes the torch CID, block CID, IPNS address, and the equivalent
manifest TOML entry for reference.

### YubiKey signing

```
torchtools namespace register --namespace us-fda-hfp --yubikey --slot 9c
torchtools manifest add ./torches/us-fda-hfp/seqsero2/2.0.0.torch --yubikey --slot 9c
```

### Inspecting a namespace chain

```
# Print reconstructed manifest (walks chain, verifies all sigs)
torchtools manifest show us-fda-hfp

# More detail including block count and genesis public key
torchtools namespace show us-fda-hfp
```

Both commands fail loudly if any block signature is invalid or if the chain has
a broken previous-CID link.

### Client resolution

Once the IPNS address is in a collaborator's `~/.torchbase/config.toml`,
`torchbase pull` resolves from the chain automatically:

```
torchbase pull us-fda-hfp/seqsero2
```

Resolution order:
1. Chain-based namespace registries (`[namespaces]` in config) — checked first
2. HTTP/IPNS flat-manifest registries (`[registries]`)
3. Error if not found in any source

### CT-style log operators

Log operators record genesis CIDs permissionlessly — they don't approve
namespaces, they only log them. Multiple independent log operators can exist;
clients don't need to trust any particular one.

Configure log operators in `~/.torchbase/config.toml`:

```toml
log_operators = [
    "https://log.torchbase.org/submit",
]
```

Or pass `--submit-to` on any publish command:

```
torchtools namespace register --namespace us-fda-hfp --submit-to https://log.example.org/submit
torchtools manifest add ./torches/us-fda-hfp/seqsero2/2.0.0.torch --submit-to https://log.example.org/submit
```

---

## Quality analysis

Every converter runs k-mer quality analysis on the allele FASTA files it
receives and writes a `quality.json` alongside the torch. The report flags:

- **Suspect overlapping pairs** — allele pairs with similarity ≥ `--overlap-threshold`
- **Duplicate pairs** — allele pairs with similarity ≥ `--duplicate-threshold`

These flags can later be used at typing time with `torchbase run
--exclude-suspect-alleles`, `--exclude-suspect-loci`, and
`--exclude-suspect-profiles`.

Adjust thresholds for highly conserved gene families:

```
torchtools convert seqsero2 \
    --overlap-threshold 0.95 \
    --duplicate-threshold 0.99 \
    --profiles serotypes.tsv \
    wzx_wzy.fasta fliC.fasta fljB.fasta
```
