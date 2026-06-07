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

### `torchtools convert seqsero2` — SeqSero2 (CDCgov)

Serotypes *Salmonella* using wzx/wzy (O-antigen), fliC (H1), and fljB (H2)
antigen loci. Corresponds to the [CDCgov/SeqSero2](https://github.com/CDCgov/SeqSero2)
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
torchtools keygen --namespace cdc
```

Writes `~/.torchbase/keys/cdc.key` (mode 0600) and `~/.torchbase/keys/cdc.pub`.
The public key is printed at the end — add it to your key registry TOML.

**YubiKey PIV:**

```
torchtools keygen --namespace cdc --yubikey [--slot 9c] [--pin <pin>]
```

Generates a key on the PIV slot (default `9c`) and exports the public key to
`~/.torchbase/keys/cdc.pub`. The private key never leaves the device.

### Signing a torch

```
torchtools sign ./torches/cdc/seqsero2/1.0.0.torch
torchtools sign ./torches/cdc/seqsero2/1.0.0.torch --yubikey [--slot 9c] [--pin <pin>]
```

Writes `signature.toml` inside the torch directory. Signing is idempotent —
re-running produces a new signature over the same content.

### Publishing (sign + IPFS add + sign CID)

```
torchtools publish ./torches/cdc/seqsero2/1.0.0.torch
torchtools publish ./torches/cdc/seqsero2/1.0.0.torch --yubikey
```

Chains three steps: signs the torch content, uploads to IPFS, then signs the
resulting CID. Prints a manifest snippet ready to paste into your registry:

```toml
["cdc/seqsero2"]
"1.0.0" = "QmABC..."
latest   = "QmABC..."

["cdc/seqsero2".signatures]
"1.0.0" = "base64url-encoded-signature"
```

### Verifying a torch

```
torchbase verify ./torches/cdc/seqsero2/1.0.0.torch
torchbase verify ./torches/cdc/seqsero2/1.0.0.torch --public-key <base64url-key>
torchbase verify ./torches/cdc/seqsero2/1.0.0.torch --require-signature
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
cdc      = "base64url-pubkey..."
pubmlst  = "base64url-pubkey..."
```

Reference it in `~/.torchbase/config.toml`:

```toml
key_registries = ["https://example.com/torchbase-keys.toml"]
# or
key_registries = ["/ipns/k51q...yourname"]

key_cache_ttl_hours = 24   # optional, default 24

[trusted_keys]
cdc = "base64url-pubkey..."  # inline fallback, no network needed
```

The registry is fetched once and cached at `~/.torchbase/key_cache.toml`;
re-fetched after `key_cache_ttl_hours` hours.

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
