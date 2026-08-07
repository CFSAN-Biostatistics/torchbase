# Plan: operon typing model for Torchbase

**Status:** implemented (v1, Phase 0 scaffold — see §11)
**Date:** 2026-08-07
**Motivation:** Generalize NCBI StxTyper's protein-space, synteny-aware operon typing into a first-class Torchbase capability, so that multi-subunit virulence systems become *configuration* rather than bespoke tools.

---

## 1. What we're building

Support for typing schemes whose call unit is a **multi-subunit operon** rather than an independent set of loci. The call depends on:

- which subunits are present, in what order, on the same contig and strand
- how far apart they are (intergenic distance)
- their **combined** identity against a reference set, thresholded per class
- optionally, specific residues at reference-anchored coordinates
- whether the operon is intact (frameshift / internal stop / truncation / extension)

Validation target for v1: reproduce StxTyper's calls on stx from a Torchbase torch.

Background analysis (StxTyper algorithm, virulence-system survey, data-source landscape) is in the session handoff at `/tmp/torchbase-operon-strategy-handoff.md`. Not duplicated here.

---

## 2. The core architectural decision

**`operon` must not be a fourth `--strategy` value.**

The existing strategies are *speed/accuracy tiers over one typing model*: `fast` (MinHash), `balanced` (MinHash + conditional alignment), `sensitive` (always align). They are interchangeable — any allelic torch runs under any of them and produces the same kind of answer at different cost.

`operon` is not interchangeable. An operon torch cannot be typed by `fast`; an MLST torch cannot be typed by `operon`. It is a different **typing model**, not a different speed tier. Adding it to `strategy_to_workflow` (`torchbase/cli.py:625`) would conflate two orthogonal axes:

| Axis | Question | Belongs to |
|---|---|---|
| **typing model** | what kind of scheme is this? | the **torch** (`metadata.toml`) |
| **strategy** | how much compute do I spend? | the **invocation** (`--strategy`) |

### Decision

Introduce `typing_model` in `metadata.toml`, defaulting to `"allelic"` when absent (backward compatible — every existing torch keeps working untouched). The CLI dispatches on the pair `(typing_model, strategy)`.

For v1 the operon model has exactly one implementation, so `--strategy` is rejected with a clear error for operon torches — the same way it is already rejected for embedded-workflow torches (`cli.py:575`). But the axis is separated correctly, so `fast`/`sensitive` operon variants (e.g. diamond vs tblastn, or Phraya protein mode) can be added later without a breaking change.

This also settles a question raised in the handoff: **operon torches bypass `Profile.__eq__` matching entirely.** `Profile` is an order-free exact/wildcard matcher with no notion of adjacency, spacing, or thresholds, and it should not be contorted to acquire one.

---

## 3. Torch layout

```
<namespace>/<name>/<version>.torch/
├── metadata.toml           # typing_model = "operon" + [operon] config block
├── profiles.tsv            # declarative list of valid subtypes (see §3.2)
├── _resources/
│   └── subunits.faa        # protein reference set, labelled headers
```

No `main.wdl` — the whole point is that this routes to a built-in workflow.

### 3.1 `metadata.toml` operon block

Straw-man schema, using stx values:

```toml
typing_model = "operon"

[operon]
subunit_order = ["A", "B"]        # required transcription order
intergenic_max = 36               # bp
intergenic_relax_factor = 2       # for partial-operon recovery pass
require_same_strand = true
require_same_contig = true

[operon.reference]
file = "_resources/subunits.faa"
# header format (resolved, see §10 Resolutions #2): >accession|subunit_role|reference_subtype|class
# e.g. AAS07582.1|A|stxA2c|2 — subunit_role is structural ("A"/"B", matches
# subunit_order); class is the already-collapsed operon-level class token
# used for threshold lookup and cross-subunit class agreement.
header_format = "accession|subunit_role|reference_subtype|class"

[operon.identity_thresholds]
# combined identity = (nident_A + nident_B) / (len_A + len_B)
default = 0.98
"1a" = 0.983
"1c" = 0.983
"1d" = 0.983
"1e" = 0.983
"2k" = 0.985
"2l" = 0.985

[operon.generalized_classes]
# classes that identity alone cannot separate; resolved by residue table
"2" = ["2a", "2c", "2d"]

[[operon.residue_rules]]
# reference-anchored coordinates, 0-based into the reference-length
# projection. StxTyper documents these positions as 1-based (A312, A318,
# B34); the stxtyper.py converter shifts by -1 at transcription time.
class = "2"
positions = [
  { subunit = "A", index = 311 },  # A312 (1-based)
  { subunit = "A", index = 317 },  # A318 (1-based)
  { subunit = "B", index = 33  },  # B34  (1-based)
]
table = [
  { call = "2a", residues = [["F","S"], ["K","E"], ["D"]] },
  { call = "2c", residues = [["F"],     ["K","E"], ["N"]] },
  { call = "2d", residues = [["S"],     ["E"],     ["N"]] },
]
fallback = "2"   # unresolved -> generalized class, flagged novel
```

**Design notes:**

- Thresholds are keyed by *class*, not subtype, matching StxTyper (`stxtyper.cpp:894-911`).
- `residue_rules` is a list so a scheme can have several generalized classes.
- Residue alternatives are lists (`["F","S"]`) because stx2a genuinely accepts either at A312.
- Coordinates are in **reference** space, requiring an alignment-to-reference projection (StxTyper's `qMap()`). This is the piece most likely to harbour off-by-one bugs; see §7.

### 3.2 What `profiles.tsv` is for

It is **not** the matching mechanism. It is a declarative manifest of valid subtypes so that `torchbase info` and torch validation work uniformly across typing models, and so that a scheme's vocabulary is inspectable without parsing FASTA headers:

```
subtype	class	subunit_A	subunit_B
stx1a	1a	stxA1a	stxB1a
stx2a	2	stxA2	stxB2a
...
```

Validation at load time: every `class` must have a threshold (explicit or default); every generalized class must have a residue rule.

---

## 4. Workflow

New file: `torchbase/workflows/builtin/operon_typing.wdl`, with tasks under `tasks/`:

- `tasks/protein_search.wdl` — `makeblastdb -dbtype nucl` on the query assembly, then `tblastn` with the subunit reference as query. Mirrors `stxtyper.cpp:887`, including `-gapextend 2`.
- `tasks/operon_assembly.wdl` — pair HSPs into candidate operons under the syntenic constraints; greedy claim with a `reported` flag; three passes at relaxing stringency (strict threshold → perfect-only → any).
- `tasks/operon_call.wdl` — combined identity, threshold comparison, residue-table resolution, status ladder, output JSON.

**Dependency note.** v1 uses BLAST+ in a pinned container inside the WDL task. This is a deliberate interim: per the Phraya discussion, protein alignment there is greenfield (no affine gaps, no substitution matrix — the crate is unit-cost edit distance throughout). When Phraya gains affine gaps plus a matrix-scored engine, `protein_search.wdl` swaps out and nothing else changes. Keeping the search behind a single task boundary is what makes that swap cheap, so **do not leak BLAST tabular format into the downstream tasks** — normalize to an internal HSP JSON at the task boundary.

---

## 5. Output schema

The current vocabulary (`known | novel_profile | novel_allele`) is too coarse. Extend rather than overload:

```json
{
  "profile_id": "stx2a",
  "profile_type": "operon_subtype",
  "scheme": "stx",
  "status": "known",
  "operon_status": "COMPLETE",
  "confidence": 0.994,
  "operon": {
    "contig": "NODE_3", "start": 12044, "stop": 13455, "strand": "+",
    "intergenic_bp": 12,
    "combined_identity": 0.994,
    "threshold_applied": 0.98,
    "subunits": {
      "A": {"reference": "AAS07596.1", "reference_subtype": "stxA2", "identity": 0.996, "coverage": 1.0},
      "B": {"reference": "...",        "reference_subtype": "stxB2a", "identity": 0.989, "coverage": 1.0}
    },
    "residue_evidence": {"A312": "F", "A318": "K", "B34": "D"}
  },
  "method": {"typing_model": "operon", "tools": ["tblastn"]}
}
```

- `status` stays for coarse cross-model compatibility: `COMPLETE` → `known`, `COMPLETE_NOVEL`/`AMBIGUOUS` → `novel_profile`, everything else → a new `incomplete`.
- `operon_status` carries the eight-value ladder: `COMPLETE > COMPLETE_NOVEL > AMBIGUOUS > PARTIAL > PARTIAL_CONTIG_END > EXTENDED > INTERNAL_STOP > FRAMESHIFT`, used both as the reported value and as the priority order when candidate operons overlap.
- `residue_evidence` is emitted **always**, not just on success — it is what makes a novel call actionable, and it is the field a curator needs when proposing a new subtype.
- Multiple operons per assembly are normal (many STEC carry both stx1 and stx2). Output must be a **list**; this differs from the allelic model's single-profile assumption and needs checking against downstream consumers.

---

## 6. Code changes

| File | Change |
|---|---|
| `torchbase/torchfs.py` | Parse `typing_model` (default `"allelic"`) and the `[operon]` block in `Torch.load()`; validate operon config at load (§3.2). Add `Torch.typing_model` and `Torch.operon_config` fields. |
| `torchbase/cli.py:606-658` | Replace the flat `strategy_to_workflow` dict with dispatch on `(typing_model, strategy)`. Reject `--strategy` for operon torches with a clear message, mirroring the existing embedded-workflow guard at `:575`. |
| `torchbase/cli.py:670-681` | The input builder currently hardcodes `allele_fasta`/`profiles_table`. Make it typing-model-aware — a small per-model input-builder function. **This is the fiddliest change**; miniwdl needs explicit `File` inputs, so "just pass the torch dir" is not available. |
| `torchbase/workflows/builtin/operon_typing.wdl` + 3 tasks | New. §4. |
| `torchbase/conversions/stxtyper.py` | New converter: fetch `stx.prot` + version from ncbi/stxtyper, emit an operon torch with the thresholds and residue table transcribed from `METHODS.md` / source. |
| `CLAUDE.md` | Document the typing-model axis and the operon torch format. |

---

## 7. Phasing

**Phase 0 — parity.** Build the stx torch, reproduce StxTyper. Success criterion is not "high concordance" but **exact agreement on COMPLETE operons**, with every discordance individually explained. Ground truth is free and large: MicroBIGG-E already contains millions of pre-called stx subtypes (AMRFinderPlus bundles StxTyper), bulk-accessible via BigQuery/GCS.

**Phase 1 — second scheme, no new code.** ETEC LT (*eltA*/*eltB*). Structurally isomorphic: two subunits, combined identity, residue-defined variants (LT2 = A:S190L/G196D/K213E/S224T + B:T75A). If this requires *only* a new torch and zero changes to `operon_typing.wdl`, the generalization is demonstrated. If it requires code changes, the config schema in §3.1 is wrong and should be revised before going further.

**Phase 2 — `torchtools derive-operon`.** The scheme-derivation pipeline (see handoff §4). This is the genuinely novel product and the thing that makes Torchbase more than a distribution path. Deliberately *after* Phases 0–1, because it should be built against two working schemes rather than zero.

**Explicitly out of scope for v1:** tripartite+ operons (*cdtABC*, *nheABC*) — the config assumes a pair in several places and generalizing to N subunits should wait for a real driver. BoNT-style schemes where cluster *arrangement* is the type are the hard boundary and are out of scope entirely.

---

## 8. Testing

- **Unit:** operon pairing under synthetic HSP layouts — correct pair accepted; gap of 37 bp rejected; wrong strand rejected; wrong order rejected; overlapping candidates resolved by status priority.
- **Unit:** residue-table resolution, including the fallback path and every row of the stx2acd table.
- **Unit:** threshold selection, including the generalized-class lookup and the `default` fallback.
- **Golden:** the `test/` fixtures from the stxtyper repo, which cover frameshift, internal stop, truncation, extension, ambiguous bases, and novel operons — a ready-made edge-case suite, and the fastest way to find where our reimplementation diverges.
- **Property:** coordinate projection round-trips (see risk below).
- **Integration:** end-to-end `torchbase run` on an assembly containing both stx1 and stx2 operons, asserting a two-element output list.

---

## 9. Risks and open questions

1. **Coordinate projection is the sharpest edge.** Reading A[312] requires projecting the alignment back into reference coordinates through gaps. An off-by-one silently emits the wrong subtype with full confidence — a wrong answer, not an error. Mitigation: property-test the projection independently of the typing logic, and assert against StxTyper's own residue output in verbose mode (`stxtyper.cpp:554` emits `"2 " + a[312] + a[318] + b[34]` under `verboseP`, which is a directly comparable ground truth).
2. **Gap placement affects which residue you read.** Related and worse: under unit-cost/linear gaps there are many co-optimal alignments, so the residue at a fixed reference coordinate is tie-break-dependent. BLAST's affine gaps largely avoid this, which is another reason not to swap in Phraya until it has affine support.
3. **StxTyper is a moving target** — roughly monthly revisions per its changelog. Parity is a snapshot, not a standing guarantee. Pin the version we claim parity against and state it in the torch metadata.
4. **Is `profiles.tsv` pulling its weight?** (§3.2) It may be pure ceremony for operon torches. Revisit after Phase 1 — if it is never read except by `info`, consider making it optional.
5. **Multi-operon output** breaks the single-profile assumption in the allelic model. Needs an audit of downstream consumers before Phase 0 lands.
6. **`--strategy` rejection may be too strict.** If someone wants `auto` to route sensibly for a mixed torch collection, the current guard will annoy. Deferred until someone hits it.

---

## 10. Decisions needed before implementation starts

- [x] Confirm `typing_model` as the metadata key and `"allelic"` as the default (§2).
- [x] Confirm the `[operon]` config schema (§3.1), or revise after a paper exercise encoding ETEC LT in it — cheap, and it de-risks Phase 1.
- [x] Decide whether `profiles.tsv` is required or optional for operon torches (§3.2, risk 4).
- [x] Confirm interim BLAST+ dependency is acceptable given the zero-binary-dependency goal stated in Phraya's README — Torchbase already shells to miniwdl/containers, so this is likely fine, but it is a stated principle worth checking against.

### Resolutions

1. **`typing_model` metadata key, `"allelic"` default — confirmed as written in §2.** Implemented in `torchfs.Torch.load()`; every existing torch keeps working untouched (no `typing_model` key → `"allelic"`).
2. **`[operon]` config schema — confirmed as the §3.1 strawman, with one concretization the strawman left ambiguous:** the reference-header format is fixed at `accession|subunit_role|reference_subtype|class` (4 fields, not 3). The original 3-field example (`accession|subunit_tag|class_label`) conflated two different things — the structural role ("A"/"B", matched against `subunit_order`) and the specific reference subtype label (e.g. "stxA2c") — which StxTyper's own `famId` encoding only disambiguates via scheme-specific string surgery (parsing the subunit letter out of a fixed-length prefix). Making the role an explicit field is scheme-agnostic and is what the ETEC LT paper exercise (Phase 1) needs to not special-case stx's naming convention. `[[operon.residue_rules]].positions[].index` is **0-based** into the reference protein sequence — StxTyper's own A312/A318/B34 documentation is 1-based; the `stxtyper.py` converter shifts by −1 at transcription time (see its `RESIDUE_RULES` constant).
3. **`profiles.tsv` required for operon torches — confirmed for v1.** `Torch.load()` validates `[operon]` metadata against it (every `class` referenced by a row has a threshold; every generalized class has a residue rule), which needs the rows to exist. Revisit per risk 4 after Phase 1 if it turns out to be pure ceremony.
4. **Interim BLAST+ dependency — confirmed acceptable.** `protein_search.wdl` runs `tblastn` in a `docker: "ncbi/blast:2.16.0"` task; Torchbase already shells to miniwdl/containers for every built-in workflow, so this doesn't introduce a new dependency *class*, only a new container image.

---

## 11. Implementation status

**Shipped:**
- `torchbase/operon.py` — reference algorithm (threshold selection, residue-table resolution, HSP synteny pairing, status ladder), unit-tested (`tests/test_operon.py`).
- `torchbase/torchfs.py` — `typing_model`/`operon_config`/`operon_profiles` on `Torch`, `[operon]` validation at load time (`tests/test_operon_torch_loading.py`).
- `torchbase/cli.py` — dispatch on `(typing_model, strategy)`; `--strategy` rejected for operon torches; routes to `operon_typing.wdl` with explicit `subunit_reference`/`profiles_table`/`operon_config_json` File inputs (`tests/test_operon_cli_routing.py`).
- `torchbase/workflows/builtin/operon_typing.wdl` + `tasks/{protein_search,operon_assembly,operon_call}.wdl` — passes `miniwdl check`; the embedded Python in each task was smoke-tested standalone against synthetic HSP/candidate fixtures (not via `miniwdl run`, which needs Docker + BLAST+ neither available in the authoring sandbox).
- `torchbase/conversions/stxtyper.py` — `torchtools convert stxtyper [--download]`, transcribes `stx.prot` into an operon torch with the real StxTyper thresholds and stx2acd residue table (`tests/test_stxtyper_conversion.py`, synthetic fixture — no network in-sandbox).
- `CLAUDE.md` — documents the typing-model axis and operon torch/workflow layout.

**Known gaps, honestly scoped rather than silently shipped:**
- **Frameshift detection is a stub** (`protein_search.wdl` HSPs always report `frameshift: false`). Real detection needs stitching co-linear HSPs split by a frame change across a subunit/contig pair — flagged in §7 risk 1-2 as the sharpest edge, and correctly belongs in `operon_assembly.wdl` once it exists, but wasn't implemented because it can't be validated without real frameshifted assemblies.
- **Phase 0 parity (§7) is unstarted.** "Exact agreement with StxTyper on COMPLETE operons" requires running the workflow against MicroBIGG-E ground truth via BLAST+ — infrastructure unavailable in the authoring sandbox (no Docker, no BLAST+, no bulk data access). The algorithm itself (pairing, thresholds, residue resolution) is unit-tested and was smoke-tested end-to-end on synthetic data producing the correct §5 output shape, but this is not the same claim as StxTyper parity.
- **Phase 1 (ETEC LT) and Phase 2 (`derive-operon`) are unstarted**, as scoped — Phase 0 needs to land first per §7's own ordering.
