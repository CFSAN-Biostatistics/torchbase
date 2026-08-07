# Handoff: operon typing model, cluster session 2026-08-07 (machine 2)

**Repo:** `CFSAN-Biostatistics/torchbase`, branch `main`
**Supersedes:** `docs/handoff-2026-08-07.md` (every item in its §6 "next steps" is
now closed or explicitly reassigned below).
**State:** committed locally, **not pushed** — review and push when you're happy.

This session ran on Windows + the Reedling2 SLURM cluster (login node only, no
batch jobs submitted, nothing left running). That gave it what the previous
sandbox lacked: BLAST+, apptainer, a working `miniwdl run`, and network.

---

## 1. Headline

**Phase 0 (StxTyper parity) and Phase 1 (second scheme) of
`docs/operon-strategy-plan.md` are both complete.**

- Parity with StxTyper 1.0.45 on its entire golden test suite: **182/182
  contigs** across `test/{basic,cases,synthetics,virulence_ecoli}.fa` agree on
  *every* reported field — stx type, operon status, combined identity,
  start/stop/strand, and each subunit's reference accession, reference subtype,
  identity and coverage.
- A second scheme, ETEC's heat-labile toxin (*eltA*/*eltB*), types real data
  through the same workflow with **no algorithm change** — one new config key.

Both are pinned as offline regression tests
(`tests/test_operon_parity.py`, `tests/test_operon_etec_lt.py`), so neither
needs BLAST, Docker or network in CI.

Test suite: **698 passed / 25 skipped on Linux**, **697 passed / 26 skipped on
Windows** (the extra skip is a POSIX-file-mode assertion). Both previously-known
failures are fixed; there are no known failures on either platform.

## 2. What the workflow proving exercise found (old §5.1)

`torchbase run <stx-torch> -c <assembly>` had never actually run. It did not
work, and the failures were exactly where the old handoff predicted plus a few
worse ones. In the order they surfaced:

1. `-c/--contigs` reached miniwdl as a Python repr (`<FileReaderWithPath
   object at 0x…>`) — `str()` of the click reader, not a path. Affected the
   operon *and* both other built-in-workflow paths. Fixed with `_input_path()`
   in `cli.py`, now used by every branch.
2. `Float evalue = 1e-10` in `protein_search.wdl` rendered as `0.000000`; BLAST
   rejects it. WDL formats `Float` with `%f`. It is a `String` now.
3. tblastn was missing StxTyper's own search parameters
   (`-comp_based_stats 0 -seg no -dbsize 10000 -word_size 5 -max_target_seqs
   10000 -db_gencode 11`). Without them, hit sets and therefore calls differ.
4. **The reference set's terminal `*` was not handled.** stx.prot proteins carry
   their stop codon; untrimmed, every full-length hit reads as an internal stop
   (540 spurious `INTERNAL_STOP` calls on 7 contigs), no operon is ever
   COMPLETE, and the reported span is 3 bp too long.
5. **No locus reduction.** 160 reference proteins hit each locus, and pairing
   treated them as independent HSPs: 896 candidate operons on an assembly that
   has 7. StxTyper collapses them (`paretoBest`); we now do too.
6. Residue-table indices were shifted by one. StxTyper's A312/A318/B34 are
   already 0-based offsets into `qMap()`; the converter's `-1` shift read the
   wrong residues and mis-resolved *every* stx2 subtype. This is §9 risk 1
   happening for real — a wrong answer at full confidence, no error.
7. The converter derived a reference's class from `stx.prot`'s subclass field.
   The class is the *collapsed* type code from `famId` (2a/2c/2d → 2); the
   subclass field is the reported reference subtype. Symptom: `stxA2a` labelled
   class "2" and "2a" simultaneously.
8. The status ladder had `EXTENDED` one rung above `PARTIAL_CONTIG_END`.
9. Greedy partner choice by smallest intergenic gap picked wrong-class partners
   over better-scoring ones. Pairing now enumerates all syntenic chains and
   `select_operons` decides on identity, mirroring StxTyper's
   generate-then-select shape.
10. Loci were clustered on coordinates alone, so same-offset loci on *different
    contigs* merged and calls silently vanished (5 contigs lost their only call).

## 3. Architecture change worth knowing about

The three WDL tasks used to carry inline Python copies of `torchbase/operon.py`
"kept in lockstep". They diverged, which is how several of the bugs above
survived unit tests. The tasks now take `operon.py` as a `File` input
(`operon_module`) and import it inside their containers: **one implementation,
covered by the unit tests, executed by the workflow.** Do not reintroduce
inline copies.

## 4. Phase 1: ETEC LT (old §5.4)

`examples/etec_lt/1.0.0.torch` — hand-encoded from GenBank S60731 (LT1,
H10407) and EU113242-EU113255 (LT3-LT16; Lasaro et al. 2008, PMID 18223074).
Verified end to end on the cluster:

| Assembly | Call |
|---|---|
| ETEC H10407 complete genome (FN649414-FN649418, 5.4 Mb) | `LT1 COMPLETE`, 99.74% identity, on plasmid FN649417 |
| Real LT2-signature *elt* operon (EU113255) padded into a synthetic contig | `LT2 COMPLETE`, 100% identity |

Findings:

- **The generalization holds.** One new config key was needed:
  `intergenic_min` (default 0). *eltA* and *eltB* overlap, and a floor of zero
  rejects every real LT operon. StxTyper hard-codes non-overlap; that is an stx
  fact, not an operon fact. No pairing/scoring/calling logic changed.
- **The LT literature's residue numbering is mature-protein 1-based**, verified
  against the sequences: only that convention puts S/G/K/S at A190/196/213/224
  and T at B75 in LT1, and only it makes the four LT2-producing strains' records
  read L/D/E/T + A. In precursor offsets (what the config uses) that is
  A207/A213/A230/A241 and B95.
- The torch discriminates LT1 vs LT2 only; the other 14 types differ at
  positions its table does not read and come out as class-level `LTI`. It is a
  proof and an example, not a curated production scheme.

## 5. Also closed from the old handoff

- **§2, the two pre-existing failures.** Root cause was ordering, not
  `get_unified_files()` itself: `cli.py` materialized concatenated allele and
  profile temp files *before* validating the invocation, so tests that mocked
  only the guard-relevant attributes blew up first. Generation now happens after
  workflow resolution and is skipped entirely for operon torches (which never
  used those files). The operon config temp dir it leaked is cleaned up too.
- **§4, `uv.lock`.** The file did not transfer to this machine, and nothing in
  the repo is uv-managed (`pyproject.toml` + pip, per `CLAUDE.md`). Decision:
  `.gitignore`d as a local side effect of `uv pip install`. Revisit only if the
  project actually adopts uv.
- **§3, the environment quirk.** Not a factor here: Python 3.12, `miniwdl` 1.8.0
  works. On the cluster, `miniwdl` runs against apptainer 1.5.2 with
  `MINIWDL__SCHEDULER__CONTAINER_BACKEND=singularity`; `ncbi/blast:2.16.0` pulls
  and runs fine.
- **§5.2, frameshift detection.** Implemented, not a stub. Co-linear HSPs of one
  reference in different reading frames are stitched into a single alignment
  flagged `frameshift=True`; the stitched alignment is offered *alongside* its
  parts so locus reduction decides whether the frameshift reading wins — which
  is what keeps StxTyper's terminal-frame-change cases reported as `EXTENDED`
  rather than `FRAMESHIFT`. Validated on `stx2_fs` and the `cases`/`synthetics`
  terminal cases.

## 6. Windows portability (new, incidental)

The suite could not run on this machine at all: 74 failures. Causes were
`csv` writers/readers opened without `newline=""` (CRLF → blank rows →
`Profile.parse` crash) plus two genuine bugs where a Windows path leaked into a
platform-neutral string: WDL input JSON (`workflow_params.py`) and an IPFS
unixfs upload filename (`chain.py`). All fixed; one test asserting POSIX mode
bits is skipped on Windows (`os.chmod` in `signing.py` is unchanged and remains
the right call). `zstandard` was pinned `~=0.20.0`, which has no wheels for
Python ≥3.12 and needs MSVC to build; relaxed to `>=0.20,<1.0`.

## 7. Suggested next steps

1. **Push.** Commits are local only.
2. **Phase 2, `torchtools derive-operon`** — §7 gated it on two working schemes;
   both now exist, so it is unblocked and is the next real product step.
3. **Broader stx concordance against MicroBIGG-E.** Parity is proven on curated
   edge cases; the bulk ground truth is the wider check, and it needs a decision
   about BigQuery/GCS access from wherever it runs.
4. Optional hardening, in rough value order:
   - Reference the frameshift junction rule against StxTyper's `Hsp::Merge`
     rather than its output. The current rule (drop the codon spanning the frame
     change) reproduces its alignment length and identity exactly on the
     fixtures, but it was inferred.
   - Decide whether `profiles.tsv` still earns its place for operon torches
     (plan §9 risk 4). It is now load-time validation *and* the source of the
     subtype-label prefix, so it does earn it — worth confirming that is the
     intent rather than an accident.
   - `EXTENDED` is only reachable when the reference carries a terminal `*`;
     schemes whose references lack one can never report it. Fine for stx and LT,
     worth a note if a third scheme appears.

## 8. Reproducing the cluster verification

Nothing is left on the cluster (no jobs, no files). To redo it:

```bash
ssh -K Justin.Payne@Reedling2.fda.gov
module load python/3.12.5 blast/2.16.0
git clone <torchbase> ~/torchbase && cd ~/torchbase
python3 -m venv .venv && .venv/bin/pip install -e . pytest hypothesis
export MINIWDL__SCHEDULER__CONTAINER_BACKEND=singularity
export MINIWDL__SINGULARITY__IMAGE_CACHE=~/sif
git clone --depth 1 https://github.com/ncbi/stxtyper ~/stxtyper
.venv/bin/torchtools convert stxtyper --download --output ~/torches
.venv/bin/torchbase run ~/torches/ncbi/stxtyper/*.torch -c ~/stxtyper/test/basic.fa
```

Compare the `operon_calls.json` against `~/stxtyper/test/basic.expected`; the
same comparison for `basic.fa` runs offline as `tests/test_operon_parity.py`.
