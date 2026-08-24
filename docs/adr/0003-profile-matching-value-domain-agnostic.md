# ADR 0003: Profile Matching Is Value-Domain-Agnostic by Design

**Status:** Accepted (retroactive — documents rationale for shipped design)
**Date:** 2026-08-12
**Context:** Exploring whether ShigaTyper's serotyping algorithm generalizes into a
new torchbase typing model, or requires a torch-embedded override workflow.

## Decision

`profile_match.match_profile()` (and the `Profile` class it descends from) matches
a table row against a query by **string equality per column, with a wildcard
escape** — nothing in the matcher inspects, parses, or assumes numeric structure
in a locus's value. A locus's value is an opaque token: `"4"` at `adk` and
`"present"` at `wzx` are the same kind of thing to the matcher.

This was a deliberate choice, not an accident of using `csv`/`str` throughout: it
means **the matcher's domain is columns-of-tokens against rows-of-tokens**, not
"columns of allele IDs." Anything that reduces to that shape — including
presence/absence marker panels — is matchable with the code as it already
stands, no new mechanism, no new typing model.

## Rationale

The question that produced this ADR: does ShigaTyper's typing algorithm (detect
a panel of marker genes by presence/absence, then look up the combination
against a table, with priority-ordered exceptions and a few gating/exclusion
rules) need a new "presence/absence" typing model alongside `allelic` and
`operon`? Checked against a second real upstream tool (LisSero, whose entire
algorithm is BLAST hit/no-hit at a coverage+identity threshold against five
genes, then an `if/elif` cascade over set membership to a serotype), the answer
is no — with two conditions:

1. Every branch of a presence/absence cascade of the form `if X present and Y
   absent and Z absent: call = ...` reduces to one `profiles.tsv` row: the
   columns that are checked get an explicit `"present"`/`"absent"` value, the
   columns that are not checked get the wildcard. `match_profile`'s "first
   consistent row wins" is already the cascade's priority order — the
   *transcription* is mechanical (the same category of work `operon.py`'s
   residue tables do for StxTyper), not fabrication, precisely because the
   matcher never cared whether a token came from an allele-identity call or a
   presence/absence call.
2. What *isn't* already generalized is upstream of the matcher, in **calling**:
   `allele_calls.py` only knows how to pick "best similarity among
   alternatives." A presence/absence locus needs a calling mode that reports
   hit/no-hit against a threshold, and — because a marker panel can have
   mutually exclusive members (only one `wzx` variant should ever be present)
   — an explicit ambiguous/no-call status when that invariant is violated or
   when there is no evidence at all. Both are additions to the calling step,
   not to the matcher, and not a new typing model.

The candidate alternative — a third typing model generalizing "presence/absence
marker panel with rule table" — was rejected once the matcher's actual contract
was re-examined: it would have duplicated `match_profile`'s row-matching logic
for no reason, because that logic was never allele-ID-specific in the first
place.

## Consequences

- `hfp/shigatyper`'s reference set needs re-converting: it was split into 95
  single-sequence pseudo-loci (one file per marker *variant*, e.g.
  `Sb11_wzx.fasta`, `Sb12_wzx.fasta`, ...) instead of the correct shape — `wzx`
  as one locus with ~30 alleles, `wzy` likewise, and each accessory gene as its
  own small locus. Calling then picks the best allele per locus exactly as MLST
  does; presence/absence loci (single-reference genes like `cadA`, `ipaB`) call
  hit/no-hit instead.
- `profiles.tsv`'s ID column is hardcoded to `ST` (any casing) in
  `load_profiles`/`st_column`. ECTyper and SeqSero2 (already-present converters)
  use `Serotype`, so their torches report `novel_profile` unconditionally today.
  Making the ID column configurable is a prerequisite this ADR's reasoning
  shares with those two converters, not a ShigaTyper-specific fix.
- Any future typing scheme whose calling primitive is "detect discrete
  observed states, then match against a rule table with fallback" — not just
  presence/absence, but e.g. a coverage/depth tier, or a small enum of variant
  classes — is a calling-step addition, evaluated the same way: does it produce
  a token per locus, and does the existing table-matching contract already
  consume it? Reach for a new typing model only when the answer to the second
  question is genuinely no (as it was for `operon`, which needs synteny,
  intergenic distance, and combined-identity thresholds the matcher has no
  concept of).
