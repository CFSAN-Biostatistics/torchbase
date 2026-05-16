# ADR 0002: Migration Wizard Design

**Status:** Accepted  
**Date:** 2026-05-16  
**Context:** Creating reusable skill for migrating arbitrary typing tools to torchbase

## Decision

Build an interactive migration wizard (SK-C approach) that automates mechanical extraction while asking domain questions for biological/logical interpretation.

## Wizard Phases

### Phase 1: Source Identification
- Repository URL or local path
- Published papers for validation

### Phase 2: Database Discovery
- Scan for common patterns (*.fasta, *_db/, *_profiles.tsv)
- Present findings, ask for confirmation
- Identify additional/missing files

### Phase 3: Schema Interpretation (most critical)
- Parse sequence organization (per-locus / merged / by-marker)
- Attempt to understand non-obvious files (pickles, custom formats)
- **Read typing logic code** and present interpretation for confirmation
- Map allele combinations to types/serovars
- Explode branching/conditional logic into profile tables with wildcards

### Phase 4: Workflow Decisions
- Identify external tools (BWA, BLAST, etc.)
- Propose replacements with torch defaults (MinHash+alignment)
- Determine if assembly needed or if depth-filtered k-mers sufficient
- Flag domain-specific steps for preservation

### Phase 5: Quality Parameters
- Suggest k-mer sizes based on gene lengths
- Propose match thresholds from code analysis
- Define ambiguity handling

### Phase 6: Metadata
- Namespace (author/org)
- Versioning strategy
- Citations and provenance

## Question Format

- Multiple choice where possible
- Free text for explanations
- "Skip for now" with ability to revisit
- "Show me the file/code" action to inspect source
- Present skill's interpretation before asking for confirmation

## Emphasis

Phase 3 (schema interpretation) receives most effort. Skill must:
1. Parse typing logic from source code
2. Present clear explanation of discovered logic
3. Ask for confirmation/correction
4. Generate profile table or typing rules

Other phases more automated with human checkpoints.

## Rationale

**Why interactive over full automation:** Typing logic often contains domain knowledge not obvious from code alone. Human expertise necessary for interpreting biological significance.

**Why code analysis first:** Skill should do the work of reading/parsing, not burden user with "go read the code." Present findings, let user validate.

**Why Phase 3 emphasis:** Database extraction and tool replacement are mechanical. Understanding "how does allele X + allele Y → serotype Z" requires domain interpretation.

## Deliverables

Wizard produces:
```
migration_output/
├── myorg/toolname-torch/1.0.0.torch/
│   ├── metadata.toml
│   ├── MIGRATION.md  # conversion decisions, interpretations, warnings
│   ├── main.wdl (if custom workflow needed)
│   ├── schemes/
│   └── quality.json
```

MIGRATION.md documents:
- Source tool URL and version
- Wizard's interpretation of typing logic
- Tool replacements (BWA→MinHash, etc.)
- Known limitations or deviations
- Suggested validation approach

## Scope

**Initial implementation:** One-shot conversion producing working torch.

**Future enhancements:**
- Update tracking (re-run on source changes)
- Automated validation dataset discovery
- Test scaffold generation
- Direct registry submission

## Consequences

- Skill must parse Python, shell scripts, WDL, and other common bioinformatics languages
- Need UI for multi-step wizard (CLI with rich formatting? TUI? Web interface?)
- Generated torches require validation testing before publication
- Documentation should include skill's interpretation as comments in profile tables
- Validation focuses on biological correctness, not replication of source tool behavior
