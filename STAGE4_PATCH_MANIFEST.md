# CML-RepairNet Stage 4 Quick-Subset Fix v0.5.2

This patch fixes the Stage 4 quick-mode failure:

```text
KeyError: 'uniprot_accession'
```

## Root cause

`prepare_candidate_table()` returns `uniprot_accession_raw`. The canonical
`uniprot_accession` column is created only by the later UniProt annotation step.
Quick-mode selection was incorrectly requiring the post-annotation column during
prepare.

## Changes

- Supports `uniprot_accession_raw` before annotation.
- Continues to support `uniprot_accession` after annotation.
- Normalizes isoform suffixes only for protein counting and selection.
- Preserves original accession columns in outputs.
- Retains deterministic evidence ranking and species balancing.
- Replaces an opaque `KeyError` with explicit schema validation.
- Adds regression tests for both pipeline schemas.

Apply over Stage 4 v0.5.1.
