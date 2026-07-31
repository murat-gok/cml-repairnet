# Stage 4 v0.5.3 — UniProt Accession Validation Fix

## Root cause

The Stage 4 parser used an over-permissive 10-character expression (`[A-Z0-9]{10}`).
This allowed uppercase protein-name words such as `ACTIVATING`, `ASSOCIATED`, and
`CYTOCHROME` to enter the UniProt REST query as if they were accessions.

## Changes

- Implements the official UniProtKB 6/10-character accession pattern.
- Rejects ordinary uppercase words and malformed identifiers.
- Selects the accession column primarily from valid identifier content rather than a
  broad `Protein` header match.
- Filters invalid identifiers before every UniProt API request.
- Writes rejected rows to `results/stage4/stage4_excluded_sequence.csv` with reason
  `invalid_uniprot_accession_format`.
- Adds clearer UniProt HTTP errors containing the validated request chunk.
- Adds three regression tests.

## Required action

After applying the patch, delete the previous candidate table and rerun:

```bash
python scripts/run_stage4.py prepare --mode quick
python scripts/run_stage4.py annotate
```
