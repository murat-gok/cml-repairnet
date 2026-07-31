# CML-RepairNet Stage 4 Parser & Evidence Fix v0.5.4

Changed files:

- `src/cml_repairnet/atlas/sources.py`
- `src/cml_repairnet/atlas/pipeline.py`
- `tests/test_stage4_parser_v054.py`
- `pyproject.toml`

The patch does not include or overwrite user data, cached UniProt/AlphaFold files,
or previous Stage 1–3 results.

## Main fixes

1. Strict site extraction rejects digits embedded in UniProt accessions.
2. Header detection prefers dense real headers over narrative preambles.
3. Cellular human/mouse overlap rows are emitted for both species.
4. Old-vs-young sheets are parsed as CML-site evidence; FDR <= 0.05 defines
   `age_associated=True`.
5. Proteostasis strings such as `Heart, Kidney (Kidney324)` produce site K324.
6. PRM peptides are stored separately and mapped to canonical UniProt sequences
   during annotation.
7. Protein names containing `PRM` no longer create false PRM validation.
8. Accession columns can no longer be reused as peptide or modification columns.
9. PRM mapping and exclusion audit files are added to the Stage 4 pipeline.
