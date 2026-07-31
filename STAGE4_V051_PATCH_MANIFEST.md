# CML-RepairNet Stage 4 Parser Correction v0.5.1

This patch replaces the Stage 4 supplementary-workbook parser.

Changes:
- Detects real Excel headers below narrative preambles.
- Supports one-row and multi-row headers.
- Uses specific alias scoring instead of first fuzzy match.
- Adds content-based fallback for UniProt accession and site columns.
- Handles additional position/accession header variants.
- Avoids misclassifying position columns as modification columns.
- Writes parse diagnostics even if zero candidates are parsed.
- Adds header row/depth and column previews to diagnostics.
- Adds three regression tests.
