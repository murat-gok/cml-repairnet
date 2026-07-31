# Scaffold validation report

Validation date: 2026-07-27

## Checks completed

- Unit tests: **4 passed**.
- Article source-data parser: **33/33** reported BSA sites parsed.
- Numbering audit: paper positions plus the 24-residue precursor offset mapped to lysine at **33/33** UniProt P02769 positions.
- Sequence feature extraction: **33 sites × 104 columns** (identifiers/window plus 101 numeric sequence-derived features).
- Structure feature extraction: passed a synthetic-PDB unit test, including residue identity, B-factor and neighborhood counts.
- Third-party raw files are excluded from the distributable scaffold.

## Important source-data observation

The published Figure 3C repair fraction and a simple `1 - mean(active)/mean(control)` reconstruction are not identical for every site. After clipping negative reconstructed values to zero, the mean absolute discrepancy across 33 sites was approximately **0.030**, with a maximum of approximately **0.196**. Therefore:

1. `repair_published` remains the primary benchmark target.
2. Replicate-derived reconstruction is retained as a diagnostic, not silently substituted for the published target.
3. The PRIDE processed/raw workflow must be rerun before confirmatory modelling to determine whether the discrepancy arises from normalization, peptide aggregation or plotting logic.

## Not completed in this runtime

The real 4F5S coordinate download could not be materialized in this restricted validation runtime. The public-data downloader is configured for RCSB and the structure parser itself passed the synthetic unit test. A full end-to-end run with 4F5S and PRIDE data remains the first local execution checkpoint.

## Scope revision

As of 2026-07-27, the project is fully in silico. Earlier wet-lab access, peptide-panel and
full-protein validation concepts have been removed. Future validation means public-data
reconstruction, external consistency, blocked generalization, uncertainty/OOD analysis and
multi-fidelity computational simulation only.
