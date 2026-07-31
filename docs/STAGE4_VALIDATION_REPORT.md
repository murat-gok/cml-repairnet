# Stage 4 Validation Report

## Package validation

- Full project assembly: base v0.2.0 + Stage 2 robustness + Stage 2.5 + Stage 3 v0.4.1 + Stage 4 v0.5.0
- Unit tests: **25 passed**
- Python compile check: PASS
- Synthetic end-to-end transfer-safe model training: PASS
- Synthetic candidate scoring: PASS
- Conformal interval generation: PASS
- Structural KNN OOD scoring: PASS
- Pareto sorting and diversity-constrained portfolio: PASS

## Explicitly tested safeguards

- AlphaFold pLDDT is not included in the BSA B-factor model feature vector.
- `bfactor_ca` and `neighbor_mean_bfactor_8a` are absent from `TRANSFER_FEATURES`.
- CML site coordinates are validated against the canonical UniProt sequence.
- AlphaFold fragments are selected only when the UniProt site falls inside the fragment range.
- Site residue identity must be lysine in both UniProt sequence and parsed structure.
- Candidate exclusions are retained in auditable CSV files.
- Pareto ranking is performed separately by species.

## Runtime limitation

Network-dependent downloads were not executed in the isolated build environment.
The downloader uses official public endpoints, retry/backoff, cache files, and a
manual CSV/workbook fallback. The user's WSL environment will perform the live API test.
