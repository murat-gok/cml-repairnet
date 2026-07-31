# CML-RepairNet Project Status

## Current release

- Version: `0.2.0`
- Scope: fully in silico
- Primary development environment: VS Code + WSL + Python 3.11

## Completed

- Modern `src/` package and `pyproject.toml` layout
- VS Code tasks, debugger profiles and extension recommendations
- Linux/WSL and Windows bootstrap scripts
- Official-data downloader with retries, atomic writes and SHA-256 manifest
- Replicate-aware extraction of 33 BSA repair sites
- UniProt precursor mapping validation
- Correct separation of mature BSA/PDB and UniProt precursor coordinates
- Transparent sequence descriptors
- Static structure descriptors
- Stage 1 scientific data audit
- Leakage-controlled baseline framework
- Multi-fidelity virtual candidate prioritization
- Unit tests and GitHub Actions CI

## Validation completed in the development environment

- Seven unit tests pass.
- Python compilation passes for `src`, `scripts` and `tests`.
- The official article source workbook yields 33 labelled sites.
- All 33 mapped UniProt residues are lysines.
- Sequence feature generation yields 101 derived descriptors per site.
- Published Figure 3C and replicate-reconstructed repair values are retained separately.

## Immediate milestone

Run the complete Stage 1 pipeline on the user's machine after official public files are downloaded:

```bash
cml-repairnet download
cml-repairnet stage1 --no-download
```

Completion criterion: `results/stage1_data_audit.md` reports `PASS` with zero missing structure
sites.

## Next development milestone

After Stage 1 passes:

1. freeze a benchmark data manifest,
2. add feature-family ablations,
3. implement repeated spatial GroupKFold,
4. add conformal prediction intervals,
5. establish a no-skill baseline and permutation significance analysis.

ESM-2 embeddings are intentionally postponed until the transparent benchmark is stable.
