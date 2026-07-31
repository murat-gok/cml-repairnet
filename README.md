# CML-RepairNet

**A fully in silico, structure-aware and uncertainty-calibrated framework for predicting site-specific enzymatic CML repairability.**

CML-RepairNet models how strongly a CML-modified lysine site may respond to CMLase under the
conditions represented by the published BSA benchmark. The repository is designed for VS Code,
Python 3.11 and reproducible command-line execution.

## Non-negotiable scope

This project is **100% in silico**. It contains no wet-lab, animal, clinical or prospective
experimental work. Published measurements are reused only as public training, benchmarking and
external-consistency data. Proteome-wide outputs are predictions and must always include
uncertainty and out-of-distribution flags.

## Start here: VS Code + WSL

For best filesystem and scientific-computing performance, keep the repository in the WSL Linux
filesystem, not under `/mnt/c`.

```bash
mkdir -p ~/projects
cd ~/projects
unzip /mnt/c/Users/user/Downloads/CML-RepairNet_VSCode.zip
cd CML-RepairNet
code .
```

In the VS Code integrated WSL terminal:

```bash
bash scripts/bootstrap.sh
source .venv/bin/activate
cml-repairnet doctor
```

Then build the public-data foundation:

```bash
cml-repairnet download
cml-repairnet stage1 --no-download
```

Expected Stage 1 outputs:

- `data/processed/bsa_labels.csv`
- `data/processed/bsa_sequence_features.csv`
- `data/processed/bsa_structure_features.csv`
- `results/stage1_data_audit.json`
- `results/stage1_data_audit.md`

The same operations are available through **Terminal → Run Task**:

1. `CML: 0 - Bootstrap environment`
2. `CML: 1 - Download public data`
3. `CML: 2 - Run Stage 1 data foundation`

Detailed setup: `docs/VS_CODE_SETUP_TR.md`.

## Why Stage 1 comes before model training

The initial benchmark contains only 33 site-resolved BSA observations. The first milestone is
therefore not a deep model. It is a traceable, audited dataset with correct numbering across:

- mature BSA positions used in the CMLase paper,
- precursor positions in UniProt P02769,
- mature 583-residue coordinates in PDB 4F5S.

The project keeps the published Figure 3C repair value as the primary target and independently
reconstructs a sensitivity target from Figure 3B replicates. These values are never silently
substituted for one another.

## Repository structure

```text
CML-RepairNet/
├── .vscode/                    # Tasks, debugger and recommended extensions
├── configs/base.yaml           # Single source of pipeline configuration
├── data/                       # Downloaded/generated data; ignored by Git
├── docs/                       # Protocol, setup and data documentation
├── results/                    # Generated audit/model outputs; ignored by Git
├── scripts/                    # Thin executable wrappers
├── src/cml_repairnet/          # Reusable Python package
└── tests/                      # Unit and scientific-integrity tests
```

## Main commands

```bash
cml-repairnet doctor
cml-repairnet download
cml-repairnet prepare-labels
cml-repairnet stage1 --no-download
python scripts/train_baselines.py
python -m pytest -q
ruff check src scripts tests
```

## Current computational design

1. Reconstruct replicate-aware labels from the article source data.
2. Audit mature/precursor/PDB residue mapping.
3. Build transparent sequence and static-structure descriptors.
4. Establish leakage-controlled small-data baselines.
5. Add frozen protein-language-model embeddings only after the core benchmark passes.
6. Quantify uncertainty and out-of-distribution risk.
7. Expand to published human/mouse CML candidate sites for virtual prioritization.
8. Allocate higher-cost docking/MD calculations through multi-objective optimization.

## Public inputs

- Article/source data: DOI `10.1038/s41467-026-75141-2`
- ProteomeXchange/PRIDE: `PXD070124`
- Bovine serum albumin: UniProt `P02769`
- Experimental BSA structure: PDB `4F5S`

Third-party files are downloaded locally and excluded from Git. See
`docs/DATA_AND_LICENSES.md`.

## Python support

Python `3.10–3.12` is supported; Python `3.11` is the reference environment. The project uses a
modern `pyproject.toml` package layout and an editable installation for development.

## License

Code: Apache-2.0. Third-party datasets retain their original terms and are not covered by the
software license.


## Stage 2.5: regularized boosting and structural robustness

Install the optional experiment dependencies and run the prespecified quick analysis:

```bash
python -m pip install -e ".[dev,experiment]"
python -u scripts/run_stage25.py --mode quick
```

The full Stage 2.5 matrix adds both targets, 3/4/5 spatial groups, clipped/square-root label weights, nested correlation pruning, XGBoost, LightGBM, a leakage-controlled Ridge–Extra Trees ensemble and held-out permutation importance. See `docs/STAGE25_GUIDE_TR.md`.


## Stage 3: generalization, mechanism and uncertainty

Stage 3 freezes the Stage 2.5 model portfolio and evaluates partition sensitivity,
leave-one-domain-out extrapolation, boundary calibration, conformal intervals, OOD risk and
the conditional contribution of the local cysteine neighbourhood.

```bash
python -u scripts/run_stage3.py --mode quick --n-jobs 1
python -u scripts/run_stage3.py --mode full --n-jobs 1
```

See `docs/STAGE3_GUIDE_TR.md`.
