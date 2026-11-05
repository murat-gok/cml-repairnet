# CML-RepairNet

**A reanalysis of enzymatic carboxymethyllysine repair data, asking how far
site-resolved PTM response measurements can be reused as model labels.**

This repository accompanies the manuscript

> *Condition identity and label weighting govern structure-response inference
> from reused site-resolved PTM data: a reanalysis of enzymatic
> carboxymethyllysine repair.*

The work is a reanalysis, not a predictor. It reuses the site-level
quantification deposited with the first report of an enzyme that reverses
N&epsilon;-carboxymethyllysine (ProteomeXchange `PXD070124`) and asks two
questions about using such data as supervised labels:

1. **Does the reporting condition matter?** The archive contains four
   substrate loads; the published site-level figure reports one without
   stating which. We identify it as the 500 ng condition and show that the
   structure-response conclusion is not interchangeable across loads.
2. **Is inverse-variance weighting safe here?** It is the textbook response
   to heteroscedastic measurements, and it turns out to be counterproductive:
   the weights span up to five orders of magnitude, the Kish effective sample
   size falls to 0.23-0.32 of the nominal site count, and the
   permutation-controlled signal is recovered in three of four conditions
   weighted against four of four unweighted. Clipping does not repair it.

The claims are bounded to what the deposited data support. There is no
wet-lab work here, and nothing in this repository should be read as a
validated prediction about any site not measured in the source experiment.

## Reproducing the manuscript

Every number, table and figure in the paper comes out of one command:

```bash
python -m pip install -e ".[dev,experiment]"
python scripts/run_multicondition_reanalysis.py --check
```

`--check` compares the regenerated values against the values printed in the
manuscript and exits non-zero if any of them has moved, so the paper itself
is covered by a regression test. The run fits a few thousand small forests
and takes several minutes.

Inputs, both tracked in this repository:

| Path | What it is |
| --- | --- |
| `data/external/pxd070124_parsed.json` | the four PXD070124 quantification tables, parsed; records the source URL, the column schema and a SHA-256 prefix per source file |
| `data/processed/bsa_structure_features_union44.csv` | nine transfer-safe structural descriptors from PDB 4F5S, for the union of confidently quantified sites |

Outputs land in `results/multicondition/`, which is tracked so the published
tables can be compared against a fresh run without re-executing anything.

Figures and LaTeX for the submission:

```bash
python scripts/build_jpr_tables.py     # the four manuscript tables
python scripts/make_jpr_figures.py     # figures 1-3, vector PDF + 600 dpi TIFF
python scripts/make_jpr_figureS1.py    # supporting figure S1
python scripts/make_jpr_toc.py         # table-of-contents graphic
```

The submission package, including the compiled PDF, is in
`docs/submission_JPR/`. The manuscript text is generated from
`docs/submission_JPR/MANUSCRIPT_JPR_TR.md` and
`docs/revision_A1/A3_METHODS_PROTEOMICS_TR.md`; the `.tex` is built, not
hand-edited.

## What the reanalysis does

The analysis reproduces the published pipeline rather than improving on it,
because the point is to test that pipeline:

- **Response.** The truncated relative reduction
  `clip(1 - mean(active) / mean(inactive), 0, 1)`.
- **Label variance.** Delta method for a ratio of two means from three
  replicates each: `var(A)/C² + A²·var(C)/C⁴`.
- **Weights.** Inverse-variance with a floor, median-normalised within a
  load. Compared against unweighted, square-root and clipped schemes.
- **Partitions.** Spatially blocked: Ward linkage on Cα coordinates, so a
  held-out block is a contiguous region of the fold.
- **Null and significance.** The weighted training-fold mean as the null;
  P values from 200 site-structure permutations, for which 0.005 is the
  resolution limit.
- **Composition control.** The twenty lysines quantified confidently at
  every load are re-run as a fixed set, and read against 200 random subsets
  of the same size, which separates composition from loss of power.

Implementation: `src/cml_repairnet/revision/multicondition.py`. Every
parameter above is a field on `MulticonditionConfig`, so a reviewer can vary
one choice and see what it costs.

## Repository layout

```text
CML-RepairNet/
├── configs/                      # Pipeline configuration
├── data/external/                # Parsed source quantification (tracked)
├── data/processed/               # Derived label and descriptor tables (tracked)
├── docs/revision_A1/             # Methods source text
├── docs/submission_JPR/          # Submission package: LaTeX, PDF, figures, cover letter
├── results/multicondition/       # Reanalysis outputs (tracked)
├── scripts/                      # Executable entry points
├── src/cml_repairnet/revision/   # The reanalysis
└── tests/                        # Unit and scientific-integrity tests
```

## Public inputs

- Source article: DOI `10.1038/s41467-026-75141-2`
- ProteomeXchange / PRIDE: `PXD070124`
- Bovine serum albumin: UniProt `P02769`
- Experimental BSA structure: PDB `4F5S`

Residue numbering is handled explicitly throughout: `site_pdb =
site_uniprot − 24` maps precursor positions in P02769 to mature coordinates
in 4F5S chain A. Large third-party downloads are excluded from Git; see
`docs/DATA_AND_LICENSES.md`.

## Environment

Python `3.10`–`3.12`, with `3.11` as the reference. Editable install via
`pyproject.toml`.

```bash
bash scripts/bootstrap.sh
source .venv/bin/activate
python -m pytest -q
ruff check src scripts tests
```

## Earlier exploratory work, not part of the manuscript

This repository began as a broader attempt to build a proteome-scale
repairability predictor, and that history is still here: the staged atlas
pipeline, proteome-wide candidate prioritization, functional enrichment,
STRING-based network analysis and orthology mapping, under
`src/cml_repairnet/atlas/`, `src/cml_repairnet/biology/` and the
`scripts/run_stage*.py` entry points, with their guides in `docs/`.

**None of it supports a claim in the manuscript, and it is retained for
provenance only.** Those analyses did not survive the controls reported in
the paper and the accompanying revision notes; the proteome-scale
prioritization and the human enrichment results in particular were dropped
for that reason. Treat anything under those paths as a record of what was
tried, not as a result.

The revision record is in `docs/REVIZYON_YOL_HARITASI_TR.md` (Turkish).

## License

Code: MIT, per the `LICENSE` file at the repository root. Third-party
datasets retain their original terms and are not covered by the software
license.
