# CML-RepairNet Stage 2.5 Patch v0.3.0

This patch adds:

- Regularized XGBoost and LightGBM regressors
- Structure-only nested spatial GroupKFold evaluation
- Training-fold correlation pruning
- Clipped and square-root assay-weight sensitivity
- Leakage-controlled Ridge–Extra Trees convex ensemble
- Held-out permutation importance
- Optional descriptive SHAP export
- Dummy and reference paired bootstrap/sign-flip comparisons with BH-FDR
- VS Code tasks, Make targets, tests and Turkish usage guide

Install with:

```bash
python -m pip install -e ".[dev,experiment]"
python -m pytest -q
python -u scripts/run_stage25.py --mode quick
```
