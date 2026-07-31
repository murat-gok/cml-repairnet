# Analysis protocol — fully in silico

## Confirmatory endpoint

Out-of-fold mean absolute error (MAE) for `repair_published` under spatially blocked nested
cross-validation.

## Secondary endpoints

- RMSE
- Spearman rank correlation
- R² with an explicit small-n warning
- probabilistic calibration and interval coverage
- resistant-site (<0.10) and highly-repairable-site (>0.90) discrimination
- stability across fold definitions and label reconstructions

## Required model hierarchy

1. Mean and median predictors
2. Nearest-neighbor and simple SASA-only controls
3. Sequence-only Ridge/Elastic Net
4. Structure-only Gaussian Process
5. Combined sequence + structure model
6. Frozen PLM embedding + low-dimensional structural features
7. Conformational-ensemble model
8. Optional multi-fidelity model using docking/MD-derived auxiliary scores

A complex model is retained only if it improves blocked MAE, calibration or robustness over the
best simpler model.

## Validation hierarchy

1. Spatial GroupKFold
2. Sequence-motif cluster holdout
3. Structural-environment cluster holdout
4. Label-definition sensitivity
5. PDB-versus-AlphaFold feature sensitivity
6. Conformational-ensemble sensitivity
7. Protein-level external-consistency analysis using published non-BSA measurements, when numeric
   data are compatible

## Leakage controls

- Scaling, imputation, feature selection, PCA and hyperparameter tuning occur inside training folds.
- Spatial coordinates define groups but are excluded from predictive features in the initial benchmark.
- Sites from overlapping peptides or highly similar local windows remain in the same fold.
- Protein-level splitting is mandatory once multiple proteins carry site-resolved labels.
- Proteome-wide candidate data never tune the BSA benchmark.
- Biological annotations used in target prioritization are excluded from the repairability model unless
  explicitly evaluated as a separate hypothesis.

## Uncertainty and OOD

- Aleatoric uncertainty is approximated from replicate variation.
- Epistemic uncertainty is estimated using GPR, ensembles or bootstrap models.
- Conformal intervals are calibrated only from held-out predictions.
- OOD scores combine sequence-window distance, structural-feature distance, PLM embedding distance
  and structure-confidence penalties.
- Proteome predictions are published with evidence tiers, not as uniformly validated outputs.

## Statistical tests

- Bootstrap 95% confidence intervals for out-of-fold metrics
- Label-permutation test for the full modelling pipeline
- Paired bootstrap comparison with the best baseline
- Sensitivity analysis with inverse-variance weighting
- Sensitivity analysis excluding low-control-abundance and near-zero-repair sites
- Multiple-testing correction for feature-level exploratory analyses

## Interpretability

Prefer coefficient stability, permutation importance, accumulated local effects, feature-family
ablation and counterfactual local-window mutation. SHAP is exploratory and must not be interpreted
as causal evidence.

## Multi-fidelity rule

Docking, MD, MM/GBSA or QM-derived quantities are auxiliary computational evidence. They are not
substitutes for experimental ground truth and cannot be used to relabel candidates as confirmed.
