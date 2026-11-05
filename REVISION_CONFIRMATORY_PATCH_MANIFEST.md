# CML-RepairNet Reviewer Confirmatory Tests v0.9.0

Added/updated:

- `src/cml_repairnet/revision/feature_ablation.py`
  - structure-only proxy family
  - B-factor-only and proxy-only diagnostics
  - direct paired B-factor/proxy contrasts
  - optional targeted family selection
  - meaningless dummy-model Spearman suppressed
- `src/cml_repairnet/revision/confirmatory.py`
  - multi-partition resumable structure permutation control
  - nested spatial model selection
  - paired spatial-cluster bootstrap and sign-flip test
- `scripts/run_revision_confirmatory.py`
- `configs/revision_confirmatory.yaml`
- `tests/test_revision_ablation.py`
- `tests/test_revision_confirmatory.py`
- `docs/REVISION_CONFIRMATORY_GUIDE_TR.md`
- `docs/PRELIMINARY_REVISION_INTERPRETATION_TR.md`

No Stage 4–6 atlas, ranking or frozen result file is deleted or overwritten. New outputs default to
`results/revision_confirmatory/`.
