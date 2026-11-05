# CML-RepairNet Reviewer Revision Test Patch v0.8.0

Added:

- `scripts/backup_before_revision.sh`
- `scripts/run_revision_tests.py`
- `src/cml_repairnet/revision/__init__.py`
- `src/cml_repairnet/revision/feature_ablation.py`
- `tests/test_revision_ablation.py`
- `configs/revision_tests.yaml`
- `docs/REVISION_TESTS_GUIDE_TR.md`

The patch does not delete or overwrite frozen Stage 1–6 result files. New analysis outputs are
written under `results/revision_tests/`.
