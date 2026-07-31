# CML-RepairNet Stage 4 v0.5.6 Patch

Changed files:

- `src/cml_repairnet/atlas/prioritize.py`
- `src/cml_repairnet/atlas/modeling.py`
- `src/cml_repairnet/atlas/pipeline.py`
- `tests/test_stage4_pareto_v056.py`
- `docs/STAGE4_V056_GUIDE_TR.md`

Purpose: prevent high-dimensional Pareto-front inflation and exclude borderline/OOD candidates from the primary Pareto front while preserving them in a separate watchlist.
