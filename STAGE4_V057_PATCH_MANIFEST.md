# CML-RepairNet Stage 4 v0.5.7 Patch

This patch fixes full-atlas structure completion after a quick run.

## Changes
- Per-protein structure caches are now candidate-ID aware.
- Sites added by a later full run are processed instead of silently omitted.
- Cached `structure_processing_error:*` rows are automatically retried.
- Stable cached structure rows are retained; the full cache does not need deletion.
- AlphaFold request operations receive an additional outer retry/backoff layer.
- Stage 4 audit reports retryable structure failures and cache-incomplete rows.
- Model-card/project version updated to 0.5.7.

## Files
- `src/cml_repairnet/atlas/structures.py`
- `src/cml_repairnet/atlas/pipeline.py`
- `src/cml_repairnet/atlas/modeling.py`
- `tests/test_stage4_structure_cache_v057.py`
- `pyproject.toml`
- `docs/STAGE4_V057_GUIDE_TR.md`
