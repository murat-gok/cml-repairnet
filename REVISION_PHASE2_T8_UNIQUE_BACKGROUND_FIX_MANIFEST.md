# Patch manifest — v0.11.2

Changed files:

- `src/cml_repairnet/revision/phase2.py`
- `tests/test_revision_phase2.py`
- `README_T8_UNIQUE_BACKGROUND_FIX_TR.md`

Scientific behavior:

- collapses duplicate background gene sets by deterministic fingerprint;
- reports label count and distinct-background count separately;
- excludes ontology/database root nodes from the informative robust-term count;
- leaves model training, predictions, candidate tiers, compartment layers and STRING edge rules unchanged.

Validation:

- `PYTHONPATH=src python -m pytest -q tests/test_revision_phase2.py`
- 6 tests passed.
