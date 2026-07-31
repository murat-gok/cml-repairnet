# CML-RepairNet Stage 5 Patch v0.6.0

## Added

- `src/cml_repairnet/biology/`
  - API clients with retry, caching and offline mode
  - UniProt annotation and family summaries
  - g:Profiler custom-background enrichment
  - STRING v12.0 network and PPI enrichment
  - human–mouse ortholog concordance
  - candidate modules
  - publication figure generation
  - audit and interpretation report
- `scripts/run_stage5.py`
- `configs/stage5.yaml`
- `tests/test_stage5.py` — 9 regression tests
- Stage 5 analysis plan, Turkish guide and validation report
- VS Code tasks and Makefile targets

## Changed

- Project version: `0.6.0`
- Added runtime dependencies: `networkx`, `tabulate`

## Scientific safeguards

- Stage 4 model remains frozen.
- Enrichment uses species-specific in-domain custom backgrounds.
- Human and mouse are tested separately.
- STRING adds no external nodes.
- External API versions/manifests are retained.
- Stage 5 results are hypothesis-generating biological context, not experimental validation.
