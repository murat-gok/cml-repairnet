# Stage 5 v0.6.0 Validation Report

## Code validation

Nine Stage 5 regression tests pass:

1. Stage 4 atlas schema and gene-symbol extraction.
2. Atlas loading and derived site labels.
3. Species-specific custom enrichment background.
4. g:Profiler response flattening and fold-enrichment calculation.
5. UniProt family parsing.
6. STRING edge deduplication and network metrics.
7. Human-to-mouse ortholog concordance.
8. Deterministic module assignment.
9. Local Stage 5 pipeline/report integration.

The Stage 5 source tree also passes Python byte-code compilation.

## Real-data smoke validation

The package was tested against the final Stage 4 v0.5.7 atlas:

- 100 priority sites loaded.
- 30 Pareto sites loaded.
- 749 in-domain sites loaded as the background universe.
- 85 unique portfolio proteins identified.
- Offline run completed without external caches.
- 100 candidate module assignments and two data-driven modules were produced.
- Local PNG/PDF figures and the biological interpretation report were generated.
- Missing external cache entries were recorded explicitly rather than silently ignored.

## External API validation boundary

The g:Profiler, STRING and UniProt clients use documented official REST
interfaces, retry/backoff, response caching and request manifests. In this build
environment, live external calls were not used for the end-to-end smoke test.
The first online run in the user's environment must therefore verify:

- nonzero UniProt annotation coverage;
- completed g:Profiler query manifests;
- recorded g:Profiler data versions;
- nonzero STRING identifier/network coverage;
- successful g:Orth mappings.

A failed external service does not corrupt the Stage 4 atlas. Its manifest is
marked `failed` or `offline_cache_missing`, and locally computable Stage 5
outputs continue.

## Expected test total in the current user project

Stage 4 v0.5.7 contained 46 tests. This patch adds 9 tests, so the expected total
is 55 passing tests after installation.

## Packaging environment note

The patch archive and source layout were validated, but an isolated editable-install build was not executed in this build environment because its package index did not expose `hatchling`. The user project already uses Hatchling and should install the patch with the existing `.venv`.
