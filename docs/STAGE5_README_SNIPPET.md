## Stage 5 — Biological interpretation

Stage 5 interprets the frozen Stage 4 atlas using UniProt annotations,
g:Profiler enrichment with species-specific in-domain backgrounds, STRING v12.0
networks, human–mouse orthology, candidate modules, and manuscript-ready figures.

```bash
python scripts/run_stage5.py doctor
python -u scripts/run_stage5.py all 2>&1 | tee results/stage5_full_console.txt
```

Results are biological-context and hypothesis-generation outputs, not independent
experimental validation of enzymatic CML repairability.
