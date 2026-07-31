# Stage 6 v0.7.2 patch manifest

Modified files:

- `src/cml_repairnet/robustness/ranking.py`
  - complete species-specific nondominated front sorting;
  - exact Stage 4 front-rank + priority portfolio order;
  - frozen-baseline reproduction diagnostics.
- `src/cml_repairnet/robustness/pipeline.py`
  - version 0.7.2 audit fields;
  - baseline reproduction count/fraction;
  - corrected methodological scope.
- `src/cml_repairnet/robustness/report.py`
  - reports baseline 100/100 reproduction and complete-front stress design.
- `src/cml_repairnet/robustness/figures.py`
  - corrected stability figure title.
- `tests/test_stage6.py`
  - complete-front rank test;
  - lower-front precedence test;
  - baseline reproduction test.
- `docs/STAGE6_V072_GUIDE_TR.md`

Unchanged analyses:

- residue-level conservation;
- enrichment driver sensitivity;
- STRING threshold, hub-removal and topology-null analyses.
