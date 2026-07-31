# CML-RepairNet Stage 6 v0.7.1 patch manifest

## Purpose

Correct the Stage 6 portfolio/Pareto stability analysis so that it respects the
frozen Stage 4 prioritization design.

## Corrected defects

1. Stage 6 v0.7.0 computed Pareto membership across the combined human/mouse
   universe. Stage 4 computes Pareto fronts separately within each species.
2. Stage 6 v0.7.0 called an unconstrained global top-100 frequency "portfolio
   retention", although the Stage 4 portfolio applies approximate 50/50 species
   balance and a maximum of two sites per protein.
3. Tier 1/2 labels used the unconstrained top-100 frequency and were therefore
   not aligned with the diversity-controlled Stage 4 portfolio.
4. Weight-scenario overlap compared unconstrained global top-100 sets rather
   than constrained portfolios.

## Files replaced

- `configs/stage6.yaml`
- `src/cml_repairnet/robustness/ranking.py`
- `src/cml_repairnet/robustness/tiers.py`
- `src/cml_repairnet/robustness/pipeline.py`
- `src/cml_repairnet/robustness/report.py`
- `src/cml_repairnet/robustness/figures.py`
- `tests/test_stage6.py`

## New output fields

- `portfolio_selection_frequency`
- `portfolio_overlap_count`
- `portfolio_jaccard`

The existing `top10/20/50/100_selection_frequency` fields remain available as
unconstrained global-rank diagnostics and are no longer used to define tiers.

## Unchanged analyses

Residue conservation, enrichment sensitivity, STRING threshold/hub/null tests,
and their existing outputs are unchanged.
