# Stage 6 prespecified analysis plan

## Objective

Stage 6 evaluates robustness and manuscript readiness of the frozen Stage 4 and
Stage 5 outputs. It must not be used to select a new predictive algorithm or to
retune the 33-site BSA benchmark.

## Prespecified modules

1. **Rank-stability stress test** — model-mixture, interval-scaled score,
   evidence, reliability and objective-weight perturbations.
2. **Pareto membership stability** — frequency under the same stress design.
3. **Residue-level ortholog alignment** — reciprocal human/mouse site mapping,
   lysine identity and ±7-residue local sequence identity.
4. **Enrichment driver sensitivity** — single-intersection-gene removal using a
   transparent hypergeometric diagnostic; Stage 5 g:SCS-adjusted p-values remain
   authoritative.
5. **STRING robustness** — 0.70/0.80/0.90 thresholds, hub removal and
   degree-preserving rewiring.
6. **Tier 1–3 integration** — deterministic rules in `configs/stage6.yaml`.

## Interpretation safeguards

- Rank frequencies are stress-test frequencies, not bootstrap confidence.
- Residue conservation is sequence-level evidence, not mechanistic validation.
- Enrichment leave-one-out raw p-values do not reproduce g:SCS correction.
- Degree-preserving network rewiring is not a random-gene-set comparison.
- Tier labels are computational prioritization categories.
