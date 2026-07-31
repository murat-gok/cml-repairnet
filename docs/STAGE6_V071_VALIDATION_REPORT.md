# Stage 6 v0.7.1 validation report

## Tests

- Stage 6 regression tests: 9/9 passed.
- Python byte-code compilation: passed.
- Full real-atlas run: passed with exit code 0.

## Real-atlas validation

- In-domain candidates: 749
- Monte Carlo iterations: 2,000
- Pareto frequency iterations: 500
- Constrained portfolio slots per iteration:
  - human: 50
  - mouse: 50
- Baseline Stage 4 portfolio candidates with constrained selection frequency
  >=0.80: 21/100
- Residue-level conserved lysine pairs: 8/10
- Significant enrichment terms assessed: 65
- Fragile terms: 2
- Tier 1: 4
- Tier 2: 33
- Tier 3: 712
- Watchlist: 519

## Interpretation

v0.7.0 enrichment, network and residue-conservation results remain valid. Only
rank/portfolio stability, weight-scenario overlap, candidate tiers, the rank
figure, audit and report require regeneration.
