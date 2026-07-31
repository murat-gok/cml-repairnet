# Stage 6 v0.7.0 validation report

## Test environment

The patch was applied over the frozen Stage 4 v0.5.7 and Stage 5 v0.6.1
result packages and executed with Python 3.12.

## Regression tests

- Stage 6 tests: 7/7 passed.
- Stage 5 + Stage 6 tests available in the validation fixture: 18/18 passed.
- Expected complete user-project suite after the prior 57-test freeze: 64 tests.

## Real-data end-to-end validation

The `all` command completed with:

- 749 in-domain candidates;
- 2,000 rank-stress iterations;
- 500 Pareto-stress iterations;
- 4 objective-weight scenarios;
- 10 human/mouse ortholog-gene pairs aligned;
- 8 reciprocal residue-level conserved lysine pairs;
- 65 significant Stage 5 enrichment terms assessed;
- 2 leave-one-driver-fragile terms;
- 12 STRING threshold rows, 12 hub-removal rows, and 4 topology-null rows;
- 1,268 final candidate labels: 2 Tier 1, 31 Tier 2, 716 Tier 3, and 519 Watchlist;
- 5 PNG/PDF figure sets;
- non-empty CSV, JSON, Markdown, PNG, and PDF outputs.

## Key safeguards validated

- The analysis uses only frozen local Stage 4/5 files.
- Human/mouse residue mapping uses the complete UniProt sequences already
  embedded in the Stage 5 pair table.
- TUBA1B and STAT1 are correctly rejected as reciprocal residue-level matches,
  despite high local sequence identity, because their declared CML lysines do
  not map to each other in the global alignment.
- STRING robustness includes the complete query-gene universe from the Stage 5
  network manifest, including isolated nodes.
- Network rewiring is labelled as a degree-preserving topology null and not as
  a random-gene-set analysis.
- Stage 5 adjusted enrichment p-values remain authoritative; Stage 6 reports
  transparent raw hypergeometric driver sensitivity only.
