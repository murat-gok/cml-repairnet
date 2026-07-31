# Stage 4 v0.5.4 Validation Report

Validation was performed against the two original supplementary workbooks supplied
by the user.

## Test result

- Full available local project test suite: 36/36 passed.
- New v0.5.4 regression tests: 6.

## Real-workbook parser result

- Aggregated candidate sites: 1681
- Unique UniProt accessions: 778
- Human sites: 298
- Mouse sites: 1383
- PRM CML-marked peptides retained for sequence mapping: 17
- Unique age-tested candidate sites: 287
- Unique FDR-significant age-associated sites: 81
- Candidates supported by ortholog-overlap evidence: 302

## Regression safeguards

- `P01942` is not parsed as site 1942.
- `P62806` is not parsed as site 62806.
- PRMT-containing protein names do not trigger `prm_validated=True`.
- Accession identifiers are not stored as peptide examples.
- The overlap table produces separate HUVEC and MEF evidence rows.
- Heart/kidney/liver age tables are parsed with adjusted P value <= 0.05 as the
  prespecified age-association rule.

Actual PRM sequence mapping requires UniProt network/cache access and is performed
in the user's connected Stage 4 annotation run. Synthetic and cached-sequence
mapping regression tests passed.
