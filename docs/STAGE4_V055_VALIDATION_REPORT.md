# Stage 4 v0.5.5 Validation Report

## Regression tests

Four new tests cover:

1. Canonical model preference over a shorter, higher-confidence isoform.
2. Rejection of an isoform whose target residue conflicts with UniProt.
3. Parsing of post-sunset AlphaFold API fields.
4. Correct fragment-local mapping when a direct coordinate number is a decoy.

The cumulative reconstructed project test suite passed.

## User-result diagnosis

The submitted quick structures run contained 66 sequence-valid candidates:

- 56 quality-eligible structures
- 10 excluded structures
- 9 exclusions used a hyphenated isoform AlphaFold entry
- 1 exclusion used the canonical model and failed the pLDDT threshold

The nine isoform-selected records were:

- O60701 K220
- P02545 K155 and K311
- P04406 K194, K254 and K259
- P06733 K89
- P06748 K236 and K239

Seven failed because the structure residue was not lysine. Two P06748 sites
reached the quality stage but had low pLDDT in the incorrectly selected isoform.

The 56 currently eligible rows had complete transfer features and median site
pLDDT 96.88. The existing scored table was stale: it contained 46 rows and only
16 candidate IDs overlapped the new structure-feature table.
