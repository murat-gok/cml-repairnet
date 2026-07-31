# Stage 5 Prespecified Biological Analysis Plan

## Frozen inputs

Stage 5 uses only the frozen Stage 4 v0.5.7 outputs:

- 100-site priority portfolio
- 30-site in-domain Pareto front
- 749-site in-domain atlas as the enrichment universe
- domain-shift watchlist for sensitivity reporting only

No Stage 5 result is allowed to alter or retrain the Stage 4 repairability model.

## Prespecified gene sets

For human and mouse separately:

1. Priority portfolio proteins.
2. Pareto-front proteins.
3. Evidence-rich portfolio proteins: PRM validation, age association, or ortholog support.
4. High-repairability in-domain proteins: top quartile of Extra Trees predictions.

## Enrichment universe

The statistical background is the set of unique genes represented in the
species-specific Stage 4 in-domain atlas. Whole-genome background is not the
primary analysis because the candidate-generation and structure-quality filters
create a restricted observable universe.

## Multiple testing

- g:Profiler significance method: FDR.
- Threshold: adjusted p <= 0.05.
- Results are retained even when nonsignificant, with a `significant` flag.
- Human and mouse are analysed independently.

## Network analysis

- STRING version: 12.0 stable endpoint.
- Network type: functional.
- Required confidence: 700/1000.
- Added nodes: zero.
- Portfolio and Pareto networks are analysed separately.
- Centrality is descriptive and not a causal measure.

## Orthology

Human portfolio genes are mapped to mouse using g:Orth. Cross-species comparison
is performed at gene level using the highest-priority site per species. Site
number equivalence is not assumed unless separately supported by Stage 4
ortholog evidence.

## Candidate modules

The 100 portfolio sites are clustered using standardized:

- Extra Trees repairability
- experimental evidence score
- transfer reliability
- model disagreement
- site pLDDT
- relative lysine SASA
- 8 A cysteine count
- local hydrophobic fraction

The number of modules is selected by silhouette score over k=2..6.

## Interpretation boundary

Enrichment, PPI connectivity, protein family membership, and ortholog agreement
are biological-context signals. They do not constitute experimental evidence
that a site is enzymatically repairable.
