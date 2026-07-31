from __future__ import annotations

import json
from pathlib import Path

import pandas as pd


def _table(frame: pd.DataFrame, columns: list[str], n: int = 15) -> str:
    existing = [column for column in columns if column in frame]
    if not existing or frame.empty:
        return "No rows available."
    return frame[existing].head(n).to_markdown(index=False)


def build_report(
    stability: pd.DataFrame,
    weight_sensitivity: pd.DataFrame,
    conservation: pd.DataFrame,
    enrichment: pd.DataFrame,
    network_thresholds: pd.DataFrame,
    hub_removal: pd.DataFrame,
    network_null: pd.DataFrame,
    tiers: pd.DataFrame,
    audit: dict,
    output_path: Path,
) -> None:
    conserved = int(conservation.get("residue_level_conserved_lysine", pd.Series(dtype=bool)).fillna(False).sum())
    portfolio_stability = stability.loc[stability["baseline_portfolio_member"] == 1]
    stable_portfolio = int((portfolio_stability.get("portfolio_selection_frequency", pd.Series(dtype=float)) >= 0.8).sum())
    fragile_terms = int(enrichment.get("driver_sensitivity_class", pd.Series(dtype=str)).isin(["single_driver", "leave_one_driver_fragile"]).sum())
    tier_counts = tiers["candidate_tier"].value_counts().to_dict()
    text = f"""# CML-RepairNet Stage 6 — Robustness and manuscript validation

## Scope

Stage 6 stress-tests the frozen Stage 4–5 outputs. It does not search for a new
model, tune Stage 3 performance, or create experimental validation.

## Completion summary

- In-domain candidates stress-tested: {audit['in_domain_candidates']}
- Monte Carlo rank-stress iterations: {audit['rank_iterations']}
- Frozen Stage 4 portfolio reproduced without perturbation: {audit.get('baseline_portfolio_reproduction_count', 0)}/{len(portfolio_stability)}
- Baseline portfolio candidates retained in the diversity-constrained portfolio with frequency >=0.80: {stable_portfolio}/{len(portfolio_stability)}
- Residue-level conserved lysine pairs: {conserved}/{len(conservation)}
- Significant enrichment terms assessed: {len(enrichment)}
- Potentially single-driver/leave-one-driver-fragile terms: {fragile_terms}
- Tier 1 sites: {tier_counts.get('Tier 1', 0)}
- Tier 2 sites: {tier_counts.get('Tier 2', 0)}
- Tier 3 sites: {tier_counts.get('Tier 3', 0)}
- Watchlist sites: {tier_counts.get('Watchlist', 0)}

## Rank and priority robustness

The Monte Carlo analysis is an uncertainty and weighting stress test, not a
sampling-distribution bootstrap. It perturbs model mixtures, interval-scaled
scores, evidence, reliability, and objective weights under the prespecified
Stage 6 configuration.

### Most stable candidates

{_table(stability.sort_values(['portfolio_selection_frequency','pareto_membership_frequency'], ascending=False), ['site_label','species','baseline_rank','mean_stress_rank','rank_iqr','portfolio_selection_frequency','pareto_membership_frequency'])}

The constrained portfolio stress recomputes complete species-specific nondominated fronts, orders candidates by weighted priority within each front, and reapplies the Stage 4 species-balance and per-protein safeguards. Global top-k frequencies remain available only as unconstrained ranking diagnostics.

### Weight scenarios

{_table(weight_sensitivity, list(weight_sensitivity.columns), n=20)}

## Residue-level cross-species conservation

A global protein-sequence alignment maps each human candidate lysine to the
mouse sequence and requires reciprocal position mapping plus lysine identity.
This supports residue-level conservation only for pairs that pass those checks.

{_table(conservation, ['human_gene','mouse_gene','human_site','mouse_site','mapped_mouse_site_from_human','residue_level_conserved_lysine','local_sequence_identity','global_sequence_identity','repairability_abs_difference'], n=20)}

## Enrichment driver sensitivity

Leave-one-intersection-gene raw hypergeometric sensitivity is used to flag
terms whose signal can be lost by removing a single contributing gene. Stage 5
adjusted p-values remain the authoritative multiple-testing results; Stage 6
does not attempt to reproduce g:SCS correction.

{_table(enrichment.sort_values(['driver_sensitivity_class','adjusted_p_value_stage5']), ['species','query_name','source','term_name','intersection_size','adjusted_p_value_stage5','worst_intersection_gene_removal_raw_p','driver_sensitivity_class','redundancy_group_size'])}

## STRING robustness

Threshold and hub-removal analyses quantify dependence on confidence cutoffs
and highly connected proteins. Degree-preserving rewiring tests topology given
the observed degree sequence; it is explicitly not a random-gene-set test.

{_table(network_thresholds, ['species','query_name','threshold','nodes','edges','largest_component_fraction','isolated_nodes','transitivity'], n=20)}

{_table(hub_removal, ['species','query_name','removed_hub_count','removed_hubs','edge_retention_fraction','largest_component_retention_fraction'], n=20)}

{_table(network_null, list(network_null.columns), n=20)}

## Final candidate tiers

Tier labels are deterministic evidence-integration categories defined in
`configs/stage6.yaml`. They are prioritization labels, not clinical or
experimental grades.

{_table(tiers, ['candidate_tier','site_label','species','repairability_extra_trees','weighted_priority_score','portfolio_selection_frequency','pareto_membership_frequency','key_evidence_count','residue_level_conservation_support','tier_reason'], n=30)}

## Interpretation boundary

- Stage 6 does not increase the 33-site BSA training sample.
- Rank frequencies depend on the prespecified stress distributions.
- Residue conservation is sequence-based and does not prove conserved repair mechanism.
- Enrichment and STRING results remain computational hypothesis generation.
- Tier 1–3 labels must not be described as experimentally validated classes.
"""
    output_path.write_text(text, encoding="utf-8")
