from __future__ import annotations

import numpy as np
import pandas as pd

from .schema import add_labels


def run_candidate_tiers(
    atlas_all: pd.DataFrame,
    stability: pd.DataFrame,
    conservation: pd.DataFrame,
    modules: pd.DataFrame,
    config: dict,
) -> pd.DataFrame:
    atlas = add_labels(atlas_all)
    stability_columns = [
        "candidate_id", "mean_stress_rank", "rank_iqr", "top50_selection_frequency",
        "top100_selection_frequency", "portfolio_selection_frequency", "pareto_membership_frequency",
    ]
    merged = atlas.merge(stability[stability_columns], on="candidate_id", how="left")
    if "module_label" in modules.columns:
        merged = merged.merge(modules[["candidate_id", "module_label"]], on="candidate_id", how="left")

    conserved_ids: set[str] = set()
    if not conservation.empty:
        good = conservation.loc[conservation["residue_level_conserved_lysine"].fillna(False)]
        for row in good.itertuples(index=False):
            human = atlas.loc[(atlas["species"] == "human") & (atlas["uniprot_accession"] == row.human_accession) & (atlas["site_uniprot"] == row.human_site), "candidate_id"]
            mouse = atlas.loc[(atlas["species"] == "mouse") & (atlas["uniprot_accession"] == row.mouse_accession) & (atlas["site_uniprot"] == row.mouse_site), "candidate_id"]
            conserved_ids.update(human.astype(str)); conserved_ids.update(mouse.astype(str))
    merged["residue_level_conservation_support"] = merged["candidate_id"].astype(str).isin(conserved_ids)
    evidence = (
        merged[[column for column in ["prm_validated", "age_associated", "ortholog_overlap"] if column in merged]]
        .fillna(False).astype(bool).sum(axis=1)
    )
    merged["key_evidence_count"] = evidence

    tier1 = (
        (merged["applicability_status"] == "in_domain")
        & (merged["portfolio_selection_frequency"].fillna(0) >= float(config.get("tier1_portfolio_frequency", config.get("tier1_top100_frequency", 0.85))))
        & (merged["pareto_membership_frequency"].fillna(0) >= float(config["tier1_pareto_frequency"]))
        & (pd.to_numeric(merged["transfer_reliability_score"], errors="coerce").fillna(0) >= float(config["tier1_transfer_reliability"]))
        & (pd.to_numeric(merged["model_disagreement"], errors="coerce").fillna(np.inf) <= float(config["tier1_max_model_disagreement"]))
        & ((merged["key_evidence_count"] >= 1) | merged["residue_level_conservation_support"])
    )
    tier2 = (
        (merged["applicability_status"] == "in_domain")
        & (merged["portfolio_selection_frequency"].fillna(0) >= float(config.get("tier2_portfolio_frequency", config.get("tier2_top100_frequency", 0.60))))
        & (pd.to_numeric(merged["transfer_reliability_score"], errors="coerce").fillna(0) >= float(config["tier2_transfer_reliability"]))
    )
    tier3 = merged["applicability_status"] == "in_domain"
    merged["candidate_tier"] = np.select(
        [tier1, tier2, tier3],
        ["Tier 1", "Tier 2", "Tier 3"],
        default="Watchlist",
    )
    merged["tier_reason"] = np.select(
        [tier1, tier2, tier3],
        [
            "stable Pareto/diversity-constrained portfolio candidate with evidence or residue conservation",
            "stable diversity-constrained portfolio candidate with acceptable transfer reliability",
            "in-domain candidate with limited stability or supporting evidence",
        ],
        default="borderline/out-of-domain transfer candidate",
    )
    tier_order = pd.Categorical(merged["candidate_tier"], ["Tier 1", "Tier 2", "Tier 3", "Watchlist"], ordered=True)
    merged["_tier_order"] = tier_order
    return merged.sort_values(
        ["_tier_order", "mean_stress_rank", "weighted_priority_score"],
        ascending=[True, True, False],
        na_position="last",
    ).drop(columns="_tier_order").reset_index(drop=True)
