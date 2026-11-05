from __future__ import annotations

import numpy as np
import pandas as pd

from cml_repairnet.robustness.conservation import run_residue_conservation
from cml_repairnet.robustness.enrichment_sensitivity import run_enrichment_sensitivity
from cml_repairnet.robustness.network_robustness import run_network_robustness
from cml_repairnet.robustness.ranking import (
    nondominated_front_ranks,
    pareto_mask,
    run_rank_stability,
    select_diversity_constrained_indices,
    specieswise_pareto_mask,
)
from cml_repairnet.robustness.tiers import run_candidate_tiers


def _atlas(n: int = 12) -> pd.DataFrame:
    values = np.linspace(0.2, 0.9, n)
    return pd.DataFrame(
        {
            "candidate_id": [f"c{i}" for i in range(n)],
            "species": ["human"] * n,
            "uniprot_accession": [f"P{i:05d}" for i in range(n)],
            "site_uniprot": np.arange(1, n + 1),
            "gene_names": [f"G{i}" for i in range(n)],
            "repairability_extra_trees": values,
            "repairability_lightgbm": np.clip(values + 0.02, 0, 1),
            "repairability_ridge": np.clip(values - 0.03, 0, 1),
            "repairability_score_norm": values,
            "evidence_score_norm": np.linspace(0.1, 0.8, n),
            "evidence_score_raw": np.arange(n),
            "transfer_reliability_score": np.linspace(0.4, 0.9, n),
            "prediction_interval_width": np.repeat(0.2, n),
            "weighted_priority_score": values,
            "model_disagreement": np.repeat(0.05, n),
            "applicability_status": ["in_domain"] * n,
            "prm_validated": [False] * (n - 1) + [True],
            "age_associated": [False] * n,
            "ortholog_overlap": [False] * n,
        }
    )


def test_pareto_mask_known_front():
    values = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1], [0.2, 0.2, 0.2], [0.1, 0.1, 0.1]])
    mask = pareto_mask(values)
    assert mask.tolist() == [True, True, True, True, False]


def test_rank_stability_is_deterministic():
    atlas = _atlas()
    config = {
        "iterations": 40,
        "pareto_iterations": 20,
        "random_seed": 3,
        "model_weight_center": [0.6, 0.25, 0.15],
        "model_weight_concentration": 50,
        "objective_weight_center": [0.5, 0.25, 0.25],
        "objective_weight_concentration": 50,
        "interval_noise_fraction": 0.1,
        "evidence_jitter_sd": 0.01,
        "reliability_jitter_sd": 0.01,
        "top_k": [3, 5, 10],
        "weight_scenarios": [
            {"name": "balanced", "repairability": 0.5, "evidence": 0.25, "reliability": 0.25}
        ],
    }
    first, _ = run_rank_stability(atlas, atlas.tail(5), atlas.tail(3), config)
    second, _ = run_rank_stability(atlas, atlas.tail(5), atlas.tail(3), config)
    pd.testing.assert_frame_equal(first, second)


def test_conservation_detects_offset_conserved_lysine():
    pairs = pd.DataFrame(
        [
            {
                "source_gene": "H",
                "target_gene": "M",
                "uniprot_accession_human": "P1",
                "uniprot_accession_mouse": "P2",
                "site_uniprot_human": 4,
                "site_uniprot_mouse": 5,
                "sequence_human": "AAAKBBB",
                "sequence_mouse": "XAAAKBBB",
                "repairability_extra_trees_human": 0.8,
                "repairability_extra_trees_mouse": 0.79,
            }
        ]
    )
    result = run_residue_conservation(pairs, {"local_window_radius": 2})
    assert bool(result.loc[0, "residue_level_conserved_lysine"])
    assert result.loc[0, "mapped_mouse_site_from_human"] == 5


def test_conservation_rejects_different_residue():
    pairs = pd.DataFrame(
        [
            {
                "source_gene": "H",
                "target_gene": "M",
                "uniprot_accession_human": "P1",
                "uniprot_accession_mouse": "P2",
                "site_uniprot_human": 4,
                "site_uniprot_mouse": 4,
                "sequence_human": "AAAKBBB",
                "sequence_mouse": "AAARBBB",
                "repairability_extra_trees_human": 0.8,
                "repairability_extra_trees_mouse": 0.79,
            }
        ]
    )
    result = run_residue_conservation(pairs, {"local_window_radius": 2})
    assert not bool(result.loc[0, "residue_level_conserved_lysine"])


def test_enrichment_single_driver_flag():
    frame = pd.DataFrame(
        [
            {
                "species": "mouse",
                "query_name": "portfolio",
                "source": "GO:BP",
                "term_id": "T1",
                "term_name": "term",
                "p_value_adjusted": 0.01,
                "significant": True,
                "query_size_mapped": 10,
                "intersection_size": 1,
                "term_size": 5,
                "effective_domain_size": 100,
                "intersection_genes": "G1",
            }
        ]
    )
    result = run_enrichment_sensitivity(
        frame, {"significant_only": True, "raw_p_threshold": 0.05, "redundancy_jaccard": 0.7}
    )
    assert result.loc[0, "driver_sensitivity_class"] == "single_driver"


def test_network_threshold_edges_decline():
    edges = pd.DataFrame(
        {
            "species": ["human"] * 3,
            "query_name": ["portfolio"] * 3,
            "node_a": ["A", "A", "B"],
            "node_b": ["B", "C", "C"],
            "combined_score": [0.71, 0.81, 0.91],
        }
    )
    nodes = pd.DataFrame(
        {"species": ["human"] * 3, "query_name": ["portfolio"] * 3, "node": ["A", "B", "C"]}
    )
    manifest = pd.DataFrame(
        {"species": ["human"], "query_name": ["portfolio"], "genes": ["A | B | C | D"]}
    )
    thresholds, _, _ = run_network_robustness(
        edges,
        nodes,
        manifest,
        {
            "thresholds": [0.7, 0.8, 0.9],
            "hub_removal_counts": [1],
            "rewiring_iterations": 2,
            "rewiring_swaps_per_edge": 1,
        },
        1,
    )
    assert thresholds["edges"].tolist() == [3, 2, 1]
    assert thresholds["nodes"].tolist() == [4, 4, 4]


def test_tier_watchlist_for_ood():
    atlas = _atlas(4)
    atlas.loc[0, "applicability_status"] = "out_of_domain"
    stability = pd.DataFrame(
        {
            "candidate_id": atlas["candidate_id"],
            "mean_stress_rank": [1, 2, 3, 4],
            "rank_iqr": [1] * 4,
            "top50_selection_frequency": [1] * 4,
            "top100_selection_frequency": [1] * 4,
            "portfolio_selection_frequency": [1] * 4,
            "pareto_membership_frequency": [1] * 4,
        }
    )
    result = run_candidate_tiers(
        atlas,
        stability,
        pd.DataFrame(),
        pd.DataFrame({"candidate_id": atlas["candidate_id"]}),
        {
            "tier1_portfolio_frequency": 0.8,
            "tier1_pareto_frequency": 0.4,
            "tier1_transfer_reliability": 0.5,
            "tier1_max_model_disagreement": 0.2,
            "tier2_portfolio_frequency": 0.6,
            "tier2_transfer_reliability": 0.4,
        },
    )
    assert result.loc[result["candidate_id"] == "c0", "candidate_tier"].item() == "Watchlist"


def test_specieswise_pareto_does_not_compare_across_species():
    values = np.array(
        [
            [0.9, 0.9, 0.9],  # human front
            [0.8, 0.8, 0.8],  # human dominated
            [0.2, 0.2, 0.2],  # mouse front despite lower absolute values
            [0.1, 0.1, 0.1],  # mouse dominated
        ]
    )
    species = np.array(["human", "human", "mouse", "mouse"])
    mask = specieswise_pareto_mask(values, species)
    assert mask.tolist() == [True, False, True, False]


def test_diversity_constrained_portfolio_enforces_species_and_protein_caps():
    priority = np.array([0.99, 0.98, 0.97, 0.96, 0.95, 0.94, 0.93, 0.92])
    first_front = np.array([True, True, False, False, True, True, False, False])
    species = np.array(["human"] * 4 + ["mouse"] * 4)
    proteins = np.array(["H1", "H1", "H1", "H2", "M1", "M1", "M1", "M2"])
    selected = select_diversity_constrained_indices(
        priority, first_front, species, proteins, top_n=4, max_sites_per_protein=1
    )
    selected_species = species[selected]
    selected_proteins = proteins[selected]
    assert (selected_species == "human").sum() == 2
    assert (selected_species == "mouse").sum() == 2
    assert len(set(zip(selected_species, selected_proteins))) == 4


def test_nondominated_front_ranks_known_layers():
    values = np.array(
        [
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [0.6, 0.6, 0.6],
            [0.5, 0.5, 0.5],
            [0.4, 0.4, 0.4],
        ]
    )
    ranks = nondominated_front_ranks(values)
    assert ranks.tolist() == [1, 1, 1, 1, 2, 3]


def test_complete_front_selection_precedes_high_priority_lower_front():
    priority = np.array([0.40, 0.99, 0.30, 0.98])
    front_ranks = np.array([1, 2, 1, 2])
    species = np.array(["human", "human", "mouse", "mouse"])
    proteins = np.array(["H1", "H2", "M1", "M2"])
    selected = select_diversity_constrained_indices(
        priority, front_ranks, species, proteins, top_n=2, max_sites_per_protein=1
    )
    assert set(selected.tolist()) == {0, 2}


def test_rank_stability_reports_exact_baseline_reproduction():
    atlas = _atlas(12)
    objectives = atlas[
        ["repairability_score_norm", "evidence_score_norm", "transfer_reliability_score"]
    ].to_numpy()
    fronts = nondominated_front_ranks(objectives)
    chosen = select_diversity_constrained_indices(
        atlas["weighted_priority_score"].to_numpy(),
        fronts,
        atlas["species"].to_numpy(),
        atlas["uniprot_accession"].to_numpy(),
        top_n=5,
        max_sites_per_protein=2,
    )
    portfolio = atlas.iloc[chosen].copy()
    config = {
        "iterations": 10,
        "pareto_iterations": 5,
        "random_seed": 3,
        "portfolio_top_n": 5,
        "max_sites_per_protein": 2,
        "model_weight_center": [0.6, 0.25, 0.15],
        "model_weight_concentration": 50,
        "objective_weight_center": [0.5, 0.25, 0.25],
        "objective_weight_concentration": 50,
        "interval_noise_fraction": 0.1,
        "evidence_jitter_sd": 0.01,
        "reliability_jitter_sd": 0.01,
        "top_k": [3, 5],
        "weight_scenarios": [
            {"name": "balanced", "repairability": 0.5, "evidence": 0.25, "reliability": 0.25}
        ],
    }
    _, weights = run_rank_stability(atlas, portfolio, atlas.loc[fronts == 1], config)
    assert weights["baseline_reproduction_count"].iloc[0] == 5
    assert weights["baseline_reproduction_fraction"].iloc[0] == 1.0
