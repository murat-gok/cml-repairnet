from __future__ import annotations

import pandas as pd

from cml_repairnet.revision.external_scaffold import (
    QualityThresholds,
    add_quality_flags,
    aggregate_predictions,
    compare_rankings,
    observed_summary,
    select_mature_range,
)


def test_select_mature_range_prefers_longest_chain() -> None:
    record = {
        "sequence": {"value": "M" * 120},
        "features": [
            {"type": "Signal", "location": {"start": {"value": 1}, "end": {"value": 15}}},
            {"type": "Chain", "description": "mature", "location": {"start": {"value": 16}, "end": {"value": 110}}},
            {"type": "Chain", "description": "small", "location": {"start": {"value": 20}, "end": {"value": 40}}},
        ],
    }
    assert select_mature_range(record)[:2] == (16, 110)


def synthetic_sites() -> pd.DataFrame:
    rows = []
    scores = {"casein": 0.9, "hemoglobin": 0.8, "collagen": 0.6}
    for scaffold, score in scores.items():
        for component_idx in range(2):
            for site in range(4):
                rows.append(
                    {
                        "scaffold": scaffold,
                        "component": f"c{component_idx}",
                        "site_uniprot": site + 1,
                        "structure_missing": False,
                        "site_plddt": 90.0,
                        "neighbor_mean_plddt_8a": 90.0,
                        "ood_percentile": 0.2,
                        "relative_sasa_lysine": 0.5 + site / 10,
                        "stoichiometric_weight": 0.5,
                        "equal_component_weight": 0.5,
                        "pred_extra_trees": score,
                        "pred_lightgbm": score - 0.01,
                        "pred_ridge": score - 0.02,
                        "pred_consensus": score - 0.01,
                    }
                )
    frame = pd.DataFrame(rows)
    return add_quality_flags(frame, QualityThresholds())


def source_labels() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "scaffold": ["casein"] * 3 + ["hemoglobin"] * 3 + ["collagen"] * 3,
            "repair_fraction": [0.92, 0.97, 1.34, 0.82, 0.87, 0.89, 0.60, 0.65, 0.80],
        }
    )


def test_observed_order_is_stable() -> None:
    summary = observed_summary(source_labels())
    assert summary["scaffold"].tolist() == ["casein", "hemoglobin", "collagen"]


def test_aggregation_and_rank_comparison() -> None:
    aggregates = aggregate_predictions(synthetic_sites())
    rankings = compare_rankings(aggregates, source_labels())
    primary = rankings[
        (rankings["quality_mode"] == "strict")
        & (rankings["domain_mode"] == "in_domain")
        & (rankings["component_scheme"] == "stoichiometric")
        & (rankings["model"] == "extra_trees")
        & (rankings["aggregation"] == "sasa_weighted_mean")
    ].iloc[0]
    assert bool(primary["exact_rank_match"])
    assert primary["pairwise_order_accuracy"] == 1.0
