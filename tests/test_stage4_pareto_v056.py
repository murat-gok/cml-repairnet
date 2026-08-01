from __future__ import annotations

import pandas as pd

from cml_repairnet.atlas.prioritize import add_pareto_ranking, add_priority_components


def test_pareto_excludes_domain_shift_and_uses_transfer_reliability() -> None:
    frame = pd.DataFrame(
        {
            "candidate_id": ["h1", "h2", "h3", "m1", "m2"],
            "species": ["human", "human", "human", "mouse", "mouse"],
            "uniprot_accession": ["P1", "P2", "P3", "Q1", "Q2"],
            "repairability_extra_trees": [0.9, 0.7, 0.99, 0.8, 0.6],
            "repairability_lightgbm": [0.85, 0.65, 0.95, 0.75, 0.55],
            "repairability_ridge": [0.8, 0.6, 0.9, 0.7, 0.5],
            "consensus_rank_percentile": [1.0, 0.5, 1.0, 1.0, 0.5],
            "evidence_score_raw": [3.0, 2.0, 4.0, 3.0, 1.0],
            "structure_quality_score": [0.95, 0.9, 0.95, 0.9, 0.8],
            "ood_percentile": [0.1, 0.2, 1.0, 0.1, 0.4],
            "prediction_interval_width": [0.4, 0.5, 0.3, 0.4, 0.6],
            "model_disagreement": [0.03, 0.05, 0.01, 0.04, 0.06],
            "applicability_status": [
                "in_domain",
                "in_domain",
                "out_of_domain",
                "in_domain",
                "in_domain",
            ],
        }
    )
    ranked = add_pareto_ranking(add_priority_components(frame))
    assert "transfer_reliability_score" in ranked.columns
    assert ranked.loc[ranked["candidate_id"] == "h3", "pareto_front"].item() == 0
    assert (ranked.loc[ranked["pareto_front"] == 1, "applicability_status"] == "in_domain").all()
