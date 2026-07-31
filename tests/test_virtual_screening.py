import pandas as pd

from cml_repairnet.virtual_screening import select_multifidelity_batch


def test_multifidelity_selection_returns_ranked_unique_candidates() -> None:
    frame = pd.DataFrame(
        {
            "candidate_id": ["a", "b", "c", "d"],
            "pred_mean": [0.9, 0.2, 0.6, 0.5],
            "pred_std": [0.1, 0.8, 0.4, 0.5],
            "biological_priority": [0.8, 0.7, 0.2, 0.9],
            "ood_score": [0.1, 0.9, 0.2, 0.3],
            "compute_cost": [0.2, 0.1, 0.8, 0.4],
            "f1": [0.0, 0.0, 1.0, 1.0],
            "f2": [0.0, 1.0, 0.0, 1.0],
        }
    )
    selected = select_multifidelity_batch(frame, ["f1", "f2"], batch_size=3)
    assert len(selected) == 3
    assert selected["candidate_id"].is_unique
    assert selected["selection_rank"].tolist() == [1.0, 2.0, 3.0]
