import pandas as pd

from cml_repairnet.active_learning.select_candidates import select_diverse_uncertain_batch


def test_selection_returns_requested_batch() -> None:
    candidates = pd.DataFrame(
        {
            "pred_std": [0.1, 0.9, 0.4, 0.7],
            "biological_priority": [0.9, 0.1, 0.5, 0.6],
            "f1": [0.0, 1.0, 0.2, 0.8],
            "f2": [0.0, 1.0, 0.9, 0.1],
        }
    )
    selected = select_diverse_uncertain_batch(candidates, ["f1", "f2"], batch_size=2)
    assert len(selected) == 2
    assert selected["selection_rank"].tolist() == [1.0, 2.0]
