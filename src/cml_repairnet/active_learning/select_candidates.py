from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import pairwise_distances
from sklearn.preprocessing import MinMaxScaler, StandardScaler


def select_diverse_uncertain_batch(
    candidates: pd.DataFrame,
    feature_columns: list[str],
    batch_size: int = 24,
    uncertainty_weight: float = 0.45,
    biological_priority_weight: float = 0.35,
    diversity_weight: float = 0.20,
) -> pd.DataFrame:
    """Backward-compatible selector for higher-fidelity computational candidates.

    Required columns: `pred_std` and `biological_priority`. This function does not select
    wet-lab experiments. New code should prefer `virtual_screening.select_multifidelity_batch`.
    """
    required = {"pred_std", "biological_priority"}
    missing = required.difference(candidates.columns)
    if missing:
        raise ValueError(f"Missing required candidate columns: {sorted(missing)}")
    if batch_size > len(candidates):
        raise ValueError("batch_size exceeds the number of candidates")

    frame = candidates.copy().reset_index(drop=True)
    frame[["uncertainty_scaled", "priority_scaled"]] = MinMaxScaler().fit_transform(
        frame[["pred_std", "biological_priority"]]
    )
    feature_matrix = StandardScaler().fit_transform(frame[feature_columns].fillna(0.0))
    distances = pairwise_distances(feature_matrix, metric="euclidean")

    selected: list[int] = []
    remaining = set(range(len(frame)))
    while len(selected) < batch_size:
        best_idx = None
        best_score = -np.inf
        for idx in remaining:
            diversity = 1.0 if not selected else float(np.min(distances[idx, selected]))
            score = (
                uncertainty_weight * frame.loc[idx, "uncertainty_scaled"]
                + biological_priority_weight * frame.loc[idx, "priority_scaled"]
                + diversity_weight * diversity
            )
            if score > best_score:
                best_score = score
                best_idx = idx
        assert best_idx is not None
        frame.loc[best_idx, "selection_score"] = best_score
        frame.loc[best_idx, "selection_rank"] = len(selected) + 1
        selected.append(best_idx)
        remaining.remove(best_idx)

    return frame.loc[selected].sort_values("selection_rank").reset_index(drop=True)


def run_selection(
    candidates_csv: Path,
    output_csv: Path,
    feature_columns: list[str],
    batch_size: int = 24,
) -> pd.DataFrame:
    candidates = pd.read_csv(candidates_csv)
    selected = select_diverse_uncertain_batch(candidates, feature_columns, batch_size=batch_size)
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    selected.to_csv(output_csv, index=False)
    return selected
