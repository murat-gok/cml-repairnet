from __future__ import annotations

import numpy as np
import pandas as pd
from sklearn.metrics import pairwise_distances
from sklearn.preprocessing import MinMaxScaler, StandardScaler


def _scaled_column(frame: pd.DataFrame, column: str, default: float) -> np.ndarray:
    if column not in frame.columns:
        return np.full(len(frame), default, dtype=float)
    values = frame[[column]].astype(float).replace([np.inf, -np.inf], np.nan)
    values = values.fillna(values.median()).fillna(default)
    if len(frame) == 1 or float(values.max().iloc[0]) == float(values.min().iloc[0]):
        return np.full(len(frame), default, dtype=float)
    return MinMaxScaler().fit_transform(values).ravel()


def select_multifidelity_batch(
    candidates: pd.DataFrame,
    feature_columns: list[str],
    batch_size: int = 24,
    repairability_weight: float = 0.20,
    uncertainty_weight: float = 0.25,
    biological_priority_weight: float = 0.20,
    diversity_weight: float = 0.20,
    ood_penalty_weight: float = 0.10,
    compute_cost_penalty_weight: float = 0.05,
) -> pd.DataFrame:
    """Select candidates for higher-fidelity *computational* evaluation.

    Expected columns are `pred_mean`, `pred_std`, `biological_priority`, `ood_score`, and
    `compute_cost`. Missing optional columns receive neutral defaults. The function does not
    select wet-lab experiments and does not treat any computational score as ground truth.
    """
    if not feature_columns:
        raise ValueError("feature_columns must not be empty")
    missing_features = [column for column in feature_columns if column not in candidates.columns]
    if missing_features:
        raise ValueError(f"Missing feature columns: {missing_features}")
    if batch_size < 1 or batch_size > len(candidates):
        raise ValueError("batch_size must be between 1 and the number of candidates")

    frame = candidates.copy().reset_index(drop=True)
    frame["repairability_scaled"] = _scaled_column(frame, "pred_mean", 0.5)
    frame["uncertainty_scaled"] = _scaled_column(frame, "pred_std", 0.5)
    frame["priority_scaled"] = _scaled_column(frame, "biological_priority", 0.5)
    frame["ood_scaled"] = _scaled_column(frame, "ood_score", 0.0)
    frame["cost_scaled"] = _scaled_column(frame, "compute_cost", 0.0)

    numeric_features = frame[feature_columns].apply(pd.to_numeric, errors="coerce")
    numeric_features = numeric_features.fillna(numeric_features.median()).fillna(0.0)
    feature_matrix = StandardScaler().fit_transform(numeric_features)
    distances = pairwise_distances(feature_matrix, metric="euclidean")
    if distances.max() > 0:
        distances = distances / distances.max()

    selected: list[int] = []
    remaining = set(range(len(frame)))
    while len(selected) < batch_size:
        best_idx: int | None = None
        best_score = -np.inf
        for idx in remaining:
            diversity = 1.0 if not selected else float(np.min(distances[idx, selected]))
            score = (
                repairability_weight * frame.loc[idx, "repairability_scaled"]
                + uncertainty_weight * frame.loc[idx, "uncertainty_scaled"]
                + biological_priority_weight * frame.loc[idx, "priority_scaled"]
                + diversity_weight * diversity
                - ood_penalty_weight * frame.loc[idx, "ood_scaled"]
                - compute_cost_penalty_weight * frame.loc[idx, "cost_scaled"]
            )
            if score > best_score:
                best_score = score
                best_idx = idx
        if best_idx is None:
            raise RuntimeError("Candidate selection failed")
        frame.loc[best_idx, "selection_score"] = best_score
        frame.loc[best_idx, "selection_rank"] = len(selected) + 1
        selected.append(best_idx)
        remaining.remove(best_idx)

    return frame.loc[selected].sort_values("selection_rank").reset_index(drop=True)
