from __future__ import annotations

import math

import numpy as np
from scipy.stats import spearmanr
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score


def regression_metrics(y_true: np.ndarray, y_pred: np.ndarray) -> dict[str, float]:
    rho = spearmanr(y_true, y_pred, nan_policy="omit").statistic
    return {
        "mae": float(mean_absolute_error(y_true, y_pred)),
        "rmse": float(math.sqrt(mean_squared_error(y_true, y_pred))),
        "spearman_rho": float(rho) if np.isfinite(rho) else float("nan"),
        "r2": float(r2_score(y_true, y_pred)),
    }


def bootstrap_metric_ci(
    y_true: np.ndarray,
    y_pred: np.ndarray,
    metric: str = "mae",
    iterations: int = 2000,
    seed: int = 20260727,
) -> tuple[float, float]:
    rng = np.random.default_rng(seed)
    values: list[float] = []
    n = len(y_true)
    for _ in range(iterations):
        idx = rng.integers(0, n, size=n)
        metrics = regression_metrics(y_true[idx], y_pred[idx])
        value = metrics[metric]
        if np.isfinite(value):
            values.append(value)
    if not values:
        return float("nan"), float("nan")
    return float(np.quantile(values, 0.025)), float(np.quantile(values, 0.975))
