import numpy as np

from cml_repairnet.evaluation.metrics import regression_metrics


def test_perfect_regression_metrics() -> None:
    y = np.array([0.0, 0.5, 1.0])
    metrics = regression_metrics(y, y)
    assert metrics["mae"] == 0.0
    assert metrics["rmse"] == 0.0
    assert metrics["spearman_rho"] == 1.0
