import numpy as np
import pandas as pd
import pytest

from cml_repairnet.models.stage25 import (
    CorrelationPruner,
    _choose_convex_weight,
    build_sample_weights,
    stage25_model_specs,
)


def test_correlation_pruner_is_fit_on_training_matrix() -> None:
    frame = pd.DataFrame(
        {
            "a": [0.0, 1.0, 2.0, 3.0],
            "a_duplicate": [0.0, 1.0, 2.0, 3.0],
            "b": [0.0, 1.0, 0.0, 1.0],
        }
    )
    pruner = CorrelationPruner(threshold=0.95).fit(frame)
    assert pruner.get_feature_names_out().tolist() == ["a", "b"]
    assert pruner.transform(frame).shape == (4, 2)


def test_sample_weights_are_clipped_and_normalized() -> None:
    frame = pd.DataFrame({"label_weight": [0.001, 1.0, 2.0, 1000.0]})
    for mode in ("none", "clipped", "sqrt"):
        weights = build_sample_weights(frame, mode)
        assert np.all(weights > 0)
        assert weights.mean() == pytest.approx(1.0)
    assert build_sample_weights(frame, "sqrt").max() < build_sample_weights(frame, "clipped").max()


def test_convex_weight_prefers_tree_when_tree_is_better() -> None:
    y = np.array([0.0, 0.5, 1.0])
    ridge = np.array([0.4, 0.4, 0.4])
    tree = np.array([0.0, 0.5, 1.0])
    assert _choose_convex_weight(y, ridge, tree) == 1.0


def test_stage25_specs_include_boosting_models() -> None:
    pytest.importorskip("xgboost")
    pytest.importorskip("lightgbm")
    specs = stage25_model_specs(seed=1, feature_mode="full")
    assert "xgboost" in specs
    assert "lightgbm" in specs
    assert specs["xgboost"].tree_model
    assert specs["lightgbm"].tree_model
