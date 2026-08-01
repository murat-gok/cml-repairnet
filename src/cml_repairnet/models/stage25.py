from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Sequence

import numpy as np
import pandas as pd
from sklearn.base import BaseEstimator, TransformerMixin, clone
from sklearn.dummy import DummyRegressor
from sklearn.ensemble import ExtraTreesRegressor, RandomForestRegressor
from sklearn.impute import SimpleImputer
from sklearn.inspection import permutation_importance
from sklearn.linear_model import ElasticNet, Ridge
from sklearn.model_selection import GridSearchCV, GroupKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

from cml_repairnet.evaluation.metrics import bootstrap_metric_ci, regression_metrics
from cml_repairnet.models.baselines import _numeric_feature_columns, _spatial_groups


@dataclass(frozen=True)
class Stage25ModelSpec:
    estimator: Any
    param_grid: dict[str, list[Any]] | list[dict[str, list[Any]]]
    supports_sample_weight: bool = True
    tree_model: bool = False


class CorrelationPruner(BaseEstimator, TransformerMixin):
    """Drop later features whose absolute training-fold correlation exceeds a threshold.

    The transformer is deliberately placed inside the sklearn Pipeline, so correlations are
    estimated on each training fold only. A threshold greater than 1 keeps every feature.
    """

    def __init__(self, threshold: float = 1.01):
        self.threshold = threshold

    def fit(self, X: Any, y: Any = None) -> "CorrelationPruner":
        frame = self._as_frame(X)
        self.feature_names_in_ = np.asarray(frame.columns, dtype=object)
        if not 0.0 < float(self.threshold) <= 1.01:
            raise ValueError("threshold must be in (0, 1.01]")

        numeric = frame.apply(pd.to_numeric, errors="coerce")
        corr = numeric.corr(method="pearson").abs().fillna(0.0)
        keep: list[int] = []
        for idx in range(corr.shape[0]):
            if all(corr.iat[idx, kept_idx] < float(self.threshold) for kept_idx in keep):
                keep.append(idx)
        if not keep:
            keep = list(range(frame.shape[1]))
        self.keep_indices_ = np.asarray(keep, dtype=int)
        self.selected_feature_names_ = self.feature_names_in_[self.keep_indices_]
        return self

    def transform(self, X: Any) -> Any:
        if not hasattr(self, "keep_indices_"):
            raise RuntimeError("CorrelationPruner must be fitted before transform")
        array = X.to_numpy() if isinstance(X, pd.DataFrame) else np.asarray(X)
        return array[:, self.keep_indices_]

    def get_feature_names_out(self, input_features: Sequence[str] | None = None) -> np.ndarray:
        if not hasattr(self, "selected_feature_names_"):
            raise RuntimeError("CorrelationPruner must be fitted before get_feature_names_out")
        return np.asarray(self.selected_feature_names_, dtype=object)

    @staticmethod
    def _as_frame(X: Any) -> pd.DataFrame:
        if isinstance(X, pd.DataFrame):
            return X.copy()
        array = np.asarray(X)
        return pd.DataFrame(array, columns=[f"x{i}" for i in range(array.shape[1])])


def _single_value_grid(**params: Any) -> dict[str, list[Any]]:
    return {key: [value] for key, value in params.items()}


def _xgboost_spec(seed: int, feature_mode: str) -> Stage25ModelSpec:
    try:
        from xgboost import XGBRegressor
    except ImportError as exc:  # pragma: no cover - exercised through CLI dependency check
        raise RuntimeError(
            "XGBoost is not installed. Run: python -m pip install -e '.[experiment]'"
        ) from exc

    estimator = XGBRegressor(
        objective="reg:squarederror",
        tree_method="hist",
        random_state=seed,
        n_jobs=1,
        verbosity=0,
    )
    pipeline = _tree_pipeline(estimator, feature_mode)
    candidates = [
        dict(n_estimators=60, max_depth=1, learning_rate=0.08, min_child_weight=3,
             subsample=1.0, colsample_bytree=1.0, reg_alpha=0.0, reg_lambda=5.0),
        dict(n_estimators=120, max_depth=1, learning_rate=0.04, min_child_weight=3,
             subsample=0.8, colsample_bytree=0.8, reg_alpha=0.25, reg_lambda=10.0),
        dict(n_estimators=200, max_depth=1, learning_rate=0.025, min_child_weight=5,
             subsample=0.8, colsample_bytree=1.0, reg_alpha=0.5, reg_lambda=15.0),
        dict(n_estimators=60, max_depth=2, learning_rate=0.06, min_child_weight=3,
             subsample=0.8, colsample_bytree=0.8, reg_alpha=0.25, reg_lambda=10.0),
        dict(n_estimators=120, max_depth=2, learning_rate=0.035, min_child_weight=5,
             subsample=0.75, colsample_bytree=1.0, reg_alpha=0.5, reg_lambda=15.0),
        dict(n_estimators=180, max_depth=2, learning_rate=0.025, min_child_weight=7,
             subsample=0.75, colsample_bytree=0.8, reg_alpha=1.0, reg_lambda=20.0),
    ]
    grids = [
        _single_value_grid(**{f"model__{key}": value for key, value in candidate.items()})
        for candidate in candidates
    ]
    _add_pruner_grid(grids, feature_mode)
    return Stage25ModelSpec(pipeline, grids, tree_model=True)


def _lightgbm_spec(seed: int, feature_mode: str) -> Stage25ModelSpec:
    try:
        from lightgbm import LGBMRegressor
    except ImportError as exc:  # pragma: no cover - exercised through CLI dependency check
        raise RuntimeError(
            "LightGBM is not installed. Run: python -m pip install -e '.[experiment]'"
        ) from exc

    estimator = LGBMRegressor(
        objective="regression_l1",
        random_state=seed,
        n_jobs=1,
        verbosity=-1,
        deterministic=True,
        force_col_wise=True,
        feature_pre_filter=False,
        subsample_freq=1,
    )
    pipeline = _tree_pipeline(estimator, feature_mode)
    candidates = [
        dict(n_estimators=60, learning_rate=0.06, num_leaves=2, max_depth=1,
             min_child_samples=4, subsample=1.0, colsample_bytree=1.0,
             reg_alpha=0.0, reg_lambda=5.0),
        dict(n_estimators=120, learning_rate=0.035, num_leaves=2, max_depth=1,
             min_child_samples=5, subsample=0.8, colsample_bytree=0.8,
             reg_alpha=0.25, reg_lambda=10.0),
        dict(n_estimators=200, learning_rate=0.02, num_leaves=2, max_depth=1,
             min_child_samples=6, subsample=0.8, colsample_bytree=1.0,
             reg_alpha=0.5, reg_lambda=15.0),
        dict(n_estimators=60, learning_rate=0.05, num_leaves=3, max_depth=2,
             min_child_samples=4, subsample=0.8, colsample_bytree=0.8,
             reg_alpha=0.25, reg_lambda=10.0),
        dict(n_estimators=120, learning_rate=0.03, num_leaves=3, max_depth=2,
             min_child_samples=5, subsample=0.75, colsample_bytree=1.0,
             reg_alpha=0.5, reg_lambda=15.0),
        dict(n_estimators=180, learning_rate=0.02, num_leaves=4, max_depth=2,
             min_child_samples=6, subsample=0.75, colsample_bytree=0.8,
             reg_alpha=1.0, reg_lambda=20.0),
    ]
    grids = [
        _single_value_grid(**{f"model__{key}": value for key, value in candidate.items()})
        for candidate in candidates
    ]
    _add_pruner_grid(grids, feature_mode)
    return Stage25ModelSpec(pipeline, grids, tree_model=True)


def _scaled_pipeline(model: Any, feature_mode: str) -> Pipeline:
    return Pipeline(
        [
            ("prune", CorrelationPruner(_fixed_pruner_threshold(feature_mode))),
            ("imputer", SimpleImputer(strategy="median")),
            ("scale", StandardScaler()),
            ("model", model),
        ]
    )


def _tree_pipeline(model: Any, feature_mode: str) -> Pipeline:
    return Pipeline(
        [
            ("prune", CorrelationPruner(_fixed_pruner_threshold(feature_mode))),
            ("imputer", SimpleImputer(strategy="median")),
            ("model", model),
        ]
    )


def _fixed_pruner_threshold(feature_mode: str) -> float:
    if feature_mode == "full":
        return 1.01
    if feature_mode == "nested_pruned":
        return 0.95
    raise ValueError("feature_mode must be 'full' or 'nested_pruned'")


def _add_pruner_grid(
    grids: list[dict[str, list[Any]]], feature_mode: str
) -> None:
    if feature_mode == "nested_pruned":
        expanded: list[dict[str, list[Any]]] = []
        for grid in grids:
            for threshold in (0.85, 0.90, 0.95):
                expanded.append({**grid, "prune__threshold": [threshold]})
        grids[:] = expanded


def stage25_model_specs(seed: int, feature_mode: str = "full") -> dict[str, Stage25ModelSpec]:
    ridge_grid: dict[str, list[Any]] = {
        "model__alpha": [0.01, 0.1, 1.0, 10.0, 100.0]
    }
    elastic_grid: dict[str, list[Any]] = {
        "model__alpha": [0.001, 0.01, 0.1, 1.0],
        "model__l1_ratio": [0.1, 0.5, 0.9],
    }
    forest_candidates = [
        dict(max_depth=2, min_samples_leaf=1, max_features=1.0),
        dict(max_depth=2, min_samples_leaf=2, max_features=0.8),
        dict(max_depth=4, min_samples_leaf=1, max_features=0.8),
        dict(max_depth=4, min_samples_leaf=2, max_features=1.0),
        dict(max_depth=None, min_samples_leaf=2, max_features=0.8),
        dict(max_depth=None, min_samples_leaf=4, max_features=1.0),
    ]
    rf_grids: list[dict[str, list[Any]]] = [
        _single_value_grid(
            model__n_estimators=500,
            **{f"model__{key}": value for key, value in candidate.items()},
        )
        for candidate in forest_candidates
    ]
    et_grids: list[dict[str, list[Any]]] = [
        _single_value_grid(
            model__n_estimators=500,
            **{f"model__{key}": value for key, value in candidate.items()},
        )
        for candidate in forest_candidates
    ]
    if feature_mode == "nested_pruned":
        for grid in (ridge_grid, elastic_grid):
            grid["prune__threshold"] = [0.85, 0.90, 0.95]
        _add_pruner_grid(rf_grids, feature_mode)
        _add_pruner_grid(et_grids, feature_mode)

    return {
        "dummy_mean": Stage25ModelSpec(
            _tree_pipeline(DummyRegressor(strategy="mean"), feature_mode), {}
        ),
        "ridge": Stage25ModelSpec(_scaled_pipeline(Ridge(), feature_mode), ridge_grid),
        "elastic_net": Stage25ModelSpec(
            _scaled_pipeline(ElasticNet(max_iter=30_000, random_state=seed), feature_mode),
            elastic_grid,
        ),
        "random_forest": Stage25ModelSpec(
            _tree_pipeline(
                RandomForestRegressor(random_state=seed, n_jobs=1), feature_mode
            ),
            rf_grids,
            tree_model=True,
        ),
        "extra_trees": Stage25ModelSpec(
            _tree_pipeline(ExtraTreesRegressor(random_state=seed, n_jobs=1), feature_mode),
            et_grids,
            tree_model=True,
        ),
        "xgboost": _xgboost_spec(seed, feature_mode),
        "lightgbm": _lightgbm_spec(seed, feature_mode),
    }


def build_sample_weights(frame: pd.DataFrame, mode: str) -> np.ndarray:
    """Return non-negative weights normalized to mean one.

    label_weight is inverse-variance-like and can be extremely skewed. The two weighted modes
    therefore clip it within the current outer-training subset before normalization. Outer-test
    observations are never used to set clipping bounds.
    """

    mode = mode.lower()
    if mode == "none":
        return np.ones(len(frame), dtype=float)
    if "label_weight" not in frame.columns:
        raise ValueError("label_weight is required for weighted learning")

    raw = frame["label_weight"].to_numpy(dtype=float)
    finite_positive = raw[np.isfinite(raw) & (raw > 0)]
    if finite_positive.size == 0:
        return np.ones(len(frame), dtype=float)
    replacement = float(np.median(finite_positive))
    raw = np.where(np.isfinite(raw) & (raw > 0), raw, replacement)
    lower, upper = np.quantile(raw, [0.10, 0.90])
    clipped = np.clip(raw, lower, upper)
    if mode == "clipped":
        weights = clipped
    elif mode == "sqrt":
        weights = np.sqrt(clipped)
    else:
        raise ValueError("weight mode must be one of: none, clipped, sqrt")
    mean = float(np.mean(weights))
    return weights / mean if mean > 0 else np.ones(len(frame), dtype=float)


def _fit_params(spec: Stage25ModelSpec, weights: np.ndarray) -> dict[str, np.ndarray]:
    if not spec.supports_sample_weight:
        return {}
    return {"model__sample_weight": weights}


def _selected_feature_names(fitted: Pipeline, original_features: Sequence[str]) -> list[str]:
    pruner = fitted.named_steps.get("prune")
    if pruner is None:
        return list(original_features)
    try:
        return [str(value) for value in pruner.get_feature_names_out(original_features)]
    except Exception:
        return list(original_features)


def _tune_estimator(
    spec: Stage25ModelSpec,
    X: pd.DataFrame,
    y: np.ndarray,
    groups: np.ndarray,
    weights: np.ndarray,
    n_jobs: int,
) -> tuple[Any, dict[str, Any]]:
    unique_groups = np.unique(groups)
    if spec.param_grid and len(unique_groups) >= 2:
        inner_splits = min(3, len(unique_groups))
        search = GridSearchCV(
            estimator=spec.estimator,
            param_grid=spec.param_grid,
            scoring="neg_mean_absolute_error",
            cv=GroupKFold(n_splits=inner_splits),
            n_jobs=n_jobs,
            refit=True,
            error_score="raise",
        )
        search.fit(X, y, groups=groups, **_fit_params(spec, weights))
        return search.best_estimator_, dict(search.best_params_)
    fitted = clone(spec.estimator).fit(X, y, **_fit_params(spec, weights))
    return fitted, {}


def _inner_oof_predictions(
    estimator: Any,
    X: pd.DataFrame,
    y: np.ndarray,
    groups: np.ndarray,
    weights: np.ndarray,
) -> np.ndarray:
    unique_groups = np.unique(groups)
    inner_splits = min(3, len(unique_groups))
    cv = GroupKFold(n_splits=inner_splits)
    predictions = np.full(len(y), np.nan, dtype=float)
    for train_idx, valid_idx in cv.split(X, y, groups):
        model = clone(estimator)
        model.fit(
            X.iloc[train_idx],
            y[train_idx],
            model__sample_weight=weights[train_idx],
        )
        predictions[valid_idx] = model.predict(X.iloc[valid_idx])
    return np.clip(predictions, 0.0, 1.0)


def _choose_convex_weight(y: np.ndarray, ridge_pred: np.ndarray, tree_pred: np.ndarray) -> float:
    candidates = np.linspace(0.0, 1.0, 21)
    errors = [
        np.mean(np.abs(y - (alpha * tree_pred + (1.0 - alpha) * ridge_pred)))
        for alpha in candidates
    ]
    return float(candidates[int(np.argmin(errors))])


def run_stage25_nested_cv(
    merged: pd.DataFrame,
    target: str,
    output_dir: Path,
    *,
    n_spatial_groups: int = 4,
    weight_mode: str = "none",
    feature_mode: str = "full",
    bootstrap_iterations: int = 5000,
    permutation_repeats: int = 50,
    n_jobs: int = 1,
    seed: int = 20260727,
    model_names: Iterable[str] | None = None,
    include_convex_ensemble: bool = True,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    feature_columns = _numeric_feature_columns(merged)
    if not feature_columns:
        raise ValueError("No numeric feature columns were found")
    X = merged[feature_columns].copy()
    y = merged[target].to_numpy(dtype=float)
    groups = _spatial_groups(merged, n_spatial_groups)
    outer_cv = GroupKFold(n_splits=n_spatial_groups)
    specs = stage25_model_specs(seed, feature_mode)
    selected_names = list(model_names) if model_names is not None else list(specs)
    unknown = sorted(set(selected_names).difference(specs))
    if unknown:
        raise ValueError(f"Unknown Stage 2.5 models: {unknown}")

    predictions_by_model = {
        name: np.full(len(y), np.nan, dtype=float) for name in selected_names
    }
    fold_by_site = np.full(len(y), -1, dtype=int)
    params_by_model: dict[str, list[dict[str, Any]]] = {name: [] for name in selected_names}
    selected_features_by_model: dict[str, list[list[str]]] = {
        name: [] for name in selected_names
    }
    importance_rows: list[dict[str, Any]] = []

    ensemble_predictions = np.full(len(y), np.nan, dtype=float)
    ensemble_alphas: list[float] = []
    ensemble_enabled = include_convex_ensemble and {"ridge", "extra_trees"}.issubset(
        selected_names
    )

    for fold, (train_idx, test_idx) in enumerate(outer_cv.split(X, y, groups), start=1):
        fold_by_site[test_idx] = fold
        train_frame = merged.iloc[train_idx]
        train_weights = build_sample_weights(train_frame, weight_mode)
        fitted_models: dict[str, Any] = {}

        for model_name in selected_names:
            spec = specs[model_name]
            fitted, best_params = _tune_estimator(
                spec,
                X.iloc[train_idx],
                y[train_idx],
                groups[train_idx],
                train_weights,
                n_jobs,
            )
            fitted_models[model_name] = fitted
            params_by_model[model_name].append(best_params)
            selected_features_by_model[model_name].append(
                _selected_feature_names(fitted, feature_columns)
            )
            predictions_by_model[model_name][test_idx] = np.clip(
                fitted.predict(X.iloc[test_idx]), 0.0, 1.0
            )

            if model_name != "dummy_mean" and permutation_repeats > 0:
                importance = permutation_importance(
                    fitted,
                    X.iloc[test_idx],
                    y[test_idx],
                    scoring="neg_mean_absolute_error",
                    n_repeats=permutation_repeats,
                    random_state=seed + fold,
                    n_jobs=1,
                )
                for feature, mean, std in zip(
                    feature_columns, importance.importances_mean, importance.importances_std
                ):
                    importance_rows.append(
                        {
                            "fold": fold,
                            "model": model_name,
                            "feature": feature,
                            "importance_mean": float(mean),
                            "importance_std": float(std),
                            "target": target,
                            "spatial_groups": n_spatial_groups,
                            "weight_mode": weight_mode,
                            "feature_mode": feature_mode,
                        }
                    )

        if ensemble_enabled:
            ridge_spec = specs["ridge"]
            tree_spec = specs["extra_trees"]
            ridge_best = fitted_models["ridge"]
            tree_best = fitted_models["extra_trees"]
            ridge_inner = _inner_oof_predictions(
                ridge_best,
                X.iloc[train_idx].reset_index(drop=True),
                y[train_idx],
                groups[train_idx],
                train_weights,
            )
            tree_inner = _inner_oof_predictions(
                tree_best,
                X.iloc[train_idx].reset_index(drop=True),
                y[train_idx],
                groups[train_idx],
                train_weights,
            )
            alpha = _choose_convex_weight(y[train_idx], ridge_inner, tree_inner)
            ensemble_alphas.append(alpha)
            ridge_test = predictions_by_model["ridge"][test_idx]
            tree_test = predictions_by_model["extra_trees"][test_idx]
            ensemble_predictions[test_idx] = np.clip(
                alpha * tree_test + (1.0 - alpha) * ridge_test, 0.0, 1.0
            )

    metric_rows: list[dict[str, Any]] = []
    prediction_frames: list[pd.DataFrame] = []
    for model_name, predictions in predictions_by_model.items():
        metrics = regression_metrics(y, predictions)
        ci_low, ci_high = bootstrap_metric_ci(
            y,
            predictions,
            metric="mae",
            iterations=bootstrap_iterations,
            seed=seed,
        )
        selected_counts = [len(values) for values in selected_features_by_model[model_name]]
        metric_rows.append(
            {
                "model": model_name,
                **metrics,
                "mae_ci_2.5": ci_low,
                "mae_ci_97.5": ci_high,
                "n_sites": len(y),
                "n_input_features": len(feature_columns),
                "mean_selected_features": float(np.mean(selected_counts)),
                "target": target,
                "spatial_groups": n_spatial_groups,
                "weight_mode": weight_mode,
                "feature_mode": feature_mode,
                "validation": f"{n_spatial_groups}-fold spatial GroupKFold with nested tuning",
            }
        )
        prediction_frames.append(
            pd.DataFrame(
                {
                    "protein_id": merged["protein_id"],
                    "site_uniprot": merged["site_uniprot"],
                    "site_paper": merged["site_paper"],
                    "spatial_group": groups,
                    "outer_fold": fold_by_site,
                    "model": model_name,
                    "y_true": y,
                    "y_pred": predictions,
                    "target": target,
                    "weight_mode": weight_mode,
                    "feature_mode": feature_mode,
                }
            )
        )

    if ensemble_enabled:
        model_name = "ridge_extra_trees_convex"
        metrics = regression_metrics(y, ensemble_predictions)
        ci_low, ci_high = bootstrap_metric_ci(
            y,
            ensemble_predictions,
            metric="mae",
            iterations=bootstrap_iterations,
            seed=seed,
        )
        metric_rows.append(
            {
                "model": model_name,
                **metrics,
                "mae_ci_2.5": ci_low,
                "mae_ci_97.5": ci_high,
                "n_sites": len(y),
                "n_input_features": len(feature_columns),
                "mean_selected_features": float("nan"),
                "target": target,
                "spatial_groups": n_spatial_groups,
                "weight_mode": weight_mode,
                "feature_mode": feature_mode,
                "validation": (
                    f"{n_spatial_groups}-fold spatial GroupKFold; convex weight selected "
                    "inside each outer training fold"
                ),
                "mean_extra_trees_weight": float(np.mean(ensemble_alphas)),
            }
        )
        prediction_frames.append(
            pd.DataFrame(
                {
                    "protein_id": merged["protein_id"],
                    "site_uniprot": merged["site_uniprot"],
                    "site_paper": merged["site_paper"],
                    "spatial_group": groups,
                    "outer_fold": fold_by_site,
                    "model": model_name,
                    "y_true": y,
                    "y_pred": ensemble_predictions,
                    "target": target,
                    "weight_mode": weight_mode,
                    "feature_mode": feature_mode,
                }
            )
        )

    metrics_frame = pd.DataFrame(metric_rows).sort_values("mae").reset_index(drop=True)
    predictions_frame = pd.concat(prediction_frames, ignore_index=True)
    importance_frame = pd.DataFrame(importance_rows)

    output_dir.mkdir(parents=True, exist_ok=True)
    metrics_frame.to_csv(output_dir / "stage25_metrics.csv", index=False)
    predictions_frame.to_csv(output_dir / "stage25_oof_predictions.csv", index=False)
    importance_frame.to_csv(output_dir / "stage25_permutation_importance.csv", index=False)
    (output_dir / "stage25_input_features.txt").write_text(
        "\n".join(feature_columns), encoding="utf-8"
    )
    (output_dir / "stage25_fold_details.json").write_text(
        json.dumps(
            {
                "best_params": params_by_model,
                "selected_features": selected_features_by_model,
                "ensemble_extra_trees_weights": ensemble_alphas,
            },
            indent=2,
            sort_keys=True,
        ),
        encoding="utf-8",
    )
    return metrics_frame, predictions_frame, importance_frame


def aggregate_permutation_importance(frame: pd.DataFrame) -> pd.DataFrame:
    if frame.empty:
        return frame.copy()
    return (
        frame.groupby(["model", "feature"], as_index=False)
        .agg(
            importance_mean=("importance_mean", "mean"),
            importance_sd_across_folds=("importance_mean", "std"),
            mean_within_fold_sd=("importance_std", "mean"),
            folds=("fold", "nunique"),
        )
        .sort_values(["model", "importance_mean"], ascending=[True, False])
        .reset_index(drop=True)
    )
