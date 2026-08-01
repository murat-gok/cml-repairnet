from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.cluster import AgglomerativeClustering
from sklearn.compose import TransformedTargetRegressor
from sklearn.dummy import DummyRegressor
from sklearn.ensemble import ExtraTreesRegressor, RandomForestRegressor
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import ConstantKernel, Matern, WhiteKernel
from sklearn.impute import SimpleImputer
from sklearn.linear_model import ElasticNet, Ridge
from sklearn.model_selection import GridSearchCV, GroupKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVR

from cml_repairnet.evaluation.metrics import bootstrap_metric_ci, regression_metrics


@dataclass
class ModelSpec:
    estimator: Any
    param_grid: dict[str, list[Any]]


def _model_specs(seed: int) -> dict[str, ModelSpec]:
    scaled = lambda model: Pipeline(
        [("imputer", SimpleImputer(strategy="median")), ("scale", StandardScaler()), ("model", model)]
    )
    tree = lambda model: Pipeline(
        [("imputer", SimpleImputer(strategy="median")), ("model", model)]
    )
    return {
        "dummy_mean": ModelSpec(Pipeline([("model", DummyRegressor(strategy="mean"))]), {}),
        "ridge": ModelSpec(scaled(Ridge()), {"model__alpha": [0.01, 0.1, 1.0, 10.0, 100.0]}),
        "elastic_net": ModelSpec(
            scaled(ElasticNet(max_iter=20000, random_state=seed)),
            {"model__alpha": [0.001, 0.01, 0.1, 1.0], "model__l1_ratio": [0.1, 0.5, 0.9]},
        ),
        "svr_rbf": ModelSpec(
            scaled(SVR(kernel="rbf")),
            {"model__C": [0.1, 1.0, 10.0, 100.0], "model__gamma": ["scale", 0.01, 0.1, 1.0]},
        ),
        "random_forest": ModelSpec(
            tree(RandomForestRegressor(random_state=seed, n_jobs=-1)),
            {"model__n_estimators": [300], "model__max_depth": [2, 4, None], "model__min_samples_leaf": [1, 2, 4]},
        ),
        "extra_trees": ModelSpec(
            tree(ExtraTreesRegressor(random_state=seed, n_jobs=-1)),
            {"model__n_estimators": [300], "model__max_depth": [2, 4, None], "model__min_samples_leaf": [1, 2, 4]},
        ),
        "gpr_matern": ModelSpec(
            scaled(
                GaussianProcessRegressor(
                    kernel=ConstantKernel(1.0) * Matern(length_scale=1.0, nu=1.5) + WhiteKernel(0.1),
                    normalize_y=True,
                    random_state=seed,
                    n_restarts_optimizer=5,
                )
            ),
            {},
        ),
    }


def _numeric_feature_columns(frame: pd.DataFrame) -> list[str]:
    excluded = {
        "site_paper", "site_pdb", "site_uniprot", "repair_published", "repair_from_means_raw",
        "repair_from_means_clipped", "repair_variance_approx", "label_weight",
        "active_mean", "active_sd", "control_mean", "control_sd", "bfactor_source",
        "ca_x", "ca_y", "ca_z", "structure_missing",
    }
    return [
        column for column in frame.select_dtypes(include=[np.number]).columns
        if column not in excluded and not column.startswith("active_rep") and not column.startswith("control_rep")
    ]


def _spatial_groups(frame: pd.DataFrame, n_groups: int) -> np.ndarray:
    coordinates = frame[["ca_x", "ca_y", "ca_z"]].to_numpy(dtype=float)
    if np.isnan(coordinates).any():
        raise ValueError("Missing CA coordinates; spatially blocked evaluation cannot proceed.")
    return AgglomerativeClustering(n_clusters=n_groups, linkage="ward").fit_predict(coordinates)


def run_nested_group_cv(
    merged: pd.DataFrame,
    target: str,
    output_dir: Path,
    n_spatial_groups: int = 4,
    bootstrap_iterations: int = 2000,
    seed: int = 20260727,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    feature_columns = _numeric_feature_columns(merged)
    if not feature_columns:
        raise ValueError("No numeric feature columns were found.")
    X = merged[feature_columns].to_numpy(dtype=float)
    y = merged[target].to_numpy(dtype=float)
    groups = _spatial_groups(merged, n_spatial_groups)
    outer_cv = GroupKFold(n_splits=n_spatial_groups)

    all_predictions: list[pd.DataFrame] = []
    metric_rows: list[dict[str, object]] = []

    for model_name, spec in _model_specs(seed).items():
        predictions = np.full(len(y), np.nan)
        uncertainties = np.full(len(y), np.nan)

        for fold, (train_idx, test_idx) in enumerate(outer_cv.split(X, y, groups), start=1):
            train_groups = groups[train_idx]
            unique_train_groups = np.unique(train_groups)
            if spec.param_grid and len(unique_train_groups) >= 2:
                inner_splits = min(3, len(unique_train_groups))
                search = GridSearchCV(
                    estimator=spec.estimator,
                    param_grid=spec.param_grid,
                    scoring="neg_mean_absolute_error",
                    cv=GroupKFold(n_splits=inner_splits),
                    n_jobs=-1,
                    refit=True,
                )
                search.fit(X[train_idx], y[train_idx], groups=train_groups)
                fitted = search.best_estimator_
            else:
                fitted = spec.estimator.fit(X[train_idx], y[train_idx])

            if model_name == "gpr_matern":
                # Pipeline output does not expose return_std, so transform before final model call.
                transformed = fitted[:-1].transform(X[test_idx])
                mean, std = fitted[-1].predict(transformed, return_std=True)
                predictions[test_idx] = mean
                uncertainties[test_idx] = std
            else:
                predictions[test_idx] = fitted.predict(X[test_idx])

        predictions = np.clip(predictions, 0.0, 1.0)
        metrics = regression_metrics(y, predictions)
        ci_low, ci_high = bootstrap_metric_ci(
            y, predictions, metric="mae", iterations=bootstrap_iterations, seed=seed
        )
        metric_rows.append(
            {
                "model": model_name,
                **metrics,
                "mae_ci_2.5": ci_low,
                "mae_ci_97.5": ci_high,
                "n_sites": len(y),
                "n_features": len(feature_columns),
                "validation": f"{n_spatial_groups}-fold spatial GroupKFold with nested tuning",
            }
        )
        all_predictions.append(
            pd.DataFrame(
                {
                    "protein_id": merged["protein_id"],
                    "site_uniprot": merged["site_uniprot"],
                    "site_paper": merged["site_paper"],
                    "spatial_group": groups,
                    "model": model_name,
                    "y_true": y,
                    "y_pred": predictions,
                    "y_std": uncertainties,
                }
            )
        )

    metrics_frame = pd.DataFrame(metric_rows).sort_values("mae").reset_index(drop=True)
    predictions_frame = pd.concat(all_predictions, ignore_index=True)
    output_dir.mkdir(parents=True, exist_ok=True)
    metrics_frame.to_csv(output_dir / "baseline_metrics.csv", index=False)
    predictions_frame.to_csv(output_dir / "out_of_fold_predictions.csv", index=False)
    (output_dir / "feature_columns.txt").write_text("\n".join(feature_columns), encoding="utf-8")
    return metrics_frame, predictions_frame


def load_merged_features(
    labels_csv: Path,
    sequence_features_csv: Path,
    structure_features_csv: Path,
) -> pd.DataFrame:
    labels = pd.read_csv(labels_csv)
    sequence = pd.read_csv(sequence_features_csv)
    structure = pd.read_csv(structure_features_csv)
    keys = ["protein_id", "site_uniprot"]
    return labels.merge(sequence, on=keys, validate="one_to_one").merge(
        structure, on=keys, validate="one_to_one"
    )
