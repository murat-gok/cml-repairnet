from __future__ import annotations

import hashlib
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import numpy as np
import pandas as pd
from sklearn.base import clone
from sklearn.cluster import AgglomerativeClustering, KMeans
from sklearn.ensemble import RandomForestRegressor
from sklearn.isotonic import IsotonicRegression
from sklearn.linear_model import Ridge
from sklearn.metrics import adjusted_rand_score, mean_absolute_error
from sklearn.model_selection import GroupKFold
from sklearn.neighbors import NearestNeighbors
from sklearn.preprocessing import StandardScaler

from cml_repairnet.evaluation.metrics import bootstrap_metric_ci, regression_metrics
from cml_repairnet.models.baselines import _numeric_feature_columns, _spatial_groups
from cml_repairnet.models.stage25 import (
    _inner_oof_predictions,
    _selected_feature_names,
    build_sample_weights,
    stage25_model_specs,
)


PRIMARY_MODELS = ("extra_trees", "lightgbm", "ridge")
DEFAULT_WEIGHT_MODES: Mapping[str, str] = {
    "extra_trees": "none",
    "lightgbm": "clipped",
    "ridge": "none",
}


def stage3_frozen_estimators(seed: int = 20260727) -> dict[str, Any]:
    """Return the Stage 2.5-prespecified model portfolio with frozen settings.

    Stage 3 is validation, not another hyperparameter search. Settings are frozen from the
    most frequently selected Stage 2.5 configurations. Correlation pruning is fixed at 0.85.
    """

    specs = stage25_model_specs(seed, feature_mode="nested_pruned")
    extra_trees = clone(specs["extra_trees"].estimator).set_params(
        prune__threshold=0.85,
        model__n_estimators=500,
        model__max_depth=4,
        model__min_samples_leaf=2,
        model__max_features=1.0,
    )
    ridge = clone(specs["ridge"].estimator).set_params(
        prune__threshold=0.85,
        model__alpha=10.0,
    )
    lightgbm = clone(specs["lightgbm"].estimator).set_params(
        prune__threshold=0.85,
        model__n_estimators=60,
        model__learning_rate=0.05,
        model__num_leaves=3,
        model__max_depth=2,
        model__min_child_samples=4,
        model__subsample=0.8,
        model__colsample_bytree=0.8,
        model__reg_alpha=0.25,
        model__reg_lambda=10.0,
    )
    return {"extra_trees": extra_trees, "lightgbm": lightgbm, "ridge": ridge}


@dataclass(frozen=True)
class DomainBoundary:
    name: str
    start: int
    end: int


DEFAULT_BSA_DOMAINS = (
    DomainBoundary("domain_I", 1, 195),
    DomainBoundary("domain_II", 196, 383),
    DomainBoundary("domain_III", 384, 583),
)


def assign_bsa_domains(
    frame: pd.DataFrame,
    boundaries: Sequence[DomainBoundary] = DEFAULT_BSA_DOMAINS,
    *,
    site_column: str = "site_paper",
) -> np.ndarray:
    """Assign mature-BSA sites to coarse three-domain labels.

    Boundaries are explicit and configurable because albumin domain annotations can differ
    slightly across resources. The defaults use mature 4F5S/BSA numbering.
    """

    if site_column not in frame.columns:
        raise ValueError(f"Missing domain coordinate column: {site_column}")
    sites = frame[site_column].to_numpy(dtype=int)
    labels = np.full(len(frame), "unassigned", dtype=object)
    for boundary in boundaries:
        mask = (sites >= boundary.start) & (sites <= boundary.end)
        labels[mask] = boundary.name
    if np.any(labels == "unassigned"):
        missing = sorted(set(sites[labels == "unassigned"].tolist()))
        raise ValueError(f"Sites outside configured BSA domains: {missing}")
    return labels.astype(str)


def partition_fingerprint(groups: Sequence[object]) -> str:
    """Return a cluster-label-invariant fingerprint for one spatial partition."""

    labels = np.asarray(groups)
    if labels.ndim != 1:
        raise ValueError("groups must be one-dimensional")
    co_membership = labels[:, None] == labels[None, :]
    upper = co_membership[np.triu_indices(len(labels), k=1)].astype(np.uint8)
    return hashlib.sha256(upper.tobytes()).hexdigest()


def seeded_spatial_groups(
    frame: pd.DataFrame,
    n_groups: int,
    seed: int,
    *,
    method: str = "kmeans_random",
    min_group_size: int = 2,
    max_attempts: int = 25,
) -> np.ndarray:
    """Generate reproducible, spatially coherent blocks with genuine seed sensitivity.

    The original Stage 3 implementation used K-means with ``n_init=50``. On the BSA
    coordinates this converged to the same optimum for every seed, so the apparent repeated
    validation did not probe partition sensitivity. The corrected default uses one random
    initialization per seed and retries only when a cluster is too small. ``method='ward'``
    remains the deterministic Stage 2 anchor.
    """

    coordinates = frame[["ca_x", "ca_y", "ca_z"]].to_numpy(dtype=float)
    if np.isnan(coordinates).any():
        raise ValueError("Missing CA coordinates; spatial validation cannot proceed")
    if not 2 <= int(n_groups) < len(frame):
        raise ValueError("n_groups must be between 2 and n_sites - 1")
    if min_group_size < 1:
        raise ValueError("min_group_size must be at least 1")
    scaled = StandardScaler().fit_transform(coordinates)
    method = method.lower()
    if method in {"kmeans", "kmeans_random"}:
        for attempt in range(max_attempts):
            attempt_seed = int(seed) + attempt * 1009
            labels = KMeans(
                n_clusters=int(n_groups),
                random_state=attempt_seed,
                n_init=1,
                init="random",
                algorithm="lloyd",
            ).fit_predict(scaled)
            if np.bincount(labels).min() >= int(min_group_size):
                return labels
        raise RuntimeError(
            f"Unable to generate a spatial partition with min_group_size={min_group_size}"
        )
    if method == "kmeans_pp":
        return KMeans(
            n_clusters=int(n_groups),
            random_state=int(seed),
            n_init=1,
            init="k-means++",
            algorithm="lloyd",
        ).fit_predict(scaled)
    if method == "ward":
        return AgglomerativeClustering(n_clusters=int(n_groups), linkage="ward").fit_predict(
            scaled
        )
    raise ValueError("method must be 'kmeans_random', 'kmeans_pp', or 'ward'")


def generate_unique_spatial_partitions(
    frame: pd.DataFrame,
    n_groups: int,
    *,
    start_seed: int,
    target_count: int,
    max_candidates: int | None = None,
    min_group_size: int = 2,
) -> list[tuple[int, np.ndarray, str]]:
    """Generate a requested number of unique, label-invariant spatial partitions."""

    if target_count < 1:
        raise ValueError("target_count must be at least 1")
    max_candidates = max_candidates or max(100, target_count * 20)
    unique: dict[str, tuple[int, np.ndarray, str]] = {}
    for offset in range(max_candidates):
        candidate_seed = int(start_seed) + offset
        groups = seeded_spatial_groups(
            frame,
            n_groups,
            candidate_seed,
            method="kmeans_random",
            min_group_size=min_group_size,
        )
        fingerprint = partition_fingerprint(groups)
        unique.setdefault(fingerprint, (candidate_seed, groups, fingerprint))
        if len(unique) >= target_count:
            break
    if len(unique) < target_count:
        raise RuntimeError(
            f"Requested {target_count} unique partitions but found {len(unique)} "
            f"after {max_candidates} candidates"
        )
    return list(unique.values())[:target_count]


def _conformal_quantile(residuals: np.ndarray, alpha: float) -> float:
    values = np.asarray(residuals, dtype=float)
    values = values[np.isfinite(values)]
    if values.size == 0:
        return float("nan")
    if not 0.0 < alpha < 1.0:
        raise ValueError("alpha must be in (0, 1)")
    level = min(1.0, math.ceil((values.size + 1) * (1.0 - alpha)) / values.size)
    try:
        return float(np.quantile(values, level, method="higher"))
    except TypeError:  # numpy < 1.22 compatibility
        return float(np.quantile(values, level, interpolation="higher"))


def _fit_affine_calibrator(predictions: np.ndarray, y: np.ndarray) -> Ridge:
    model = Ridge(alpha=1.0)
    model.fit(np.asarray(predictions).reshape(-1, 1), np.asarray(y, dtype=float))
    return model


def _fit_isotonic_calibrator(predictions: np.ndarray, y: np.ndarray) -> IsotonicRegression | None:
    predictions = np.asarray(predictions, dtype=float)
    if np.unique(predictions).size < 4:
        return None
    model = IsotonicRegression(y_min=0.0, y_max=1.0, out_of_bounds="clip")
    model.fit(predictions, np.asarray(y, dtype=float))
    return model


def _ood_scores(
    X_train: pd.DataFrame,
    X_test: pd.DataFrame,
    *,
    n_neighbors: int = 5,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return test distance, empirical training-percentile and OOD flag.

    Distances are calculated in a standardized feature space. The reference distribution is
    the leave-one-out mean k-nearest-neighbour distance among training observations.
    """

    train = X_train.apply(pd.to_numeric, errors="coerce").copy()
    test = X_test.apply(pd.to_numeric, errors="coerce").copy()
    medians = train.median(axis=0).fillna(0.0)
    train = train.fillna(medians)
    test = test.fillna(medians)
    scaler = StandardScaler().fit(train)
    train_z = scaler.transform(train)
    test_z = scaler.transform(test)
    k = max(1, min(int(n_neighbors), len(train_z) - 1))
    reference_nn = NearestNeighbors(n_neighbors=k + 1).fit(train_z)
    train_dist = reference_nn.kneighbors(train_z, return_distance=True)[0][:, 1:].mean(axis=1)
    test_nn = NearestNeighbors(n_neighbors=k).fit(train_z)
    test_dist = test_nn.kneighbors(test_z, return_distance=True)[0].mean(axis=1)
    percentiles = np.asarray([(train_dist <= value).mean() for value in test_dist], dtype=float)
    flags = percentiles >= 0.95
    return test_dist, percentiles, flags


def _feature_subset(frame: pd.DataFrame, excluded: Iterable[str] = ()) -> list[str]:
    excluded_set = set(excluded)
    return [column for column in _numeric_feature_columns(frame) if column not in excluded_set]


def _evaluate_group_partition(
    merged: pd.DataFrame,
    target: str,
    groups: np.ndarray,
    *,
    partition_id: str,
    model_names: Sequence[str] = PRIMARY_MODELS,
    feature_mode: str = "nested_pruned",
    model_weight_modes: Mapping[str, str] = DEFAULT_WEIGHT_MODES,
    excluded_features: Iterable[str] = (),
    n_jobs: int = 1,
    seed: int = 20260727,
    conformal_alpha: float = 0.10,
    include_calibration: bool = True,
    include_uncertainty: bool = True,
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    feature_columns = _feature_subset(merged, excluded_features)
    if not feature_columns:
        raise ValueError("No model features remain after exclusion")
    X = merged[feature_columns].copy()
    y = merged[target].to_numpy(dtype=float)
    groups = np.asarray(groups)
    unique_groups = np.unique(groups)
    if unique_groups.size < 2:
        raise ValueError("At least two outer groups are required")
    outer_cv = GroupKFold(n_splits=unique_groups.size)
    frozen_estimators = stage3_frozen_estimators(seed)

    predictions: dict[str, np.ndarray] = {
        model: np.full(len(y), np.nan, dtype=float) for model in model_names
    }
    fold_by_site = np.full(len(y), -1, dtype=int)
    selected_features: dict[str, list[list[str]]] = {model: [] for model in model_names}
    best_params: dict[str, list[dict[str, Any]]] = {model: [] for model in model_names}

    calibrated_affine = np.full(len(y), np.nan, dtype=float)
    calibrated_isotonic = np.full(len(y), np.nan, dtype=float)
    interval_low = np.full(len(y), np.nan, dtype=float)
    interval_high = np.full(len(y), np.nan, dtype=float)
    conformal_q = np.full(len(y), np.nan, dtype=float)
    ood_distance = np.full(len(y), np.nan, dtype=float)
    ood_percentile = np.full(len(y), np.nan, dtype=float)
    ood_flag = np.full(len(y), False, dtype=bool)

    for fold, (train_idx, test_idx) in enumerate(outer_cv.split(X, y, groups), start=1):
        fold_by_site[test_idx] = fold
        fitted_models: dict[str, Any] = {}
        train_weights_by_model: dict[str, np.ndarray] = {}

        for model_name in model_names:
            if model_name not in frozen_estimators:
                raise ValueError(f"Unknown Stage 3 model: {model_name}")
            mode = model_weight_modes.get(model_name, "none")
            train_weights = build_sample_weights(merged.iloc[train_idx], mode)
            train_weights_by_model[model_name] = train_weights
            fitted = clone(frozen_estimators[model_name]).fit(
                X.iloc[train_idx],
                y[train_idx],
                model__sample_weight=train_weights,
            )
            params = fitted.get_params(deep=True)
            fitted_models[model_name] = fitted
            best_params[model_name].append(
                {key: value for key, value in params.items() if key.startswith(("model__", "prune__"))}
            )
            selected_features[model_name].append(
                _selected_feature_names(fitted, feature_columns)
            )
            predictions[model_name][test_idx] = np.clip(
                fitted.predict(X.iloc[test_idx]), 0.0, 1.0
            )

        if "extra_trees" in fitted_models and (include_calibration or include_uncertainty):
            main_model = fitted_models["extra_trees"]
            main_weights = train_weights_by_model["extra_trees"]
            inner_pred = _inner_oof_predictions(
                main_model,
                X.iloc[train_idx].reset_index(drop=True),
                y[train_idx],
                groups[train_idx],
                main_weights,
            )
            raw_test = predictions["extra_trees"][test_idx]

            if include_calibration:
                affine = _fit_affine_calibrator(inner_pred, y[train_idx])
                calibrated_affine[test_idx] = np.clip(
                    affine.predict(raw_test.reshape(-1, 1)), 0.0, 1.0
                )
                isotonic = _fit_isotonic_calibrator(inner_pred, y[train_idx])
                calibrated_isotonic[test_idx] = (
                    raw_test if isotonic is None else isotonic.predict(raw_test)
                )

            if include_uncertainty:
                q = _conformal_quantile(np.abs(y[train_idx] - inner_pred), conformal_alpha)
                conformal_q[test_idx] = q
                interval_low[test_idx] = np.clip(raw_test - q, 0.0, 1.0)
                interval_high[test_idx] = np.clip(raw_test + q, 0.0, 1.0)
                distance, percentile, flag = _ood_scores(
                    X.iloc[train_idx], X.iloc[test_idx]
                )
                ood_distance[test_idx] = distance
                ood_percentile[test_idx] = percentile
                ood_flag[test_idx] = flag

    metric_rows: list[dict[str, Any]] = []
    prediction_frames: list[pd.DataFrame] = []
    for model_name, model_pred in predictions.items():
        metrics = regression_metrics(y, model_pred)
        low, high = bootstrap_metric_ci(
            y, model_pred, metric="mae", iterations=2000, seed=seed
        )
        selected_count = [len(items) for items in selected_features[model_name]]
        metric_rows.append(
            {
                "partition_id": partition_id,
                "model": model_name,
                **metrics,
                "mae_ci_2.5": low,
                "mae_ci_97.5": high,
                "target": target,
                "n_groups": int(unique_groups.size),
                "n_sites": len(y),
                "n_input_features": len(feature_columns),
                "mean_selected_features": float(np.mean(selected_count)),
                "feature_mode": feature_mode,
                "weight_mode": model_weight_modes.get(model_name, "none"),
            }
        )
        prediction_frames.append(
            pd.DataFrame(
                {
                    "protein_id": merged["protein_id"],
                    "site_uniprot": merged["site_uniprot"],
                    "site_paper": merged["site_paper"],
                    "partition_id": partition_id,
                    "outer_group": groups.astype(str),
                    "outer_fold": fold_by_site,
                    "model": model_name,
                    "y_true": y,
                    "y_pred": model_pred,
                    "target": target,
                }
            )
        )

    if include_calibration and np.isfinite(calibrated_affine).all():
        for name, values in (
            ("extra_trees_affine_calibrated", calibrated_affine),
            ("extra_trees_isotonic_calibrated", calibrated_isotonic),
        ):
            metrics = regression_metrics(y, values)
            low, high = bootstrap_metric_ci(
                y, values, metric="mae", iterations=2000, seed=seed
            )
            low_mask = y <= 0.10
            high_mask = y >= 0.90
            metric_rows.append(
                {
                    "partition_id": partition_id,
                    "model": name,
                    **metrics,
                    "mae_ci_2.5": low,
                    "mae_ci_97.5": high,
                    "target": target,
                    "n_groups": int(unique_groups.size),
                    "n_sites": len(y),
                    "n_input_features": len(feature_columns),
                    "mean_selected_features": float("nan"),
                    "feature_mode": feature_mode,
                    "weight_mode": "none",
                    "low_boundary_mae": (
                        float(mean_absolute_error(y[low_mask], values[low_mask]))
                        if np.any(low_mask) else float("nan")
                    ),
                    "high_boundary_mae": (
                        float(mean_absolute_error(y[high_mask], values[high_mask]))
                        if np.any(high_mask) else float("nan")
                    ),
                }
            )
            prediction_frames.append(
                pd.DataFrame(
                    {
                        "protein_id": merged["protein_id"],
                        "site_uniprot": merged["site_uniprot"],
                        "site_paper": merged["site_paper"],
                        "partition_id": partition_id,
                        "outer_group": groups.astype(str),
                        "outer_fold": fold_by_site,
                        "model": name,
                        "y_true": y,
                        "y_pred": values,
                        "target": target,
                    }
                )
            )

    uncertainty_frame = pd.DataFrame()
    if include_uncertainty and np.isfinite(interval_low).all():
        raw = predictions["extra_trees"]
        covered = (y >= interval_low) & (y <= interval_high)
        uncertainty_frame = pd.DataFrame(
            {
                "protein_id": merged["protein_id"],
                "site_uniprot": merged["site_uniprot"],
                "site_paper": merged["site_paper"],
                "partition_id": partition_id,
                "outer_group": groups.astype(str),
                "outer_fold": fold_by_site,
                "y_true": y,
                "y_pred": raw,
                "interval_low": interval_low,
                "interval_high": interval_high,
                "interval_width": interval_high - interval_low,
                "covered": covered,
                "conformal_q": conformal_q,
                "ood_distance": ood_distance,
                "ood_percentile": ood_percentile,
                "ood_flag": ood_flag,
                "absolute_error": np.abs(y - raw),
            }
        )

    details = {
        "partition_id": partition_id,
        "target": target,
        "feature_columns": feature_columns,
        "best_params": best_params,
        "selected_features": selected_features,
        "conformal_alpha": conformal_alpha,
    }
    return (
        pd.DataFrame(metric_rows).sort_values("mae").reset_index(drop=True),
        pd.concat(prediction_frames, ignore_index=True),
        {"details": details, "uncertainty": uncertainty_frame},
    )


def run_repeated_spatial_validation(
    merged: pd.DataFrame,
    target: str,
    output_dir: Path,
    *,
    seeds: Sequence[int],
    group_counts: Sequence[int] = (4,),
    model_names: Sequence[str] = PRIMARY_MODELS,
    n_jobs: int = 1,
    feature_mode: str = "nested_pruned",
    include_anchor: bool = True,
    model_seed: int = 20260727,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Evaluate frozen models over genuinely distinct spatial partitions.

    Estimator randomness is held fixed with ``model_seed``. Therefore variability across
    rows reflects partition choice rather than a mixture of partition and model randomness.
    """

    output_dir.mkdir(parents=True, exist_ok=True)
    metric_frames: list[pd.DataFrame] = []
    prediction_frames: list[pd.DataFrame] = []
    manifest_rows: list[dict[str, Any]] = []

    requested_count = len(seeds)
    if requested_count < 1:
        raise ValueError("At least one seed/partition is required")
    start_seed = int(min(seeds))

    for n_groups in group_counts:
        partition_groups: list[tuple[str, str, int, np.ndarray, str]] = []
        if include_anchor:
            groups = _spatial_groups(merged, int(n_groups))
            partition_groups.append(
                (
                    f"ward_g{n_groups}",
                    "ward",
                    -1,
                    groups,
                    partition_fingerprint(groups),
                )
            )

        generated = generate_unique_spatial_partitions(
            merged,
            int(n_groups),
            start_seed=start_seed,
            target_count=requested_count,
            min_group_size=2,
        )
        for actual_seed, groups, fingerprint in generated:
            partition_groups.append(
                (
                    f"kmeans_random_seed{actual_seed}_g{n_groups}",
                    "kmeans_random",
                    int(actual_seed),
                    groups,
                    fingerprint,
                )
            )

        for partition_id, method, clustering_seed, groups, fingerprint in partition_groups:
            metrics, predictions, _ = _evaluate_group_partition(
                merged,
                target,
                groups,
                partition_id=partition_id,
                model_names=model_names,
                feature_mode=feature_mode,
                n_jobs=n_jobs,
                seed=model_seed,
                include_calibration=False,
                include_uncertainty=False,
            )
            group_sizes = np.bincount(pd.factorize(groups)[0])
            metrics["clustering_method"] = method
            metrics["clustering_seed"] = clustering_seed
            metrics["partition_fingerprint"] = fingerprint
            predictions["clustering_method"] = method
            predictions["clustering_seed"] = clustering_seed
            predictions["partition_fingerprint"] = fingerprint
            metric_frames.append(metrics)
            prediction_frames.append(predictions)
            manifest_rows.append(
                {
                    "partition_id": partition_id,
                    "n_groups": int(n_groups),
                    "clustering_method": method,
                    "clustering_seed": clustering_seed,
                    "partition_fingerprint": fingerprint,
                    "min_group_size": int(group_sizes.min()),
                    "max_group_size": int(group_sizes.max()),
                }
            )

    all_metrics = pd.concat(metric_frames, ignore_index=True)
    all_predictions = pd.concat(prediction_frames, ignore_index=True)
    manifest = pd.DataFrame(manifest_rows)

    diversity_rows: list[dict[str, Any]] = []
    for n_groups, group_manifest in manifest.groupby("n_groups"):
        partition_ids = group_manifest["partition_id"].tolist()
        assignment = (
            all_predictions[all_predictions["model"] == model_names[0]]
            .pivot(index="site_uniprot", columns="partition_id", values="outer_group")
        )
        aris: list[float] = []
        for i, first in enumerate(partition_ids):
            for second in partition_ids[i + 1 :]:
                aris.append(
                    float(adjusted_rand_score(assignment[first], assignment[second]))
                )
        diversity_rows.append(
            {
                "n_groups": int(n_groups),
                "partition_count_diagnostic": len(partition_ids),
                "unique_fingerprints": int(group_manifest["partition_fingerprint"].nunique()),
                "pairwise_ari_min": float(np.min(aris)) if aris else float("nan"),
                "pairwise_ari_median": float(np.median(aris)) if aris else float("nan"),
                "pairwise_ari_max": float(np.max(aris)) if aris else float("nan"),
            }
        )
    diversity = pd.DataFrame(diversity_rows)

    summary = (
        all_metrics.groupby(["model", "n_groups"], as_index=False)
        .agg(
            partitions=("partition_id", "nunique"),
            mean_mae=("mae", "mean"),
            sd_mae=("mae", "std"),
            median_mae=("mae", "median"),
            mae_q05=("mae", lambda values: float(np.quantile(values, 0.05))),
            mae_q95=("mae", lambda values: float(np.quantile(values, 0.95))),
            mean_spearman=("spearman_rho", "mean"),
            spearman_q05=("spearman_rho", lambda values: float(np.quantile(values, 0.05))),
            spearman_q95=("spearman_rho", lambda values: float(np.quantile(values, 0.95))),
            mean_r2=("r2", "mean"),
        )
        .sort_values(["n_groups", "mean_mae"])
        .reset_index(drop=True)
    )
    winner_rows = (
        all_metrics.loc[
            all_metrics.groupby(["partition_id", "n_groups"])["mae"].idxmin(),
            ["partition_id", "n_groups", "model"],
        ]
        .groupby(["model", "n_groups"])
        .size()
        .rename("partition_wins")
        .reset_index()
    )
    summary = summary.merge(winner_rows, on=["model", "n_groups"], how="left")
    summary["partition_wins"] = summary["partition_wins"].fillna(0).astype(int)
    summary = summary.merge(diversity, on="n_groups", how="left")

    all_metrics.to_csv(output_dir / "repeated_spatial_metrics.csv", index=False)
    all_predictions.to_csv(output_dir / "repeated_spatial_predictions.csv", index=False)
    summary.to_csv(output_dir / "repeated_spatial_summary.csv", index=False)
    manifest.to_csv(output_dir / "repeated_spatial_partition_manifest.csv", index=False)
    diversity.to_csv(output_dir / "repeated_spatial_partition_diversity.csv", index=False)
    return all_metrics, all_predictions, summary


def run_leave_one_domain_out(
    merged: pd.DataFrame,
    target: str,
    output_dir: Path,
    *,
    model_names: Sequence[str] = PRIMARY_MODELS,
    n_jobs: int = 1,
    feature_mode: str = "nested_pruned",
    conformal_alpha: float = 0.10,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    output_dir.mkdir(parents=True, exist_ok=True)
    domains = assign_bsa_domains(merged)
    metrics, predictions, extras = _evaluate_group_partition(
        merged,
        target,
        domains,
        partition_id="leave_one_bsa_domain_out",
        model_names=model_names,
        feature_mode=feature_mode,
        n_jobs=n_jobs,
        seed=20260727,
        conformal_alpha=conformal_alpha,
        include_calibration=True,
        include_uncertainty=True,
    )
    uncertainty = extras["uncertainty"]
    metrics.to_csv(output_dir / "domain_metrics.csv", index=False)
    predictions.to_csv(output_dir / "domain_predictions.csv", index=False)
    uncertainty.to_csv(output_dir / "domain_uncertainty.csv", index=False)
    (output_dir / "domain_fold_details.json").write_text(
        json.dumps(extras["details"], indent=2, sort_keys=True), encoding="utf-8"
    )
    return metrics, predictions, uncertainty


def _conditional_permutation_fold(
    fitted: Any,
    X_train: pd.DataFrame,
    X_test: pd.DataFrame,
    y_test: np.ndarray,
    feature: str,
    *,
    repeats: int,
    seed: int,
) -> tuple[np.ndarray, np.ndarray]:
    if feature not in X_train.columns:
        raise ValueError(f"Conditional feature not found: {feature}")
    other_columns = [column for column in X_train.columns if column != feature]
    conditional_model = RandomForestRegressor(
        n_estimators=200,
        min_samples_leaf=3,
        max_features=0.8,
        random_state=seed,
        n_jobs=1,
    )
    conditional_model.fit(X_train[other_columns], X_train[feature])
    expected_test = conditional_model.predict(X_test[other_columns])
    residuals = X_test[feature].to_numpy(dtype=float) - expected_test
    baseline = mean_absolute_error(y_test, np.clip(fitted.predict(X_test), 0.0, 1.0))
    rng = np.random.default_rng(seed)
    unconditional = np.empty(repeats, dtype=float)
    conditional = np.empty(repeats, dtype=float)
    train_min = float(X_train[feature].min())
    train_max = float(X_train[feature].max())

    for i in range(repeats):
        permuted = X_test.copy()
        permuted[feature] = rng.permutation(permuted[feature].to_numpy())
        unconditional[i] = (
            mean_absolute_error(y_test, np.clip(fitted.predict(permuted), 0.0, 1.0))
            - baseline
        )

        conditional_values = expected_test + rng.permutation(residuals)
        conditional_values = np.clip(conditional_values, train_min, train_max)
        permuted_conditional = X_test.copy()
        permuted_conditional[feature] = conditional_values
        conditional[i] = (
            mean_absolute_error(
                y_test, np.clip(fitted.predict(permuted_conditional), 0.0, 1.0)
            )
            - baseline
        )
    return unconditional, conditional


def run_cysteine_mechanism_analysis(
    merged: pd.DataFrame,
    target: str,
    output_dir: Path,
    *,
    feature: str = "neighbor_cysteine_count_8a",
    n_spatial_groups: int = 4,
    permutation_repeats: int = 500,
    n_jobs: int = 1,
    seed: int = 20260727,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Run leakage-controlled drop-column and conditional permutation analyses."""

    output_dir.mkdir(parents=True, exist_ok=True)
    groups = _spatial_groups(merged, n_spatial_groups)
    full_metrics, full_predictions, _ = _evaluate_group_partition(
        merged,
        target,
        groups,
        partition_id="cysteine_full",
        model_names=("extra_trees", "ridge"),
        feature_mode="nested_pruned",
        n_jobs=n_jobs,
        seed=seed,
        include_calibration=False,
        include_uncertainty=False,
    )
    drop_metrics, drop_predictions, _ = _evaluate_group_partition(
        merged,
        target,
        groups,
        partition_id="cysteine_drop_column",
        model_names=("extra_trees", "ridge"),
        feature_mode="nested_pruned",
        excluded_features=(feature,),
        n_jobs=n_jobs,
        seed=seed,
        include_calibration=False,
        include_uncertainty=False,
    )

    comparison_rows: list[dict[str, Any]] = []
    for model_name in ("extra_trees", "ridge"):
        full_row = full_metrics.loc[full_metrics["model"] == model_name].iloc[0]
        drop_row = drop_metrics.loc[drop_metrics["model"] == model_name].iloc[0]
        comparison_rows.append(
            {
                "model": model_name,
                "feature": feature,
                "full_mae": full_row["mae"],
                "drop_column_mae": drop_row["mae"],
                "drop_column_mae_increase": drop_row["mae"] - full_row["mae"],
                "full_spearman": full_row["spearman_rho"],
                "drop_column_spearman": drop_row["spearman_rho"],
            }
        )

    feature_columns = _feature_subset(merged)
    X = merged[feature_columns].copy()
    y = merged[target].to_numpy(dtype=float)
    outer_cv = GroupKFold(n_splits=n_spatial_groups)
    estimator = stage3_frozen_estimators(seed)["extra_trees"]
    unconditional_values: list[float] = []
    conditional_values: list[float] = []
    per_fold_rows: list[dict[str, Any]] = []

    for fold, (train_idx, test_idx) in enumerate(outer_cv.split(X, y, groups), start=1):
        weights = build_sample_weights(merged.iloc[train_idx], "none")
        fitted = clone(estimator).fit(
            X.iloc[train_idx],
            y[train_idx],
            model__sample_weight=weights,
        )
        params = {
            key: value for key, value in fitted.get_params(deep=True).items()
            if key.startswith(("model__", "prune__"))
        }
        uncond, cond = _conditional_permutation_fold(
            fitted,
            X.iloc[train_idx],
            X.iloc[test_idx],
            y[test_idx],
            feature,
            repeats=permutation_repeats,
            seed=seed + fold,
        )
        unconditional_values.extend(uncond.tolist())
        conditional_values.extend(cond.tolist())
        per_fold_rows.append(
            {
                "fold": fold,
                "feature": feature,
                "unconditional_mae_increase_mean": float(np.mean(uncond)),
                "conditional_mae_increase_mean": float(np.mean(cond)),
                "best_params": json.dumps(params, sort_keys=True),
            }
        )

    unconditional_array = np.asarray(unconditional_values)
    conditional_array = np.asarray(conditional_values)
    comparison_rows.append(
        {
            "model": "extra_trees",
            "feature": feature,
            "analysis": "permutation",
            "unconditional_mae_increase_mean": float(np.mean(unconditional_array)),
            "unconditional_ci_2.5": float(np.quantile(unconditional_array, 0.025)),
            "unconditional_ci_97.5": float(np.quantile(unconditional_array, 0.975)),
            "conditional_mae_increase_mean": float(np.mean(conditional_array)),
            "conditional_ci_2.5": float(np.quantile(conditional_array, 0.025)),
            "conditional_ci_97.5": float(np.quantile(conditional_array, 0.975)),
        }
    )

    comparison = pd.DataFrame(comparison_rows)
    per_fold = pd.DataFrame(per_fold_rows)
    comparison.to_csv(output_dir / "cysteine_mechanism_summary.csv", index=False)
    per_fold.to_csv(output_dir / "cysteine_conditional_permutation_by_fold.csv", index=False)
    pd.concat([full_predictions, drop_predictions], ignore_index=True).to_csv(
        output_dir / "cysteine_drop_column_predictions.csv", index=False
    )
    return comparison, per_fold


def summarize_uncertainty(frame: pd.DataFrame) -> pd.DataFrame:
    if frame.empty:
        return pd.DataFrame()
    return pd.DataFrame(
        [
            {
                "n_sites": len(frame),
                "coverage": float(frame["covered"].mean()),
                "mean_interval_width": float(frame["interval_width"].mean()),
                "median_interval_width": float(frame["interval_width"].median()),
                "ood_count": int(frame["ood_flag"].sum()),
                "ood_mae": (
                    float(frame.loc[frame["ood_flag"], "absolute_error"].mean())
                    if frame["ood_flag"].any() else float("nan")
                ),
                "in_domain_mae": (
                    float(frame.loc[~frame["ood_flag"], "absolute_error"].mean())
                    if (~frame["ood_flag"]).any() else float("nan")
                ),
                "error_ood_spearman": float(
                    pd.Series(frame["absolute_error"]).corr(
                        pd.Series(frame["ood_percentile"]), method="spearman"
                    )
                ),
            }
        ]
    )
