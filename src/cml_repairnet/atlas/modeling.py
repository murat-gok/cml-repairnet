from __future__ import annotations

import json
from collections.abc import Mapping
from pathlib import Path
from typing import Any

import joblib
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.base import clone
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from sklearn.model_selection import GroupKFold
from sklearn.neighbors import NearestNeighbors
from sklearn.preprocessing import StandardScaler

from cml_repairnet.atlas.structures import TRANSFER_FEATURES
from cml_repairnet.models.stage3 import (
    _conformal_quantile,
    seeded_spatial_groups,
    stage3_frozen_estimators,
)
from cml_repairnet.models.stage25 import _selected_feature_names, build_sample_weights

MODEL_NAMES = ("extra_trees", "lightgbm", "ridge")
MODEL_WEIGHT_MODES: Mapping[str, str] = {
    "extra_trees": "none",
    "lightgbm": "clipped",
    "ridge": "none",
}


def _metrics(y_true: np.ndarray, y_pred: np.ndarray) -> dict[str, float]:
    rho = spearmanr(y_true, y_pred).statistic
    return {
        "mae": float(mean_absolute_error(y_true, y_pred)),
        "rmse": float(mean_squared_error(y_true, y_pred) ** 0.5),
        "spearman_rho": float(rho) if np.isfinite(rho) else float("nan"),
        "r2": float(r2_score(y_true, y_pred)),
    }


def _fit_with_weights(model: Any, X: pd.DataFrame, y: np.ndarray, weights: np.ndarray) -> Any:
    fitted = clone(model)
    try:
        fitted.fit(X, y, model__sample_weight=weights)
    except TypeError:
        fitted.fit(X, y)
    return fitted


def _build_ood_reference(X: pd.DataFrame, n_neighbors: int = 5) -> dict[str, Any]:
    numeric = X.apply(pd.to_numeric, errors="coerce").copy()
    medians = numeric.median(axis=0).fillna(0.0)
    filled = numeric.fillna(medians)
    scaler = StandardScaler().fit(filled)
    values = scaler.transform(filled)
    k = max(1, min(int(n_neighbors), len(values) - 1))
    reference_nn = NearestNeighbors(n_neighbors=k + 1).fit(values)
    reference_distances = reference_nn.kneighbors(values, return_distance=True)[0][:, 1:].mean(
        axis=1
    )
    neighbor_model = NearestNeighbors(n_neighbors=k).fit(values)
    return {
        "feature_columns": list(X.columns),
        "medians": medians.to_dict(),
        "scaler": scaler,
        "neighbor_model": neighbor_model,
        "reference_distances": reference_distances,
        "n_neighbors": k,
    }


def _score_ood(reference: Mapping[str, Any], X: pd.DataFrame) -> tuple[np.ndarray, np.ndarray]:
    columns = list(reference["feature_columns"])
    numeric = X.reindex(columns=columns).apply(pd.to_numeric, errors="coerce")
    medians = pd.Series(reference["medians"], dtype=float)
    numeric = numeric.fillna(medians)
    transformed = reference["scaler"].transform(numeric)
    distances = (
        reference["neighbor_model"].kneighbors(transformed, return_distance=True)[0].mean(axis=1)
    )
    reference_distances = np.asarray(reference["reference_distances"], dtype=float)
    percentiles = np.asarray(
        [(reference_distances <= value).mean() for value in distances], dtype=float
    )
    return distances, percentiles


def fit_transfer_safe_bundle(
    labels_csv: Path,
    bsa_structure_features_csv: Path,
    output_dir: Path,
    *,
    target: str = "repair_published",
    n_groups: int = 4,
    seed: int = 20260727,
    conformal_alpha: float = 0.10,
) -> dict[str, Any]:
    labels = pd.read_csv(labels_csv)
    structures = pd.read_csv(bsa_structure_features_csv)
    merge_keys = [
        column
        for column in ("protein_id", "site_uniprot")
        if column in labels.columns and column in structures.columns
    ]
    if not merge_keys:
        merge_keys = ["site_uniprot"]
    frame = labels.merge(structures, on=merge_keys, how="inner", suffixes=("", "_structure"))
    missing = [feature for feature in TRANSFER_FEATURES if feature not in frame.columns]
    if missing:
        raise ValueError(f"BSA structure table is missing transfer-safe features: {missing}")
    if target not in frame.columns:
        raise ValueError(f"Target column not found: {target}")
    X = frame[list(TRANSFER_FEATURES)].copy()
    y = frame[target].to_numpy(dtype=float)
    groups = seeded_spatial_groups(frame, n_groups, seed, method="ward", min_group_size=2)
    estimators = stage3_frozen_estimators(seed)
    outer_cv = GroupKFold(n_splits=len(np.unique(groups)))
    oof = {name: np.full(len(frame), np.nan, dtype=float) for name in MODEL_NAMES}
    fold_index = np.full(len(frame), -1, dtype=int)
    selected_by_model: dict[str, list[list[str]]] = {name: [] for name in MODEL_NAMES}

    for fold, (train_idx, test_idx) in enumerate(outer_cv.split(X, y, groups), start=1):
        fold_index[test_idx] = fold
        for name in MODEL_NAMES:
            weights = build_sample_weights(frame.iloc[train_idx], MODEL_WEIGHT_MODES[name])
            fitted = _fit_with_weights(estimators[name], X.iloc[train_idx], y[train_idx], weights)
            oof[name][test_idx] = np.clip(fitted.predict(X.iloc[test_idx]), 0.0, 1.0)
            selected_by_model[name].append(_selected_feature_names(fitted, list(TRANSFER_FEATURES)))

    residuals = np.abs(y - oof["extra_trees"])
    conformal_q = _conformal_quantile(residuals, conformal_alpha)
    fitted_models: dict[str, Any] = {}
    selected_final: dict[str, list[str]] = {}
    for name in MODEL_NAMES:
        weights = build_sample_weights(frame, MODEL_WEIGHT_MODES[name])
        fitted = _fit_with_weights(estimators[name], X, y, weights)
        fitted_models[name] = fitted
        selected_final[name] = _selected_feature_names(fitted, list(TRANSFER_FEATURES))

    model_metrics = {name: _metrics(y, oof[name]) for name in MODEL_NAMES}
    oof_frame = frame[
        [
            column
            for column in ("protein_id", "site_uniprot", "site_paper")
            if column in frame.columns
        ]
    ].copy()
    oof_frame["y_true"] = y
    oof_frame["spatial_group"] = groups
    oof_frame["fold"] = fold_index
    for name in MODEL_NAMES:
        oof_frame[f"prediction_{name}"] = oof[name]
    output_dir.mkdir(parents=True, exist_ok=True)
    oof_frame.to_csv(output_dir / "stage4_transfer_safe_oof_predictions.csv", index=False)

    bundle = {
        "version": "0.5.7",
        "target": target,
        "feature_columns": list(TRANSFER_FEATURES),
        "models": fitted_models,
        "model_weight_modes": dict(MODEL_WEIGHT_MODES),
        "selected_features": selected_final,
        "selected_features_by_fold": selected_by_model,
        "model_metrics": model_metrics,
        "conformal_alpha": conformal_alpha,
        "conformal_q": conformal_q,
        "ood_reference": _build_ood_reference(X),
        "training_site_count": len(frame),
        "training_groups": n_groups,
        "seed": seed,
        "transfer_note": (
            "B-factor features are excluded because AlphaFold pLDDT is not physically "
            "equivalent to crystallographic B-factor. pLDDT is used only for quality/OOD."
        ),
    }
    joblib.dump(bundle, output_dir / "stage4_transfer_safe_model_bundle.joblib")
    card = {
        key: value
        for key, value in bundle.items()
        if key not in {"models", "ood_reference", "selected_features_by_fold"}
    }
    card["model_metrics"] = model_metrics
    card["selected_features"] = selected_final
    (output_dir / "stage4_transfer_safe_model_card.json").write_text(
        json.dumps(card, indent=2, default=str), encoding="utf-8"
    )
    return bundle


def load_or_fit_transfer_bundle(
    model_dir: Path,
    labels_csv: Path,
    bsa_structure_features_csv: Path,
    *,
    refit: bool = False,
) -> dict[str, Any]:
    path = model_dir / "stage4_transfer_safe_model_bundle.joblib"
    if path.exists() and not refit:
        return joblib.load(path)
    return fit_transfer_safe_bundle(labels_csv, bsa_structure_features_csv, model_dir)


def score_candidate_features(
    features_csv: Path,
    output_csv: Path,
    model_dir: Path,
    labels_csv: Path,
    bsa_structure_features_csv: Path,
    *,
    refit_model: bool = False,
) -> pd.DataFrame:
    frame = pd.read_csv(features_csv)
    bundle = load_or_fit_transfer_bundle(
        model_dir,
        labels_csv,
        bsa_structure_features_csv,
        refit=refit_model,
    )
    features = list(bundle["feature_columns"])
    missing = [feature for feature in features if feature not in frame.columns]
    if missing:
        raise ValueError(f"Candidate structure table is missing model features: {missing}")
    X = frame[features].copy()
    for name, model in bundle["models"].items():
        frame[f"repairability_{name}"] = np.clip(model.predict(X), 0.0, 1.0)

    ood_distance, ood_percentile = _score_ood(bundle["ood_reference"], X)
    frame["ood_distance"] = ood_distance
    frame["ood_percentile"] = ood_percentile
    frame["ood_flag"] = frame["ood_percentile"] >= 0.95
    q = float(bundle["conformal_q"])
    frame["prediction_interval_low"] = np.clip(frame["repairability_extra_trees"] - q, 0.0, 1.0)
    frame["prediction_interval_high"] = np.clip(frame["repairability_extra_trees"] + q, 0.0, 1.0)
    frame["prediction_interval_width"] = (
        frame["prediction_interval_high"] - frame["prediction_interval_low"]
    )
    site_confidence = pd.to_numeric(frame.get("site_plddt"), errors="coerce").clip(0, 100) / 100
    neighbor_confidence = (
        pd.to_numeric(frame.get("neighbor_mean_plddt_8a"), errors="coerce").clip(0, 100) / 100
    )
    low_confidence = pd.to_numeric(
        frame.get("low_confidence_neighbor_fraction_8a"), errors="coerce"
    ).clip(0, 1)
    frame["structure_quality_score"] = (
        0.5 * site_confidence.fillna(0)
        + 0.3 * neighbor_confidence.fillna(0)
        + 0.2 * (1.0 - low_confidence.fillna(1.0))
    )
    frame["applicability_status"] = np.select(
        [
            (frame["ood_percentile"] >= 0.99) | (site_confidence < 0.70),
            (frame["ood_percentile"] >= 0.95) | (site_confidence < 0.80),
        ],
        ["out_of_domain", "borderline"],
        default="in_domain",
    )

    rank_columns = []
    for name in MODEL_NAMES:
        rank_column = f"rank_percentile_{name}"
        frame[rank_column] = frame.groupby("species")[f"repairability_{name}"].rank(
            pct=True, method="average"
        )
        rank_columns.append(rank_column)
    frame["consensus_rank_percentile"] = frame[rank_columns].mean(axis=1)
    frame["model_disagreement"] = frame[[f"repairability_{name}" for name in MODEL_NAMES]].std(
        axis=1
    )
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(output_csv, index=False)
    return frame
