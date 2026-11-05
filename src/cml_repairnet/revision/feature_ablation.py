from __future__ import annotations

import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import numpy as np
import pandas as pd
from sklearn.base import clone
from sklearn.dummy import DummyRegressor
from sklearn.ensemble import ExtraTreesRegressor
from sklearn.impute import SimpleImputer
from sklearn.linear_model import Ridge
from sklearn.model_selection import GroupKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

from cml_repairnet.evaluation.metrics import regression_metrics
from cml_repairnet.evaluation.robustness import build_feature_view
from cml_repairnet.models.baselines import _numeric_feature_columns, _spatial_groups
from cml_repairnet.models.stage25 import CorrelationPruner, build_sample_weights
from cml_repairnet.models.stage3 import (
    generate_unique_spatial_partitions,
    partition_fingerprint,
)

KEY_COLUMNS = ("protein_id", "site_uniprot")
COORDINATE_COLUMNS = ("ca_x", "ca_y", "ca_z")
BFACTOR_COLUMNS = ("bfactor_ca", "neighbor_mean_bfactor_8a")
DEFAULT_TARGET = "repair_published"
DEFAULT_SEED = 20260727
REVISION_METADATA_COLUMNS = {
    "spatial_group", "outer_group", "outer_fold", "fold", "fold_id",
    "partition_id", "partition_seed", "cluster", "cluster_id",
}
DEFAULT_MODEL_WEIGHT_MODES: Mapping[str, str] = {
    "dummy_mean": "none",
    "ridge": "none",
    "extra_trees": "none",
    "lightgbm": "clipped",
}


@dataclass(frozen=True)
class SpatialPartition:
    partition_id: str
    groups: np.ndarray
    method: str
    seed: int
    fingerprint: str
    n_groups: int


@dataclass(frozen=True)
class FeatureFamily:
    name: str
    frame: pd.DataFrame
    feature_columns: tuple[str, ...]
    structure_columns: tuple[str, ...]
    description: str


def _validate_unique(frame: pd.DataFrame, name: str) -> None:
    missing = set(KEY_COLUMNS).difference(frame.columns)
    if missing:
        raise ValueError(f"{name} is missing key columns: {sorted(missing)}")
    duplicated = frame.duplicated(list(KEY_COLUMNS), keep=False)
    if duplicated.any():
        examples = frame.loc[duplicated, list(KEY_COLUMNS)].head(5).to_dict("records")
        raise ValueError(f"{name} has duplicated site keys; examples={examples}")


def load_bsa_tables(
    labels_path: Path,
    sequence_path: Path,
    structure_path: Path,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    for path in (labels_path, sequence_path, structure_path):
        if not path.exists():
            raise FileNotFoundError(f"Required input was not found: {path}")
    labels = pd.read_csv(labels_path)
    sequence = pd.read_csv(sequence_path)
    structure = pd.read_csv(structure_path)
    _validate_unique(labels, "labels")
    _validate_unique(sequence, "sequence features")
    _validate_unique(structure, "structure features")
    return labels, sequence, structure


def add_contact_flexibility_proxies(frame: pd.DataFrame) -> pd.DataFrame:
    """Add transfer-safe, contact-based flexibility *proxies*.

    These descriptors are deliberately named proxies: they do not claim to be experimental
    B-factors, molecular-dynamics fluctuations or AlphaFold pLDDT. They quantify loose packing,
    shell growth and exposure from static coordinates and are evaluated only through ablation.
    """

    out = frame.copy()
    required = {
        "ca_contact_count_8a",
        "ca_contact_count_10a",
        "ca_contact_count_12a",
        "relative_sasa_lysine",
    }
    missing = required.difference(out.columns)
    if missing:
        raise ValueError(
            "Cannot construct contact-based flexibility proxies; missing columns: "
            f"{sorted(missing)}"
        )

    c8 = pd.to_numeric(out["ca_contact_count_8a"], errors="coerce")
    c10 = pd.to_numeric(out["ca_contact_count_10a"], errors="coerce")
    c12 = pd.to_numeric(out["ca_contact_count_12a"], errors="coerce")
    rsa = pd.to_numeric(out["relative_sasa_lysine"], errors="coerce").clip(0.0, 1.5)

    out["flex_proxy_inverse_contacts_8a"] = 1.0 / (1.0 + c8)
    out["flex_proxy_inverse_contacts_10a"] = 1.0 / (1.0 + c10)
    out["flex_proxy_shell_growth_8_to_12a"] = (c12 - c8).clip(lower=0.0)
    out["flex_proxy_contact_growth_per_angstrom"] = (c12 - c8) / 4.0
    out["flex_proxy_exposure_per_contact"] = rsa / (1.0 + c8)
    out["flex_proxy_burial_packing"] = (1.0 - rsa.clip(upper=1.0)) * np.log1p(c8)
    return out


def _feature_columns(
    frame: pd.DataFrame,
    *,
    excluded: Iterable[str] = (),
) -> tuple[str, ...]:
    excluded_set = set(excluded).union(REVISION_METADATA_COLUMNS)
    columns = tuple(
        column
        for column in _numeric_feature_columns(frame)
        if column not in excluded_set
    )
    if not columns:
        raise ValueError("No numeric model features remain after exclusions")
    return columns


def _structure_numeric_columns(
    structure: pd.DataFrame,
    *,
    excluded: Iterable[str] = (),
) -> tuple[str, ...]:
    excluded_set = set(excluded).union(COORDINATE_COLUMNS).union(REVISION_METADATA_COLUMNS)
    return tuple(
        column
        for column in _numeric_feature_columns(structure)
        if column not in excluded_set
    )


def build_feature_families(
    labels: pd.DataFrame,
    sequence: pd.DataFrame,
    structure: pd.DataFrame,
) -> dict[str, FeatureFamily]:
    """Build prespecified sequence/structure/flexibility ablation families."""

    sequence_frame = build_feature_view(labels, sequence, structure, "sequence")
    structure_frame = build_feature_view(labels, sequence, structure, "structure")
    combined_frame = build_feature_view(labels, sequence, structure, "combined")
    structure_proxy_frame = add_contact_flexibility_proxies(structure_frame)
    combined_proxy_frame = add_contact_flexibility_proxies(combined_frame)

    sequence_columns = set(_numeric_feature_columns(sequence))
    structure_full_columns = set(_structure_numeric_columns(structure))
    structure_transfer_columns = set(
        _structure_numeric_columns(structure, excluded=BFACTOR_COLUMNS)
    )
    proxy_columns = {
        column for column in combined_proxy_frame.columns if column.startswith("flex_proxy_")
    }

    families = {
        "sequence": FeatureFamily(
            name="sequence",
            frame=sequence_frame,
            feature_columns=_feature_columns(sequence_frame),
            structure_columns=tuple(),
            description="Local sequence descriptors only; coordinates retained solely for blocking.",
        ),
        "structure_full": FeatureFamily(
            name="structure_full",
            frame=structure_frame,
            feature_columns=_feature_columns(structure_frame),
            structure_columns=tuple(sorted(structure_full_columns)),
            description="Experimental BSA structure descriptors including crystallographic B-factors.",
        ),
        "structure_transfer_safe": FeatureFamily(
            name="structure_transfer_safe",
            frame=structure_frame,
            feature_columns=_feature_columns(structure_frame, excluded=BFACTOR_COLUMNS),
            structure_columns=tuple(sorted(structure_transfer_columns)),
            description="Structure descriptors excluding B-factor variables unavailable in AlphaFold transfer.",
        ),
        "structure_transfer_safe_plus_proxy": FeatureFamily(
            name="structure_transfer_safe_plus_proxy",
            frame=structure_proxy_frame,
            feature_columns=_feature_columns(structure_proxy_frame, excluded=BFACTOR_COLUMNS),
            structure_columns=tuple(sorted(structure_transfer_columns.union(proxy_columns))),
            description=(
                "Transfer-safe structure-only descriptors plus static contact/exposure flexibility "
                "proxies; no sequence descriptors and no pLDDT-as-flexibility substitution."
            ),
        ),
        "structure_full_plus_proxy": FeatureFamily(
            name="structure_full_plus_proxy",
            frame=structure_proxy_frame,
            feature_columns=_feature_columns(structure_proxy_frame),
            structure_columns=tuple(sorted(structure_full_columns.union(proxy_columns))),
            description="Experimental structure descriptors including B-factors plus contact proxies.",
        ),
        "structure_bfactor_only": FeatureFamily(
            name="structure_bfactor_only",
            frame=structure_frame,
            feature_columns=tuple(column for column in BFACTOR_COLUMNS if column in structure_frame.columns),
            structure_columns=tuple(column for column in BFACTOR_COLUMNS if column in structure_frame.columns),
            description="Crystallographic B-factor variables only; mechanistic diagnostic, never transferred.",
        ),
        "structure_proxy_only": FeatureFamily(
            name="structure_proxy_only",
            frame=structure_proxy_frame,
            feature_columns=tuple(sorted(proxy_columns)),
            structure_columns=tuple(sorted(proxy_columns)),
            description="Static contact/exposure flexibility proxies only; diagnostic ablation.",
        ),
        "combined_full": FeatureFamily(
            name="combined_full",
            frame=combined_frame,
            feature_columns=_feature_columns(combined_frame),
            structure_columns=tuple(sorted(structure_full_columns)),
            description="Sequence plus all experimental BSA structure descriptors.",
        ),
        "combined_transfer_safe": FeatureFamily(
            name="combined_transfer_safe",
            frame=combined_frame,
            feature_columns=_feature_columns(combined_frame, excluded=BFACTOR_COLUMNS),
            structure_columns=tuple(sorted(structure_transfer_columns)),
            description="Sequence plus transfer-safe structure descriptors; B-factors excluded.",
        ),
        "combined_transfer_safe_plus_proxy": FeatureFamily(
            name="combined_transfer_safe_plus_proxy",
            frame=combined_proxy_frame,
            feature_columns=_feature_columns(combined_proxy_frame, excluded=BFACTOR_COLUMNS),
            structure_columns=tuple(sorted(structure_transfer_columns.union(proxy_columns))),
            description=(
                "Transfer-safe combined features plus static contact/exposure flexibility proxies; "
                "no pLDDT-as-flexibility substitution."
            ),
        ),
    }

    # Guard against accidental family contamination.
    if sequence_columns:
        sequence_selected = set(families["sequence"].feature_columns)
        unexpected = sequence_selected.difference(sequence_columns).difference(REVISION_METADATA_COLUMNS)
        if unexpected:
            raise ValueError(f"Sequence family contains unexpected features: {sorted(unexpected)}")
    for family_name in (
        "structure_transfer_safe",
        "structure_transfer_safe_plus_proxy",
        "structure_proxy_only",
        "combined_transfer_safe",
        "combined_transfer_safe_plus_proxy",
    ):
        leaked = set(BFACTOR_COLUMNS).intersection(families[family_name].feature_columns)
        if leaked:
            raise ValueError(f"Transfer-safe family {family_name} leaked B-factor columns: {sorted(leaked)}")
    return families


def make_frozen_estimator(
    model_name: str,
    *,
    seed: int = DEFAULT_SEED,
    n_jobs: int = 1,
) -> Pipeline:
    """Construct Stage-3-compatible frozen estimators without reopening model search."""

    name = model_name.lower()
    if name == "dummy_mean":
        return Pipeline([("model", DummyRegressor(strategy="mean"))])
    if name == "ridge":
        return Pipeline(
            [
                ("prune", CorrelationPruner(0.85)),
                ("imputer", SimpleImputer(strategy="median")),
                ("scale", StandardScaler()),
                ("model", Ridge(alpha=10.0)),
            ]
        )
    if name == "extra_trees":
        return Pipeline(
            [
                ("prune", CorrelationPruner(0.85)),
                ("imputer", SimpleImputer(strategy="median")),
                (
                    "model",
                    ExtraTreesRegressor(
                        n_estimators=500,
                        max_depth=4,
                        min_samples_leaf=2,
                        max_features=1.0,
                        random_state=int(seed),
                        n_jobs=int(n_jobs),
                    ),
                ),
            ]
        )
    if name == "lightgbm":
        try:
            from lightgbm import LGBMRegressor
        except ImportError as exc:  # pragma: no cover - environment dependent
            raise RuntimeError(
                "LightGBM is required for model='lightgbm'. Install project experiment extras: "
                "python -m pip install -e '.[dev,experiment]'"
            ) from exc
        return Pipeline(
            [
                ("prune", CorrelationPruner(0.85)),
                ("imputer", SimpleImputer(strategy="median")),
                (
                    "model",
                    LGBMRegressor(
                        objective="regression_l1",
                        n_estimators=60,
                        learning_rate=0.05,
                        num_leaves=3,
                        max_depth=2,
                        min_child_samples=4,
                        subsample=0.8,
                        colsample_bytree=0.8,
                        reg_alpha=0.25,
                        reg_lambda=10.0,
                        random_state=int(seed),
                        n_jobs=int(n_jobs),
                        verbosity=-1,
                        deterministic=True,
                        force_col_wise=True,
                        feature_pre_filter=False,
                        subsample_freq=1,
                    ),
                ),
            ]
        )
    raise ValueError("model_name must be one of: dummy_mean, ridge, extra_trees, lightgbm")


def build_spatial_partitions(
    frame: pd.DataFrame,
    *,
    group_counts: Sequence[int] = (4,),
    random_partition_count: int = 5,
    start_seed: int = 20260729,
    include_ward_anchor: bool = True,
) -> list[SpatialPartition]:
    partitions: list[SpatialPartition] = []
    for n_groups in group_counts:
        if include_ward_anchor:
            groups = _spatial_groups(frame, int(n_groups))
            partitions.append(
                SpatialPartition(
                    partition_id=f"ward_g{n_groups}",
                    groups=np.asarray(groups),
                    method="ward",
                    seed=-1,
                    fingerprint=partition_fingerprint(groups),
                    n_groups=int(n_groups),
                )
            )
        if int(random_partition_count) > 0:
            generated = generate_unique_spatial_partitions(
                frame,
                int(n_groups),
                start_seed=int(start_seed),
                target_count=int(random_partition_count),
                min_group_size=2,
            )
            for seed, groups, fingerprint in generated:
                partitions.append(
                    SpatialPartition(
                        partition_id=f"kmeans_random_seed{seed}_g{n_groups}",
                        groups=np.asarray(groups),
                        method="kmeans_random",
                        seed=int(seed),
                        fingerprint=fingerprint,
                        n_groups=int(n_groups),
                    )
                )
    return partitions


def _fit_pipeline(
    estimator: Pipeline,
    X: pd.DataFrame,
    y: np.ndarray,
    weights: np.ndarray,
    model_name: str,
) -> Pipeline:
    if model_name == "dummy_mean":
        return estimator.fit(X, y)
    return estimator.fit(X, y, model__sample_weight=weights)


def _mean_selected_features(fitted_models: Sequence[Pipeline], n_input: int) -> float:
    counts: list[int] = []
    for fitted in fitted_models:
        pruner = fitted.named_steps.get("prune")
        if pruner is None:
            counts.append(n_input)
        else:
            selected = getattr(pruner, "selected_feature_names_", None)
            counts.append(len(selected) if selected is not None else n_input)
    return float(np.mean(counts)) if counts else float("nan")


def evaluate_feature_family_partition(
    family: FeatureFamily,
    partition: SpatialPartition,
    *,
    target: str,
    model_names: Sequence[str],
    seed: int = DEFAULT_SEED,
    n_jobs: int = 1,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    frame = family.frame
    if target not in frame.columns:
        raise ValueError(f"Target {target!r} is absent from feature family {family.name!r}")
    y = pd.to_numeric(frame[target], errors="raise").to_numpy(dtype=float)
    groups = np.asarray(partition.groups)
    if len(groups) != len(frame):
        raise ValueError("Partition length does not match feature frame")
    X = frame.loc[:, list(family.feature_columns)].copy()
    outer_cv = GroupKFold(n_splits=np.unique(groups).size)

    metric_rows: list[dict[str, Any]] = []
    prediction_rows: list[pd.DataFrame] = []
    for model_name in model_names:
        predictions = np.full(len(frame), np.nan, dtype=float)
        fold_ids = np.full(len(frame), -1, dtype=int)
        fitted_models: list[Pipeline] = []
        for fold, (train_idx, test_idx) in enumerate(outer_cv.split(X, y, groups), start=1):
            estimator = make_frozen_estimator(model_name, seed=seed, n_jobs=n_jobs)
            mode = DEFAULT_MODEL_WEIGHT_MODES.get(model_name, "none")
            weights = build_sample_weights(frame.iloc[train_idx], mode)
            fitted = _fit_pipeline(
                estimator,
                X.iloc[train_idx],
                y[train_idx],
                weights,
                model_name,
            )
            fitted_models.append(fitted)
            predictions[test_idx] = np.clip(fitted.predict(X.iloc[test_idx]), 0.0, 1.0)
            fold_ids[test_idx] = fold
        if not np.isfinite(predictions).all():
            raise RuntimeError(f"Incomplete OOF predictions for {family.name}/{model_name}")
        metrics = regression_metrics(y, predictions)
        metric_rows.append(
            {
                "partition_id": partition.partition_id,
                "partition_method": partition.method,
                "partition_seed": partition.seed,
                "partition_fingerprint": partition.fingerprint,
                "n_groups": partition.n_groups,
                "feature_set": family.name,
                "model": model_name,
                "target": target,
                "n_sites": len(frame),
                "n_input_features": len(family.feature_columns),
                "mean_selected_features": _mean_selected_features(
                    fitted_models, len(family.feature_columns)
                ),
                "weight_mode": DEFAULT_MODEL_WEIGHT_MODES.get(model_name, "none"),
                **metrics,
            }
        )
        prediction_rows.append(
            pd.DataFrame(
                {
                    "protein_id": frame["protein_id"].astype(str),
                    "site_uniprot": frame["site_uniprot"],
                    "site_paper": frame.get("site_paper", pd.Series(np.nan, index=frame.index)),
                    "partition_id": partition.partition_id,
                    "partition_method": partition.method,
                    "outer_group": groups.astype(str),
                    "outer_fold": fold_ids,
                    "feature_set": family.name,
                    "model": model_name,
                    "target": target,
                    "y_true": y,
                    "y_pred": predictions,
                    "absolute_error": np.abs(y - predictions),
                }
            )
        )
    return pd.DataFrame(metric_rows), pd.concat(prediction_rows, ignore_index=True)


def _bootstrap_mean_ci(
    values: np.ndarray,
    *,
    iterations: int,
    seed: int,
) -> tuple[float, float]:
    finite = np.asarray(values, dtype=float)
    finite = finite[np.isfinite(finite)]
    if finite.size == 0:
        return float("nan"), float("nan")
    if finite.size == 1:
        return float(finite[0]), float(finite[0])
    rng = np.random.default_rng(seed)
    draws = np.empty(iterations, dtype=float)
    for i in range(iterations):
        draws[i] = rng.choice(finite, size=finite.size, replace=True).mean()
    low, high = np.quantile(draws, [0.025, 0.975])
    return float(low), float(high)


def summarize_ablation_metrics(
    metrics: pd.DataFrame,
    *,
    bootstrap_iterations: int = 10_000,
    seed: int = DEFAULT_SEED,
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for (feature_set, model), group in metrics.groupby(["feature_set", "model"], sort=True):
        mae_values = group["mae"].to_numpy(dtype=float)
        mae_low, mae_high = _bootstrap_mean_ci(
            mae_values, iterations=bootstrap_iterations, seed=seed
        )
        rho_low, rho_high = _bootstrap_mean_ci(
            group["spearman_rho"].to_numpy(dtype=float),
            iterations=bootstrap_iterations,
            seed=seed + 1,
        )
        spearman_mean = float(np.nanmean(group["spearman_rho"]))
        if model == "dummy_mean":
            # Fold-wise training means can vary across held-out groups, producing a numerical
            # rank correlation that has no meaningful null-model interpretation.
            spearman_mean = float("nan")
            rho_low = float("nan")
            rho_high = float("nan")
        rows.append(
            {
                "feature_set": feature_set,
                "model": model,
                "n_partitions": len(group),
                "mae_mean": float(np.mean(mae_values)),
                "mae_sd": float(np.std(mae_values, ddof=1)) if len(group) > 1 else 0.0,
                "mae_ci_2.5": mae_low,
                "mae_ci_97.5": mae_high,
                "spearman_mean": spearman_mean,
                "spearman_ci_2.5": rho_low,
                "spearman_ci_97.5": rho_high,
                "r2_mean": float(np.nanmean(group["r2"])),
                "mean_selected_features": float(np.nanmean(group["mean_selected_features"])),
            }
        )
    return pd.DataFrame(rows).sort_values(["model", "mae_mean"]).reset_index(drop=True)


def build_pairwise_feature_deltas(
    metrics: pd.DataFrame,
    *,
    bootstrap_iterations: int = 10_000,
    seed: int = DEFAULT_SEED,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    comparisons = (
        ("sequence", "combined_transfer_safe", "structure_increment_over_sequence"),
        ("structure_transfer_safe", "combined_transfer_safe", "sequence_increment_over_structure"),
        (
            "structure_transfer_safe",
            "structure_transfer_safe_plus_proxy",
            "contact_flexibility_proxy_increment_structure_only",
        ),
        (
            "combined_transfer_safe",
            "combined_transfer_safe_plus_proxy",
            "contact_flexibility_proxy_increment_combined",
        ),
        (
            "structure_transfer_safe",
            "structure_full",
            "crystallographic_bfactor_increment_structure_only",
        ),
        (
            "combined_transfer_safe",
            "combined_full",
            "crystallographic_bfactor_increment_combined",
        ),
        (
            "structure_transfer_safe_plus_proxy",
            "structure_full_plus_proxy",
            "crystallographic_bfactor_increment_with_proxy",
        ),
    )
    detail_rows: list[dict[str, Any]] = []
    summary_rows: list[dict[str, Any]] = []
    for baseline, candidate, comparison_name in comparisons:
        left = metrics.loc[metrics["feature_set"] == baseline].copy()
        right = metrics.loc[metrics["feature_set"] == candidate].copy()
        paired = left.merge(
            right,
            on=["partition_id", "model", "target", "n_groups"],
            suffixes=("_baseline", "_candidate"),
            validate="one_to_one",
        )
        if paired.empty:
            continue
        paired["delta_mae"] = paired["mae_baseline"] - paired["mae_candidate"]
        paired["comparison"] = comparison_name
        paired["baseline_feature_set"] = baseline
        paired["candidate_feature_set"] = candidate
        for row in paired.itertuples(index=False):
            detail_rows.append(
                {
                    "comparison": comparison_name,
                    "baseline_feature_set": baseline,
                    "candidate_feature_set": candidate,
                    "partition_id": row.partition_id,
                    "model": row.model,
                    "target": row.target,
                    "n_groups": row.n_groups,
                    "baseline_mae": row.mae_baseline,
                    "candidate_mae": row.mae_candidate,
                    "delta_mae": row.delta_mae,
                }
            )
        for model, model_group in paired.groupby("model"):
            values = model_group["delta_mae"].to_numpy(dtype=float)
            low, high = _bootstrap_mean_ci(
                values, iterations=bootstrap_iterations, seed=seed
            )
            summary_rows.append(
                {
                    "comparison": comparison_name,
                    "baseline_feature_set": baseline,
                    "candidate_feature_set": candidate,
                    "model": model,
                    "n_partitions": len(values),
                    "mean_delta_mae": float(np.mean(values)),
                    "median_delta_mae": float(np.median(values)),
                    "delta_mae_ci_2.5": low,
                    "delta_mae_ci_97.5": high,
                    "candidate_win_fraction": float(np.mean(values > 0.0)),
                    "interpretation": "positive_delta_favors_candidate",
                }
            )
    return pd.DataFrame(detail_rows), pd.DataFrame(summary_rows)


def _partition_manifest(partitions: Sequence[SpatialPartition]) -> pd.DataFrame:
    rows = []
    for partition in partitions:
        _, counts = np.unique(partition.groups, return_counts=True)
        rows.append(
            {
                "partition_id": partition.partition_id,
                "method": partition.method,
                "seed": partition.seed,
                "fingerprint": partition.fingerprint,
                "n_groups": partition.n_groups,
                "min_group_size": int(counts.min()),
                "max_group_size": int(counts.max()),
            }
        )
    return pd.DataFrame(rows)


def run_feature_family_ablation(
    labels: pd.DataFrame,
    sequence: pd.DataFrame,
    structure: pd.DataFrame,
    output_dir: Path,
    *,
    target: str = DEFAULT_TARGET,
    model_names: Sequence[str] = ("extra_trees", "ridge", "lightgbm"),
    family_names: Sequence[str] | None = None,
    group_counts: Sequence[int] = (4,),
    random_partition_count: int = 5,
    start_seed: int = 20260729,
    n_jobs: int = 1,
    bootstrap_iterations: int = 10_000,
    seed: int = DEFAULT_SEED,
) -> dict[str, pd.DataFrame]:
    output_dir.mkdir(parents=True, exist_ok=True)
    families = build_feature_families(labels, sequence, structure)
    if family_names is not None:
        unknown = set(family_names).difference(families)
        if unknown:
            raise ValueError(f"Unknown feature families: {sorted(unknown)}")
        families = {name: families[name] for name in family_names}
    master_family = (
        families.get("combined_transfer_safe")
        or families.get("structure_transfer_safe")
        or next(iter(families.values()))
    )
    master = master_family.frame
    partitions = build_spatial_partitions(
        master,
        group_counts=group_counts,
        random_partition_count=random_partition_count,
        start_seed=start_seed,
        include_ward_anchor=True,
    )

    metric_frames: list[pd.DataFrame] = []
    prediction_frames: list[pd.DataFrame] = []
    for partition in partitions:
        # Dummy is evaluated once per partition and is independent of feature family.
        dummy_family = FeatureFamily(
            name="null_mean",
            frame=master.assign(_constant_feature=0.0),
            feature_columns=("_constant_feature",),
            structure_columns=tuple(),
            description="Training-fold mean null model.",
        )
        metrics, predictions = evaluate_feature_family_partition(
            dummy_family,
            partition,
            target=target,
            model_names=("dummy_mean",),
            seed=seed,
            n_jobs=n_jobs,
        )
        metric_frames.append(metrics)
        prediction_frames.append(predictions)

        for family in families.values():
            metrics, predictions = evaluate_feature_family_partition(
                family,
                partition,
                target=target,
                model_names=model_names,
                seed=seed,
                n_jobs=n_jobs,
            )
            metric_frames.append(metrics)
            prediction_frames.append(predictions)

    all_metrics = pd.concat(metric_frames, ignore_index=True)
    all_predictions = pd.concat(prediction_frames, ignore_index=True)
    summary = summarize_ablation_metrics(
        all_metrics,
        bootstrap_iterations=bootstrap_iterations,
        seed=seed,
    )
    pairwise_detail, pairwise_summary = build_pairwise_feature_deltas(
        all_metrics,
        bootstrap_iterations=bootstrap_iterations,
        seed=seed,
    )
    manifest = _partition_manifest(partitions)
    family_manifest = pd.DataFrame(
        [
            {
                "feature_set": family.name,
                "n_input_features": len(family.feature_columns),
                "n_structure_features": len(family.structure_columns),
                "description": family.description,
                "features": "|".join(family.feature_columns),
            }
            for family in families.values()
        ]
    )

    outputs = {
        "metrics": all_metrics,
        "predictions": all_predictions,
        "summary": summary,
        "pairwise_detail": pairwise_detail,
        "pairwise_summary": pairwise_summary,
        "partition_manifest": manifest,
        "feature_manifest": family_manifest,
    }
    filenames = {
        "metrics": "revision_t1_ablation_partition_metrics.csv",
        "predictions": "revision_t1_ablation_oof_predictions.csv",
        "summary": "revision_t1_ablation_summary.csv",
        "pairwise_detail": "revision_t1_pairwise_deltas.csv",
        "pairwise_summary": "revision_t1_pairwise_summary.csv",
        "partition_manifest": "revision_t1_partition_manifest.csv",
        "feature_manifest": "revision_t1_feature_manifest.csv",
    }
    for key, frame in outputs.items():
        frame.to_csv(output_dir / filenames[key], index=False)

    config = {
        "analysis": "reviewer_revision_feature_family_ablation",
        "target": target,
        "models": list(model_names),
        "feature_families": list(families),
        "group_counts": list(group_counts),
        "random_partition_count_per_group_count": random_partition_count,
        "ward_anchor_included": True,
        "start_seed": start_seed,
        "model_seed": seed,
        "bootstrap_iterations": bootstrap_iterations,
        "b_factor_columns_excluded_from_transfer_safe": list(BFACTOR_COLUMNS),
        "n_sites": len(labels),
    }
    (output_dir / "revision_t1_config.json").write_text(
        json.dumps(config, indent=2, sort_keys=True), encoding="utf-8"
    )
    return outputs


def _permuted_oof_mae(
    family: FeatureFamily,
    partition: SpatialPartition,
    *,
    target: str,
    model_name: str,
    iteration: int,
    seed: int,
    n_jobs: int,
) -> float:
    frame = family.frame
    X = frame.loc[:, list(family.feature_columns)].copy()
    y = pd.to_numeric(frame[target], errors="raise").to_numpy(dtype=float)
    groups = np.asarray(partition.groups)
    structure_columns = [column for column in family.structure_columns if column in X.columns]
    if not structure_columns:
        raise ValueError(f"Feature family {family.name!r} has no structure columns to permute")

    predictions = np.full(len(frame), np.nan, dtype=float)
    cv = GroupKFold(n_splits=np.unique(groups).size)
    for fold, (train_idx, test_idx) in enumerate(cv.split(X, y, groups), start=1):
        X_train = X.iloc[train_idx].copy()
        X_test = X.iloc[test_idx]
        rng = np.random.default_rng(int(seed) + int(iteration) * 1009 + fold * 7919)
        row_permutation = rng.permutation(len(train_idx))
        # One row permutation is shared by all structure columns, preserving their joint
        # covariance while destroying site-to-structure correspondence in the training fold.
        X_train.loc[:, structure_columns] = (
            X_train.loc[:, structure_columns].iloc[row_permutation].to_numpy()
        )
        estimator = make_frozen_estimator(model_name, seed=seed, n_jobs=n_jobs)
        mode = DEFAULT_MODEL_WEIGHT_MODES.get(model_name, "none")
        weights = build_sample_weights(frame.iloc[train_idx], mode)
        fitted = _fit_pipeline(estimator, X_train, y[train_idx], weights, model_name)
        predictions[test_idx] = np.clip(fitted.predict(X_test), 0.0, 1.0)
    return float(np.mean(np.abs(y - predictions)))


def run_structure_permutation_control(
    labels: pd.DataFrame,
    sequence: pd.DataFrame,
    structure: pd.DataFrame,
    output_dir: Path,
    *,
    target: str = DEFAULT_TARGET,
    feature_set: str = "combined_transfer_safe",
    model_name: str = "extra_trees",
    iterations: int = 200,
    n_groups: int = 4,
    partition_method: str = "ward",
    partition_seed: int = 20260729,
    n_jobs: int = 1,
    seed: int = DEFAULT_SEED,
) -> dict[str, pd.DataFrame]:
    if iterations < 1:
        raise ValueError("iterations must be at least 1")
    output_dir.mkdir(parents=True, exist_ok=True)
    families = build_feature_families(labels, sequence, structure)
    if feature_set not in families:
        raise ValueError(f"Unknown feature_set {feature_set!r}; choose from {sorted(families)}")
    family = families[feature_set]

    if partition_method == "ward":
        groups = _spatial_groups(family.frame, int(n_groups))
        partition = SpatialPartition(
            partition_id=f"ward_g{n_groups}",
            groups=np.asarray(groups),
            method="ward",
            seed=-1,
            fingerprint=partition_fingerprint(groups),
            n_groups=int(n_groups),
        )
    elif partition_method == "kmeans_random":
        generated = generate_unique_spatial_partitions(
            family.frame,
            int(n_groups),
            start_seed=int(partition_seed),
            target_count=1,
            min_group_size=2,
        )[0]
        actual_seed, groups, fingerprint = generated
        partition = SpatialPartition(
            partition_id=f"kmeans_random_seed{actual_seed}_g{n_groups}",
            groups=np.asarray(groups),
            method="kmeans_random",
            seed=int(actual_seed),
            fingerprint=fingerprint,
            n_groups=int(n_groups),
        )
    else:
        raise ValueError("partition_method must be 'ward' or 'kmeans_random'")

    real_metrics, real_predictions = evaluate_feature_family_partition(
        family,
        partition,
        target=target,
        model_names=(model_name,),
        seed=seed,
        n_jobs=n_jobs,
    )
    real_mae = float(real_metrics.iloc[0]["mae"])

    rows: list[dict[str, Any]] = []
    for iteration in range(iterations):
        permuted_mae = _permuted_oof_mae(
            family,
            partition,
            target=target,
            model_name=model_name,
            iteration=iteration,
            seed=seed,
            n_jobs=n_jobs,
        )
        rows.append(
            {
                "iteration": iteration,
                "feature_set": feature_set,
                "model": model_name,
                "partition_id": partition.partition_id,
                "real_mae": real_mae,
                "permuted_structure_mae": permuted_mae,
                "delta_mae_permuted_minus_real": permuted_mae - real_mae,
            }
        )
    iterations_frame = pd.DataFrame(rows)
    deltas = iterations_frame["delta_mae_permuted_minus_real"].to_numpy(dtype=float)
    p_value = float((1 + np.sum(deltas <= 0.0)) / (iterations + 1))
    delta_low, delta_high = np.quantile(deltas, [0.025, 0.975])
    summary = pd.DataFrame(
        [
            {
                "feature_set": feature_set,
                "model": model_name,
                "partition_id": partition.partition_id,
                "partition_fingerprint": partition.fingerprint,
                "n_groups": n_groups,
                "iterations": iterations,
                "real_mae": real_mae,
                "permuted_mae_mean": float(np.mean(iterations_frame["permuted_structure_mae"])),
                "permuted_mae_median": float(np.median(iterations_frame["permuted_structure_mae"])),
                "mean_delta_mae": float(np.mean(deltas)),
                "delta_mae_ci_2.5": float(delta_low),
                "delta_mae_ci_97.5": float(delta_high),
                "empirical_p_one_sided": p_value,
                "null_definition": (
                    "training-fold structure rows jointly permuted; sequence and test folds unchanged"
                ),
                "interpretation": (
                    "positive delta means correctly mapped structure outperforms the permuted control"
                ),
            }
        ]
    )

    iterations_frame.to_csv(
        output_dir / "revision_t2_structure_permutation_iterations.csv", index=False
    )
    summary.to_csv(
        output_dir / "revision_t2_structure_permutation_summary.csv", index=False
    )
    real_predictions.to_csv(
        output_dir / "revision_t2_real_structure_oof_predictions.csv", index=False
    )
    config = {
        "analysis": "reviewer_revision_permuted_structure_negative_control",
        "target": target,
        "feature_set": feature_set,
        "model": model_name,
        "iterations": iterations,
        "n_groups": n_groups,
        "partition_method": partition_method,
        "partition_seed": partition.seed,
        "model_seed": seed,
        "permutation_unit": "training rows, shared across all structure columns within fold",
        "test_structure_features": "left unmodified",
    }
    (output_dir / "revision_t2_config.json").write_text(
        json.dumps(config, indent=2, sort_keys=True), encoding="utf-8"
    )
    return {
        "iterations": iterations_frame,
        "summary": summary,
        "real_predictions": real_predictions,
    }


def write_revision_report(
    output_dir: Path,
    *,
    ablation_summary: pd.DataFrame | None = None,
    pairwise_summary: pd.DataFrame | None = None,
    permutation_summary: pd.DataFrame | None = None,
) -> Path:
    lines = [
        "# CML-RepairNet Reviewer Revision Tests",
        "",
        "This report is generated from frozen, leakage-controlled spatial validation. Positive",
        "pairwise delta-MAE values favor the candidate feature family. Contact-based flexibility",
        "variables are proxies from static geometry and are not interpreted as pLDDT or measured dynamics.",
        "",
    ]
    if ablation_summary is not None and not ablation_summary.empty:
        lines.extend(
            [
                "## T1 — Feature-family ablation",
                "",
                "```text",
                ablation_summary.to_string(index=False),
                "```",
                "",
            ]
        )
    if pairwise_summary is not None and not pairwise_summary.empty:
        lines.extend(
            [
                "## T1 — Prespecified paired contrasts",
                "",
                "```text",
                pairwise_summary.to_string(index=False),
                "```",
                "",
            ]
        )
    if permutation_summary is not None and not permutation_summary.empty:
        lines.extend(
            [
                "## T2 — Permuted-structure negative control",
                "",
                "```text",
                permutation_summary.to_string(index=False),
                "```",
                "",
            ]
        )
    lines.extend(
        [
            "## Decision rules",
            "",
            "- Structure adds evidence only when the sequence → combined contrast has positive",
            "  delta-MAE across most partitions and its bootstrap interval is predominantly above zero.",
            "- The contact-flexibility proxy is retained only if it improves transfer-safe combined",
            "  features without relying on crystallographic B-factors.",
            "- The structure-aware claim is weakened if permuted structure performs as well as the",
            "  correctly mapped structure (one-sided empirical P not supportive).",
            "- All conclusions remain BSA/CrGO-897 assay-specific; these tests do not establish",
            "  cross-protein experimental calibration.",
            "",
        ]
    )
    path = output_dir / "revision_tests_report.md"
    path.write_text("\n".join(lines), encoding="utf-8")
    return path
