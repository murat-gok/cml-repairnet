from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Sequence

import numpy as np
import pandas as pd
from scipy.stats import binomtest, spearmanr
from sklearn.model_selection import GroupKFold

from cml_repairnet.evaluation.metrics import regression_metrics
from cml_repairnet.models.stage25 import build_sample_weights
from cml_repairnet.revision.feature_ablation import (
    DEFAULT_MODEL_WEIGHT_MODES,
    DEFAULT_SEED,
    FeatureFamily,
    SpatialPartition,
    _fit_pipeline,
    _permuted_oof_mae,
    _bootstrap_mean_ci,
    build_feature_families,
    build_spatial_partitions,
    evaluate_feature_family_partition,
    make_frozen_estimator,
)


def run_multi_partition_permutation_control(
    labels: pd.DataFrame,
    sequence: pd.DataFrame,
    structure: pd.DataFrame,
    output_dir: Path,
    *,
    target: str = "repair_published",
    feature_set: str = "structure_transfer_safe",
    model_name: str = "extra_trees",
    iterations: int = 2_000,
    group_counts: Sequence[int] = (4,),
    random_partition_count: int = 5,
    start_seed: int = 20260729,
    n_jobs: int = 1,
    seed: int = DEFAULT_SEED,
    resume: bool = True,
    checkpoint_every: int = 25,
) -> dict[str, pd.DataFrame]:
    """Run the structure-row permutation control over every prespecified partition.

    Each training fold receives one joint row permutation across all structure columns.
    Test-fold features, labels and sequence descriptors are never permuted. The aggregate
    null statistic is the mean delta-MAE across partitions within each permutation index.
    """

    if iterations < 1:
        raise ValueError("iterations must be at least 1")
    output_dir.mkdir(parents=True, exist_ok=True)
    families = build_feature_families(labels, sequence, structure)
    if feature_set not in families:
        raise ValueError(f"Unknown feature_set {feature_set!r}; choose from {sorted(families)}")
    family = families[feature_set]
    partitions = build_spatial_partitions(
        family.frame,
        group_counts=group_counts,
        random_partition_count=random_partition_count,
        start_seed=start_seed,
        include_ward_anchor=True,
    )

    checkpoint_path = output_dir / "revision_t2_full_structure_permutation_iterations.csv"
    existing = pd.DataFrame()
    if resume and checkpoint_path.exists():
        existing = pd.read_csv(checkpoint_path)
        required = {"iteration", "partition_id", "feature_set", "model"}
        if not required.issubset(existing.columns):
            existing = pd.DataFrame()
        else:
            allowed_partitions = {partition.partition_id for partition in partitions}
            existing = existing.loc[
                existing["partition_id"].isin(allowed_partitions)
                & existing["feature_set"].eq(feature_set)
                & existing["model"].eq(model_name)
                & existing["iteration"].lt(iterations)
            ].copy()

    iteration_rows: list[dict[str, Any]] = existing.to_dict("records") if not existing.empty else []
    partition_rows: list[dict[str, Any]] = []
    prediction_frames: list[pd.DataFrame] = []
    for partition_index, partition in enumerate(partitions):
        real_metrics, real_predictions = evaluate_feature_family_partition(
            family,
            partition,
            target=target,
            model_names=(model_name,),
            seed=seed,
            n_jobs=n_jobs,
        )
        real_mae = float(real_metrics.iloc[0]["mae"])
        prediction_frames.append(real_predictions)
        existing_partition = existing.loc[existing["partition_id"].eq(partition.partition_id)] if not existing.empty else pd.DataFrame()
        completed = set(existing_partition["iteration"].astype(int).tolist()) if not existing_partition.empty else set()
        deltas_by_iteration = {
            int(row.iteration): float(row.delta_mae_permuted_minus_real)
            for row in existing_partition.itertuples(index=False)
        }
        since_checkpoint = 0
        for iteration in range(iterations):
            if iteration in completed:
                continue
            permuted_mae = _permuted_oof_mae(
                family,
                partition,
                target=target,
                model_name=model_name,
                iteration=iteration,
                seed=seed + partition_index * 1_000_003,
                n_jobs=n_jobs,
            )
            delta = float(permuted_mae - real_mae)
            deltas_by_iteration[iteration] = delta
            iteration_rows.append(
                {
                    "iteration": iteration,
                    "partition_id": partition.partition_id,
                    "partition_method": partition.method,
                    "partition_seed": partition.seed,
                    "partition_fingerprint": partition.fingerprint,
                    "n_groups": partition.n_groups,
                    "feature_set": feature_set,
                    "model": model_name,
                    "real_mae": real_mae,
                    "permuted_structure_mae": permuted_mae,
                    "delta_mae_permuted_minus_real": delta,
                }
            )
            since_checkpoint += 1
            if checkpoint_every > 0 and since_checkpoint >= checkpoint_every:
                pd.DataFrame(iteration_rows).drop_duplicates(
                    ["partition_id", "iteration", "feature_set", "model"], keep="last"
                ).sort_values(["partition_id", "iteration"]).to_csv(checkpoint_path, index=False)
                since_checkpoint = 0
        delta_array = np.asarray(
            [deltas_by_iteration[i] for i in range(iterations)], dtype=float
        )
        low, high = np.quantile(delta_array, [0.025, 0.975])
        partition_rows.append(
            {
                "partition_id": partition.partition_id,
                "partition_method": partition.method,
                "partition_seed": partition.seed,
                "partition_fingerprint": partition.fingerprint,
                "n_groups": partition.n_groups,
                "feature_set": feature_set,
                "model": model_name,
                "iterations": iterations,
                "real_mae": real_mae,
                "permuted_mae_mean": float(real_mae + np.mean(delta_array)),
                "mean_delta_mae": float(np.mean(delta_array)),
                "median_delta_mae": float(np.median(delta_array)),
                "delta_mae_ci_2.5": float(low),
                "delta_mae_ci_97.5": float(high),
                "empirical_p_one_sided": float((1 + np.sum(delta_array <= 0.0)) / (iterations + 1)),
                "correct_mapping_win_fraction": float(np.mean(delta_array > 0.0)),
            }
        )

    iterations_frame = pd.DataFrame(iteration_rows).drop_duplicates(
        ["partition_id", "iteration", "feature_set", "model"], keep="last"
    ).sort_values(["partition_id", "iteration"]).reset_index(drop=True)
    iterations_frame.to_csv(checkpoint_path, index=False)
    partition_summary = pd.DataFrame(partition_rows)
    aggregate_draws = (
        iterations_frame.groupby("iteration", as_index=False)["delta_mae_permuted_minus_real"]
        .mean()
        .rename(columns={"delta_mae_permuted_minus_real": "mean_delta_across_partitions"})
    )
    aggregate_values = aggregate_draws["mean_delta_across_partitions"].to_numpy(dtype=float)
    agg_low, agg_high = np.quantile(aggregate_values, [0.025, 0.975])
    positive_partitions = int((partition_summary["mean_delta_mae"] > 0.0).sum())
    n_partitions = len(partition_summary)
    sign_p = float(binomtest(positive_partitions, n_partitions, p=0.5, alternative="greater").pvalue)
    aggregate_summary = pd.DataFrame(
        [
            {
                "feature_set": feature_set,
                "model": model_name,
                "n_partitions": n_partitions,
                "iterations_per_partition": iterations,
                "mean_real_mae": float(partition_summary["real_mae"].mean()),
                "mean_permuted_mae": float(partition_summary["permuted_mae_mean"].mean()),
                "mean_delta_mae": float(aggregate_values.mean()),
                "delta_mae_ci_2.5": float(agg_low),
                "delta_mae_ci_97.5": float(agg_high),
                "aggregate_empirical_p_one_sided": float(
                    (1 + np.sum(aggregate_values <= 0.0)) / (iterations + 1)
                ),
                "positive_partition_count": positive_partitions,
                "partition_win_fraction": float(positive_partitions / n_partitions),
                "partition_sign_test_p_one_sided": sign_p,
                "null_definition": (
                    "training-fold structure rows jointly permuted within every spatial partition; "
                    "test folds and labels unchanged"
                ),
            }
        ]
    )

    real_predictions = pd.concat(prediction_frames, ignore_index=True)
    outputs = {
        "iterations": iterations_frame,
        "partition_summary": partition_summary,
        "aggregate_draws": aggregate_draws,
        "aggregate_summary": aggregate_summary,
        "real_predictions": real_predictions,
    }
    filenames = {
        "iterations": "revision_t2_full_structure_permutation_iterations.csv",
        "partition_summary": "revision_t2_full_partition_summary.csv",
        "aggregate_draws": "revision_t2_full_aggregate_draws.csv",
        "aggregate_summary": "revision_t2_full_aggregate_summary.csv",
        "real_predictions": "revision_t2_full_real_oof_predictions.csv",
    }
    for key, frame in outputs.items():
        frame.to_csv(output_dir / filenames[key], index=False)
    (output_dir / "revision_t2_full_config.json").write_text(
        json.dumps(
            {
                "analysis": "multi_partition_permuted_structure_control",
                "target": target,
                "feature_set": feature_set,
                "model": model_name,
                "iterations_per_partition": iterations,
                "group_counts": list(group_counts),
                "random_partition_count_per_group_count": random_partition_count,
                "ward_anchor_included": True,
                "start_seed": start_seed,
                "model_seed": seed,
                "resume": resume,
                "checkpoint_every": checkpoint_every,
            },
            indent=2,
            sort_keys=True,
        ),
        encoding="utf-8",
    )
    return outputs


def _inner_oof_mae(
    family: FeatureFamily,
    outer_train_idx: np.ndarray,
    groups: np.ndarray,
    *,
    target: str,
    model_name: str,
    seed: int,
    n_jobs: int,
) -> float:
    frame = family.frame
    X = frame.loc[:, list(family.feature_columns)]
    y = pd.to_numeric(frame[target], errors="raise").to_numpy(dtype=float)
    train_groups = groups[outer_train_idx]
    unique = np.unique(train_groups)
    if len(unique) < 2:
        raise ValueError("Nested model selection needs at least two spatial groups in outer training")
    cv = GroupKFold(n_splits=len(unique))
    predictions = np.full(len(outer_train_idx), np.nan, dtype=float)
    for inner_fold, (inner_train_local, inner_test_local) in enumerate(
        cv.split(X.iloc[outer_train_idx], y[outer_train_idx], train_groups), start=1
    ):
        inner_train_global = outer_train_idx[inner_train_local]
        inner_test_global = outer_train_idx[inner_test_local]
        estimator = make_frozen_estimator(model_name, seed=seed, n_jobs=n_jobs)
        weights = build_sample_weights(
            frame.iloc[inner_train_global], DEFAULT_MODEL_WEIGHT_MODES.get(model_name, "none")
        )
        fitted = _fit_pipeline(
            estimator,
            X.iloc[inner_train_global],
            y[inner_train_global],
            weights,
            model_name,
        )
        predictions[inner_test_local] = np.clip(
            fitted.predict(X.iloc[inner_test_global]), 0.0, 1.0
        )
    if not np.isfinite(predictions).all():
        raise RuntimeError("Incomplete inner OOF predictions")
    return float(np.mean(np.abs(y[outer_train_idx] - predictions)))


def run_nested_model_selection(
    labels: pd.DataFrame,
    sequence: pd.DataFrame,
    structure: pd.DataFrame,
    output_dir: Path,
    *,
    target: str = "repair_published",
    feature_set: str = "structure_transfer_safe",
    model_names: Sequence[str] = ("extra_trees", "lightgbm", "ridge"),
    fixed_model: str = "extra_trees",
    group_counts: Sequence[int] = (4,),
    random_partition_count: int = 5,
    start_seed: int = 20260729,
    n_jobs: int = 1,
    bootstrap_iterations: int = 10_000,
    seed: int = DEFAULT_SEED,
) -> dict[str, pd.DataFrame]:
    """Nested spatial model selection versus naive same-data selection and fixed Extra Trees."""

    output_dir.mkdir(parents=True, exist_ok=True)
    families = build_feature_families(labels, sequence, structure)
    if feature_set not in families:
        raise ValueError(f"Unknown feature_set {feature_set!r}")
    family = families[feature_set]
    frame = family.frame
    X = frame.loc[:, list(family.feature_columns)]
    y = pd.to_numeric(frame[target], errors="raise").to_numpy(dtype=float)
    partitions = build_spatial_partitions(
        frame,
        group_counts=group_counts,
        random_partition_count=random_partition_count,
        start_seed=start_seed,
        include_ward_anchor=True,
    )

    fold_rows: list[dict[str, Any]] = []
    prediction_rows: list[pd.DataFrame] = []
    partition_rows: list[dict[str, Any]] = []
    for partition in partitions:
        groups = np.asarray(partition.groups)
        cv = GroupKFold(n_splits=np.unique(groups).size)
        nested_pred = np.full(len(frame), np.nan, dtype=float)
        fixed_pred = np.full(len(frame), np.nan, dtype=float)
        nested_model_by_site = np.empty(len(frame), dtype=object)
        outer_fold_ids = np.full(len(frame), -1, dtype=int)

        for outer_fold, (train_idx, test_idx) in enumerate(cv.split(X, y, groups), start=1):
            inner_scores = {
                model_name: _inner_oof_mae(
                    family,
                    np.asarray(train_idx),
                    groups,
                    target=target,
                    model_name=model_name,
                    seed=seed,
                    n_jobs=n_jobs,
                )
                for model_name in model_names
            }
            selected_model = min(model_names, key=lambda name: (inner_scores[name], model_names.index(name)))
            estimator = make_frozen_estimator(
                selected_model, seed=seed, n_jobs=n_jobs
            )
            weights = build_sample_weights(
                frame.iloc[train_idx], DEFAULT_MODEL_WEIGHT_MODES.get(selected_model, "none")
            )
            fitted = _fit_pipeline(
                estimator, X.iloc[train_idx], y[train_idx], weights, selected_model
            )
            nested_pred[test_idx] = np.clip(fitted.predict(X.iloc[test_idx]), 0.0, 1.0)
            nested_model_by_site[test_idx] = selected_model

            fixed_estimator = make_frozen_estimator(
                fixed_model, seed=seed, n_jobs=n_jobs
            )
            fixed_weights = build_sample_weights(
                frame.iloc[train_idx], DEFAULT_MODEL_WEIGHT_MODES.get(fixed_model, "none")
            )
            fixed_fitted = _fit_pipeline(
                fixed_estimator,
                X.iloc[train_idx],
                y[train_idx],
                fixed_weights,
                fixed_model,
            )
            fixed_pred[test_idx] = np.clip(
                fixed_fitted.predict(X.iloc[test_idx]), 0.0, 1.0
            )
            outer_fold_ids[test_idx] = outer_fold
            fold_rows.append(
                {
                    "partition_id": partition.partition_id,
                    "partition_method": partition.method,
                    "partition_seed": partition.seed,
                    "n_groups": partition.n_groups,
                    "outer_fold": outer_fold,
                    "outer_test_group": "|".join(map(str, np.unique(groups[test_idx]))),
                    "n_train": len(train_idx),
                    "n_test": len(test_idx),
                    "selected_model": selected_model,
                    **{f"inner_mae_{name}": score for name, score in inner_scores.items()},
                }
            )

        nested_metrics = regression_metrics(y, nested_pred)
        fixed_metrics = regression_metrics(y, fixed_pred)
        all_model_metrics, _ = evaluate_feature_family_partition(
            family,
            partition,
            target=target,
            model_names=model_names,
            seed=seed,
            n_jobs=n_jobs,
        )
        naive_row = all_model_metrics.sort_values(["mae", "model"]).iloc[0]
        partition_rows.append(
            {
                "partition_id": partition.partition_id,
                "partition_method": partition.method,
                "partition_seed": partition.seed,
                "partition_fingerprint": partition.fingerprint,
                "n_groups": partition.n_groups,
                "feature_set": feature_set,
                "nested_mae": nested_metrics["mae"],
                "nested_spearman": nested_metrics["spearman_rho"],
                "nested_r2": nested_metrics["r2"],
                "fixed_model": fixed_model,
                "fixed_mae": fixed_metrics["mae"],
                "fixed_spearman": fixed_metrics["spearman_rho"],
                "fixed_r2": fixed_metrics["r2"],
                "naive_best_model": str(naive_row["model"]),
                "naive_best_mae": float(naive_row["mae"]),
                "selection_optimism_penalty_nested_minus_naive": float(
                    nested_metrics["mae"] - float(naive_row["mae"])
                ),
                "delta_mae_fixed_minus_nested": float(fixed_metrics["mae"] - nested_metrics["mae"]),
            }
        )
        prediction_rows.append(
            pd.DataFrame(
                {
                    "protein_id": frame["protein_id"].astype(str),
                    "site_uniprot": frame["site_uniprot"],
                    "partition_id": partition.partition_id,
                    "outer_group": groups.astype(str),
                    "outer_fold": outer_fold_ids,
                    "y_true": y,
                    "nested_y_pred": nested_pred,
                    "nested_selected_model": nested_model_by_site,
                    "fixed_model": fixed_model,
                    "fixed_y_pred": fixed_pred,
                }
            )
        )

    fold_selection = pd.DataFrame(fold_rows)
    predictions = pd.concat(prediction_rows, ignore_index=True)
    partition_metrics = pd.DataFrame(partition_rows)
    selection_frequency = (
        fold_selection.groupby("selected_model", as_index=False)
        .size()
        .rename(columns={"size": "selected_outer_folds"})
    )
    selection_frequency["selection_fraction"] = (
        selection_frequency["selected_outer_folds"] / len(fold_selection)
    )

    summary_rows: list[dict[str, Any]] = []
    for metric in (
        "nested_mae",
        "fixed_mae",
        "naive_best_mae",
        "selection_optimism_penalty_nested_minus_naive",
        "delta_mae_fixed_minus_nested",
    ):
        values = partition_metrics[metric].to_numpy(dtype=float)
        low, high = _bootstrap_mean_ci(values, iterations=bootstrap_iterations, seed=seed)
        summary_rows.append(
            {
                "metric": metric,
                "n_partitions": len(values),
                "mean": float(np.mean(values)),
                "median": float(np.median(values)),
                "ci_2.5": low,
                "ci_97.5": high,
                "positive_fraction": float(np.mean(values > 0.0)),
            }
        )
    summary = pd.DataFrame(summary_rows)

    outputs = {
        "fold_selection": fold_selection,
        "predictions": predictions,
        "partition_metrics": partition_metrics,
        "selection_frequency": selection_frequency,
        "summary": summary,
    }
    filenames = {
        "fold_selection": "revision_t3_nested_fold_selection.csv",
        "predictions": "revision_t3_nested_oof_predictions.csv",
        "partition_metrics": "revision_t3_nested_partition_metrics.csv",
        "selection_frequency": "revision_t3_nested_selection_frequency.csv",
        "summary": "revision_t3_nested_summary.csv",
    }
    for key, frame_out in outputs.items():
        frame_out.to_csv(output_dir / filenames[key], index=False)
    (output_dir / "revision_t3_config.json").write_text(
        json.dumps(
            {
                "analysis": "nested_spatial_model_selection",
                "target": target,
                "feature_set": feature_set,
                "candidate_models": list(model_names),
                "fixed_model": fixed_model,
                "group_counts": list(group_counts),
                "random_partition_count_per_group_count": random_partition_count,
                "ward_anchor_included": True,
                "start_seed": start_seed,
                "bootstrap_iterations": bootstrap_iterations,
            },
            indent=2,
            sort_keys=True,
        ),
        encoding="utf-8",
    )
    return outputs


def _metric_bundle(y: np.ndarray, primary: np.ndarray, null: np.ndarray) -> dict[str, float]:
    primary_metrics = regression_metrics(y, primary)
    null_mae = float(np.mean(np.abs(y - null)))
    delta = float(null_mae - primary_metrics["mae"])
    return {
        "primary_mae": primary_metrics["mae"],
        "primary_spearman": primary_metrics["spearman_rho"],
        "primary_r2": primary_metrics["r2"],
        "null_mae": null_mae,
        "delta_mae_null_minus_primary": delta,
        "relative_mae_improvement": float(delta / null_mae) if null_mae > 0 else float("nan"),
    }


def _sign_flip_pvalue(values: np.ndarray, *, iterations: int, seed: int) -> float:
    values = np.asarray(values, dtype=float)
    observed = float(np.mean(values))
    if not np.isfinite(observed):
        return float("nan")
    rng = np.random.default_rng(seed)
    exceed = 0
    done = 0
    chunk = 10_000
    while done < iterations:
        n = min(chunk, iterations - done)
        signs = rng.choice(np.array([-1.0, 1.0]), size=(n, len(values)))
        draws = (signs * values).mean(axis=1)
        exceed += int(np.sum(draws >= observed))
        done += n
    return float((1 + exceed) / (iterations + 1))


def run_primary_cluster_bootstrap(
    predictions_path: Path,
    output_dir: Path,
    *,
    partition_id: str = "ward_g4",
    primary_feature_set: str = "structure_transfer_safe",
    primary_model: str = "extra_trees",
    iterations: int = 20_000,
    sign_flip_iterations: int = 100_000,
    seed: int = DEFAULT_SEED,
) -> dict[str, pd.DataFrame]:
    """Paired uncertainty analysis for the frozen primary OOF model.

    Spatial-group cluster bootstrap is primary. Ordinary site bootstrap is retained only as a
    sensitivity analysis because sites within a structural block are not independent.
    """

    if not predictions_path.exists():
        raise FileNotFoundError(predictions_path)
    output_dir.mkdir(parents=True, exist_ok=True)
    pred = pd.read_csv(predictions_path)
    primary = pred.loc[
        (pred["partition_id"] == partition_id)
        & (pred["feature_set"] == primary_feature_set)
        & (pred["model"] == primary_model)
    ].copy()
    null = pred.loc[
        (pred["partition_id"] == partition_id)
        & (pred["feature_set"] == "null_mean")
        & (pred["model"] == "dummy_mean")
    ].copy()
    if primary.empty or null.empty:
        raise ValueError("Requested primary/null OOF predictions were not found")
    keys = ["protein_id", "site_uniprot", "partition_id"]
    merged = primary.merge(
        null[keys + ["y_pred"]].rename(columns={"y_pred": "null_y_pred"}),
        on=keys,
        how="inner",
        validate="one_to_one",
    )
    if len(merged) != len(primary):
        raise ValueError("Primary and null OOF rows do not align one-to-one")
    y = merged["y_true"].to_numpy(dtype=float)
    p = merged["y_pred"].to_numpy(dtype=float)
    n = merged["null_y_pred"].to_numpy(dtype=float)
    observed = _metric_bundle(y, p, n)

    draw_rows: list[dict[str, Any]] = []
    rng = np.random.default_rng(seed)
    group_values = merged["outer_group"].astype(str).to_numpy()
    unique_groups = np.unique(group_values)
    for method in ("spatial_group_cluster", "site"):
        for iteration in range(iterations):
            if method == "spatial_group_cluster":
                sampled_groups = rng.choice(unique_groups, size=len(unique_groups), replace=True)
                pieces = [np.flatnonzero(group_values == group) for group in sampled_groups]
                idx = np.concatenate(pieces)
            else:
                idx = rng.integers(0, len(merged), size=len(merged))
            metrics = _metric_bundle(y[idx], p[idx], n[idx])
            draw_rows.append({"method": method, "iteration": iteration, **metrics})
    draws = pd.DataFrame(draw_rows)
    summary_rows: list[dict[str, Any]] = []
    for method, method_group in draws.groupby("method"):
        for metric, observed_value in observed.items():
            values = method_group[metric].to_numpy(dtype=float)
            values = values[np.isfinite(values)]
            low, high = (float("nan"), float("nan"))
            if len(values):
                low, high = map(float, np.quantile(values, [0.025, 0.975]))
            summary_rows.append(
                {
                    "method": method,
                    "metric": metric,
                    "observed": observed_value,
                    "bootstrap_mean": float(np.mean(values)) if len(values) else float("nan"),
                    "ci_2.5": low,
                    "ci_97.5": high,
                    "valid_draws": len(values),
                }
            )
    summary = pd.DataFrame(summary_rows)
    error_advantage = np.abs(y - n) - np.abs(y - p)
    tests = pd.DataFrame(
        [
            {
                "test": "paired_sign_flip_mean_absolute_error_advantage",
                "observed_mean_advantage": float(error_advantage.mean()),
                "iterations": sign_flip_iterations,
                "p_one_sided": _sign_flip_pvalue(
                    error_advantage, iterations=sign_flip_iterations, seed=seed + 17
                ),
                "positive_site_fraction": float(np.mean(error_advantage > 0.0)),
            }
        ]
    )
    observed_frame = pd.DataFrame(
        [
            {
                "partition_id": partition_id,
                "primary_feature_set": primary_feature_set,
                "primary_model": primary_model,
                "n_sites": len(merged),
                "n_spatial_groups": len(unique_groups),
                **observed,
            }
        ]
    )

    outputs = {
        "observed": observed_frame,
        "draws": draws,
        "summary": summary,
        "tests": tests,
    }
    filenames = {
        "observed": "revision_t4_primary_observed_metrics.csv",
        "draws": "revision_t4_bootstrap_draws.csv",
        "summary": "revision_t4_bootstrap_summary.csv",
        "tests": "revision_t4_paired_tests.csv",
    }
    for key, frame_out in outputs.items():
        frame_out.to_csv(output_dir / filenames[key], index=False)
    (output_dir / "revision_t4_config.json").write_text(
        json.dumps(
            {
                "analysis": "paired_cluster_bootstrap_primary_model",
                "predictions_path": str(predictions_path),
                "partition_id": partition_id,
                "primary_feature_set": primary_feature_set,
                "primary_model": primary_model,
                "bootstrap_iterations_per_method": iterations,
                "sign_flip_iterations": sign_flip_iterations,
                "primary_bootstrap_unit": "spatial group",
                "sensitivity_bootstrap_unit": "site",
            },
            indent=2,
            sort_keys=True,
        ),
        encoding="utf-8",
    )
    return outputs


def write_confirmatory_report(
    output_dir: Path,
    *,
    t2_aggregate: pd.DataFrame | None = None,
    t3_summary: pd.DataFrame | None = None,
    t3_frequency: pd.DataFrame | None = None,
    t4_observed: pd.DataFrame | None = None,
    t4_summary: pd.DataFrame | None = None,
    t4_tests: pd.DataFrame | None = None,
) -> Path:
    lines = [
        "# CML-RepairNet Confirmatory Reviewer Tests",
        "",
        "All analyses remain BSA/CrGO-897-assay specific and fully in silico.",
        "Structure-aware evidence does not establish cross-protein calibration.",
        "",
    ]
    sections = [
        ("T2 — Multi-partition permuted-structure control", t2_aggregate),
        ("T3 — Nested spatial model-selection summary", t3_summary),
        ("T3 — Nested selection frequencies", t3_frequency),
        ("T4 — Frozen primary observed metrics", t4_observed),
        ("T4 — Paired bootstrap intervals", t4_summary),
        ("T4 — Paired tests", t4_tests),
    ]
    for title, frame in sections:
        if frame is not None and not frame.empty:
            lines.extend([f"## {title}", "", "```text", frame.to_string(index=False), "```", ""])
    path = output_dir / "revision_confirmatory_report.md"
    path.write_text("\n".join(lines), encoding="utf-8")
    return path
