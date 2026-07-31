from __future__ import annotations

import argparse
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

from cml_repairnet.evaluation.robustness import build_feature_view, paired_error_comparison
from cml_repairnet.models.stage25 import (
    aggregate_permutation_importance,
    run_stage25_nested_cv,
)
from cml_repairnet.utils.config import load_config
from cml_repairnet.utils.reproducibility import set_global_seed


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Run CML-RepairNet Stage 2.5: regularized XGBoost/LightGBM, "
            "weight sensitivity, nested correlation pruning and convex ensembling."
        )
    )
    parser.add_argument("--mode", choices=["quick", "full"], default="quick")
    parser.add_argument("--config", default="configs/base.yaml")
    parser.add_argument("--n-jobs", type=int, default=1)
    parser.add_argument("--bootstrap", type=int, default=5000)
    parser.add_argument("--permutations", type=int, default=20_000)
    parser.add_argument("--permutation-repeats", type=int, default=50)
    parser.add_argument(
        "--models",
        nargs="+",
        default=[
            "dummy_mean",
            "ridge",
            "extra_trees",
            "xgboost",
            "lightgbm",
        ],
    )
    return parser.parse_args()


def _analysis_plan(mode: str) -> list[tuple[str, int, str, str]]:
    if mode == "quick":
        return [("repair_published", 4, "none", "full")]

    plan: list[tuple[str, int, str, str]] = []
    # Primary robustness matrix: both labels and all spatial group counts, unweighted/full.
    for target in ("repair_published", "repair_from_means_clipped"):
        for groups in (3, 4, 5):
            plan.append((target, groups, "none", "full"))

    # Sensitivity analyses are intentionally anchored to the prespecified primary condition.
    for weight_mode in ("clipped", "sqrt"):
        plan.append(("repair_published", 4, weight_mode, "full"))
    plan.append(("repair_published", 4, "none", "nested_pruned"))
    return plan


def add_bh_fdr(frame: pd.DataFrame, p_column: str = "permutation_p_one_sided") -> pd.DataFrame:
    result = frame.copy()
    if result.empty:
        result["fdr_q_global"] = pd.Series(dtype=float)
        return result
    p = result[p_column].to_numpy(dtype=float)
    order = np.argsort(p)
    ranked = p[order]
    adjusted = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    q = np.empty_like(adjusted)
    q[order] = np.clip(adjusted, 0.0, 1.0)
    result["fdr_q_global"] = q
    return result


def write_report(
    metrics: pd.DataFrame,
    comparisons: pd.DataFrame,
    reference_comparisons: pd.DataFrame,
    importance: pd.DataFrame,
    output_path: Path,
) -> None:
    best = (
        metrics.sort_values("mae")
        .groupby(["target", "spatial_groups", "weight_mode", "feature_mode"], as_index=False)
        .first()
    )
    primary = metrics[
        (metrics["target"] == "repair_published")
        & (metrics["spatial_groups"] == 4)
        & (metrics["weight_mode"] == "none")
        & (metrics["feature_mode"] == "full")
    ].sort_values("mae")

    top_importance = (
        importance.sort_values(["model", "importance_mean"], ascending=[True, False])
        .groupby("model", as_index=False)
        .head(8)
        if not importance.empty
        else importance.copy()
    )
    primary_reference = (
        reference_comparisons[
            (reference_comparisons["target"] == "repair_published")
            & (reference_comparisons["spatial_groups"] == 4)
            & (reference_comparisons["weight_mode"] == "none")
            & (reference_comparisons["feature_mode"] == "full")
        ]
        if not reference_comparisons.empty
        else reference_comparisons.copy()
    )
    reference_text = (
        primary_reference[
            [
                "model", "baseline", "model_mae", "baseline_mae",
                "relative_improvement_pct", "improvement_ci_2.5",
                "improvement_ci_97.5", "permutation_p_one_sided", "fdr_q_global",
            ]
        ].to_string(index=False)
        if not primary_reference.empty
        else "No Ridge/Extra Trees reference comparisons were requested."
    )
    importance_text = (
        top_importance[
            ["model", "feature", "importance_mean", "importance_sd_across_folds", "folds"]
        ].to_string(index=False)
        if not top_importance.empty
        else "Permutation importance was disabled or unavailable."
    )
    lines = [
        "# CML-RepairNet Stage 2.5 Report",
        "",
        "## Prespecified primary condition",
        "",
        "Target: `repair_published`; spatial groups: 4; unweighted; full structural panel.",
        "",
        "```text",
        primary[
            [
                "model", "mae", "rmse", "spearman_rho", "r2",
                "mae_ci_2.5", "mae_ci_97.5", "mean_selected_features",
            ]
        ].to_string(index=False),
        "```",
        "",
        "## Best model per Stage 2.5 analysis",
        "",
        "```text",
        best[
            [
                "target", "spatial_groups", "weight_mode", "feature_mode",
                "model", "mae", "spearman_rho", "r2",
            ]
        ].to_string(index=False),
        "```",
        "",
        "## Paired OOF comparisons against dummy_mean",
        "",
        "```text",
        comparisons[
            [
                "target", "spatial_groups", "weight_mode", "feature_mode", "model",
                "relative_improvement_pct", "improvement_ci_2.5",
                "improvement_ci_97.5", "permutation_p_one_sided",
            ]
        ].to_string(index=False),
        "```",
        "",
        "## Primary comparisons against Ridge and Extra Trees",
        "",
        "```text",
        reference_text,
        "```",
        "",
        "## Held-out permutation importance (top eight per model)",
        "",
        "```text",
        importance_text,
        "```",
        "",
        "## Safeguards",
        "",
        "- All reported performance is outer-fold out-of-fold performance.",
        "- XGBoost and LightGBM use small, prespecified, strongly regularized candidate sets.",
        "- Correlation pruning is fitted inside training folds only.",
        "- Sample weights are clipped and normalized within each outer training fold.",
        "- Convex ensemble weights are selected without using outer-test labels.",
        "- Permutation importance is calculated on held-out outer-fold samples.",
        "",
    ]
    output_path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    warnings.filterwarnings(
        "ignore",
        message="X does not have valid feature names, but LGBMRegressor was fitted with feature names",
        category=UserWarning,
    )
    args = parse_args()
    config = load_config(args.config)
    seed = int(config["project"]["seed"])
    set_global_seed(seed)
    paths = config["paths"]

    labels = pd.read_csv(paths["labels"])
    sequence = pd.read_csv(paths["sequence_features"])
    structure = pd.read_csv(paths["structure_features"])
    merged = build_feature_view(labels, sequence, structure, "structure")

    result_root = Path(paths["results_dir"]) / "stage25"
    all_metrics: list[pd.DataFrame] = []
    all_predictions: list[pd.DataFrame] = []
    all_importances: list[pd.DataFrame] = []
    comparison_rows: list[dict[str, object]] = []
    reference_comparison_rows: list[dict[str, object]] = []

    for target, groups, weight_mode, feature_mode in _analysis_plan(args.mode):
        run_dir = result_root / target / f"g{groups}" / weight_mode / feature_mode
        print(
            f"Running target={target}, groups={groups}, weight={weight_mode}, "
            f"features={feature_mode}"
        )
        metrics, predictions, importance = run_stage25_nested_cv(
            merged,
            target,
            run_dir,
            n_spatial_groups=groups,
            weight_mode=weight_mode,
            feature_mode=feature_mode,
            bootstrap_iterations=args.bootstrap,
            permutation_repeats=args.permutation_repeats,
            n_jobs=args.n_jobs,
            seed=seed,
            model_names=args.models,
            include_convex_ensemble=True,
        )
        all_metrics.append(metrics)
        all_predictions.append(predictions)
        if not importance.empty:
            importance = importance.copy()
            all_importances.append(importance)

        for model in sorted(set(predictions["model"]) - {"dummy_mean"}):
            row = paired_error_comparison(
                predictions,
                model=model,
                bootstrap_iterations=args.bootstrap,
                permutation_iterations=args.permutations,
                seed=seed + groups,
            ).as_dict()
            row.update(
                target=target,
                spatial_groups=groups,
                weight_mode=weight_mode,
                feature_mode=feature_mode,
            )
            comparison_rows.append(row)

        available_models = set(predictions["model"])
        for baseline in ("ridge", "extra_trees"):
            if baseline not in available_models:
                continue
            for model in ("xgboost", "lightgbm", "ridge_extra_trees_convex"):
                if model not in available_models or model == baseline:
                    continue
                row = paired_error_comparison(
                    predictions,
                    model=model,
                    baseline=baseline,
                    bootstrap_iterations=args.bootstrap,
                    permutation_iterations=args.permutations,
                    seed=seed + groups + 101,
                ).as_dict()
                row.update(
                    target=target,
                    spatial_groups=groups,
                    weight_mode=weight_mode,
                    feature_mode=feature_mode,
                )
                reference_comparison_rows.append(row)

    metrics_frame = pd.concat(all_metrics, ignore_index=True).sort_values(
        ["target", "spatial_groups", "weight_mode", "feature_mode", "mae"]
    )
    predictions_frame = pd.concat(all_predictions, ignore_index=True)
    comparisons_frame = pd.DataFrame(comparison_rows)
    if not comparisons_frame.empty:
        comparisons_frame = comparisons_frame.sort_values(
            ["target", "spatial_groups", "weight_mode", "feature_mode", "model"]
        )
    comparisons_frame = add_bh_fdr(comparisons_frame)

    reference_comparisons_frame = pd.DataFrame(reference_comparison_rows)
    if not reference_comparisons_frame.empty:
        reference_comparisons_frame = reference_comparisons_frame.sort_values(
            ["target", "spatial_groups", "weight_mode", "feature_mode", "baseline", "model"]
        )
    reference_comparisons_frame = add_bh_fdr(reference_comparisons_frame)
    raw_importance = (
        pd.concat(all_importances, ignore_index=True) if all_importances else pd.DataFrame()
    )
    if raw_importance.empty:
        importance_frame = pd.DataFrame(
            columns=[
                "model", "feature", "importance_mean",
                "importance_sd_across_folds", "mean_within_fold_sd", "folds",
            ]
        )
    else:
        primary_importance = raw_importance[
            (raw_importance["target"] == "repair_published")
            & (raw_importance["spatial_groups"] == 4)
            & (raw_importance["weight_mode"] == "none")
            & (raw_importance["feature_mode"] == "full")
        ]
        importance_frame = aggregate_permutation_importance(primary_importance)

    result_root.mkdir(parents=True, exist_ok=True)
    metrics_path = result_root / "stage25_all_metrics.csv"
    predictions_path = result_root / "stage25_all_oof_predictions.csv"
    comparisons_path = result_root / "stage25_paired_comparisons.csv"
    reference_comparisons_path = result_root / "stage25_reference_comparisons.csv"
    importance_path = result_root / "stage25_primary_permutation_importance.csv"
    report_path = result_root / "stage25_report.md"
    metrics_frame.to_csv(metrics_path, index=False)
    predictions_frame.to_csv(predictions_path, index=False)
    comparisons_frame.to_csv(comparisons_path, index=False)
    reference_comparisons_frame.to_csv(reference_comparisons_path, index=False)
    importance_frame.to_csv(importance_path, index=False)
    write_report(
        metrics_frame, comparisons_frame, reference_comparisons_frame, importance_frame, report_path
    )

    primary = metrics_frame[
        (metrics_frame["target"] == "repair_published")
        & (metrics_frame["spatial_groups"] == 4)
        & (metrics_frame["weight_mode"] == "none")
        & (metrics_frame["feature_mode"] == "full")
    ].sort_values("mae")
    print("\nPrimary Stage 2.5 ranking:")
    print(
        primary[["model", "mae", "spearman_rho", "r2", "mae_ci_2.5", "mae_ci_97.5"]]
        .to_string(index=False)
    )
    print(f"\nMetrics: {metrics_path}")
    print(f"Predictions: {predictions_path}")
    print(f"Comparisons: {comparisons_path}")
    print(f"Reference comparisons: {reference_comparisons_path}")
    print(f"Importance: {importance_path}")
    print(f"Report: {report_path}")


if __name__ == "__main__":
    main()
