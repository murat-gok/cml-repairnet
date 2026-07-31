from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from cml_repairnet.evaluation.robustness import (
    build_feature_view,
    paired_error_comparison,
    write_markdown_summary,
)
from cml_repairnet.models.baselines import run_nested_group_cv
from cml_repairnet.utils.config import load_config
from cml_repairnet.utils.reproducibility import set_global_seed


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run Stage 2 robustness analyses.")
    parser.add_argument("--mode", choices=["quick", "full"], default="quick")
    parser.add_argument("--config", default="configs/base.yaml")
    parser.add_argument("--bootstrap", type=int, default=10_000)
    parser.add_argument("--permutations", type=int, default=20_000)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    config = load_config(args.config)
    seed = int(config["project"]["seed"])
    set_global_seed(seed)
    paths = config["paths"]

    labels = pd.read_csv(paths["labels"])
    sequence = pd.read_csv(paths["sequence_features"])
    structure = pd.read_csv(paths["structure_features"])

    feature_sets = ["sequence", "structure", "combined"]
    targets = ["repair_published"] if args.mode == "quick" else [
        "repair_published", "repair_from_means_clipped"
    ]
    group_counts = [4] if args.mode == "quick" else [3, 4, 5]

    result_root = Path(paths["results_dir"]) / "stage2_robustness"
    metric_frames: list[pd.DataFrame] = []
    comparison_rows: list[dict[str, object]] = []

    for feature_set in feature_sets:
        merged = build_feature_view(labels, sequence, structure, feature_set)
        for target in targets:
            for n_groups in group_counts:
                run_dir = result_root / feature_set / target / f"g{n_groups}"
                print(f"Running feature_set={feature_set}, target={target}, groups={n_groups}")
                metrics, predictions = run_nested_group_cv(
                    merged,
                    target=target,
                    output_dir=run_dir,
                    n_spatial_groups=n_groups,
                    bootstrap_iterations=int(config["validation"]["bootstrap_iterations"]),
                    seed=seed,
                )
                metrics = metrics.copy()
                metrics.insert(0, "feature_set", feature_set)
                metrics.insert(1, "target", target)
                metrics.insert(2, "spatial_groups", n_groups)
                metric_frames.append(metrics)

                for model in sorted(set(predictions["model"]) - {"dummy_mean"}):
                    row = paired_error_comparison(
                        predictions,
                        model=model,
                        bootstrap_iterations=args.bootstrap,
                        permutation_iterations=args.permutations,
                        seed=seed + n_groups,
                    ).as_dict()
                    row.update(
                        feature_set=feature_set,
                        target=target,
                        spatial_groups=n_groups,
                    )
                    comparison_rows.append(row)

    all_metrics = pd.concat(metric_frames, ignore_index=True).sort_values(
        ["feature_set", "target", "spatial_groups", "mae"]
    )
    all_comparisons = pd.DataFrame(comparison_rows).sort_values(
        ["feature_set", "target", "spatial_groups", "model"]
    )
    result_root.mkdir(parents=True, exist_ok=True)
    metrics_path = result_root / "stage2_robustness_metrics.csv"
    comparisons_path = result_root / "stage2_paired_comparisons.csv"
    report_path = result_root / "stage2_robustness_report.md"
    all_metrics.to_csv(metrics_path, index=False)
    all_comparisons.to_csv(comparisons_path, index=False)
    write_markdown_summary(all_metrics, all_comparisons, report_path)

    best = (
        all_metrics.sort_values("mae")
        .groupby(["feature_set", "target", "spatial_groups"], as_index=False)
        .first()
    )
    print("\nBest model per analysis:")
    print(best[["feature_set", "target", "spatial_groups", "model", "mae", "spearman_rho", "r2", "n_features"]].to_string(index=False))
    print(f"\nMetrics: {metrics_path}")
    print(f"Paired comparisons: {comparisons_path}")
    print(f"Report: {report_path}")


if __name__ == "__main__":
    main()
