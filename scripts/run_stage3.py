from __future__ import annotations

import argparse
import json
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

from cml_repairnet.evaluation.robustness import build_feature_view
from cml_repairnet.models.stage3 import (
    run_cysteine_mechanism_analysis,
    run_leave_one_domain_out,
    run_repeated_spatial_validation,
    summarize_uncertainty,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run CML-RepairNet Stage 3 validation")
    parser.add_argument(
        "--mode",
        choices=("quick", "full", "exhaustive"),
        default="quick",
        help=(
            "quick: 5 unique random-init g4 partitions; full: 50 unique g4 partitions; "
            "exhaustive: 50 unique partitions for g3/g4/g5"
        ),
    )
    parser.add_argument(
        "--target",
        choices=("repair_published", "repair_from_means_clipped"),
        default="repair_published",
    )
    parser.add_argument("--seeds", type=int, default=None)
    parser.add_argument("--groups", type=int, nargs="+", default=None)
    parser.add_argument("--n-jobs", type=int, default=1)
    parser.add_argument("--conditional-repeats", type=int, default=500)
    parser.add_argument("--conformal-alpha", type=float, default=0.10)
    parser.add_argument("--skip-repeated", action="store_true")
    parser.add_argument("--skip-domain", action="store_true")
    parser.add_argument("--skip-cysteine", action="store_true")
    parser.add_argument("--output", type=Path, default=Path("results/stage3"))
    return parser.parse_args()


def _mode_plan(args: argparse.Namespace) -> tuple[list[int], list[int]]:
    if args.seeds is not None:
        n_seeds = args.seeds
    elif args.mode == "quick":
        n_seeds = 5
    else:
        n_seeds = 50
    if n_seeds < 1:
        raise ValueError("--seeds must be at least 1")

    if args.groups is not None:
        groups = args.groups
    elif args.mode == "exhaustive":
        groups = [3, 4, 5]
    else:
        groups = [4]
    return list(range(20260729, 20260729 + n_seeds)), groups


def _load_structure_view() -> pd.DataFrame:
    labels = pd.read_csv("data/processed/bsa_labels.csv")
    sequence = pd.read_csv("data/processed/bsa_sequence_features.csv")
    structure = pd.read_csv("data/processed/bsa_structure_features.csv")
    return build_feature_view(labels, sequence, structure, "structure")


def _write_report(
    output: Path,
    repeated_summary: pd.DataFrame,
    domain_metrics: pd.DataFrame,
    uncertainty_summary: pd.DataFrame,
    cysteine_summary: pd.DataFrame,
    config: dict[str, object],
) -> None:
    lines = [
        "# CML-RepairNet Stage 3 Validation Report",
        "",
        "## Configuration",
        "",
        "```json",
        json.dumps(config, indent=2, sort_keys=True),
        "```",
        "",
    ]
    if not repeated_summary.empty:
        lines += [
            "## Repeated spatial validation",
            "",
            "```text",
            repeated_summary.to_string(index=False),
            "```",
            "",
        ]
    if not domain_metrics.empty:
        display_columns = [
            "model", "mae", "rmse", "spearman_rho", "r2",
            "low_boundary_mae", "high_boundary_mae",
        ]
        display_columns = [column for column in display_columns if column in domain_metrics]
        lines += [
            "## Leave-one-domain-out validation",
            "",
            "```text",
            domain_metrics[display_columns].to_string(index=False),
            "```",
            "",
        ]
    if not uncertainty_summary.empty:
        lines += [
            "## Conformal and OOD summary",
            "",
            "```text",
            uncertainty_summary.to_string(index=False),
            "```",
            "",
        ]
    if not cysteine_summary.empty:
        lines += [
            "## Cysteine-neighbourhood mechanism checks",
            "",
            "```text",
            cysteine_summary.to_string(index=False),
            "```",
            "",
        ]
    lines += [
        "## Prespecified interpretation",
        "",
        "- Extra Trees is the primary continuous predictor.",
        "- LightGBM is the uncertainty-aware ranking model.",
        "- Ridge is the linear interpretation model.",
        "- Repeated spatial partitions, domain holdout and conditional feature tests are",
        "  treated as validation, not as opportunities for selecting a new winner.",
        "- Human-proteome generalization must not be claimed from BSA-only labels.",
        "",
    ]
    (output / "stage3_report.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    warnings.filterwarnings(
        "ignore",
        message="X does not have valid feature names, but LGBMRegressor was fitted with feature names",
    )
    args = parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    merged = _load_structure_view()
    seeds, group_counts = _mode_plan(args)

    repeated_summary = pd.DataFrame()
    domain_metrics = pd.DataFrame()
    uncertainty_summary = pd.DataFrame()
    cysteine_summary = pd.DataFrame()

    if not args.skip_repeated:
        print(
            f"Repeated spatial validation: requested_unique_partitions={len(seeds)}, "
            f"groups={group_counts}, target={args.target}"
        )
        _, _, repeated_summary = run_repeated_spatial_validation(
            merged,
            args.target,
            args.output / "repeated_spatial",
            seeds=seeds,
            group_counts=group_counts,
            n_jobs=args.n_jobs,
            feature_mode="nested_pruned",
            include_anchor=True,
        )
        print(repeated_summary.to_string(index=False))

    if not args.skip_domain:
        print("Running leave-one-domain-out validation")
        domain_metrics, domain_predictions, uncertainty = run_leave_one_domain_out(
            merged,
            args.target,
            args.output / "domain_holdout",
            n_jobs=args.n_jobs,
            feature_mode="nested_pruned",
            conformal_alpha=args.conformal_alpha,
        )
        uncertainty_summary = summarize_uncertainty(uncertainty)
        uncertainty_summary.to_csv(
            args.output / "domain_holdout" / "domain_uncertainty_summary.csv", index=False
        )
        domain_rows = []
        for (model_name, domain_name), group in domain_predictions.groupby(
            ["model", "outer_group"]
        ):
            from cml_repairnet.evaluation.metrics import regression_metrics
            domain_rows.append(
                {
                    "model": model_name,
                    "domain": domain_name,
                    "n_sites": len(group),
                    **regression_metrics(
                        group["y_true"].to_numpy(dtype=float),
                        group["y_pred"].to_numpy(dtype=float),
                    ),
                }
            )
        pd.DataFrame(domain_rows).to_csv(
            args.output / "domain_holdout" / "domain_metrics_by_domain.csv", index=False
        )
        if not uncertainty.empty:
            uncertainty_by_domain = (
                uncertainty.groupby("outer_group", as_index=False)
                .agg(
                    n_sites=("covered", "size"),
                    coverage=("covered", "mean"),
                    mean_interval_width=("interval_width", "mean"),
                    median_interval_width=("interval_width", "median"),
                    mae=("absolute_error", "mean"),
                    ood_count=("ood_flag", "sum"),
                )
            )
            uncertainty_by_domain.to_csv(
                args.output / "domain_holdout" / "domain_uncertainty_by_domain.csv",
                index=False,
            )
        print(domain_metrics[["model", "mae", "spearman_rho", "r2"]].to_string(index=False))
        print(uncertainty_summary.to_string(index=False))

    if not args.skip_cysteine:
        print("Running cysteine-neighbourhood mechanism analysis")
        cysteine_summary, _ = run_cysteine_mechanism_analysis(
            merged,
            args.target,
            args.output / "cysteine_mechanism",
            n_spatial_groups=4,
            permutation_repeats=args.conditional_repeats,
            n_jobs=args.n_jobs,
        )
        print(cysteine_summary.to_string(index=False))

    config = {
        "stage3_patch_version": "0.4.1",
        "mode": args.mode,
        "target": args.target,
        "seed_count": len(seeds),
        "group_counts": group_counts,
        "n_jobs": args.n_jobs,
        "conditional_repeats": args.conditional_repeats,
        "conformal_alpha": args.conformal_alpha,
        "primary_model": "extra_trees",
        "ranking_model": "lightgbm",
        "interpretation_model": "ridge",
    }
    (args.output / "stage3_config.json").write_text(
        json.dumps(config, indent=2, sort_keys=True), encoding="utf-8"
    )
    _write_report(
        args.output,
        repeated_summary,
        domain_metrics,
        uncertainty_summary,
        cysteine_summary,
        config,
    )
    print(f"Report: {args.output / 'stage3_report.md'}")


if __name__ == "__main__":
    main()
