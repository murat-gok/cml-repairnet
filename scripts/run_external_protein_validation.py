from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.metrics import mean_absolute_error


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Aggregate site-level predictions and compare them with published protein-level "
            "repair measurements. This is an external-consistency test, not site-level validation."
        )
    )
    parser.add_argument("--site-predictions", type=Path, required=True)
    parser.add_argument("--protein-labels", type=Path, required=True)
    parser.add_argument(
        "--output", type=Path, default=Path("results/protein_external_validation.csv")
    )
    parser.add_argument(
        "--aggregation", choices=["mean", "median", "confidence_weighted"], default="mean"
    )
    return parser.parse_args()


def aggregate(frame: pd.DataFrame, method: str) -> pd.DataFrame:
    required = {"protein_id", "pred_mean"}
    missing = required.difference(frame.columns)
    if missing:
        raise ValueError(f"Missing columns in site predictions: {sorted(missing)}")

    if method == "mean":
        grouped = frame.groupby("protein_id", as_index=False)["pred_mean"].mean()
    elif method == "median":
        grouped = frame.groupby("protein_id", as_index=False)["pred_mean"].median()
    else:
        if "pred_std" not in frame.columns:
            raise ValueError("confidence_weighted aggregation requires pred_std")
        work = frame.copy()
        work["weight"] = 1.0 / np.maximum(work["pred_std"].astype(float) ** 2, 1e-6)
        work["weighted"] = work["pred_mean"] * work["weight"]
        grouped = work.groupby("protein_id", as_index=False).agg(
            weighted=("weighted", "sum"), weight=("weight", "sum")
        )
        grouped["pred_mean"] = grouped["weighted"] / grouped["weight"]
        grouped = grouped[["protein_id", "pred_mean"]]
    return grouped.rename(columns={"pred_mean": "predicted_protein_repair"})


def main() -> None:
    args = parse_args()
    sites = pd.read_csv(args.site_predictions)
    labels = pd.read_csv(args.protein_labels)
    if not {"protein_id", "observed_protein_repair"}.issubset(labels.columns):
        raise ValueError("protein labels require protein_id and observed_protein_repair")

    merged = aggregate(sites, args.aggregation).merge(
        labels, on="protein_id", validate="one_to_one"
    )
    if len(merged) < 3:
        raise ValueError("At least three proteins are required for an informative consistency test")
    rho = spearmanr(merged["observed_protein_repair"], merged["predicted_protein_repair"]).statistic
    mae = mean_absolute_error(merged["observed_protein_repair"], merged["predicted_protein_repair"])
    merged["global_spearman"] = rho
    merged["global_mae"] = mae
    merged["interpretation"] = "external consistency only; not site-level validation"
    args.output.parent.mkdir(parents=True, exist_ok=True)
    merged.to_csv(args.output, index=False)
    print(merged.to_string(index=False))


if __name__ == "__main__":
    main()
