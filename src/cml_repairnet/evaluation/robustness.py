from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

KEY_COLUMNS = ["protein_id", "site_uniprot", "site_paper", "spatial_group"]


@dataclass(frozen=True)
class PairedComparison:
    model: str
    baseline: str
    model_mae: float
    baseline_mae: float
    absolute_improvement: float
    relative_improvement_pct: float
    improvement_ci_2_5: float
    improvement_ci_97_5: float
    permutation_p_one_sided: float
    n_sites: int

    def as_dict(self) -> dict[str, object]:
        return {
            "model": self.model,
            "baseline": self.baseline,
            "model_mae": self.model_mae,
            "baseline_mae": self.baseline_mae,
            "absolute_improvement": self.absolute_improvement,
            "relative_improvement_pct": self.relative_improvement_pct,
            "improvement_ci_2.5": self.improvement_ci_2_5,
            "improvement_ci_97.5": self.improvement_ci_97_5,
            "permutation_p_one_sided": self.permutation_p_one_sided,
            "n_sites": self.n_sites,
        }


def build_feature_view(
    labels: pd.DataFrame,
    sequence: pd.DataFrame,
    structure: pd.DataFrame,
    feature_set: str,
) -> pd.DataFrame:
    """Construct an ablation table while retaining coordinates only for spatial grouping."""
    keys = ["protein_id", "site_uniprot"]
    feature_set = feature_set.lower()
    if feature_set == "combined":
        return labels.merge(sequence, on=keys, validate="one_to_one").merge(
            structure, on=keys, validate="one_to_one"
        )
    if feature_set == "structure":
        return labels.merge(structure, on=keys, validate="one_to_one")
    if feature_set == "sequence":
        coordinate_columns = keys + ["ca_x", "ca_y", "ca_z", "structure_missing"]
        coordinates = structure[coordinate_columns].copy()
        return labels.merge(sequence, on=keys, validate="one_to_one").merge(
            coordinates, on=keys, validate="one_to_one"
        )
    raise ValueError("feature_set must be one of: sequence, structure, combined")


def paired_error_comparison(
    predictions: pd.DataFrame,
    model: str,
    baseline: str = "dummy_mean",
    *,
    bootstrap_iterations: int = 10_000,
    permutation_iterations: int = 20_000,
    seed: int = 20260727,
) -> PairedComparison:
    """Compare two OOF models on identical sites using paired uncertainty estimates."""
    required = set(KEY_COLUMNS + ["model", "y_true", "y_pred"])
    missing = required.difference(predictions.columns)
    if missing:
        raise ValueError(f"Missing prediction columns: {sorted(missing)}")

    columns = KEY_COLUMNS + ["y_true", "y_pred"]
    left = predictions.loc[predictions["model"] == model, columns].copy()
    right = predictions.loc[predictions["model"] == baseline, columns].copy()
    if left.empty:
        raise ValueError(f"Model {model!r} was not found in predictions")
    if right.empty:
        raise ValueError(f"Baseline {baseline!r} was not found in predictions")

    paired = left.merge(
        right,
        on=KEY_COLUMNS,
        suffixes=("_model", "_baseline"),
        validate="one_to_one",
    )
    if not np.allclose(paired["y_true_model"], paired["y_true_baseline"]):
        raise ValueError("Model and baseline y_true values do not match")

    y = paired["y_true_model"].to_numpy(dtype=float)
    model_error = np.abs(y - paired["y_pred_model"].to_numpy(dtype=float))
    baseline_error = np.abs(y - paired["y_pred_baseline"].to_numpy(dtype=float))
    site_improvement = baseline_error - model_error

    model_mae = float(model_error.mean())
    baseline_mae = float(baseline_error.mean())
    absolute_improvement = float(site_improvement.mean())
    relative_improvement = (
        100.0 * absolute_improvement / baseline_mae if baseline_mae > 0 else float("nan")
    )

    rng = np.random.default_rng(seed)
    n = len(site_improvement)
    bootstrap_values = np.empty(bootstrap_iterations, dtype=float)
    for i in range(bootstrap_iterations):
        idx = rng.integers(0, n, size=n)
        bootstrap_values[i] = site_improvement[idx].mean()
    ci_low, ci_high = np.quantile(bootstrap_values, [0.025, 0.975])

    exceedances = 0
    for _ in range(permutation_iterations):
        signs = rng.choice(np.array([-1.0, 1.0]), size=n)
        exceedances += bool(np.mean(signs * site_improvement) >= absolute_improvement)
    p_value = float((1 + exceedances) / (permutation_iterations + 1))

    return PairedComparison(
        model=model,
        baseline=baseline,
        model_mae=model_mae,
        baseline_mae=baseline_mae,
        absolute_improvement=absolute_improvement,
        relative_improvement_pct=float(relative_improvement),
        improvement_ci_2_5=float(ci_low),
        improvement_ci_97_5=float(ci_high),
        permutation_p_one_sided=p_value,
        n_sites=n,
    )


def write_markdown_summary(
    metrics: pd.DataFrame,
    comparisons: pd.DataFrame,
    output_path: Path,
) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    best = (
        metrics.sort_values("mae")
        .groupby(["feature_set", "target", "spatial_groups"], as_index=False)
        .first()
    )
    best_columns = [
        "feature_set", "target", "spatial_groups", "model", "mae",
        "spearman_rho", "r2", "n_features",
    ]
    comparison_columns = [
        "feature_set", "target", "spatial_groups", "model",
        "relative_improvement_pct", "improvement_ci_2.5",
        "improvement_ci_97.5", "permutation_p_one_sided",
    ]
    lines = [
        "# CML-RepairNet Stage 2 Robustness Report",
        "",
        "## Best model per analysis",
        "",
        "```text",
        best[best_columns].to_string(index=False),
        "```",
        "",
        "## Paired comparisons against dummy_mean",
        "",
        "```text",
        comparisons[comparison_columns].to_string(index=False),
        "```",
        "",
        "## Interpretation rule",
        "",
        "A model is considered robust only when its point estimate improves on the dummy "
        "baseline, the paired bootstrap interval is predominantly above zero, and the "
        "one-sided sign-flip permutation result is supportive. Feature-set and target "
        "sensitivity must also be reported; no single split is treated as definitive.",
        "",
    ]
    output_path.write_text("\n".join(lines), encoding="utf-8")
