from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.model_selection import GridSearchCV, GroupKFold

from cml_repairnet.evaluation.robustness import build_feature_view
from cml_repairnet.models.baselines import _numeric_feature_columns, _spatial_groups
from cml_repairnet.models.stage25 import build_sample_weights, stage25_model_specs
from cml_repairnet.utils.config import load_config
from cml_repairnet.utils.reproducibility import set_global_seed


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Fit a descriptive full-data tree model and export SHAP values."
    )
    parser.add_argument("--model", choices=["extra_trees", "xgboost", "lightgbm"], default="extra_trees")
    parser.add_argument("--target", default="repair_published")
    parser.add_argument("--groups", type=int, default=4)
    parser.add_argument("--weight-mode", choices=["none", "clipped", "sqrt"], default="none")
    parser.add_argument("--feature-mode", choices=["full", "nested_pruned"], default="full")
    parser.add_argument("--config", default="configs/base.yaml")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    try:
        import shap
    except ImportError as exc:
        raise SystemExit("SHAP is not installed. Run: python -m pip install -e '.[experiment]'") from exc

    config = load_config(args.config)
    seed = int(config["project"]["seed"])
    set_global_seed(seed)
    paths = config["paths"]

    labels = pd.read_csv(paths["labels"])
    sequence = pd.read_csv(paths["sequence_features"])
    structure = pd.read_csv(paths["structure_features"])
    merged = build_feature_view(labels, sequence, structure, "structure")
    features = _numeric_feature_columns(merged)
    X = merged[features].copy()
    y = merged[args.target].to_numpy(dtype=float)
    groups = _spatial_groups(merged, args.groups)
    weights = build_sample_weights(merged, args.weight_mode)

    spec = stage25_model_specs(seed, args.feature_mode)[args.model]
    search = GridSearchCV(
        spec.estimator,
        spec.param_grid,
        scoring="neg_mean_absolute_error",
        cv=GroupKFold(n_splits=args.groups),
        n_jobs=1,
        refit=True,
        error_score="raise",
    )
    search.fit(X, y, groups=groups, model__sample_weight=weights)
    pipeline = search.best_estimator_

    transformed = pipeline[:-1].transform(X)
    pruner = pipeline.named_steps["prune"]
    selected_features = [str(v) for v in pruner.get_feature_names_out(features)]
    model = pipeline.named_steps["model"]
    explainer = shap.TreeExplainer(model)
    shap_values = explainer(transformed)

    output_dir = Path(paths["results_dir"]) / "stage25" / "shap" / args.model
    output_dir.mkdir(parents=True, exist_ok=True)
    values = np.asarray(shap_values.values)
    pd.DataFrame(values, columns=selected_features).assign(
        protein_id=merged["protein_id"].to_numpy(),
        site_uniprot=merged["site_uniprot"].to_numpy(),
        site_paper=merged["site_paper"].to_numpy(),
        y_true=y,
    ).to_csv(output_dir / "shap_values.csv", index=False)
    pd.DataFrame(
        {
            "feature": selected_features,
            "mean_abs_shap": np.mean(np.abs(values), axis=0),
        }
    ).sort_values("mean_abs_shap", ascending=False).to_csv(
        output_dir / "shap_feature_importance.csv", index=False
    )

    shap.summary_plot(values, transformed, feature_names=selected_features, show=False)
    plt.tight_layout()
    plt.savefig(output_dir / "shap_beeswarm.png", dpi=240, bbox_inches="tight")
    plt.close()
    shap.summary_plot(
        values,
        transformed,
        feature_names=selected_features,
        plot_type="bar",
        show=False,
    )
    plt.tight_layout()
    plt.savefig(output_dir / "shap_bar.png", dpi=240, bbox_inches="tight")
    plt.close()

    (output_dir / "best_params.txt").write_text(str(search.best_params_), encoding="utf-8")
    print(f"SHAP outputs: {output_dir}")
    print("Important: SHAP is descriptive; it is not an additional performance estimate.")


if __name__ == "__main__":
    main()
