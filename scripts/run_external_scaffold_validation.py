from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd
import yaml

from cml_repairnet.revision.external_scaffold import (
    DEFAULT_FEATURES,
    QualityThresholds,
    add_quality_flags,
    aggregate_predictions,
    bootstrap_rank_probability,
    compare_rankings,
    data_fingerprint,
    fit_and_score_models,
    fit_ood_reference,
    prepare_external_sites,
    score_ood,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Frozen protein-scaffold external-consistency analysis for CML-RepairNet"
    )
    parser.add_argument("command", choices=["doctor", "prepare", "score", "analyze", "all"])
    parser.add_argument("--config", type=Path, default=Path("configs/external_scaffold_validation.yaml"))
    parser.add_argument("--bootstrap-iterations", type=int)
    return parser.parse_args()


def load_config(path: Path) -> dict:
    return yaml.safe_load(path.read_text(encoding="utf-8"))


def paths_from_config(config: dict) -> dict[str, Path]:
    return {key: Path(value) for key, value in config["paths"].items()}


def doctor(config: dict, paths: dict[str, Path]) -> dict:
    required = [
        paths["labels_csv"],
        paths["structure_features_csv"],
        paths["source_labels_csv"],
        paths["scaffold_panel_csv"],
    ]
    missing = [str(path) for path in required if not path.exists()]
    if missing:
        raise FileNotFoundError(f"Missing required files: {missing}")
    labels = pd.read_csv(paths["labels_csv"])
    structures = pd.read_csv(paths["structure_features_csv"])
    source = pd.read_csv(paths["source_labels_csv"])
    panel = pd.read_csv(paths["scaffold_panel_csv"])
    feature_columns = config["model"].get("transfer_safe_features", list(DEFAULT_FEATURES))
    absent = [column for column in feature_columns if column not in structures]
    if absent:
        raise ValueError(f"BSA structure table lacks features: {absent}")
    if set(source["scaffold"]) != {"casein", "collagen", "hemoglobin"}:
        raise ValueError("Source labels must contain casein, collagen and hemoglobin only")
    report = {
        "status": "PASS",
        "n_bsa_labels": len(labels),
        "n_bsa_structure_rows": len(structures),
        "n_external_conditions": len(source),
        "n_external_components": len(panel),
        "external_scaffolds": sorted(source["scaffold"].unique().tolist()),
        "feature_columns": feature_columns,
        "input_fingerprint": data_fingerprint(required),
        "claim_boundary": "external consistency only; not site-level external validation",
    }
    output_dir = paths["output_dir"]
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "external_scaffold_doctor.json").write_text(
        json.dumps(report, indent=2), encoding="utf-8"
    )
    print(json.dumps(report, indent=2))
    return report


def prepare(config: dict, paths: dict[str, Path]) -> None:
    panel = pd.read_csv(paths["scaffold_panel_csv"])
    sites, manifest = prepare_external_sites(panel, paths["cache_dir"])
    output_dir = paths["output_dir"]
    output_dir.mkdir(parents=True, exist_ok=True)
    sites.to_csv(output_dir / "external_scaffold_structure_features.csv", index=False)
    manifest.to_csv(output_dir / "external_scaffold_component_manifest.csv", index=False)
    print(manifest.to_string(index=False))
    if sites.empty:
        raise RuntimeError("No external scaffold sites were prepared")


def score(config: dict, paths: dict[str, Path]) -> None:
    output_dir = paths["output_dir"]
    sites = pd.read_csv(output_dir / "external_scaffold_structure_features.csv")
    labels = pd.read_csv(paths["labels_csv"])
    structure = pd.read_csv(paths["structure_features_csv"])
    feature_columns = config["model"].get("transfer_safe_features", list(DEFAULT_FEATURES))
    scored, context = fit_and_score_models(
        labels=labels,
        structure_features=structure,
        external_sites=sites,
        target=config["model"]["target"],
        feature_columns=feature_columns,
        model_names=config["model"]["models"],
        seed=int(config["seed"]),
        model_output=output_dir / "external_scaffold_frozen_models.joblib",
    )
    reference = fit_ood_reference(
        context["train"], feature_columns, int(config["quality"]["ood_neighbors"])
    )
    scored = score_ood(scored, reference, feature_columns)
    thresholds = QualityThresholds(
        strict_site_plddt=float(config["quality"]["strict_site_plddt"]),
        strict_neighbor_plddt=float(config["quality"]["strict_neighbor_plddt"]),
        relaxed_site_plddt=float(config["quality"]["relaxed_site_plddt"]),
        relaxed_neighbor_plddt=float(config["quality"]["relaxed_neighbor_plddt"]),
        ood_percentile_cutoff=float(config["quality"]["ood_percentile_cutoff"]),
        ood_neighbors=int(config["quality"]["ood_neighbors"]),
    )
    scored = add_quality_flags(scored, thresholds)
    scored.to_csv(output_dir / "external_scaffold_site_predictions.csv", index=False)
    print(
        scored.groupby("scaffold").agg(
            sites=("site_uniprot", "size"),
            strict=("quality_strict", "sum"),
            relaxed=("quality_relaxed", "sum"),
            in_domain=("in_domain", "sum"),
        ).to_string()
    )


def analyze(config: dict, paths: dict[str, Path], iterations: int | None = None) -> None:
    output_dir = paths["output_dir"]
    sites = pd.read_csv(output_dir / "external_scaffold_site_predictions.csv")
    labels = pd.read_csv(paths["source_labels_csv"])
    aggregates = aggregate_predictions(sites)
    rankings = compare_rankings(aggregates, labels)
    bootstrap = bootstrap_rank_probability(
        sites,
        labels,
        model=config["model"]["primary_model"],
        quality_mode=config["analysis"]["primary_quality"],
        domain_mode=config["analysis"]["primary_domain"],
        component_scheme=config["analysis"]["primary_component_scheme"],
        iterations=int(iterations or config["analysis"]["bootstrap_iterations"]),
        seed=int(config["seed"]),
    )
    aggregates.to_csv(output_dir / "external_scaffold_aggregates.csv", index=False)
    rankings.to_csv(output_dir / "external_scaffold_rank_consistency.csv", index=False)
    bootstrap.to_csv(output_dir / "external_scaffold_rank_bootstrap.csv", index=False)

    primary = rankings[
        (rankings["quality_mode"] == config["analysis"]["primary_quality"])
        & (rankings["domain_mode"] == config["analysis"]["primary_domain"])
        & (rankings["component_scheme"] == config["analysis"]["primary_component_scheme"])
        & (rankings["model"] == config["model"]["primary_model"])
        & (rankings["aggregation"] == config["analysis"]["primary_aggregation"])
    ]
    report = [
        "# CML-RepairNet external scaffold consistency report",
        "",
        "## Claim boundary",
        "This is a frozen protein-level external-consistency analysis against bulk ELISA outcomes. ",
        "It is not site-resolved external validation and it does not increase the 33-site training sample.",
        "",
        "## Prespecified observed order",
        "Casein > hemoglobin > collagen across all three reported modification levels.",
        "",
        "## Primary analysis",
    ]
    if primary.empty:
        report.append("Primary result unavailable because one or more scaffolds lacked strict in-domain coverage.")
    else:
        row = primary.iloc[0]
        report.extend(
            [
                f"- Spearman rank agreement: {row['spearman']:.3f}",
                f"- Kendall tau: {row['kendall_tau']:.3f}",
                f"- Exact order recovered: {bool(row['exact_rank_match'])}",
                f"- Pairwise ordering accuracy: {row['pairwise_order_accuracy']:.3f}",
                f"- Minimum site coverage: {row['minimum_site_coverage']:.3f}",
                f"- Descriptive calibration MAE: {row['calibration_mae_descriptive']:.3f}",
            ]
        )
    if "exact_rank_probability" in bootstrap and len(bootstrap):
        report.extend(
            [
                "",
                "## Site-bootstrap stability",
                f"Exact-order recovery probability: {bootstrap['exact_rank_probability'].iloc[0]:.3f}",
            ]
        )
    report.extend(
        [
            "",
            "## Mandatory limitations",
            "- The source publication did not identify exact commercial species/isoforms for the three purified proteins.",
            "- Bovine canonical components are a prespecified representative panel, not confirmed reagent identities.",
            "- AlphaFold monomers do not represent casein micelles, hemoglobin quaternary interfaces, or the collagen triple helix.",
            "- Bulk ELISA weights modified sites non-uniformly and does not provide site-resolved labels.",
            "- With three scaffolds, rank statistics are descriptive and no conventional inferential claim is warranted.",
        ]
    )
    (output_dir / "external_scaffold_report.md").write_text("\n".join(report), encoding="utf-8")
    print("\n".join(report))


def main() -> None:
    args = parse_args()
    config = load_config(args.config)
    paths = paths_from_config(config)
    if args.command in {"doctor", "all"}:
        doctor(config, paths)
    if args.command in {"prepare", "all"}:
        prepare(config, paths)
    if args.command in {"score", "all"}:
        score(config, paths)
    if args.command in {"analyze", "all"}:
        analyze(config, paths, args.bootstrap_iterations)


if __name__ == "__main__":
    main()
