from __future__ import annotations

import hashlib
import json
import math
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import networkx as nx
import numpy as np
import pandas as pd
import requests
from scipy.stats import spearmanr
from sklearn.base import clone
from sklearn.isotonic import IsotonicRegression
from sklearn.linear_model import LinearRegression, LogisticRegression
from sklearn.metrics import (
    mean_absolute_error,
    r2_score,
    roc_auc_score,
)
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler


TRANSFER_SAFE_FEATURES: tuple[str, ...] = (
    "residue_sasa",
    "relative_sasa_lysine",
    "distance_to_protein_centroid",
    "neighbor_hydrophobic_fraction_8a",
    "neighbor_cysteine_count_8a",
    "neighbor_aromatic_count_8a",
    "ca_contact_count_8a",
    "ca_contact_count_10a",
    "ca_contact_count_12a",
)

DEFAULT_SEED = 20260801
GPROFILER_URL = "https://biit.cs.ut.ee/gprofiler/api/gost/profile/"


@dataclass(frozen=True)
class Phase2Paths:
    labels: Path
    structure: Path
    predictions: Path
    atlas: Path
    portfolio: Path | None
    pareto: Path | None
    tiers: Path | None
    enrichment: Path | None
    string_edges: Path | None


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def stable_gene_set_fingerprint(genes: Iterable[str]) -> str:
    """Return a deterministic fingerprint for a background gene set."""
    normalized = sorted({str(gene).strip() for gene in genes if str(gene).strip()})
    return hashlib.sha256("\n".join(normalized).encode("utf-8")).hexdigest()


def is_ontology_root_term(source: object, term_id: object, term_name: object) -> bool:
    """Flag non-informative ontology/database root nodes."""
    source_text = str(source or "").strip().upper()
    term_id_text = str(term_id or "").strip().upper()
    term_name_text = str(term_name or "").strip().lower()
    explicit_roots = {
        ("CORUM", "CORUM:0000000"),
        ("KEGG", "KEGG:00000"),
        ("REAC", "REAC:0000000"),
        ("WP", "WP:000000"),
    }
    if (source_text, term_id_text) in explicit_roots:
        return True
    return (
        term_name_text.endswith(" root")
        or term_name_text.endswith(" root term")
        or term_name_text in {"wikipathways", "corum root"}
    )


def first_existing(candidates: Iterable[Path]) -> Path | None:
    for path in candidates:
        if path.exists():
            return path
    return None


def discover_paths(root: Path) -> Phase2Paths:
    root = root.resolve()
    labels = root / "data/processed/bsa_labels.csv"
    structure = root / "data/processed/bsa_structure_features.csv"
    predictions = first_existing(
        (
            root / "results/revision_confirmatory/revision_t1_ablation_oof_predictions.csv",
            root / "results/revision_tests/revision_t1_ablation_oof_predictions.csv",
        )
    ) or root / "results/revision_confirmatory/revision_t1_ablation_oof_predictions.csv"
    atlas = first_existing(
        (
            root / "results/stage4/atlas/cml_repairability_atlas_all_sites.csv",
            root / "data/processed/stage4/stage4_scored_candidates.csv",
            root / "supplementary_data/Table_S2_stage4_atlas_all_sites.csv",
        )
    ) or root / "results/stage4/atlas/cml_repairability_atlas_all_sites.csv"
    portfolio = first_existing(
        (
            root / "results/stage4/atlas/cml_repairability_priority_portfolio.csv",
            root / "supplementary_data/Table_S3_stage4_priority_portfolio.csv",
        )
    )
    pareto = first_existing(
        (
            root / "results/stage4/atlas/cml_repairability_pareto_front.csv",
            root / "supplementary_data/Table_S4_stage4_pareto_front.csv",
        )
    )
    tiers = first_existing(
        (
            root / "data/processed/stage6/stage6_candidate_tiers.csv",
            root / "results/stage6/stage6_candidate_tiers.csv",
            root / "supplementary_data/Table_S10_stage6_candidate_tiers.csv",
        )
    )
    enrichment = first_existing(
        (
            root / "data/processed/stage5/stage5_gprofiler_enrichment_all.csv",
            root / "supplementary_data/Table_S6_stage5_enrichment_all.csv",
        )
    )
    string_edges = first_existing(
        (
            root / "data/processed/stage5/stage5_string_edges.csv",
            root / "supplementary_data/Table_S7_stage5_string_edges.csv",
        )
    )
    return Phase2Paths(
        labels=labels,
        structure=structure,
        predictions=predictions,
        atlas=atlas,
        portfolio=portfolio,
        pareto=pareto,
        tiers=tiers,
        enrichment=enrichment,
        string_edges=string_edges,
    )


def _read_csv(path: Path, required: Sequence[str] = ()) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(path)
    frame = pd.read_csv(path)
    missing = [column for column in required if column not in frame.columns]
    if missing:
        raise ValueError(f"{path} is missing required columns: {missing}")
    return frame


def _safe_spearman(y_true: np.ndarray, y_pred: np.ndarray) -> float:
    result = spearmanr(y_true, y_pred, nan_policy="omit")
    return float(result.statistic) if np.isfinite(result.statistic) else float("nan")


def weighted_quantile(values: np.ndarray, weights: np.ndarray, quantile: float) -> float:
    values = np.asarray(values, dtype=float)
    weights = np.asarray(weights, dtype=float)
    valid = np.isfinite(values) & np.isfinite(weights) & (weights >= 0)
    values = values[valid]
    weights = weights[valid]
    if values.size == 0 or weights.sum() <= 0:
        return float("nan")
    if not 0 <= quantile <= 1:
        raise ValueError("quantile must be in [0, 1]")
    order = np.argsort(values, kind="mergesort")
    values = values[order]
    weights = weights[order]
    cumulative = np.cumsum(weights)
    threshold = quantile * cumulative[-1]
    index = int(np.searchsorted(cumulative, threshold, side="left"))
    return float(values[min(index, len(values) - 1)])


def conformal_level(n: int, alpha: float) -> float:
    if n < 1:
        return float("nan")
    return min(1.0, math.ceil((n + 1) * (1.0 - alpha)) / n)


def _matrix(frame: pd.DataFrame, columns: Sequence[str]) -> pd.DataFrame:
    missing = [column for column in columns if column not in frame]
    if missing:
        raise ValueError(f"Missing transfer-safe columns: {missing}")
    matrix = frame.loc[:, list(columns)].apply(pd.to_numeric, errors="coerce")
    matrix = matrix.replace([np.inf, -np.inf], np.nan)
    medians = matrix.median(axis=0).fillna(0.0)
    return matrix.fillna(medians)


def _primary_oof(
    predictions: pd.DataFrame,
    *,
    partition_id: str = "ward_g4",
    feature_set: str = "structure_transfer_safe",
    model: str = "extra_trees",
) -> pd.DataFrame:
    required = {
        "site_uniprot",
        "partition_id",
        "feature_set",
        "model",
        "y_true",
        "y_pred",
    }
    missing = sorted(required - set(predictions.columns))
    if missing:
        raise ValueError(f"OOF table missing columns: {missing}")
    selected = predictions.loc[
        (predictions["partition_id"].astype(str) == partition_id)
        & (predictions["feature_set"].astype(str) == feature_set)
        & (predictions["model"].astype(str) == model)
    ].copy()
    if selected.empty:
        raise ValueError(
            f"No OOF rows for partition={partition_id}, feature_set={feature_set}, model={model}"
        )
    selected = selected.drop_duplicates("site_uniprot").sort_values("site_uniprot")
    selected["absolute_error"] = np.abs(
        pd.to_numeric(selected["y_true"], errors="raise")
        - pd.to_numeric(selected["y_pred"], errors="raise")
    )
    return selected.reset_index(drop=True)


def _density_ratio_weights(
    source: pd.DataFrame,
    target: pd.DataFrame,
    *,
    seed: int,
    clip_low: float,
    clip_high: float,
) -> tuple[np.ndarray, dict[str, float]]:
    combined = pd.concat([source, target], ignore_index=True)
    labels = np.concatenate(
        [np.zeros(len(source), dtype=int), np.ones(len(target), dtype=int)]
    )
    pipeline = make_pipeline(
        StandardScaler(),
        LogisticRegression(
            C=1.0,
            class_weight="balanced",
            max_iter=5000,
            random_state=seed,
        ),
    )
    n_splits = min(5, len(source), len(target))
    if n_splits < 2:
        raise ValueError("Density-ratio estimation needs at least two source and target rows")
    cv = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=seed)
    probabilities = cross_val_predict(
        pipeline,
        combined,
        labels,
        cv=cv,
        method="predict_proba",
    )[:, 1]
    auc = roc_auc_score(labels, probabilities)
    source_probability = np.clip(probabilities[: len(source)], 1e-6, 1 - 1e-6)
    # With class-balanced logistic regression, posterior odds are proportional to p_target/p_source.
    weights = source_probability / (1.0 - source_probability)
    weights = np.clip(weights, clip_low, clip_high)
    effective_n = float(weights.sum() ** 2 / np.square(weights).sum())
    diagnostics = {
        "shift_classifier_cv_auc": float(auc),
        "weight_min": float(weights.min()),
        "weight_median": float(np.median(weights)),
        "weight_max": float(weights.max()),
        "weight_effective_sample_size": effective_n,
    }
    return weights, diagnostics


def _localized_q(
    source_z: np.ndarray,
    residuals: np.ndarray,
    target_z: np.ndarray,
    *,
    alpha: float,
    bandwidth: float,
) -> np.ndarray:
    output = np.empty(len(target_z), dtype=float)
    level = conformal_level(len(residuals), alpha)
    for index, row in enumerate(target_z):
        distances = np.linalg.norm(source_z - row, axis=1)
        weights = np.exp(-0.5 * np.square(distances / max(bandwidth, 1e-8)))
        weights = np.maximum(weights, 1e-8)
        output[index] = weighted_quantile(residuals, weights, level)
    return output


def _interval_columns(center: np.ndarray, q: np.ndarray | float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    q_array = np.broadcast_to(np.asarray(q, dtype=float), center.shape)
    low = np.clip(center - q_array, 0.0, 1.0)
    high = np.clip(center + q_array, 0.0, 1.0)
    return low, high, high - low


def run_t5_conformal_shift_sensitivity(
    labels_path: Path,
    structure_path: Path,
    predictions_path: Path,
    atlas_path: Path,
    output_dir: Path,
    *,
    alpha: float = 0.10,
    partition_id: str = "ward_g4",
    feature_set: str = "structure_transfer_safe",
    model: str = "extra_trees",
    density_ratio_clip: tuple[float, float] = (0.05, 20.0),
    ood_inflation_max: float = 2.0,
    seed: int = DEFAULT_SEED,
) -> dict[str, pd.DataFrame]:
    output_dir.mkdir(parents=True, exist_ok=True)
    labels = _read_csv(labels_path, required=("site_uniprot",))
    structure = _read_csv(structure_path, required=("site_uniprot", *TRANSFER_SAFE_FEATURES))
    predictions = _read_csv(predictions_path)
    atlas = _read_csv(
        atlas_path,
        required=("candidate_id", "repairability_extra_trees", *TRANSFER_SAFE_FEATURES),
    )
    primary = _primary_oof(
        predictions,
        partition_id=partition_id,
        feature_set=feature_set,
        model=model,
    )
    bsa = primary.merge(
        labels.drop_duplicates("site_uniprot"),
        on="site_uniprot",
        how="left",
        suffixes=("", "_label"),
    ).merge(
        structure.drop_duplicates("site_uniprot"),
        on="site_uniprot",
        how="left",
        suffixes=("", "_structure"),
    )
    if bsa[list(TRANSFER_SAFE_FEATURES)].isna().all(axis=1).any():
        missing_sites = bsa.loc[
            bsa[list(TRANSFER_SAFE_FEATURES)].isna().all(axis=1), "site_uniprot"
        ].tolist()
        raise ValueError(f"BSA features missing for sites: {missing_sites}")

    source_x = _matrix(bsa, TRANSFER_SAFE_FEATURES)
    target_x = _matrix(atlas, TRANSFER_SAFE_FEATURES)
    y = bsa["y_true"].to_numpy(dtype=float)
    center_bsa = bsa["y_pred"].to_numpy(dtype=float)
    residuals = np.abs(y - center_bsa)
    center_atlas = pd.to_numeric(atlas["repairability_extra_trees"], errors="raise").to_numpy()
    level = conformal_level(len(residuals), alpha)
    standard_q = weighted_quantile(residuals, np.ones(len(residuals)), level)

    ratio_weights, ratio_diagnostics = _density_ratio_weights(
        source_x,
        target_x,
        seed=seed,
        clip_low=density_ratio_clip[0],
        clip_high=density_ratio_clip[1],
    )
    weighted_q = weighted_quantile(residuals, ratio_weights, level)

    scaler = StandardScaler().fit(source_x)
    source_z = scaler.transform(source_x)
    target_z = scaler.transform(target_x)
    pairwise = np.linalg.norm(source_z[:, None, :] - source_z[None, :, :], axis=2)
    positive_pairwise = pairwise[np.triu_indices_from(pairwise, k=1)]
    bandwidth = float(np.median(positive_pairwise[positive_pairwise > 0]))
    localized_atlas_q = _localized_q(
        source_z,
        residuals,
        target_z,
        alpha=alpha,
        bandwidth=bandwidth,
    )

    atlas_ood = pd.to_numeric(
        atlas.get("ood_percentile", pd.Series(0.5, index=atlas.index)), errors="coerce"
    ).fillna(1.0).clip(0.0, 1.0).to_numpy()
    inflation = 1.0 + (ood_inflation_max - 1.0) * np.maximum(0.0, atlas_ood - 0.5) / 0.5
    inflated_q = standard_q * inflation

    methods: dict[str, np.ndarray] = {
        "standard_residual": np.full(len(atlas), standard_q),
        "density_ratio_weighted": np.full(len(atlas), weighted_q),
        "localized_distance_weighted": localized_atlas_q,
        "descriptive_ood_inflated": inflated_q,
    }
    atlas_intervals = atlas.copy()
    atlas_intervals["cmlase_response_score"] = center_atlas
    for method_name, q_values in methods.items():
        low, high, width = _interval_columns(center_atlas, q_values)
        atlas_intervals[f"{method_name}_q"] = q_values
        atlas_intervals[f"{method_name}_low"] = low
        atlas_intervals[f"{method_name}_high"] = high
        atlas_intervals[f"{method_name}_width"] = width

    # Leave-one-out coverage diagnostics on BSA; these are internal benchmark diagnostics only.
    bsa_rows: list[dict[str, Any]] = []
    bsa_detail_rows: list[dict[str, Any]] = []
    for method_name in ("standard_residual", "density_ratio_weighted", "localized_distance_weighted"):
        lows = np.empty(len(bsa), dtype=float)
        highs = np.empty(len(bsa), dtype=float)
        qs = np.empty(len(bsa), dtype=float)
        for index in range(len(bsa)):
            keep = np.arange(len(bsa)) != index
            local_level = conformal_level(int(keep.sum()), alpha)
            if method_name == "standard_residual":
                q = weighted_quantile(residuals[keep], np.ones(keep.sum()), local_level)
            elif method_name == "density_ratio_weighted":
                q = weighted_quantile(residuals[keep], ratio_weights[keep], local_level)
            else:
                distances = np.linalg.norm(source_z[keep] - source_z[index], axis=1)
                weights = np.exp(-0.5 * np.square(distances / max(bandwidth, 1e-8)))
                q = weighted_quantile(residuals[keep], np.maximum(weights, 1e-8), local_level)
            qs[index] = q
            lows[index] = max(0.0, center_bsa[index] - q)
            highs[index] = min(1.0, center_bsa[index] + q)
        covered = (y >= lows) & (y <= highs)
        bsa_rows.append(
            {
                "method": method_name,
                "alpha": alpha,
                "nominal_coverage": 1.0 - alpha,
                "observed_loo_coverage": float(covered.mean()),
                "mean_interval_width": float(np.mean(highs - lows)),
                "median_interval_width": float(np.median(highs - lows)),
                "mean_q": float(np.mean(qs)),
                "n_sites": len(bsa),
                "claim_boundary": "within-BSA diagnostic; no proteome coverage guarantee",
            }
        )
        for index, row in bsa.iterrows():
            bsa_detail_rows.append(
                {
                    "site_uniprot": row["site_uniprot"],
                    "site_paper": row.get("site_paper", np.nan),
                    "method": method_name,
                    "y_true": y[index],
                    "y_pred": center_bsa[index],
                    "q": qs[index],
                    "interval_low": lows[index],
                    "interval_high": highs[index],
                    "covered": bool(covered[index]),
                }
            )

    coverage = pd.DataFrame(bsa_rows)
    bsa_details = pd.DataFrame(bsa_detail_rows)
    applicability = atlas_intervals.get(
        "applicability_status", pd.Series("unknown", index=atlas_intervals.index)
    ).fillna("unknown").astype(str)
    width_rows: list[dict[str, Any]] = []
    for method_name in methods:
        temp = pd.DataFrame(
            {
                "applicability_status": applicability,
                "width": atlas_intervals[f"{method_name}_width"],
            }
        )
        for status, group in temp.groupby("applicability_status", dropna=False):
            width_rows.append(
                {
                    "method": method_name,
                    "applicability_status": status,
                    "n_candidates": len(group),
                    "mean_width": float(group["width"].mean()),
                    "median_width": float(group["width"].median()),
                    "q05_width": float(group["width"].quantile(0.05)),
                    "q95_width": float(group["width"].quantile(0.95)),
                }
            )
    widths = pd.DataFrame(width_rows)

    method_summary = pd.DataFrame(
        [
            {
                "method": "standard_residual",
                "global_q": standard_q,
                "is_formal_finite_sample_claimed": False,
                "interpretation": "internal OOF residual conformal sensitivity",
            },
            {
                "method": "density_ratio_weighted",
                "global_q": weighted_q,
                "is_formal_finite_sample_claimed": False,
                "interpretation": (
                    "estimated covariate-shift importance weighting; approximate sensitivity only"
                ),
            },
            {
                "method": "localized_distance_weighted",
                "global_q": np.nan,
                "is_formal_finite_sample_claimed": False,
                "interpretation": "candidate-local structural-distance sensitivity",
            },
            {
                "method": "descriptive_ood_inflated",
                "global_q": standard_q,
                "is_formal_finite_sample_claimed": False,
                "interpretation": "descriptive OOD inflation, not conformal coverage",
            },
        ]
    )
    diagnostics: dict[str, Any] = {
        "analysis": "T5 covariate-shift conformal sensitivity",
        "alpha": alpha,
        "partition_id": partition_id,
        "feature_set": feature_set,
        "model": model,
        "n_bsa_sites": len(bsa),
        "n_atlas_candidates": len(atlas),
        "standard_q": standard_q,
        "density_ratio_weighted_q": weighted_q,
        "localized_bandwidth": bandwidth,
        "density_ratio_clip": list(density_ratio_clip),
        "ood_inflation_max": ood_inflation_max,
        "claim_boundary": (
            "all atlas intervals are sensitivity/descriptive intervals; no cross-protein coverage guarantee"
        ),
        **ratio_diagnostics,
        "input_sha256": {
            "labels": sha256_file(labels_path),
            "structure": sha256_file(structure_path),
            "predictions": sha256_file(predictions_path),
            "atlas": sha256_file(atlas_path),
        },
    }

    outputs = {
        "coverage": coverage,
        "bsa_details": bsa_details,
        "atlas_intervals": atlas_intervals,
        "widths": widths,
        "method_summary": method_summary,
    }
    filenames = {
        "coverage": "phase2_t5_bsa_coverage.csv",
        "bsa_details": "phase2_t5_bsa_interval_details.csv",
        "atlas_intervals": "phase2_t5_atlas_interval_sensitivity.csv",
        "widths": "phase2_t5_width_by_applicability.csv",
        "method_summary": "phase2_t5_interval_methods.csv",
    }
    for key, frame in outputs.items():
        frame.to_csv(output_dir / filenames[key], index=False)
    (output_dir / "phase2_t5_diagnostics.json").write_text(
        json.dumps(diagnostics, indent=2, sort_keys=True), encoding="utf-8"
    )
    return outputs


def _pairwise_concordance(y_true: np.ndarray, y_pred: np.ndarray) -> float:
    concordant = 0
    comparable = 0
    for i in range(len(y_true)):
        for j in range(i + 1, len(y_true)):
            true_delta = y_true[i] - y_true[j]
            pred_delta = y_pred[i] - y_pred[j]
            if true_delta == 0 or pred_delta == 0:
                continue
            comparable += 1
            concordant += int(np.sign(true_delta) == np.sign(pred_delta))
    return float(concordant / comparable) if comparable else float("nan")


def _calibration_fit(y_true: np.ndarray, y_pred: np.ndarray) -> tuple[float, float]:
    model = LinearRegression().fit(y_pred.reshape(-1, 1), y_true)
    return float(model.intercept_), float(model.coef_[0])


def _compression_fit(y_true: np.ndarray, y_pred: np.ndarray) -> tuple[float, float]:
    model = LinearRegression().fit(y_true.reshape(-1, 1), y_pred)
    return float(model.intercept_), float(model.coef_[0])


def _metric_row(name: str, y_true: np.ndarray, y_pred: np.ndarray) -> dict[str, Any]:
    intercept, slope = _calibration_fit(y_true, y_pred)
    compression_intercept, compression_slope = _compression_fit(y_true, y_pred)
    return {
        "prediction_variant": name,
        "mae": float(mean_absolute_error(y_true, y_pred)),
        "spearman_rho": _safe_spearman(y_true, y_pred),
        "r2": float(r2_score(y_true, y_pred)),
        "calibration_intercept_y_on_pred": intercept,
        "calibration_slope_y_on_pred": slope,
        "prediction_intercept_pred_on_y": compression_intercept,
        "prediction_compression_slope_pred_on_y": compression_slope,
        "pairwise_concordance": _pairwise_concordance(y_true, y_pred),
    }


def run_t6_boundary_shrinkage(
    predictions_path: Path,
    output_dir: Path,
    *,
    feature_set: str = "structure_transfer_safe",
    model: str = "extra_trees",
    low_threshold: float = 0.10,
    high_threshold: float = 0.90,
    bootstrap_iterations: int = 20_000,
    seed: int = DEFAULT_SEED,
) -> dict[str, pd.DataFrame]:
    output_dir.mkdir(parents=True, exist_ok=True)
    predictions = _read_csv(
        predictions_path,
        required=(
            "site_uniprot",
            "partition_id",
            "feature_set",
            "model",
            "y_true",
            "y_pred",
        ),
    )
    selected = predictions.loc[
        (predictions["feature_set"].astype(str) == feature_set)
        & (predictions["model"].astype(str) == model)
    ].copy()
    if selected.empty:
        raise ValueError("No rows for requested T6 model/feature set")
    site = (
        selected.groupby("site_uniprot", as_index=False)
        .agg(
            site_paper=("site_paper", "first") if "site_paper" in selected else ("site_uniprot", "first"),
            y_true=("y_true", "first"),
            y_pred_mean=("y_pred", "mean"),
            y_pred_sd=("y_pred", "std"),
            y_pred_min=("y_pred", "min"),
            y_pred_max=("y_pred", "max"),
            n_partitions=("partition_id", "nunique"),
        )
        .sort_values("site_uniprot")
        .reset_index(drop=True)
    )
    y = site["y_true"].to_numpy(dtype=float)
    raw = site["y_pred_mean"].to_numpy(dtype=float)
    loo_affine = np.empty(len(site), dtype=float)
    loo_isotonic = np.empty(len(site), dtype=float)
    for index in range(len(site)):
        keep = np.arange(len(site)) != index
        affine = LinearRegression().fit(raw[keep].reshape(-1, 1), y[keep])
        loo_affine[index] = np.clip(affine.predict(raw[index].reshape(1, -1))[0], 0.0, 1.0)
        if np.unique(raw[keep]).size >= 4:
            isotonic = IsotonicRegression(y_min=0.0, y_max=1.0, out_of_bounds="clip")
            isotonic.fit(raw[keep], y[keep])
            loo_isotonic[index] = isotonic.predict([raw[index]])[0]
        else:
            loo_isotonic[index] = raw[index]
    site["raw_absolute_error"] = np.abs(y - raw)
    site["loo_affine_prediction"] = loo_affine
    site["loo_affine_absolute_error"] = np.abs(y - loo_affine)
    site["loo_isotonic_prediction"] = loo_isotonic
    site["loo_isotonic_absolute_error"] = np.abs(y - loo_isotonic)
    site["response_stratum"] = np.select(
        [y <= low_threshold, y >= high_threshold],
        ["low_boundary", "high_boundary"],
        default="middle",
    )
    site["shrinkage_signed"] = raw - y

    summary = pd.DataFrame(
        [
            _metric_row("raw_partition_mean", y, raw),
            _metric_row("loo_affine_diagnostic", y, loo_affine),
            _metric_row("loo_isotonic_diagnostic", y, loo_isotonic),
        ]
    )
    strata = (
        site.groupby("response_stratum", as_index=False)
        .agg(
            n_sites=("site_uniprot", "size"),
            mean_y_true=("y_true", "mean"),
            mean_y_pred=("y_pred_mean", "mean"),
            mae=("raw_absolute_error", "mean"),
            mean_signed_shrinkage=("shrinkage_signed", "mean"),
            median_signed_shrinkage=("shrinkage_signed", "median"),
        )
    )
    site["prediction_bin"] = pd.qcut(
        site["y_pred_mean"], q=min(5, len(site)), duplicates="drop"
    ).astype(str)
    calibration_bins = (
        site.groupby("prediction_bin", as_index=False, observed=True)
        .agg(
            n_sites=("site_uniprot", "size"),
            mean_prediction=("y_pred_mean", "mean"),
            mean_observed=("y_true", "mean"),
            mae=("raw_absolute_error", "mean"),
        )
        .sort_values("mean_prediction")
    )

    rng = np.random.default_rng(seed)
    bootstrap_rows: list[dict[str, float]] = []
    indices = np.arange(len(site))
    for iteration in range(int(bootstrap_iterations)):
        sample = rng.choice(indices, size=len(indices), replace=True)
        sample_y = y[sample]
        sample_pred = raw[sample]
        if np.unique(sample_pred).size < 2:
            continue
        intercept, slope = _calibration_fit(sample_y, sample_pred)
        compression_intercept, compression_slope = _compression_fit(sample_y, sample_pred)
        low_mask = sample_y <= low_threshold
        high_mask = sample_y >= high_threshold
        bootstrap_rows.append(
            {
                "iteration": iteration,
                "calibration_intercept_y_on_pred": intercept,
                "calibration_slope_y_on_pred": slope,
                "prediction_intercept_pred_on_y": compression_intercept,
                "prediction_compression_slope_pred_on_y": compression_slope,
                "mae": float(mean_absolute_error(sample_y, sample_pred)),
                "spearman_rho": _safe_spearman(sample_y, sample_pred),
                "low_boundary_signed_shrinkage": (
                    float(np.mean(sample_pred[low_mask] - sample_y[low_mask]))
                    if low_mask.any()
                    else np.nan
                ),
                "high_boundary_signed_shrinkage": (
                    float(np.mean(sample_pred[high_mask] - sample_y[high_mask]))
                    if high_mask.any()
                    else np.nan
                ),
            }
        )
    bootstrap = pd.DataFrame(bootstrap_rows)
    ci_rows: list[dict[str, Any]] = []
    for column in (
        "calibration_intercept_y_on_pred",
        "calibration_slope_y_on_pred",
        "prediction_intercept_pred_on_y",
        "prediction_compression_slope_pred_on_y",
        "mae",
        "spearman_rho",
        "low_boundary_signed_shrinkage",
        "high_boundary_signed_shrinkage",
    ):
        values = pd.to_numeric(bootstrap[column], errors="coerce").dropna()
        ci_rows.append(
            {
                "metric": column,
                "bootstrap_mean": float(values.mean()),
                "ci_2.5": float(values.quantile(0.025)),
                "ci_97.5": float(values.quantile(0.975)),
                "valid_draws": len(values),
            }
        )
    bootstrap_summary = pd.DataFrame(ci_rows)
    recommendation = {
        "analysis": "T6 boundary shrinkage",
        "feature_set": feature_set,
        "model": model,
        "low_threshold": low_threshold,
        "high_threshold": high_threshold,
        "recommended_atlas_name": "cmlase_response_score",
        "do_not_call": "experimentally calibrated repair fraction outside BSA",
        "claim_boundary": (
            "absolute response estimates shrink toward the center; atlas use should emphasize relative ranking"
        ),
        "input_sha256": sha256_file(predictions_path),
    }

    outputs = {
        "summary": summary,
        "strata": strata,
        "site": site,
        "calibration_bins": calibration_bins,
        "bootstrap": bootstrap,
        "bootstrap_summary": bootstrap_summary,
    }
    filenames = {
        "summary": "phase2_t6_boundary_summary.csv",
        "strata": "phase2_t6_response_strata.csv",
        "site": "phase2_t6_site_diagnostics.csv",
        "calibration_bins": "phase2_t6_calibration_bins.csv",
        "bootstrap": "phase2_t6_bootstrap_draws.csv",
        "bootstrap_summary": "phase2_t6_bootstrap_summary.csv",
    }
    for key, frame in outputs.items():
        frame.to_csv(output_dir / filenames[key], index=False)
    (output_dir / "phase2_t6_recommendation.json").write_text(
        json.dumps(recommendation, indent=2, sort_keys=True), encoding="utf-8"
    )
    return outputs


def classify_subcellular_location(value: object) -> dict[str, Any]:
    text = "" if value is None or pd.isna(value) else str(value).lower()
    flags = {
        "flag_extracellular_secreted": any(
            token in text
            for token in (
                "secreted",
                "extracellular space",
                "extracellular matrix",
                "basement membrane",
            )
        ),
        "flag_cell_surface": any(
            token in text
            for token in (
                "cell surface",
                "extracellular side",
                "external side",
                "extracellular domain",
            )
        ),
        "flag_secretory_lumen": any(
            token in text
            for token in (
                "endoplasmic reticulum lumen",
                "golgi lumen",
                "secretory vesicle lumen",
                "lysosome lumen",
            )
        ),
        "flag_mitochondrial": "mitochond" in text,
        "flag_nuclear": any(token in text for token in ("nucleus", "nuclear")),
        "flag_cytosolic": any(token in text for token in ("cytoplasm", "cytosol")),
        "flag_other_organellar": any(
            token in text
            for token in (
                "endosome",
                "lysosome",
                "peroxisome",
                "golgi",
                "endoplasmic reticulum",
                "melanosome",
                "vesicle",
            )
        ),
    }
    major = [name for name, state in flags.items() if state]
    accessible = flags["flag_extracellular_secreted"] or flags["flag_cell_surface"]
    intracellular = (
        flags["flag_cytosolic"]
        or flags["flag_nuclear"]
        or flags["flag_mitochondrial"]
        or flags["flag_other_organellar"]
    )
    if accessible and intracellular:
        primary = "mixed_accessible_intracellular"
        layer = "U_mixed_or_unknown"
    elif flags["flag_extracellular_secreted"]:
        primary = "extracellular_secreted"
        layer = "A_extracellular_or_surface"
    elif flags["flag_cell_surface"]:
        primary = "cell_surface"
        layer = "A_extracellular_or_surface"
    elif flags["flag_secretory_lumen"]:
        primary = "secretory_lumen"
        layer = "B_secretory_lumen"
    elif flags["flag_mitochondrial"]:
        primary = "mitochondrial"
        layer = "D_mitochondrial"
    elif flags["flag_nuclear"]:
        primary = "nuclear"
        layer = "C_intracellular"
    elif flags["flag_cytosolic"]:
        primary = "cytosolic"
        layer = "C_intracellular"
    elif flags["flag_other_organellar"]:
        primary = "other_organellar"
        layer = "C_intracellular"
    else:
        primary = "unknown"
        layer = "U_mixed_or_unknown"
    return {
        **flags,
        "compartment_primary": primary,
        "accessibility_layer": layer,
        "compartment_flag_count": len(major),
        "compartment_ambiguous": len(major) > 1,
    }


def _gene_symbol(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    return str(value).strip().split()[0]


def run_t7_compartment_layering(
    atlas_path: Path,
    output_dir: Path,
    *,
    tiers_path: Path | None = None,
) -> dict[str, pd.DataFrame]:
    output_dir.mkdir(parents=True, exist_ok=True)
    atlas = _read_csv(
        atlas_path,
        required=("candidate_id", "species", "subcellular_location", "repairability_extra_trees"),
    )
    classifications = pd.DataFrame(
        [classify_subcellular_location(value) for value in atlas["subcellular_location"]]
    )
    layered = pd.concat([atlas.reset_index(drop=True), classifications], axis=1)
    layered["gene_symbol"] = layered.get(
        "gene_symbol", layered.get("gene_names", pd.Series("", index=layered.index))
    ).map(_gene_symbol)
    layered["cmlase_response_score"] = pd.to_numeric(
        layered["repairability_extra_trees"], errors="coerce"
    )
    if tiers_path is not None and tiers_path.exists():
        tiers = _read_csv(tiers_path, required=("candidate_id",))
        tier_columns = [
            column
            for column in (
                "candidate_id",
                "candidate_tier",
                "tier_reason",
                "mean_stress_rank",
                "top100_selection_frequency",
                "pareto_membership_frequency",
            )
            if column in tiers
        ]
        layered = layered.merge(
            tiers[tier_columns].drop_duplicates("candidate_id"),
            on="candidate_id",
            how="left",
            suffixes=("", "_stage6"),
        )
    applicability = layered.get(
        "applicability_status", pd.Series("unknown", index=layered.index)
    ).fillna("unknown").astype(str)
    layered["biomedical_interpretation_class"] = np.select(
        [
            layered["accessibility_layer"].eq("A_extracellular_or_surface")
            & applicability.eq("in_domain"),
            layered["accessibility_layer"].eq("A_extracellular_or_surface"),
            applicability.eq("in_domain"),
        ],
        [
            "higher_context_plausibility",
            "accessible_but_structurally_uncertain",
            "structurally_resembling_but_delivery_unresolved",
        ],
        default="secondary_hypothesis_only",
    )
    layered["hard_exclusion"] = False
    layered["interpretation_note"] = np.where(
        layered["accessibility_layer"].eq("A_extracellular_or_surface"),
        "extracellular or surface context; enzyme access is comparatively more plausible but unvalidated",
        "intracellular/organelle delivery has not been demonstrated; retain as secondary hypothesis",
    )

    summary_columns = ["species", "accessibility_layer", "compartment_primary"]
    summary = (
        layered.groupby(summary_columns, dropna=False, as_index=False)
        .agg(
            n_sites=("candidate_id", "size"),
            n_proteins=("uniprot_accession", "nunique")
            if "uniprot_accession" in layered
            else ("gene_symbol", "nunique"),
            mean_response_score=("cmlase_response_score", "mean"),
            median_response_score=("cmlase_response_score", "median"),
            in_domain_sites=(
                "applicability_status",
                lambda values: int((values.astype(str) == "in_domain").sum()),
            )
            if "applicability_status" in layered
            else ("candidate_id", lambda values: 0),
        )
    )
    if "candidate_tier" in layered:
        tier_summary = (
            layered.groupby(
                ["species", "accessibility_layer", "candidate_tier"],
                dropna=False,
                as_index=False,
            )
            .agg(n_sites=("candidate_id", "size"), n_proteins=("gene_symbol", "nunique"))
        )
    else:
        tier_summary = pd.DataFrame(
            columns=["species", "accessibility_layer", "candidate_tier", "n_sites", "n_proteins"]
        )
    actionability = (
        layered.groupby(
            ["species", "biomedical_interpretation_class"], dropna=False, as_index=False
        )
        .agg(
            n_sites=("candidate_id", "size"),
            n_proteins=("gene_symbol", "nunique"),
            mean_response_score=("cmlase_response_score", "mean"),
        )
    )
    outputs = {
        "layered": layered,
        "summary": summary,
        "tier_summary": tier_summary,
        "actionability": actionability,
    }
    filenames = {
        "layered": "phase2_t7_compartment_layered_atlas.csv",
        "summary": "phase2_t7_compartment_summary.csv",
        "tier_summary": "phase2_t7_tier_by_accessibility.csv",
        "actionability": "phase2_t7_biomedical_interpretation_summary.csv",
    }
    for key, frame in outputs.items():
        frame.to_csv(output_dir / filenames[key], index=False)
    manifest = {
        "analysis": "T7 compartment and accessibility layering",
        "hard_exclusion_applied": False,
        "claim_boundary": (
            "compartment layers contextualize plausibility; they do not establish delivery or efficacy"
        ),
        "input_sha256": {"atlas": sha256_file(atlas_path)},
    }
    if tiers_path is not None and tiers_path.exists():
        manifest["input_sha256"]["tiers"] = sha256_file(tiers_path)
    (output_dir / "phase2_t7_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8"
    )
    return outputs


def _protein_level(atlas: pd.DataFrame) -> pd.DataFrame:
    frame = atlas.copy()
    frame["gene_symbol"] = frame.get(
        "gene_symbol", frame.get("gene_names", pd.Series("", index=frame.index))
    ).map(_gene_symbol)
    frame["lysine_count"] = frame.get(
        "sequence", pd.Series("", index=frame.index)
    ).fillna("").astype(str).str.count("K")
    aggregations: dict[str, tuple[str, Any]] = {
        "gene_symbol": ("gene_symbol", "first"),
        "species": ("species", "first"),
        "protein_length": ("protein_length", "median")
        if "protein_length" in frame
        else ("site_uniprot", "max"),
        "lysine_count": ("lysine_count", "median"),
        "evidence_score": (
            "evidence_score_norm" if "evidence_score_norm" in frame else "evidence_score_raw",
            "max",
        ),
        "structure_quality_score": (
            "structure_quality_score",
            "median",
        )
        if "structure_quality_score" in frame
        else ("site_uniprot", lambda values: np.nan),
        "response_score": ("repairability_extra_trees", "max"),
        "candidate_sites": ("candidate_id", "size"),
    }
    if "uniprot_accession" in frame:
        group_key = "uniprot_accession"
    else:
        group_key = "gene_symbol"
    protein = frame.groupby(group_key, as_index=False).agg(**aggregations)
    if group_key != "uniprot_accession":
        protein["uniprot_accession"] = protein["gene_symbol"]
    return protein


def _query_sets(
    atlas: pd.DataFrame,
    portfolio: pd.DataFrame | None,
    pareto: pd.DataFrame | None,
) -> dict[tuple[str, str], set[str]]:
    frame = atlas.copy()
    frame["gene_symbol"] = frame.get(
        "gene_symbol", frame.get("gene_names", pd.Series("", index=frame.index))
    ).map(_gene_symbol)
    frame = frame.loc[frame["gene_symbol"].ne("")].copy()
    query_sets: dict[tuple[str, str], set[str]] = {}
    for species, sub in frame.groupby("species"):
        in_domain = sub.loc[
            sub.get("applicability_status", pd.Series("unknown", index=sub.index))
            .astype(str)
            .eq("in_domain")
        ].copy()
        if in_domain.empty:
            in_domain = sub.copy()
        threshold = pd.to_numeric(in_domain["repairability_extra_trees"], errors="coerce").quantile(0.90)
        query_sets[(str(species), "high_response")] = set(
            in_domain.loc[
                pd.to_numeric(in_domain["repairability_extra_trees"], errors="coerce") >= threshold,
                "gene_symbol",
            ]
        )
        evidence_column = (
            "evidence_score_norm" if "evidence_score_norm" in in_domain else "evidence_score_raw"
        )
        evidence_threshold = pd.to_numeric(in_domain[evidence_column], errors="coerce").quantile(0.75)
        query_sets[(str(species), "evidence_rich")] = set(
            in_domain.loc[
                pd.to_numeric(in_domain[evidence_column], errors="coerce") >= evidence_threshold,
                "gene_symbol",
            ]
        )
    if portfolio is not None and not portfolio.empty:
        portfolio = portfolio.copy()
        portfolio["gene_symbol"] = portfolio.get(
            "gene_symbol", portfolio.get("gene_names", pd.Series("", index=portfolio.index))
        ).map(_gene_symbol)
        for species, sub in portfolio.groupby("species"):
            query_sets[(str(species), "portfolio")] = set(sub["gene_symbol"]) - {""}
    if pareto is not None and not pareto.empty:
        pareto = pareto.copy()
        pareto["gene_symbol"] = pareto.get(
            "gene_symbol", pareto.get("gene_names", pd.Series("", index=pareto.index))
        ).map(_gene_symbol)
        for species, sub in pareto.groupby("species"):
            query_sets[(str(species), "pareto")] = set(sub["gene_symbol"]) - {""}
    return query_sets


def _background_sets(atlas: pd.DataFrame) -> dict[tuple[str, str], set[str]]:
    frame = atlas.copy()
    frame["gene_symbol"] = frame.get(
        "gene_symbol", frame.get("gene_names", pd.Series("", index=frame.index))
    ).map(_gene_symbol)
    frame = frame.loc[frame["gene_symbol"].ne("")].copy()
    backgrounds: dict[tuple[str, str], set[str]] = {}
    for species, sub in frame.groupby("species"):
        species = str(species)
        backgrounds[(species, "candidate_universe")] = set(sub["gene_symbol"])
        structure_status = sub.get(
            "structure_status", pd.Series("eligible", index=sub.index)
        ).astype(str)
        backgrounds[(species, "structure_eligible")] = set(
            sub.loc[structure_status.eq("eligible"), "gene_symbol"]
        )
        applicability = sub.get(
            "applicability_status", pd.Series("unknown", index=sub.index)
        ).astype(str)
        backgrounds[(species, "in_domain")] = set(
            sub.loc[applicability.eq("in_domain"), "gene_symbol"]
        )
    return backgrounds


def _matched_background(
    proteins: pd.DataFrame,
    query_genes: set[str],
    *,
    controls_per_query: int,
) -> tuple[set[str], pd.DataFrame]:
    frame = proteins.loc[proteins["gene_symbol"].ne("")].drop_duplicates("gene_symbol").copy()
    query_frame = frame.loc[frame["gene_symbol"].isin(query_genes)].copy()
    control_frame = frame.loc[~frame["gene_symbol"].isin(query_genes)].copy()
    feature_columns = ["protein_length", "lysine_count", "evidence_score", "structure_quality_score"]
    numeric = pd.concat([query_frame[feature_columns], control_frame[feature_columns]], ignore_index=True)
    numeric = numeric.apply(pd.to_numeric, errors="coerce")
    numeric["protein_length"] = np.log1p(numeric["protein_length"].clip(lower=0))
    numeric["lysine_count"] = np.log1p(numeric["lysine_count"].clip(lower=0))
    medians = numeric.median(axis=0).fillna(0.0)
    numeric = numeric.fillna(medians)
    scaler = StandardScaler().fit(numeric)
    query_z = scaler.transform(numeric.iloc[: len(query_frame)])
    control_z = scaler.transform(numeric.iloc[len(query_frame) :])
    selected_controls: set[str] = set()
    match_rows: list[dict[str, Any]] = []
    for query_index, query_row in query_frame.reset_index(drop=True).iterrows():
        distances = np.linalg.norm(control_z - query_z[query_index], axis=1)
        order = np.argsort(distances)
        picked = 0
        for control_index in order:
            gene = str(control_frame.iloc[control_index]["gene_symbol"])
            if gene in selected_controls:
                continue
            selected_controls.add(gene)
            match_rows.append(
                {
                    "query_gene": query_row["gene_symbol"],
                    "control_gene": gene,
                    "distance": float(distances[control_index]),
                }
            )
            picked += 1
            if picked >= controls_per_query:
                break
    background = set(query_genes) | selected_controls
    return background, pd.DataFrame(match_rows)


def _gprofiler_request(
    organism: str,
    query: Sequence[str],
    background: Sequence[str],
    *,
    sources: Sequence[str],
    cache_path: Path,
    timeout: float,
    retries: int,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    if cache_path.exists():
        payload = json.loads(cache_path.read_text(encoding="utf-8"))
        return payload.get("result", []), payload.get("meta", {})
    body = {
        "organism": organism,
        "query": sorted(set(query)),
        "background": sorted(set(background)),
        "domain_scope": "custom",
        "sources": list(sources),
        "user_threshold": 0.05,
        "significance_threshold_method": "fdr",
        "no_evidences": True,
        "ordered": False,
    }
    last_error: Exception | None = None
    for attempt in range(retries + 1):
        try:
            response = requests.post(
                GPROFILER_URL,
                json=body,
                headers={"User-Agent": "CML-RepairNet/phase2"},
                timeout=timeout,
            )
            response.raise_for_status()
            payload = response.json()
            cache_path.parent.mkdir(parents=True, exist_ok=True)
            cache_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
            return payload.get("result", []), payload.get("meta", {})
        except Exception as exc:  # pragma: no cover - network behavior
            last_error = exc
            if attempt < retries:
                time.sleep(2**attempt)
    raise RuntimeError(f"g:Profiler request failed after {retries + 1} attempts: {last_error}")


def _network_metrics(edges: pd.DataFrame, threshold: float) -> dict[str, float]:
    selected = edges.loc[pd.to_numeric(edges["combined_score"], errors="coerce") >= threshold]
    graph = nx.Graph()
    graph.add_edges_from(selected[["node_a", "node_b"]].itertuples(index=False, name=None))
    if graph.number_of_nodes() == 0:
        return {
            "nodes": 0,
            "edges": 0,
            "density": 0.0,
            "transitivity": 0.0,
            "largest_component": 0,
            "mean_degree": 0.0,
        }
    return {
        "nodes": graph.number_of_nodes(),
        "edges": graph.number_of_edges(),
        "density": float(nx.density(graph)),
        "transitivity": float(nx.transitivity(graph)),
        "largest_component": int(max(len(component) for component in nx.connected_components(graph))),
        "mean_degree": float(np.mean([degree for _, degree in graph.degree()])),
    }


def run_t8_background_and_network_sensitivity(
    atlas_path: Path,
    output_dir: Path,
    *,
    portfolio_path: Path | None = None,
    pareto_path: Path | None = None,
    existing_enrichment_path: Path | None = None,
    string_edges_path: Path | None = None,
    controls_per_query: int = 5,
    sources: Sequence[str] = ("GO:BP", "GO:MF", "GO:CC", "REAC", "KEGG", "WP", "CORUM", "HPA"),
    run_online: bool = True,
    timeout: float = 90.0,
    retries: int = 2,
) -> dict[str, pd.DataFrame]:
    output_dir.mkdir(parents=True, exist_ok=True)
    atlas = _read_csv(atlas_path, required=("candidate_id", "species", "repairability_extra_trees"))
    portfolio = _read_csv(portfolio_path) if portfolio_path and portfolio_path.exists() else None
    pareto = _read_csv(pareto_path) if pareto_path and pareto_path.exists() else None
    proteins = _protein_level(atlas)
    query_sets = _query_sets(atlas, portfolio, pareto)
    backgrounds = _background_sets(atlas)

    manifest_rows: list[dict[str, Any]] = []
    match_frames: list[pd.DataFrame] = []
    enrichment_rows: list[dict[str, Any]] = []
    cache_dir = output_dir / "gprofiler_cache"
    organism_map = {"human": "hsapiens", "mouse": "mmusculus"}
    for (species, query_name), query_genes in sorted(query_sets.items()):
        query_genes = {gene for gene in query_genes if gene}
        if len(query_genes) < 2:
            continue
        species_proteins = proteins.loc[proteins["species"].astype(str) == species].copy()
        matched, matches = _matched_background(
            species_proteins,
            query_genes,
            controls_per_query=controls_per_query,
        )
        matches["species"] = species
        matches["query_name"] = query_name
        match_frames.append(matches)
        test_backgrounds = {
            "candidate_universe": backgrounds.get((species, "candidate_universe"), set()),
            "structure_eligible": backgrounds.get((species, "structure_eligible"), set()),
            "in_domain": backgrounds.get((species, "in_domain"), set()),
            "matched": matched,
        }
        for background_name, background_genes in test_backgrounds.items():
            background_genes = set(background_genes) | query_genes
            background_fingerprint = stable_gene_set_fingerprint(background_genes)
            manifest_record = {
                "species": species,
                "query_name": query_name,
                "background_name": background_name,
                "background_fingerprint": background_fingerprint,
                "query_gene_count": len(query_genes),
                "background_gene_count": len(background_genes),
                "query_genes": " | ".join(sorted(query_genes)),
                "background_genes": " | ".join(sorted(background_genes)),
                "api_status": "not_run" if not run_online else "pending",
                "api_error": "",
            }
            if not run_online:
                manifest_rows.append(manifest_record)
                continue
            organism = organism_map.get(species)
            if organism is None:
                manifest_record["api_status"] = "unsupported_species"
                manifest_rows.append(manifest_record)
                continue
            cache_path = cache_dir / f"{species}__{query_name}__{background_name}.json"
            try:
                result, meta = _gprofiler_request(
                    organism,
                    sorted(query_genes),
                    sorted(background_genes),
                    sources=sources,
                    cache_path=cache_path,
                    timeout=timeout,
                    retries=retries,
                )
                manifest_record["api_status"] = "PASS"
            except Exception as exc:  # pragma: no cover - network behavior
                manifest_record["api_status"] = "FAIL"
                manifest_record["api_error"] = str(exc)
                manifest_rows.append(manifest_record)
                continue
            manifest_rows.append(manifest_record)
            for row in result:
                enrichment_rows.append(
                    {
                        "species": species,
                        "query_name": query_name,
                        "background_name": background_name,
                        "background_fingerprint": background_fingerprint,
                        "source": row.get("source"),
                        "term_id": row.get("native"),
                        "term_name": row.get("name"),
                        "p_value_adjusted": row.get("p_value"),
                        "significant": row.get("significant", False),
                        "query_size": row.get("query_size"),
                        "intersection_size": row.get("intersection_size"),
                        "term_size": row.get("term_size"),
                        "effective_domain_size": row.get("effective_domain_size"),
                        "precision": row.get("precision"),
                        "recall": row.get("recall"),
                        "data_version": json.dumps(meta.get("data_version", {}), sort_keys=True),
                    }
                )

    manifest = pd.DataFrame(manifest_rows)
    matches = pd.concat(match_frames, ignore_index=True) if match_frames else pd.DataFrame()
    enrichment = pd.DataFrame(enrichment_rows)
    if not enrichment.empty:
        enrichment["significant"] = enrichment["significant"].fillna(False).astype(bool)

        # Several named backgrounds can resolve to the same gene set (for example,
        # all structure-eligible proteins may equal the candidate universe). Count
        # unique gene-set fingerprints rather than treating labels as independent
        # sensitivity analyses.
        background_group_map: dict[tuple[str, str, str], str] = {}
        for (species_key, query_key, fingerprint), sub in manifest.groupby(
            ["species", "query_name", "background_fingerprint"], dropna=False
        ):
            labels = sorted(set(sub["background_name"].astype(str)))
            background_group_map[
                (str(species_key), str(query_key), str(fingerprint))
            ] = "=".join(labels)

        enrichment["background_group"] = [
            background_group_map.get(
                (str(species_key), str(query_key), str(fingerprint)),
                str(background_name),
            )
            for species_key, query_key, fingerprint, background_name in enrichment[
                ["species", "query_name", "background_fingerprint", "background_name"]
            ].itertuples(index=False, name=None)
        ]

        term_robustness = (
            enrichment.groupby(
                ["species", "query_name", "source", "term_id", "term_name"],
                as_index=False,
                dropna=False,
            )
            .agg(
                background_labels_tested=("background_name", "nunique"),
                distinct_backgrounds_tested=("background_fingerprint", "nunique"),
                significant_background_labels=("significant", "sum"),
                min_adjusted_p=("p_value_adjusted", "min"),
                max_adjusted_p=("p_value_adjusted", "max"),
                background_list=("background_name", lambda values: " | ".join(sorted(set(values)))),
                distinct_background_list=(
                    "background_group",
                    lambda values: " | ".join(sorted(set(values))),
                ),
            )
        )

        significant_unique = (
            enrichment.loc[enrichment["significant"]]
            .groupby(
                ["species", "query_name", "source", "term_id"],
                dropna=False,
            )
            .agg(
                significant_background_list=(
                    "background_name",
                    lambda values: " | ".join(sorted(set(values))),
                ),
                significant_distinct_backgrounds=(
                    "background_fingerprint",
                    "nunique",
                ),
                significant_distinct_background_list=(
                    "background_group",
                    lambda values: " | ".join(sorted(set(values))),
                ),
            )
            .reset_index()
        )
        term_robustness = term_robustness.drop(
            columns=[
                "significant_background_list",
                "significant_distinct_backgrounds",
                "significant_distinct_background_list",
            ],
            errors="ignore",
        ).merge(
            significant_unique,
            on=["species", "query_name", "source", "term_id"],
            how="left",
        )
        term_robustness["significant_distinct_backgrounds"] = (
            term_robustness["significant_distinct_backgrounds"].fillna(0).astype(int)
        )
        term_robustness["includes_matched_background"] = term_robustness[
            "significant_background_list"
        ].fillna("").str.contains(r"(?:^| \| )matched(?:$| \| )", regex=True)
        term_robustness["background_robust"] = (
            (term_robustness["significant_distinct_backgrounds"] >= 2)
            & term_robustness["includes_matched_background"]
        )
        term_robustness["ontology_root_term"] = [
            is_ontology_root_term(source, term_id, term_name)
            for source, term_id, term_name in term_robustness[
                ["source", "term_id", "term_name"]
            ].itertuples(index=False, name=None)
        ]
        term_robustness["background_robust_informative"] = (
            term_robustness["background_robust"]
            & ~term_robustness["ontology_root_term"]
        )
    else:
        term_robustness = pd.DataFrame(
            columns=[
                "species",
                "query_name",
                "source",
                "term_id",
                "term_name",
                "background_labels_tested",
                "distinct_backgrounds_tested",
                "significant_background_labels",
                "significant_distinct_backgrounds",
                "min_adjusted_p",
                "max_adjusted_p",
                "background_list",
                "distinct_background_list",
                "significant_background_list",
                "significant_distinct_background_list",
                "includes_matched_background",
                "background_robust",
                "ontology_root_term",
                "background_robust_informative",
            ]
        )

    network_summary_rows: list[dict[str, Any]] = []
    network_edges_corrected = pd.DataFrame()
    if string_edges_path is not None and string_edges_path.exists():
        edges = _read_csv(
            string_edges_path,
            required=("species", "query_name", "node_a", "node_b", "combined_score"),
        )
        for score_column in ("tscore", "escore", "dscore"):
            if score_column not in edges:
                edges[score_column] = 0.0
        edges["experimental_or_database_score"] = edges[["escore", "dscore"]].max(axis=1)
        edges["text_mining_dominated"] = (
            pd.to_numeric(edges["tscore"], errors="coerce").fillna(0.0) >= 0.70
        ) & (edges["experimental_or_database_score"] < 0.40)
        edges["study_bias_corrected_keep"] = ~edges["text_mining_dominated"]
        network_edges_corrected = edges
        for (species, query_name), sub in edges.groupby(["species", "query_name"]):
            for threshold in (0.70, 0.80, 0.90):
                for variant, variant_edges in (
                    ("all_evidence", sub),
                    ("exclude_text_mining_dominated", sub.loc[sub["study_bias_corrected_keep"]]),
                    (
                        "experimental_or_database_supported",
                        sub.loc[sub["experimental_or_database_score"] >= 0.40],
                    ),
                ):
                    metrics = _network_metrics(variant_edges, threshold)
                    network_summary_rows.append(
                        {
                            "species": species,
                            "query_name": query_name,
                            "variant": variant,
                            "combined_score_threshold": threshold,
                            **metrics,
                            "text_mining_dominated_fraction": float(
                                sub["text_mining_dominated"].mean()
                            ),
                        }
                    )
    network_summary = pd.DataFrame(network_summary_rows)

    if existing_enrichment_path is not None and existing_enrichment_path.exists():
        existing_enrichment = _read_csv(existing_enrichment_path)
        existing_summary = pd.DataFrame(
            [
                {
                    "existing_rows": len(existing_enrichment),
                    "existing_significant_rows": int(
                        existing_enrichment.get(
                            "significant", pd.Series(False, index=existing_enrichment.index)
                        )
                        .fillna(False)
                        .astype(bool)
                        .sum()
                    ),
                    "existing_background_interpretation": "Stage 5 in-domain custom background",
                    "input_sha256": sha256_file(existing_enrichment_path),
                }
            ]
        )
    else:
        existing_summary = pd.DataFrame(
            [
                {
                    "existing_rows": 0,
                    "existing_significant_rows": 0,
                    "existing_background_interpretation": "not available",
                    "input_sha256": "",
                }
            ]
        )

    outputs = {
        "manifest": manifest,
        "matches": matches,
        "enrichment": enrichment,
        "term_robustness": term_robustness,
        "network_edges": network_edges_corrected,
        "network_summary": network_summary,
        "existing_summary": existing_summary,
    }
    filenames = {
        "manifest": "phase2_t8_query_background_manifest.csv",
        "matches": "phase2_t8_matched_background_pairs.csv",
        "enrichment": "phase2_t8_gprofiler_multibackground.csv",
        "term_robustness": "phase2_t8_enrichment_background_robustness.csv",
        "network_edges": "phase2_t8_string_study_bias_edges.csv",
        "network_summary": "phase2_t8_string_study_bias_summary.csv",
        "existing_summary": "phase2_t8_existing_enrichment_summary.csv",
    }
    for key, frame in outputs.items():
        frame.to_csv(output_dir / filenames[key], index=False)
    manifest_json = {
        "analysis": "T8 multi-background enrichment and STRING study-bias sensitivity",
        "online_gprofiler_run": run_online,
        "gprofiler_endpoint": GPROFILER_URL,
        "gprofiler_domain_scope": "custom",
        "gprofiler_correction": "FDR",
        "controls_per_query": controls_per_query,
        "sources": list(sources),
        "background_independence_rule": (
            "background labels with identical gene sets are collapsed by SHA-256 fingerprint"
        ),
        "background_robust_rule": (
            "significant in at least two distinct background gene sets and in the matched background"
        ),
        "ontology_root_terms_excluded_from_informative_count": True,
        "network_correction": (
            "exclude edges dominated by STRING text mining unless experimental/database support is present"
        ),
        "claim_boundary": (
            "functional and network analyses contextualize hypotheses and are not evidence of CMLase access or efficacy"
        ),
        "input_sha256": {"atlas": sha256_file(atlas_path)},
    }
    for name, path in (
        ("portfolio", portfolio_path),
        ("pareto", pareto_path),
        ("existing_enrichment", existing_enrichment_path),
        ("string_edges", string_edges_path),
    ):
        if path is not None and path.exists():
            manifest_json["input_sha256"][name] = sha256_file(path)
    (output_dir / "phase2_t8_manifest.json").write_text(
        json.dumps(manifest_json, indent=2, sort_keys=True), encoding="utf-8"
    )
    return outputs


def run_phase2_doctor(root: Path, output_dir: Path) -> dict[str, Any]:
    paths = discover_paths(root)
    required = {
        "labels": paths.labels,
        "structure": paths.structure,
        "predictions": paths.predictions,
        "atlas": paths.atlas,
    }
    optional = {
        "portfolio": paths.portfolio,
        "pareto": paths.pareto,
        "tiers": paths.tiers,
        "enrichment": paths.enrichment,
        "string_edges": paths.string_edges,
    }
    status = "PASS" if all(path.exists() for path in required.values()) else "FAIL"
    report: dict[str, Any] = {
        "status": status,
        "required": {name: {"path": str(path), "exists": path.exists()} for name, path in required.items()},
        "optional": {
            name: {"path": str(path) if path is not None else None, "exists": bool(path and path.exists())}
            for name, path in optional.items()
        },
        "transfer_safe_features": list(TRANSFER_SAFE_FEATURES),
        "claim_boundary": (
            "Phase II adds uncertainty, calibration, compartment and functional sensitivity analyses; it does not add independent site-level labels"
        ),
    }
    if paths.predictions.exists():
        predictions = pd.read_csv(paths.predictions)
        report["n_oof_rows"] = len(predictions)
        report["n_oof_sites"] = int(predictions.get("site_uniprot", pd.Series(dtype=object)).nunique())
    if paths.atlas.exists():
        atlas = pd.read_csv(paths.atlas, low_memory=False)
        report["n_atlas_rows"] = len(atlas)
        report["n_atlas_species"] = int(atlas.get("species", pd.Series(dtype=object)).nunique())
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "phase2_doctor.json").write_text(
        json.dumps(report, indent=2, sort_keys=True), encoding="utf-8"
    )
    return report


def write_phase2_report(output_dir: Path) -> Path:
    sections = [
        "# CML-RepairNet Revision Phase II report",
        "",
        "## Claim boundary",
        "These analyses refine uncertainty, calibration, biological accessibility and functional robustness. They do not create an independent non-BSA site-resolved validation set.",
        "",
    ]
    t5 = output_dir / "phase2_t5_bsa_coverage.csv"
    if t5.exists():
        frame = pd.read_csv(t5)
        sections.extend(["## T5 — Covariate-shift interval sensitivity", "", frame.to_markdown(index=False), ""])
    t6 = output_dir / "phase2_t6_boundary_summary.csv"
    if t6.exists():
        frame = pd.read_csv(t6)
        sections.extend(["## T6 — Boundary shrinkage", "", frame.to_markdown(index=False), ""])
    t7 = output_dir / "phase2_t7_biomedical_interpretation_summary.csv"
    if t7.exists():
        frame = pd.read_csv(t7)
        sections.extend(["## T7 — Compartment/accessibility layers", "", frame.to_markdown(index=False), ""])
    t8 = output_dir / "phase2_t8_enrichment_background_robustness.csv"
    if t8.exists():
        frame = pd.read_csv(t8)
        robust = frame.loc[frame.get("background_robust", False).fillna(False)] if not frame.empty else frame
        informative = (
            frame.loc[frame["background_robust_informative"].fillna(False)]
            if not frame.empty and "background_robust_informative" in frame
            else robust
        )
        sections.extend(
            [
                "## T8 — Multi-background enrichment and network study-bias sensitivity",
                "",
                f"Background-robust enrichment terms (including ontology roots): {len(robust)}",
                f"Background-robust biologically informative terms: {len(informative)}",
                "",
            ]
        )
    path = output_dir / "revision_phase2_report.md"
    path.write_text("\n".join(sections), encoding="utf-8")
    return path
