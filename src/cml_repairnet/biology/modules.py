from __future__ import annotations

import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
from sklearn.impute import SimpleImputer
from sklearn.metrics import silhouette_score
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

MODULE_FEATURES = [
    "repairability_extra_trees",
    "evidence_score_raw",
    "transfer_reliability_score",
    "model_disagreement",
    "site_plddt",
    "relative_sasa_lysine",
    "neighbor_cysteine_count_8a",
    "neighbor_hydrophobic_fraction_8a",
]


def choose_k(
    values: np.ndarray, *, min_k: int = 2, max_k: int = 6, random_state: int = 42
) -> tuple[int, pd.DataFrame]:
    rows = []
    upper = min(max_k, len(values) - 1)
    if upper < min_k:
        return 1, pd.DataFrame([{"k": 1, "silhouette": np.nan}])
    best_k, best_score = min_k, -np.inf
    for k in range(min_k, upper + 1):
        labels = KMeans(n_clusters=k, random_state=random_state, n_init=50).fit_predict(values)
        score = silhouette_score(values, labels)
        rows.append({"k": k, "silhouette": score})
        if score > best_score:
            best_k, best_score = k, score
    return best_k, pd.DataFrame(rows)


def _module_label(profile: pd.Series, global_median: pd.Series) -> str:
    high_repair = profile["repairability_extra_trees"] >= global_median["repairability_extra_trees"]
    high_evidence = profile["evidence_score_raw"] >= global_median["evidence_score_raw"]
    high_reliability = (
        profile["transfer_reliability_score"] >= global_median["transfer_reliability_score"]
    )
    if high_repair and high_evidence and high_reliability:
        return "high-repair evidence-rich reliable"
    if high_repair and high_reliability:
        return "high-repair reliable"
    if high_evidence:
        return "evidence-rich moderate-repair"
    if not high_reliability:
        return "transfer-sensitive"
    return "moderate-priority structural"


def build_candidate_modules(
    frame: pd.DataFrame,
    *,
    random_state: int = 42,
    max_k: int = 6,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    missing = [column for column in MODULE_FEATURES if column not in frame]
    if missing:
        raise ValueError(f"Module features are missing: {missing}")
    matrix = frame[MODULE_FEATURES].apply(pd.to_numeric, errors="coerce")
    transform = make_pipeline(SimpleImputer(strategy="median"), StandardScaler())
    values = transform.fit_transform(matrix)
    best_k, selection = choose_k(values, max_k=max_k, random_state=random_state)
    labels = (
        np.zeros(len(frame), dtype=int)
        if best_k == 1
        else KMeans(n_clusters=best_k, random_state=random_state, n_init=100).fit_predict(values)
    )
    assignments = frame.copy()
    assignments["module_id"] = labels + 1
    profiles = assignments.groupby("module_id")[MODULE_FEATURES].mean()
    global_median = assignments[MODULE_FEATURES].median()
    profiles["module_label"] = profiles.apply(_module_label, axis=1, global_median=global_median)
    profiles["sites"] = assignments.groupby("module_id").size()
    profiles["proteins"] = assignments.groupby("module_id")["uniprot_accession"].nunique()
    profiles["human_sites"] = assignments.groupby("module_id")["species"].apply(
        lambda x: int((x == "human").sum())
    )
    profiles["mouse_sites"] = assignments.groupby("module_id")["species"].apply(
        lambda x: int((x == "mouse").sum())
    )
    profiles = profiles.reset_index()
    assignments = assignments.merge(
        profiles[["module_id", "module_label"]], on="module_id", how="left"
    )
    return assignments, profiles, selection
