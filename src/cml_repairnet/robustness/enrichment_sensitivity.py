from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import hypergeom


def _gene_set(value: object) -> frozenset[str]:
    if value is None or pd.isna(value):
        return frozenset()
    return frozenset(part.strip() for part in str(value).split("|") if part.strip())


def _redundancy_groups(frame: pd.DataFrame, threshold: float) -> list[int]:
    groups = [-1] * len(frame)
    next_group = 1
    sets = frame["_genes"].tolist()
    for i in range(len(frame)):
        if groups[i] != -1:
            continue
        groups[i] = next_group
        stack = [i]
        while stack:
            current = stack.pop()
            for j in range(len(frame)):
                if groups[j] != -1:
                    continue
                union = sets[current] | sets[j]
                jaccard = len(sets[current] & sets[j]) / len(union) if union else 0.0
                if jaccard >= threshold:
                    groups[j] = next_group
                    stack.append(j)
        next_group += 1
    return groups


def run_enrichment_sensitivity(enrichment: pd.DataFrame, config: dict) -> pd.DataFrame:
    required = {
        "species",
        "query_name",
        "source",
        "term_id",
        "term_name",
        "p_value_adjusted",
        "significant",
        "query_size_mapped",
        "intersection_size",
        "term_size",
        "effective_domain_size",
        "intersection_genes",
    }
    missing = sorted(required - set(enrichment.columns))
    if missing:
        raise ValueError(f"Enrichment table missing columns: {missing}")
    frame = enrichment.copy()
    if bool(config.get("significant_only", True)):
        frame = frame.loc[frame["significant"].fillna(False).astype(bool)].copy()
    frame["_genes"] = frame["intersection_genes"].map(_gene_set)
    raw_threshold = float(config.get("raw_p_threshold", 0.05))
    rows = []
    for row in frame.itertuples(index=False):
        N = int(row.effective_domain_size)
        M = int(row.term_size)
        n = int(row.query_size_mapped)
        k = int(row.intersection_size)
        baseline_raw_p = (
            float(hypergeom.sf(k - 1, N, M, n)) if all(v >= 0 for v in (N, M, n, k)) else np.nan
        )
        if k > 0 and n > 1:
            loo_raw_p = float(hypergeom.sf(k - 2, N, M, n - 1))
        else:
            loo_raw_p = np.nan
        if k <= 1:
            sensitivity = "single_driver"
        elif not np.isnan(loo_raw_p) and loo_raw_p >= raw_threshold:
            sensitivity = "leave_one_driver_fragile"
        elif k <= 3:
            sensitivity = "small_intersection"
        else:
            sensitivity = "driver_robust"
        rows.append(
            {
                "species": row.species,
                "query_name": row.query_name,
                "source": row.source,
                "term_id": row.term_id,
                "term_name": row.term_name,
                "adjusted_p_value_stage5": float(row.p_value_adjusted),
                "query_size": n,
                "intersection_size": k,
                "term_size": M,
                "effective_domain_size": N,
                "intersection_genes": row.intersection_genes,
                "baseline_hypergeom_raw_p": baseline_raw_p,
                "worst_intersection_gene_removal_raw_p": loo_raw_p,
                "remains_raw_significant_after_any_single_intersection_removal": bool(
                    not np.isnan(loo_raw_p) and loo_raw_p < raw_threshold
                ),
                "max_single_gene_share": 1.0 / k if k else np.nan,
                "driver_sensitivity_class": sensitivity,
            }
        )
    result = pd.DataFrame(rows)
    if result.empty:
        return result
    result["_genes"] = result["intersection_genes"].map(_gene_set)
    threshold = float(config.get("redundancy_jaccard", 0.70))
    group_values = np.zeros(len(result), dtype=int)
    group_offset = 0
    for _, indices in result.groupby(["species", "query_name", "source"]).groups.items():
        sub = result.loc[list(indices)].copy().reset_index()
        groups = _redundancy_groups(sub, threshold)
        mapping = {original: group + group_offset for original, group in zip(sub["index"], groups)}
        for original, group in mapping.items():
            group_values[result.index.get_loc(original)] = group
        group_offset += max(groups, default=0)
    result["redundancy_group"] = group_values
    group_sizes = result.groupby("redundancy_group")["term_id"].transform("size")
    result["redundancy_group_size"] = group_sizes
    return (
        result.drop(columns=["_genes"])
        .sort_values(["driver_sensitivity_class", "adjusted_p_value_stage5"])
        .reset_index(drop=True)
    )
