from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path

import numpy as np
import pandas as pd

DEFAULT_MAXIMIZE = (
    "repairability_extra_trees",
    "evidence_score_norm",
    "transfer_reliability_score",
)
DEFAULT_MINIMIZE: tuple[str, ...] = ()


def _minmax(series: pd.Series, *, invert: bool = False) -> pd.Series:
    numeric = pd.to_numeric(series, errors="coerce")
    finite = numeric[np.isfinite(numeric)]
    if finite.empty:
        result = pd.Series(np.zeros(len(series)), index=series.index, dtype=float)
    else:
        low, high = float(finite.min()), float(finite.max())
        if high <= low:
            result = pd.Series(np.ones(len(series)), index=series.index, dtype=float)
        else:
            result = (numeric.fillna(low) - low) / (high - low)
    return 1.0 - result if invert else result


def add_priority_components(frame: pd.DataFrame) -> pd.DataFrame:
    result = frame.copy()
    result["evidence_score_norm"] = result.groupby("species")["evidence_score_raw"].transform(
        _minmax
    )
    result["repairability_score_norm"] = result.groupby("species")[
        "repairability_extra_trees"
    ].transform(_minmax)
    result["ood_safety_score"] = 1.0 - result["ood_percentile"].clip(0, 1)
    max_width = result["prediction_interval_width"].max()
    result["interval_precision_score"] = 1.0 - (
        result["prediction_interval_width"] / max_width if max_width > 0 else 0.0
    )
    max_disagreement = result["model_disagreement"].max()
    result["agreement_score"] = 1.0 - (
        result["model_disagreement"] / max_disagreement if max_disagreement > 0 else 0.0
    )
    # One composite transfer-reliability objective prevents high-dimensional
    # Pareto-front inflation while retaining OOD, uncertainty, and agreement.
    result["transfer_reliability_score"] = (
        0.50 * result["ood_safety_score"]
        + 0.25 * result["interval_precision_score"]
        + 0.25 * result["agreement_score"]
    )
    result["weighted_priority_score"] = (
        0.30 * result["repairability_score_norm"]
        + 0.20 * result["consensus_rank_percentile"]
        + 0.15 * result["evidence_score_norm"]
        + 0.15 * result["structure_quality_score"]
        + 0.10 * result["ood_safety_score"]
        + 0.05 * result["interval_precision_score"]
        + 0.05 * result["agreement_score"]
    )
    return result


def _dominates(a: np.ndarray, b: np.ndarray) -> bool:
    return bool(np.all(a >= b) and np.any(a > b))


def non_dominated_sort(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    n = len(values)
    dominates: list[list[int]] = [[] for _ in range(n)]
    dominated_count = np.zeros(n, dtype=int)
    fronts: list[list[int]] = [[]]
    for i in range(n):
        for j in range(i + 1, n):
            if _dominates(values[i], values[j]):
                dominates[i].append(j)
                dominated_count[j] += 1
            elif _dominates(values[j], values[i]):
                dominates[j].append(i)
                dominated_count[i] += 1
        if dominated_count[i] == 0:
            fronts[0].append(i)
    ranks = np.full(n, -1, dtype=int)
    front_index = 0
    while front_index < len(fronts) and fronts[front_index]:
        next_front: list[int] = []
        for i in fronts[front_index]:
            ranks[i] = front_index + 1
            for j in dominates[i]:
                dominated_count[j] -= 1
                if dominated_count[j] == 0:
                    next_front.append(j)
        if next_front:
            fronts.append(next_front)
        front_index += 1
    return ranks


def crowding_distance(values: np.ndarray, ranks: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    distances = np.zeros(len(values), dtype=float)
    for rank in np.unique(ranks):
        indices = np.where(ranks == rank)[0]
        if len(indices) <= 2:
            distances[indices] = np.inf
            continue
        front = values[indices]
        for column in range(front.shape[1]):
            order = np.argsort(front[:, column])
            sorted_indices = indices[order]
            distances[sorted_indices[0]] = np.inf
            distances[sorted_indices[-1]] = np.inf
            low = front[order[0], column]
            high = front[order[-1], column]
            if high <= low:
                continue
            for position in range(1, len(order) - 1):
                previous_value = front[order[position - 1], column]
                next_value = front[order[position + 1], column]
                distances[sorted_indices[position]] += (next_value - previous_value) / (high - low)
    return distances


def add_pareto_ranking(
    frame: pd.DataFrame,
    *,
    maximize: Sequence[str] = DEFAULT_MAXIMIZE,
    minimize: Sequence[str] = DEFAULT_MINIMIZE,
    eligible_statuses: Sequence[str] = ("in_domain",),
) -> pd.DataFrame:
    result = frame.copy()
    missing = [column for column in (*maximize, *minimize) if column not in result.columns]
    if missing:
        raise ValueError(f"Pareto columns are missing: {missing}")
    if "applicability_status" not in result.columns:
        raise ValueError("Pareto ranking requires applicability_status")

    # Rank only candidates inside the prespecified applicability domain.
    # Borderline/OOD candidates remain in the atlas and a separate watchlist,
    # but cannot be promoted to the primary Pareto front.
    result["pareto_front"] = 0
    result["crowding_distance"] = np.nan
    eligible_mask = result["applicability_status"].isin(list(eligible_statuses))
    eligible = result.loc[eligible_mask]
    for species, group in eligible.groupby("species", sort=False):
        objective_parts = [_minmax(group[column]).to_numpy(dtype=float) for column in maximize] + [
            _minmax(group[column], invert=True).to_numpy(dtype=float) for column in minimize
        ]
        values = np.column_stack(objective_parts)
        ranks = non_dominated_sort(values)
        crowding = crowding_distance(values, ranks)
        result.loc[group.index, "pareto_front"] = ranks
        result.loc[group.index, "crowding_distance"] = crowding
    result["pareto_front"] = result["pareto_front"].astype(int)
    return result


def select_diverse_portfolio(
    frame: pd.DataFrame,
    *,
    top_n: int = 100,
    max_sites_per_protein: int = 2,
    include_borderline: bool = False,
) -> pd.DataFrame:
    allowed = ["in_domain"] + (["borderline"] if include_borderline else [])
    eligible = frame[frame["applicability_status"].isin(allowed)].copy()
    eligible = eligible.sort_values(
        ["pareto_front", "weighted_priority_score", "crowding_distance"],
        ascending=[True, False, False],
    )
    selected: list[int] = []
    protein_counts: dict[tuple[str, str], int] = {}
    species_counts = {species: 0 for species in eligible["species"].dropna().unique()}
    species_target = max(1, top_n // max(1, len(species_counts)))

    # First pass keeps approximate species balance.
    for index, row in eligible.iterrows():
        species = str(row["species"])
        protein_key = (species, str(row["uniprot_accession"]))
        if protein_counts.get(protein_key, 0) >= max_sites_per_protein:
            continue
        if species_counts.get(species, 0) >= species_target:
            continue
        selected.append(index)
        protein_counts[protein_key] = protein_counts.get(protein_key, 0) + 1
        species_counts[species] = species_counts.get(species, 0) + 1
        if len(selected) >= top_n:
            break

    # Second pass fills remaining slots without a species quota.
    if len(selected) < top_n:
        for index, row in eligible.iterrows():
            if index in selected:
                continue
            protein_key = (str(row["species"]), str(row["uniprot_accession"]))
            if protein_counts.get(protein_key, 0) >= max_sites_per_protein:
                continue
            selected.append(index)
            protein_counts[protein_key] = protein_counts.get(protein_key, 0) + 1
            if len(selected) >= top_n:
                break
    portfolio = eligible.loc[selected].copy()
    portfolio["portfolio_rank"] = np.arange(1, len(portfolio) + 1)
    return portfolio


def prioritize_atlas(
    scored_csv: Path,
    output_dir: Path,
    *,
    top_n: int = 100,
    max_sites_per_protein: int = 2,
) -> dict[str, pd.DataFrame]:
    scored = pd.read_csv(scored_csv)
    ranked = add_priority_components(scored)
    ranked = add_pareto_ranking(ranked)
    front = ranked[ranked["pareto_front"] == 1].copy()
    in_domain = ranked[ranked["applicability_status"] == "in_domain"].copy()
    watchlist = ranked[ranked["applicability_status"] != "in_domain"].copy()
    portfolio = select_diverse_portfolio(
        ranked,
        top_n=top_n,
        max_sites_per_protein=max_sites_per_protein,
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    ranked.sort_values(["pareto_front", "weighted_priority_score"], ascending=[True, False]).to_csv(
        output_dir / "cml_repairability_atlas_all_sites.csv", index=False
    )
    in_domain.sort_values("weighted_priority_score", ascending=False).to_csv(
        output_dir / "cml_repairability_atlas_in_domain.csv", index=False
    )
    front.sort_values("weighted_priority_score", ascending=False).to_csv(
        output_dir / "cml_repairability_pareto_front.csv", index=False
    )
    watchlist.sort_values(
        ["applicability_status", "weighted_priority_score"], ascending=[True, False]
    ).to_csv(output_dir / "cml_repairability_domain_shift_watchlist.csv", index=False)
    portfolio.to_csv(output_dir / "cml_repairability_priority_portfolio.csv", index=False)
    return {
        "all": ranked,
        "in_domain": in_domain,
        "pareto": front,
        "watchlist": watchlist,
        "portfolio": portfolio,
    }
