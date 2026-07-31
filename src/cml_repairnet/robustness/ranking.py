from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import rankdata, spearmanr

from .schema import add_labels


MODEL_COLUMNS = ["repairability_extra_trees", "repairability_lightgbm", "repairability_ridge"]
OBJECTIVE_COLUMNS = ["repairability_score_norm", "evidence_score_norm", "transfer_reliability_score"]


def percentile(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    return (rankdata(values, method="average") - 1.0) / max(len(values) - 1, 1)


def minmax(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    finite = values[np.isfinite(values)]
    if finite.size == 0:
        return np.zeros_like(values, dtype=float)
    low = float(finite.min())
    high = float(finite.max())
    if high <= low:
        return np.ones_like(values, dtype=float)
    return (np.nan_to_num(values, nan=low) - low) / (high - low)


def pareto_mask(values: np.ndarray) -> np.ndarray:
    """Return a maximize-all first-front mask without dense quadratic matrices."""
    values = np.asarray(values, dtype=float)
    if values.ndim != 2:
        raise ValueError("values must be a 2D array")
    n = len(values)
    result = np.zeros(n, dtype=bool)
    order = np.argsort(-values[:, 0], kind="mergesort")
    skyline: list[int] = []
    for idx in order:
        point = values[idx]
        dominated = False
        for kept in skyline:
            other = values[kept]
            if np.all(other >= point) and np.any(other > point):
                dominated = True
                break
        if not dominated:
            result[idx] = True
            skyline = [
                kept
                for kept in skyline
                if not (np.all(point >= values[kept]) and np.any(point > values[kept]))
            ]
            skyline.append(idx)
    return result


def specieswise_minmax(values: np.ndarray, species: np.ndarray) -> np.ndarray:
    result = np.zeros(len(values), dtype=float)
    for label in pd.unique(species):
        indices = np.flatnonzero(species == label)
        result[indices] = minmax(values[indices])
    return result


def specieswise_pareto_mask(values: np.ndarray, species: np.ndarray) -> np.ndarray:
    """Compute the first Pareto front independently within each species.

    Stage 4 computes Pareto fronts separately for human and mouse candidates.
    This helper preserves that applicability-universe design in Stage 6.
    """
    result = np.zeros(len(values), dtype=bool)
    for label in pd.unique(species):
        indices = np.flatnonzero(species == label)
        if len(indices):
            result[indices] = pareto_mask(values[indices])
    return result


def nondominated_front_ranks(values: np.ndarray) -> np.ndarray:
    """Assign complete maximize-all nondominated front ranks (1, 2, ...).

    Points are processed in descending lexicographic objective order. Every
    possible dominator therefore appears earlier, and a point's front is one
    plus the maximum front among its strict dominators. This exactly reproduces
    the frozen Stage 4 front labels without repeated dense matrix allocations.
    """
    values = np.asarray(values, dtype=float)
    if values.ndim != 2:
        raise ValueError("values must be a 2D array")
    if values.shape[1] != 3:
        raise ValueError("Stage 6 nondominated sorting expects three objectives")
    n = len(values)
    ranks = np.ones(n, dtype=int)
    order = np.lexsort((-values[:, 2], -values[:, 1], -values[:, 0]))
    for position, index in enumerate(order):
        if position == 0:
            continue
        previous = order[:position]
        previous_values = values[previous]
        point = values[index]
        dominates = (
            (previous_values[:, 0] >= point[0])
            & (previous_values[:, 1] >= point[1])
            & (previous_values[:, 2] >= point[2])
            & (
                (previous_values[:, 0] > point[0])
                | (previous_values[:, 1] > point[1])
                | (previous_values[:, 2] > point[2])
            )
        )
        if dominates.any():
            ranks[index] = 1 + int(ranks[previous[dominates]].max())
    return ranks


def specieswise_nondominated_front_ranks(values: np.ndarray, species: np.ndarray) -> np.ndarray:
    """Assign complete nondominated front ranks independently by species."""
    result = np.zeros(len(values), dtype=int)
    for label in pd.unique(species):
        indices = np.flatnonzero(species == label)
        if len(indices):
            result[indices] = nondominated_front_ranks(values[indices])
    return result


def select_diversity_constrained_indices(
    priority: np.ndarray,
    front_ranks: np.ndarray,
    species: np.ndarray,
    proteins: np.ndarray,
    *,
    top_n: int = 100,
    max_sites_per_protein: int = 2,
) -> np.ndarray:
    """Select a Stage-4-compatible diversity-constrained portfolio.

    Complete species-specific nondominated front ranks are placed first, then
    priority is sorted descending inside each front. The two-pass species quota
    and per-protein cap reproduce the frozen Stage 4 portfolio selection rule.
    """
    priority = np.asarray(priority, dtype=float)
    front_ranks = np.asarray(front_ranks)
    if front_ranks.dtype == bool:
        # Backward-compatible interpretation for callers that provide only a
        # first-front mask. Production Stage 6 uses complete front ranks.
        front_ranks = np.where(front_ranks, 1, 2)
    front_ranks = front_ranks.astype(int)
    species = np.asarray(species, dtype=str)
    proteins = np.asarray(proteins, dtype=str)
    if not (len(priority) == len(front_ranks) == len(species) == len(proteins)):
        raise ValueError("Portfolio selection arrays must have equal length")
    if np.any(front_ranks < 1):
        raise ValueError("Nondominated front ranks must start at 1")

    # Exact Stage 4 ordering: complete front rank first, then weighted
    # priority descending. Species quotas and the per-protein cap are applied
    # during selection.
    order = np.lexsort((-priority, front_ranks))
    unique_species = [str(value) for value in pd.unique(species)]
    species_target = max(1, int(top_n) // max(1, len(unique_species)))
    species_counts = {label: 0 for label in unique_species}
    protein_counts: dict[tuple[str, str], int] = {}
    selected: list[int] = []

    for idx in order:
        label = str(species[idx])
        key = (label, str(proteins[idx]))
        if protein_counts.get(key, 0) >= int(max_sites_per_protein):
            continue
        if species_counts.get(label, 0) >= species_target:
            continue
        selected.append(int(idx))
        species_counts[label] = species_counts.get(label, 0) + 1
        protein_counts[key] = protein_counts.get(key, 0) + 1
        if len(selected) >= int(top_n):
            return np.asarray(selected, dtype=int)

    if len(selected) < int(top_n):
        selected_set = set(selected)
        for idx in order:
            idx = int(idx)
            if idx in selected_set:
                continue
            label = str(species[idx])
            key = (label, str(proteins[idx]))
            if protein_counts.get(key, 0) >= int(max_sites_per_protein):
                continue
            selected.append(idx)
            selected_set.add(idx)
            protein_counts[key] = protein_counts.get(key, 0) + 1
            if len(selected) >= int(top_n):
                break
    return np.asarray(selected, dtype=int)


def _dirichlet(rng: np.random.Generator, center: list[float], concentration: float) -> np.ndarray:
    center_array = np.asarray(center, dtype=float)
    center_array = center_array / center_array.sum()
    return rng.dirichlet(np.maximum(center_array * concentration, 1e-6))


def run_rank_stability(
    in_domain: pd.DataFrame,
    portfolio: pd.DataFrame,
    pareto: pd.DataFrame,
    config: dict,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    required = set(
        MODEL_COLUMNS
        + OBJECTIVE_COLUMNS
        + [
            "candidate_id",
            "species",
            "uniprot_accession",
            "prediction_interval_width",
            "weighted_priority_score",
        ]
    )
    missing = sorted(required - set(in_domain.columns))
    if missing:
        raise ValueError(f"In-domain atlas missing rank-stability columns: {missing}")

    frame = add_labels(in_domain).reset_index(drop=True)
    n = len(frame)
    iterations = int(config["iterations"])
    pareto_iterations = int(config.get("pareto_iterations", min(iterations, 500)))
    seed = int(config.get("random_seed", 20260731))
    rng = np.random.default_rng(seed)
    portfolio_top_n = int(config.get("portfolio_top_n", len(portfolio) if len(portfolio) else 100))
    max_sites_per_protein = int(config.get("max_sites_per_protein", 2))

    model_matrix = frame[MODEL_COLUMNS].apply(pd.to_numeric, errors="coerce").to_numpy(float)
    if np.isnan(model_matrix).any():
        raise ValueError("Model predictions contain missing values")
    evidence = pd.to_numeric(frame["evidence_score_norm"], errors="coerce").fillna(0).to_numpy(float)
    reliability = pd.to_numeric(frame["transfer_reliability_score"], errors="coerce").fillna(0).to_numpy(float)
    interval = pd.to_numeric(frame["prediction_interval_width"], errors="coerce").fillna(0).to_numpy(float)
    species = frame["species"].astype(str).to_numpy()
    proteins = frame["uniprot_accession"].astype(str).to_numpy()

    rank_sum = np.zeros(n)
    rank_sq_sum = np.zeros(n)
    rank_samples = np.empty((iterations, n), dtype=np.float32)
    top_ks = sorted({int(value) for value in config.get("top_k", [10, 20, 50, 100])})
    top_counts = {k: np.zeros(n, dtype=int) for k in top_ks}
    portfolio_counts = np.zeros(n, dtype=int)
    pareto_counts = np.zeros(n, dtype=int)

    for iteration in range(iterations):
        model_weights = _dirichlet(
            rng,
            config["model_weight_center"],
            float(config["model_weight_concentration"]),
        )
        objective_weights = _dirichlet(
            rng,
            config["objective_weight_center"],
            float(config["objective_weight_concentration"]),
        )
        repairability = model_matrix @ model_weights
        stress_sd = interval * float(config.get("interval_noise_fraction", 0.20))
        repairability = np.clip(repairability + rng.normal(0.0, stress_sd), 0.0, 1.0)

        # Stage 4 normalizes repairability inside each species. Keep that design
        # for stressed objectives; global percentile ranks are retained only for
        # the unconstrained ranking diagnostic.
        repairability_norm = specieswise_minmax(repairability, species)
        global_repairability_rank = percentile(repairability)
        evidence_stress = np.clip(
            evidence + rng.normal(0.0, float(config.get("evidence_jitter_sd", 0.03)), n),
            0.0,
            1.0,
        )
        reliability_stress = np.clip(
            reliability + rng.normal(0.0, float(config.get("reliability_jitter_sd", 0.03)), n),
            0.0,
            1.0,
        )
        priority = (
            objective_weights[0] * repairability_norm
            + objective_weights[1] * evidence_stress
            + objective_weights[2] * reliability_stress
        )

        # Unconstrained global rank is a useful diagnostic but is not called
        # portfolio retention.
        global_priority = (
            objective_weights[0] * global_repairability_rank
            + objective_weights[1] * evidence_stress
            + objective_weights[2] * reliability_stress
        )
        order = np.argsort(-global_priority, kind="mergesort")
        ranks = np.empty(n, dtype=float)
        ranks[order] = np.arange(1, n + 1)
        rank_samples[iteration] = ranks
        rank_sum += ranks
        rank_sq_sum += ranks * ranks
        for k in top_ks:
            top_counts[k][order[: min(k, n)]] += 1

        objectives = np.column_stack([repairability_norm, evidence_stress, reliability_stress])
        front_ranks = specieswise_nondominated_front_ranks(objectives, species)
        first_front = front_ranks == 1
        if iteration < pareto_iterations:
            pareto_counts += first_front

        selected = select_diversity_constrained_indices(
            priority,
            front_ranks,
            species,
            proteins,
            top_n=portfolio_top_n,
            max_sites_per_protein=max_sites_per_protein,
        )
        portfolio_counts[selected] += 1

    mean_rank = rank_sum / iterations
    rank_sd = np.sqrt(np.maximum(rank_sq_sum / iterations - mean_rank**2, 0.0))
    result = frame[
        [
            "candidate_id",
            "species",
            "uniprot_accession",
            "site_uniprot",
            "gene_symbol",
            "site_label",
            "weighted_priority_score",
            "repairability_extra_trees",
            "evidence_score_raw",
            "transfer_reliability_score",
            "model_disagreement",
            "prediction_interval_width",
        ]
    ].copy()
    baseline_order = np.argsort(
        -pd.to_numeric(frame["weighted_priority_score"], errors="coerce").to_numpy(float)
    )
    baseline_rank = np.empty(n, dtype=int)
    baseline_rank[baseline_order] = np.arange(1, n + 1)
    result["baseline_rank"] = baseline_rank
    result["mean_stress_rank"] = mean_rank
    result["median_stress_rank"] = np.median(rank_samples, axis=0)
    result["rank_q25"] = np.quantile(rank_samples, 0.25, axis=0)
    result["rank_q75"] = np.quantile(rank_samples, 0.75, axis=0)
    result["rank_iqr"] = result["rank_q75"] - result["rank_q25"]
    result["rank_sd"] = rank_sd
    for k in top_ks:
        result[f"top{k}_selection_frequency"] = top_counts[k] / iterations
    result["portfolio_selection_frequency"] = portfolio_counts / iterations
    result["pareto_membership_frequency"] = pareto_counts / max(pareto_iterations, 1)
    result["baseline_portfolio_member"] = result["candidate_id"].isin(
        set(portfolio["candidate_id"])
    ).astype(int)
    result["baseline_pareto_member"] = result["candidate_id"].isin(
        set(pareto["candidate_id"])
    ).astype(int)
    result = result.sort_values(["mean_stress_rank", "baseline_rank"]).reset_index(drop=True)

    baseline_priority = pd.to_numeric(
        frame["weighted_priority_score"], errors="coerce"
    ).to_numpy(float)
    baseline_portfolio = set(portfolio["candidate_id"].astype(str))
    scenario_rows = []
    fixed_repairability = pd.to_numeric(
        frame["repairability_score_norm"], errors="coerce"
    ).fillna(0).to_numpy(float)
    fixed_objectives = np.column_stack([fixed_repairability, evidence, reliability])
    fixed_front_ranks = specieswise_nondominated_front_ranks(fixed_objectives, species)
    baseline_selected = select_diversity_constrained_indices(
        baseline_priority,
        fixed_front_ranks,
        species,
        proteins,
        top_n=portfolio_top_n,
        max_sites_per_protein=max_sites_per_protein,
    )
    reproduced_portfolio = set(frame.iloc[baseline_selected]["candidate_id"].astype(str))
    baseline_recovery_count = len(baseline_portfolio & reproduced_portfolio)
    baseline_recovery_fraction = baseline_recovery_count / max(len(baseline_portfolio), 1)
    for scenario in config.get("weight_scenarios", []):
        weights = np.array(
            [
                float(scenario["repairability"]),
                float(scenario["evidence"]),
                float(scenario["reliability"]),
            ]
        )
        weights /= weights.sum()
        priority = (
            weights[0] * fixed_repairability
            + weights[1] * evidence
            + weights[2] * reliability
        )
        selected = select_diversity_constrained_indices(
            priority,
            fixed_front_ranks,
            species,
            proteins,
            top_n=portfolio_top_n,
            max_sites_per_protein=max_sites_per_protein,
        )
        scenario_portfolio = set(frame.iloc[selected]["candidate_id"].astype(str))
        rho = spearmanr(baseline_priority, priority).statistic
        intersection = len(baseline_portfolio & scenario_portfolio)
        union = len(baseline_portfolio | scenario_portfolio)
        scenario_rows.append(
            {
                "scenario": scenario["name"],
                "repairability_weight": weights[0],
                "evidence_weight": weights[1],
                "reliability_weight": weights[2],
                "spearman_vs_stage4_priority": rho,
                "portfolio_overlap_count": intersection,
                "portfolio_jaccard": intersection / max(union, 1),
                "baseline_reproduction_count": baseline_recovery_count,
                "baseline_reproduction_fraction": baseline_recovery_fraction,
            }
        )
    return result, pd.DataFrame(scenario_rows)
