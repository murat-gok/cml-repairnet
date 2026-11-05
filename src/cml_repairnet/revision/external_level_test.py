"""Level-resolved external bulk-assay consistency test (revision item A2).

Rationale
---------
``external_scaffold.observed_summary`` collapses the Trabosh et al. (2026)
bulk CML-ELISA panel to one median per scaffold, reducing the observed vector
to n=3. A Spearman correlation on three points can attain at best p = 1/6,
so the published "stress test" was underpowered by construction.

This module is **additive**: it does not modify the frozen
``external_scaffold`` pipeline, its outputs, or its tests. It re-reads the
same source table at full resolution (3 scaffolds x 3 modification levels =
9 observations, each with an SD and replicate count) and separates two
hypotheses that the collapsed analysis conflated:

H-OBS  The *observed* scaffold ordering is reproducible across modification
       levels. Independent units = levels (correlated; see caveat).
H-MOD  The *structure-derived model* recovers that ordering. Independent
       units = scaffolds (3), so the attainable p-value floor is 1/6
       regardless of how many levels are measured, because the model emits a
       single level-invariant ordering.

Conflating these two inflates the apparent evidence for H-MOD. Everything
here is exact-permutation based (no asymptotic approximations, no scipy).
"""

from __future__ import annotations

import itertools
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np
import pandas as pd

PRESPECIFIED_ORDER: tuple[str, ...] = ("casein", "hemoglobin", "collagen")
LEVEL_ORDER: tuple[str, ...] = ("low", "mid", "high")
CONFIG_KEYS: tuple[str, ...] = (
    "quality_mode",
    "domain_mode",
    "component_scheme",
    "model",
    "aggregation",
)
PRIMARY_CONFIG: Mapping[str, str] = {
    "quality_mode": "strict",
    "domain_mode": "in_domain",
    "component_scheme": "stoichiometric",
    "model": "extra_trees",
    "aggregation": "sasa_weighted_mean",
}


# --------------------------------------------------------------------------
# labels
# --------------------------------------------------------------------------
def load_level_labels(path: str | Path, *, winsorize: bool = True) -> pd.DataFrame:
    """Load the 9-observation bulk panel and attach weights and quality flags.

    ``repair_fraction`` > 1 is physically impossible (more than complete
    repair) and signals assay noise; the median-collapse hid this. We keep the
    raw value, add a winsorized copy, and flag the row so downstream analyses
    can be run with and without it.

    Weight = replicates / sd^2 (inverse squared standard error). Rows with a
    missing or zero SD fall back to the median weight.
    """
    frame = pd.read_csv(path)
    required = {"scaffold", "condition_order", "repair_fraction"}
    missing = required - set(frame.columns)
    if missing:
        raise ValueError(f"bulk panel is missing columns: {sorted(missing)}")

    frame["repair_fraction"] = frame["repair_fraction"].astype(float)
    frame["implausible_gt_one"] = frame["repair_fraction"] > 1.0
    frame["repair_fraction_winsorized"] = frame["repair_fraction"].clip(0.0, 1.0)

    sd = frame.get("repair_sd", pd.Series(np.nan, index=frame.index)).astype(float)
    n = frame.get("replicates", pd.Series(1.0, index=frame.index)).astype(float).clip(lower=1.0)
    with np.errstate(divide="ignore", invalid="ignore"):
        raw_w = n / np.square(sd)
    raw_w = pd.Series(raw_w, index=frame.index).replace([np.inf, -np.inf], np.nan)
    fallback = float(np.nanmedian(raw_w)) if np.isfinite(raw_w).any() else 1.0
    frame["weight"] = raw_w.fillna(fallback)
    frame["weight"] = frame["weight"] / frame["weight"].median()

    frame["value_col_used"] = "repair_fraction_winsorized" if winsorize else "repair_fraction"
    return frame


def _value_column(winsorize: bool) -> str:
    return "repair_fraction_winsorized" if winsorize else "repair_fraction"


def _pivot(labels: pd.DataFrame, value_col: str) -> pd.DataFrame:
    table = labels.pivot_table(index="condition_order", columns="scaffold", values=value_col)
    levels = [lv for lv in LEVEL_ORDER if lv in table.index]
    return table.reindex(levels)


# --------------------------------------------------------------------------
# H-OBS: is the observed ordering reproducible across levels?
# --------------------------------------------------------------------------
def observed_level_orderings(labels: pd.DataFrame, *, winsorize: bool = True) -> pd.DataFrame:
    value_col = _value_column(winsorize)
    table = _pivot(labels, value_col)
    rows: list[dict[str, Any]] = []
    for level, series in table.iterrows():
        order = tuple(series.sort_values(ascending=False).index)
        rows.append(
            {
                "condition_order": level,
                "observed_order": " > ".join(order),
                "matches_prespecified": order == PRESPECIFIED_ORDER,
                **{f"value_{s}": float(series[s]) for s in series.index},
            }
        )
    return pd.DataFrame(rows)


def observed_reproducibility_test(
    labels: pd.DataFrame, *, winsorize: bool = True
) -> dict[str, Any]:
    """Exact test that all levels share one ordering.

    Null: each level's ordering is drawn uniformly from the k! permutations,
    independently across levels. Then
        P(all L levels share *some* single ordering) = k! * (1/k!)^L
        P(all L levels equal a *pre-specified* ordering) = (1/k!)^L

    Caveat: the levels are the same proteins and the same enzyme preparation,
    so they are correlated pseudo-replicates, not independent experiments.
    The p-values below are therefore optimistic and must be reported as an
    upper bound on the evidence, not as clean inference.
    """
    orderings = observed_level_orderings(labels, winsorize=winsorize)
    n_levels = len(orderings)
    n_scaffolds = labels["scaffold"].nunique()
    k_fact = math.factorial(n_scaffolds)
    unique_orders = orderings["observed_order"].unique()
    consistent = len(unique_orders) == 1
    all_prespecified = bool(orderings["matches_prespecified"].all())
    return {
        "n_levels": int(n_levels),
        "n_scaffolds": int(n_scaffolds),
        "n_permutations_per_level": k_fact,
        "orderings_observed": list(unique_orders),
        "ordering_consistent_across_levels": consistent,
        "all_levels_match_prespecified": all_prespecified,
        "p_exact_any_consistent_ordering": (
            float(k_fact * (1.0 / k_fact) ** n_levels) if consistent else None
        ),
        "p_exact_prespecified_ordering": (
            float((1.0 / k_fact) ** n_levels) if all_prespecified else None
        ),
        "correlation_caveat": (
            "Levels are the same proteins/enzyme batch; they are correlated "
            "pseudo-replicates. Report these p-values as an upper bound on "
            "evidence, not as independent-replicate inference."
        ),
    }


# --------------------------------------------------------------------------
# dose-response (descriptive, newly surfaced)
# --------------------------------------------------------------------------
def dose_response_summary(labels: pd.DataFrame, *, winsorize: bool = True) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Per-scaffold trend of repair fraction across modification level.

    Level is coded 0/1/2 for low/mid/high. Slope is a weighted least-squares
    fit using the inverse-squared-standard-error weights.
    """
    value_col = _value_column(winsorize)
    code = {lv: i for i, lv in enumerate(LEVEL_ORDER)}
    rows: list[dict[str, Any]] = []
    for scaffold, group in labels.groupby("scaffold"):
        g = group.assign(_x=group["condition_order"].map(code)).sort_values("_x")
        x = g["_x"].to_numpy(dtype=float)
        y = g[value_col].to_numpy(dtype=float)
        w = g["weight"].to_numpy(dtype=float)
        if len(x) >= 2:
            W = np.diag(w)
            X = np.column_stack([np.ones_like(x), x])
            beta = np.linalg.lstsq(X.T @ W @ X, X.T @ W @ y, rcond=None)[0]
            slope = float(beta[1])
        else:
            slope = float("nan")
        rows.append(
            {
                "scaffold": scaffold,
                "n_levels": int(len(x)),
                "monotone_decreasing": bool(np.all(np.diff(y) <= 0)),
                "weighted_slope_per_level": slope,
                **{f"value_{lv}": float(g.loc[g['_x'] == code[lv], value_col].iloc[0])
                   for lv in LEVEL_ORDER if (g["_x"] == code[lv]).any()},
            }
        )
    table = pd.DataFrame(rows)
    n_mono = int(table["monotone_decreasing"].sum())
    n_tot = int(len(table))
    summary = {
        "n_scaffolds_monotone_decreasing": n_mono,
        "n_scaffolds": n_tot,
        "sign_test_one_sided_p": float(0.5 ** n_tot) if n_mono == n_tot else None,
        "mean_weighted_slope_per_level": float(np.nanmean(table["weighted_slope_per_level"])),
        "interpretation": (
            "Descriptive only. With 3 scaffolds a one-sided sign test cannot "
            "fall below p=0.125, so a uniform trend is suggestive, not "
            "established. The structural model is level-invariant and does not "
            "predict this effect."
        ),
    }
    return table, summary


# --------------------------------------------------------------------------
# H-MOD: does the structural model recover the ordering?
# --------------------------------------------------------------------------
def model_ordering_agreement(
    aggregates: pd.DataFrame, *, prespecified: Sequence[str] = PRESPECIFIED_ORDER
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Per-configuration predicted ordering vs the observed ordering.

    The model prediction is level-invariant, so the comparison is made once
    per configuration against the (reproducible) observed ordering.
    """
    target = tuple(prespecified)
    rows: list[dict[str, Any]] = []
    for keys, group in aggregates.groupby(list(CONFIG_KEYS), dropna=False):
        scores = dict(zip(group["scaffold"], group["predicted_repair"].astype(float)))
        complete = set(scores) == set(target)
        order = tuple(sorted(scores, key=lambda s: -scores[s])) if complete else None
        row = dict(zip(CONFIG_KEYS, keys))
        row.update(
            {
                "n_scaffolds_available": len(scores),
                "complete": complete,
                "predicted_order": " > ".join(order) if order else None,
                "matches_observed": bool(order == target) if complete else None,
                "missing_scaffolds": ";".join(sorted(set(target) - set(scores))) or None,
                **{f"pred_{s}": scores.get(s) for s in target},
            }
        )
        rows.append(row)
    table = pd.DataFrame(rows)

    complete_tbl = table[table["complete"]]
    n_complete = int(len(complete_tbl))
    n_match = int(complete_tbl["matches_observed"].sum()) if n_complete else 0
    k_fact = math.factorial(len(target))
    summary = {
        "n_configurations_total": int(len(table)),
        "n_configurations_complete": n_complete,
        "n_configurations_matching_observed": n_match,
        "fraction_matching": float(n_match / n_complete) if n_complete else None,
        "chance_fraction": float(1.0 / k_fact),
        "beats_chance": bool((n_match / n_complete) > (1.0 / k_fact)) if n_complete else None,
        "predicted_order_distribution": (
            complete_tbl["predicted_order"].value_counts().to_dict() if n_complete else {}
        ),
        "interpretation": (
            "Configurations are nested variants of one pipeline and are "
            "therefore highly correlated; the matching fraction is descriptive "
            "and must not be read as a significance test. Selecting the "
            "matching minority post hoc would be a specification-search "
            "artifact."
        ),
    }
    return table, summary


def model_ordering_exact_test(
    aggregates: pd.DataFrame,
    *,
    config: Mapping[str, str] | None = None,
    prespecified: Sequence[str] = PRESPECIFIED_ORDER,
) -> dict[str, Any]:
    """Exact scaffold-permutation test for one configuration.

    Independent units are scaffolds, so the reference distribution has k!
    members and the attainable p-value floor is 1/k! (= 1/6 for 3 scaffolds).
    """
    cfg = dict(config or PRIMARY_CONFIG)
    mask = np.ones(len(aggregates), dtype=bool)
    for key, value in cfg.items():
        mask &= (aggregates[key] == value).to_numpy()
    subset = aggregates[mask]
    target = tuple(prespecified)
    scores = dict(zip(subset["scaffold"], subset["predicted_repair"].astype(float)))
    k_fact = math.factorial(len(target))

    if set(scores) != set(target):
        return {
            "config": cfg,
            "status": "incomplete_scaffold_coverage",
            "available_scaffolds": sorted(scores),
            "missing_scaffolds": sorted(set(target) - set(scores)),
            "p_exact": None,
            "p_floor": float(1.0 / k_fact),
            "interpretation": (
                "The pre-specified configuration cannot be evaluated because at "
                "least one scaffold has no eligible lysine. This is a "
                "representation-validity failure, not a null result."
            ),
        }

    observed_order = tuple(sorted(scores, key=lambda s: -scores[s]))
    hits = sum(
        1 for perm in itertools.permutations(target) if perm == target
    )  # == 1, kept explicit for clarity
    return {
        "config": cfg,
        "status": "evaluated",
        "predicted_order": " > ".join(observed_order),
        "observed_order": " > ".join(target),
        "exact_match": bool(observed_order == target),
        "p_exact": float(hits / k_fact) if observed_order == target else None,
        "p_floor": float(1.0 / k_fact),
        "interpretation": (
            "With 3 scaffolds the smallest attainable p-value for an exact "
            "ordering match is 1/6 = 0.167. Adding modification levels does "
            "not lower this floor, because the model emits a single "
            "level-invariant ordering and the levels are correlated."
        ),
    }


def weighted_continuous_agreement(
    aggregates: pd.DataFrame,
    labels: pd.DataFrame,
    *,
    config: Mapping[str, str] | None = None,
    winsorize: bool = True,
) -> dict[str, Any]:
    """Weighted agreement across all 9 observations for one configuration.

    The predicted value is level-invariant, so each scaffold's prediction is
    paired with all of its levels. The permutation reference reassigns the
    three *scaffold-level* predictions among scaffolds (k! = 6 draws), which
    respects the clustered design instead of pretending there are 9
    independent units.
    """
    cfg = dict(config or PRIMARY_CONFIG)
    value_col = _value_column(winsorize)
    mask = np.ones(len(aggregates), dtype=bool)
    for key, value in cfg.items():
        mask &= (aggregates[key] == value).to_numpy()
    subset = aggregates[mask]
    scores = dict(zip(subset["scaffold"], subset["predicted_repair"].astype(float)))
    scaffolds = sorted(labels["scaffold"].unique())

    if set(scores) != set(scaffolds):
        return {
            "config": cfg,
            "status": "incomplete_scaffold_coverage",
            "missing_scaffolds": sorted(set(scaffolds) - set(scores)),
            "n_observations": None,
        }

    obs = labels[["scaffold", "condition_order", value_col, "weight"]].copy()
    y = obs[value_col].to_numpy(dtype=float)
    w = obs["weight"].to_numpy(dtype=float)

    def _wcorr(pred_map: Mapping[str, float]) -> float:
        x = obs["scaffold"].map(pred_map).to_numpy(dtype=float)
        mx = np.average(x, weights=w)
        my = np.average(y, weights=w)
        cov = np.average((x - mx) * (y - my), weights=w)
        vx = np.average((x - mx) ** 2, weights=w)
        vy = np.average((y - my) ** 2, weights=w)
        return float(cov / np.sqrt(vx * vy)) if vx > 0 and vy > 0 else float("nan")

    observed_r = _wcorr(scores)
    values = [scores[s] for s in scaffolds]
    null_r: list[float] = []
    for perm in itertools.permutations(values):
        null_r.append(_wcorr(dict(zip(scaffolds, perm))))
    null_r_arr = np.asarray(null_r, dtype=float)
    p_exact = float(np.mean(null_r_arr >= observed_r)) if np.isfinite(observed_r) else None

    return {
        "config": cfg,
        "status": "evaluated",
        "n_observations": int(len(obs)),
        "n_independent_clusters": int(len(scaffolds)),
        "weighted_pearson_r": observed_r,
        "p_exact_cluster_permutation": p_exact,
        "p_floor": float(1.0 / len(null_r_arr)),
        "null_r_values": [round(v, 4) for v in sorted(null_r_arr.tolist(), reverse=True)],
        "interpretation": (
            "Using all 9 observations sharpens the estimate of the observed "
            "side and lets the >1 outlier be handled explicitly, but the "
            "independent units remain the 3 scaffolds, so the p-value floor is "
            "1/6. Nine observations do not create nine degrees of freedom."
        ),
    }


# --------------------------------------------------------------------------
# driver
# --------------------------------------------------------------------------
@dataclass
class LevelTestPaths:
    labels: Path
    aggregates: Path
    out_dir: Path


def run(paths: LevelTestPaths, *, winsorize: bool = True) -> dict[str, Any]:
    labels = load_level_labels(paths.labels, winsorize=winsorize)
    aggregates = pd.read_csv(paths.aggregates)
    paths.out_dir.mkdir(parents=True, exist_ok=True)

    orderings = observed_level_orderings(labels, winsorize=winsorize)
    repro = observed_reproducibility_test(labels, winsorize=winsorize)
    dose_tbl, dose_sum = dose_response_summary(labels, winsorize=winsorize)
    agree_tbl, agree_sum = model_ordering_agreement(aggregates)
    primary_exact = model_ordering_exact_test(aggregates)
    primary_cont = weighted_continuous_agreement(aggregates, labels, winsorize=winsorize)

    orderings.to_csv(paths.out_dir / "observed_level_orderings.csv", index=False)
    dose_tbl.to_csv(paths.out_dir / "dose_response_by_scaffold.csv", index=False)
    agree_tbl.to_csv(paths.out_dir / "model_ordering_by_configuration.csv", index=False)
    labels.to_csv(paths.out_dir / "level_labels_with_weights.csv", index=False)

    report = {
        "winsorized": winsorize,
        "observed_reproducibility": repro,
        "dose_response": dose_sum,
        "model_ordering_agreement": agree_sum,
        "primary_configuration_exact_test": primary_exact,
        "primary_configuration_continuous": primary_cont,
    }
    (paths.out_dir / "external_level_test_report.json").write_text(
        json.dumps(report, indent=2, ensure_ascii=False, default=str), encoding="utf-8"
    )
    return report
