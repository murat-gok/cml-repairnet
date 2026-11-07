"""Multi-condition reanalysis of the deposited carboxymethyllysine repair data.

This module reproduces every quantitative result reported in the manuscript
"Condition identity and label weighting govern structure-response inference
from reused site-resolved PTM data".

The analysis takes two inputs:

* a parsed copy of the four ProteomeXchange PXD070124 quantification tables
  (``pxd070124_parsed.json``), which records the source URL, the column
  schema and a SHA-256 prefix for each source file; and
* the nine transfer-safe structural descriptors computed from PDB 4F5S for
  the union of confidently quantified sites
  (``bsa_structure_features_union44.csv``).

Everything else -- labels, weights, cross-validation partitions, permutation
controls, resampling references -- is derived here, so the published numbers
can be regenerated from the deposited data alone.

Design notes
------------
The published analysis fixes several choices that are reproduced exactly
rather than re-tuned, because the point of the reanalysis is to test the
published pipeline, not to improve it:

* Response is the truncated relative reduction
  ``clip(1 - mean(active) / mean(inactive), 0, 1)``.
* Label variance follows the delta method for a ratio of two means estimated
  from three replicates each:
  ``var(A)/C^2 + A^2 var(C)/C^4``.
* Weights are inverse-variance with a floor, normalised by the within-load
  median so that a typical site carries unit weight.  The median is taken
  over all sites in that load, before the confidence filter, as in the
  published code.
* Partitions are spatially blocked: Ward linkage on Calpha coordinates.
* The null is the weighted mean of the training folds; significance comes
  from permuting the site-to-structure assignment.

Run via ``scripts/run_multicondition_reanalysis.py``.
"""

from __future__ import annotations

import json
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, spearmanr
from sklearn.cluster import AgglomerativeClustering
from sklearn.ensemble import ExtraTreesRegressor
from sklearn.metrics import mean_absolute_error
from sklearn.preprocessing import StandardScaler

__all__ = [
    "CONDITIONS",
    "FEATURES",
    "SUBSTRATE_NG",
    "MulticonditionConfig",
    "MulticonditionOutputs",
    "assemble_conditions",
    "blocked_cv",
    "build_label_table",
    "common_core",
    "data_coverage",
    "effective_sample_size",
    "label_diagnostics",
    "load_archive",
    "load_structure_features",
    "matched_core",
    "measurability_confound",
    "random_subset_reference",
    "replication_across_conditions",
    "response_dispersion",
    "summarise_subset_reference",
    "run",
    "weighting_sensitivity",
]

# ---------------------------------------------------------------------------
# Fixed analysis constants (these mirror the published pipeline)
# ---------------------------------------------------------------------------

#: Archive keys for the four substrate loads, in increasing order.
CONDITIONS: tuple[str, ...] = ("0p25", "0p5", "1p0", "2p0")

#: Substrate load in nanograms of CML-BSA per condition.  The enzyme
#: concentration is fixed at 5 uM in the source experiment, so these are
#: substrate loads and not enzyme-to-substrate ratios.
SUBSTRATE_NG: dict[str, int] = {"0p25": 250, "0p5": 500, "1p0": 1000, "2p0": 2000}

#: Legacy per-condition scale factor retained in the label table for
#: traceability against the published intermediate files.
_RATIO: dict[str, float] = {"0p25": 0.25, "0p5": 0.5, "1p0": 1.0, "2p0": 2.0}

#: The nine transfer-safe structural descriptors.  No sequence-derived or
#: annotation-derived feature enters the model.
FEATURES: tuple[str, ...] = (
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

_COORDS: tuple[str, ...] = ("ca_x", "ca_y", "ca_z")

#: Number of technical replicates per arm in the source tables.
_N_REPLICATES = 3


@dataclass(frozen=True)
class MulticonditionConfig:
    """Parameters of the reanalysis.

    The defaults reproduce the manuscript.  They are exposed so that a
    reviewer can vary any single choice and see what it costs.
    """

    #: Localisation-probability cutoff applied to site-level quantification.
    confidence_cutoff: float = 90.0
    #: Additive floor on the label variance before inversion, which keeps a
    #: near-zero-variance site from dominating the fit.
    variance_floor: float = 1e-4
    #: Spatial blocks for the main per-condition analysis.
    n_blocks: int = 4
    #: Spatial blocks for the matched core, which is too small for four.
    n_blocks_core: int = 3
    #: Site-structure permutations behind each reported P value.
    n_permutations: int = 200
    #: Draws in the random-subset reference.
    n_subset_draws: int = 200
    #: Size of each random subset, matched to the common core.
    subset_size: int = 20
    #: Seed for the model, the permutation control and the partitioning.
    seed: int = 0
    #: Separate seed for the subset resampling, as in the published run.
    subset_seed: int = 7
    #: Forest size.
    n_estimators: int = 500
    #: Forest depth.
    max_depth: int = 4
    #: Minimum samples in a leaf.
    min_samples_leaf: int = 2
    #: Weighting schemes evaluated in the sensitivity analysis: the published
    #: choice, no weighting, and a moderately clipped variant that tests
    #: whether the difficulty is a few extreme weights or the weighting
    #: itself.  A square-root scheme is available via `_apply_weight_mode`
    #: but is not run by default, as it is not reported.
    weight_modes: tuple[str, ...] = ("inv_var", "none", "clip5")
    #: Clipping range for the ``clip5`` scheme.
    clip_range: tuple[float, float] = (0.2, 5.0)


# ---------------------------------------------------------------------------
# Inputs
# ---------------------------------------------------------------------------


def load_archive(path: str | Path) -> dict:
    """Load the parsed PXD070124 quantification tables.

    The file carries provenance metadata under keys beginning with an
    underscore (``_source``, ``_schema``, ``_sha256_prefix``, ``_issues``)
    and one list of rows per condition.  Each row is

    ``[residue, confidence, active_r1..r3, inactive_r1..r3]``

    with abundances as percentages.
    """
    with open(path, encoding="utf-8") as handle:
        raw = json.load(handle)
    missing = [c for c in CONDITIONS if c not in raw]
    if missing:
        raise ValueError(f"archive is missing conditions: {missing}")
    return raw


def load_structure_features(path: str | Path) -> pd.DataFrame:
    """Load the structural descriptor table and check it is usable."""
    features = pd.read_csv(path)
    needed = {"site_uniprot", *FEATURES, *_COORDS}
    absent = needed - set(features.columns)
    if absent:
        raise ValueError(f"structure feature table lacks columns: {sorted(absent)}")
    if "structure_missing" in features.columns:
        features = features[features["structure_missing"] == 0].copy()
    return features.reset_index(drop=True)


# ---------------------------------------------------------------------------
# Labels and weights
# ---------------------------------------------------------------------------


def build_label_table(raw: dict, config: MulticonditionConfig | None = None) -> pd.DataFrame:
    """Derive the per-site response, its variance and its weight.

    Returns one row per (condition, site) with the untruncated response
    retained alongside the truncated one, so that truncation can be counted
    rather than silently absorbed.
    """
    cfg = config or MulticonditionConfig()
    records: list[dict] = []
    for cond in CONDITIONS:
        for row in raw[cond]:
            site = int(row[0])
            confidence = float(row[1])
            active = np.asarray(row[2 : 2 + _N_REPLICATES], dtype=float)
            inactive = np.asarray(
                row[2 + _N_REPLICATES : 2 + 2 * _N_REPLICATES], dtype=float
            )
            active_mean = active.mean()
            control_mean = inactive.mean()
            if control_mean == 0:
                # A ratio against a zero control carries no information.
                continue
            active_sd = active.std(ddof=1)
            control_sd = inactive.std(ddof=1)
            response_raw = 1.0 - active_mean / control_mean
            # Delta-method variance of a ratio of two independent means.
            variance = (active_sd**2 / _N_REPLICATES) / control_mean**2 + (
                active_mean**2
            ) * (control_sd**2 / _N_REPLICATES) / control_mean**4
            records.append(
                {
                    "cond": cond,
                    "substrate_ng": SUBSTRATE_NG[cond],
                    "ratio": _RATIO[cond],
                    "site_uniprot": site,
                    "conf": confidence,
                    "active_mean": active_mean,
                    "control_mean": control_mean,
                    "resp_raw": response_raw,
                    "resp": float(np.clip(response_raw, 0.0, 1.0)),
                    "var": variance,
                }
            )

    labels = pd.DataFrame.from_records(records)
    labels["var"] = labels["var"].fillna(1.0)
    labels["w_raw"] = 1.0 / (labels["var"] + cfg.variance_floor)
    # Median-normalise within a load, over all sites in that load, before the
    # confidence filter.  This is what the published code does.
    labels["weight"] = labels.groupby("cond")["w_raw"].transform(lambda s: s / s.median())
    return labels


def effective_sample_size(weights: Sequence[float] | np.ndarray) -> float:
    """Kish effective sample size, ``(sum w)^2 / sum w^2``."""
    w = np.asarray(weights, dtype=float)
    denominator = float(np.sum(w**2))
    if denominator == 0:
        return 0.0
    return float(np.sum(w) ** 2 / denominator)


def _apply_weight_mode(
    weights: np.ndarray, mode: str, config: MulticonditionConfig
) -> np.ndarray:
    """Transform inverse-variance weights into one of the tested schemes."""
    if mode == "inv_var":
        return weights
    if mode == "none":
        return np.ones_like(weights)
    if mode == "sqrt":
        return np.sqrt(weights)
    if mode == "clip5":
        low, high = config.clip_range
        return np.clip(weights, low, high)
    raise ValueError(f"unknown weight mode: {mode!r}")


def assemble_conditions(
    labels: pd.DataFrame,
    features: pd.DataFrame,
    config: MulticonditionConfig | None = None,
) -> dict[str, pd.DataFrame]:
    """Join confident labels to descriptors, one analysis frame per load."""
    cfg = config or MulticonditionConfig()
    confident = labels[labels["conf"] > cfg.confidence_cutoff]
    columns = ["site_uniprot", *FEATURES, *_COORDS]
    return {
        cond: confident[confident["cond"] == cond]
        .merge(features[columns], on="site_uniprot")
        .reset_index(drop=True)
        for cond in CONDITIONS
    }


# ---------------------------------------------------------------------------
# Spatially blocked cross-validation with a permutation control
# ---------------------------------------------------------------------------


def blocked_cv(
    frame: pd.DataFrame,
    config: MulticonditionConfig | None = None,
    *,
    weight_mode: str = "inv_var",
    n_blocks: int | None = None,
    n_permutations: int | None = None,
) -> dict:
    """Spatially blocked CV with a site-structure permutation control.

    Blocks come from Ward linkage on Calpha coordinates, so a held-out block
    is a contiguous region of the fold rather than a scatter of sites.  The
    null model is the weighted mean of the training folds.  The reported P
    value permutes the rows of the descriptor matrix, which breaks the
    site-to-structure correspondence while leaving the response distribution
    and the spatial partition untouched.

    Set ``n_permutations=0`` to skip the control, which is what the subset
    resampling does.
    """
    cfg = config or MulticonditionConfig()
    blocks = cfg.n_blocks if n_blocks is None else n_blocks
    n_perm = cfg.n_permutations if n_permutations is None else n_permutations

    X = frame[list(FEATURES)].to_numpy(dtype=float)
    y = frame["resp"].to_numpy(dtype=float)
    w = _apply_weight_mode(frame["weight"].to_numpy(dtype=float), weight_mode, cfg)
    coords = frame[list(_COORDS)].to_numpy(dtype=float)

    groups = AgglomerativeClustering(n_clusters=blocks, linkage="ward").fit_predict(coords)
    unique_groups = np.unique(groups)

    def fit_predict(matrix: np.ndarray) -> tuple[float, float]:
        oof = np.full(len(y), np.nan)
        for group in unique_groups:
            test = groups == group
            train = ~test
            if train.sum() < 4 or test.sum() < 1:
                continue
            scaler = StandardScaler().fit(matrix[train])
            model = ExtraTreesRegressor(
                n_estimators=cfg.n_estimators,
                max_depth=cfg.max_depth,
                min_samples_leaf=cfg.min_samples_leaf,
                random_state=cfg.seed,
                n_jobs=1,
            )
            model.fit(scaler.transform(matrix[train]), y[train], sample_weight=w[train])
            oof[test] = model.predict(scaler.transform(matrix[test]))
        scored = ~np.isnan(oof)
        return (
            mean_absolute_error(y[scored], oof[scored]),
            float(spearmanr(y[scored], oof[scored]).statistic),
        )

    mae, rho = fit_predict(X)

    # Null model: the weighted training-fold mean, scored on the held-out block.
    residuals = []
    for group in unique_groups:
        test = groups == group
        train = ~test
        if train.sum() < 4 or test.sum() < 1:
            continue
        residuals.append(np.abs(y[test] - np.average(y[train], weights=w[train])))
    baseline = float(np.concatenate(residuals).mean())

    result = {
        "n": int(len(y)),
        "groups": int(len(unique_groups)),
        "weight_mode": weight_mode,
        "mae": float(mae),
        "spearman": rho,
        "baseline": baseline,
        "rel_impr": float((baseline - mae) / baseline) if baseline else float("nan"),
        "dMAE_base": float(baseline - mae),
        "ess": effective_sample_size(w),
    }

    if n_perm <= 0:
        result.update(perm_mae_mean=float("nan"), dMAE=float("nan"), p_perm=float("nan"))
        return result

    rng = np.random.default_rng(cfg.seed)
    permuted = np.array(
        [fit_predict(X[rng.permutation(len(y))])[0] for _ in range(n_perm)]
    )
    result.update(
        perm_mae_mean=float(permuted.mean()),
        dMAE=float(permuted.mean() - mae),
        p_perm=float((np.sum(permuted <= mae) + 1) / (n_perm + 1)),
    )
    return result


# ---------------------------------------------------------------------------
# The reported analyses
# ---------------------------------------------------------------------------


def data_coverage(
    labels: pd.DataFrame,
    features: pd.DataFrame,
    config: MulticonditionConfig | None = None,
) -> pd.DataFrame:
    """How many sites each load contributes, and how many survive each filter."""
    cfg = config or MulticonditionConfig()
    confident = labels[labels["conf"] > cfg.confidence_cutoff]
    with_features = set(features["site_uniprot"])
    rows = []
    for cond in CONDITIONS:
        quantified = labels[labels["cond"] == cond]
        high = confident[confident["cond"] == cond]
        rows.append(
            {
                "cond": cond,
                "substrate_ng": SUBSTRATE_NG[cond],
                "quantified": int(len(quantified)),
                "confident": int(len(high)),
                "with_structure": int(len(set(high["site_uniprot"]) & with_features)),
            }
        )
    return pd.DataFrame(rows)


def label_diagnostics(
    conditions: dict[str, pd.DataFrame],
    labels: pd.DataFrame,
    config: MulticonditionConfig | None = None,
) -> pd.DataFrame:
    """Weight range, Kish ESS and truncation counts per load (manuscript Table 2)."""
    cfg = config or MulticonditionConfig()
    confident = labels[labels["conf"] > cfg.confidence_cutoff]
    rows = []
    for cond in CONDITIONS:
        frame = conditions[cond]
        weights = frame["weight"].to_numpy(dtype=float)
        analysed = confident[
            (confident["cond"] == cond)
            & confident["site_uniprot"].isin(frame["site_uniprot"])
        ]
        raw = analysed["resp_raw"].to_numpy(dtype=float)
        response = frame["resp"].to_numpy(dtype=float)
        ess = effective_sample_size(weights)
        rows.append(
            {
                "cond": cond,
                "substrate_ng": SUBSTRATE_NG[cond],
                "n": int(len(frame)),
                "w_min": float(weights.min()),
                "w_max": float(weights.max()),
                "kish_ess": ess,
                "ess_over_n": ess / len(frame) if len(frame) else float("nan"),
                "truncated": int(np.sum((raw < 0.0) | (raw > 1.0))),
                "response_at_zero": int(np.sum(response <= 0.0)),
                "response_at_one": int(np.sum(response >= 1.0)),
            }
        )
    return pd.DataFrame(rows)


def replication_across_conditions(
    conditions: dict[str, pd.DataFrame],
    config: MulticonditionConfig | None = None,
    *,
    weight_modes: Iterable[str] = ("none", "inv_var"),
) -> pd.DataFrame:
    """Does the published structure-response signal hold in all four loads?"""
    cfg = config or MulticonditionConfig()
    rows = []
    for mode in weight_modes:
        for cond in CONDITIONS:
            result = blocked_cv(conditions[cond], cfg, weight_mode=mode)
            rows.append({"cond": cond, "substrate_ng": SUBSTRATE_NG[cond], **result})
    return pd.DataFrame(rows)


def weighting_sensitivity(
    conditions: dict[str, pd.DataFrame], config: MulticonditionConfig | None = None
) -> pd.DataFrame:
    """Compare the four weighting schemes across all loads."""
    cfg = config or MulticonditionConfig()
    return replication_across_conditions(conditions, cfg, weight_modes=cfg.weight_modes)


def common_core(conditions: dict[str, pd.DataFrame]) -> list[int]:
    """Sites confidently quantified with descriptors in every load."""
    sets = [set(conditions[c]["site_uniprot"]) for c in CONDITIONS]
    return sorted(set.intersection(*sets))


def matched_core(
    conditions: dict[str, pd.DataFrame],
    config: MulticonditionConfig | None = None,
    *,
    weight_mode: str = "none",
) -> pd.DataFrame:
    """Re-run every load on the fixed common core, holding composition constant."""
    cfg = config or MulticonditionConfig()
    core = common_core(conditions)
    rows = []
    for cond in CONDITIONS:
        frame = conditions[cond]
        subset = frame[frame["site_uniprot"].isin(core)].reset_index(drop=True)
        result = blocked_cv(
            subset, cfg, weight_mode=weight_mode, n_blocks=cfg.n_blocks_core
        )
        rows.append(
            {"cond": cond, "substrate_ng": SUBSTRATE_NG[cond], "core_size": len(core), **result}
        )
    return pd.DataFrame(rows)


def random_subset_reference(
    conditions: dict[str, pd.DataFrame],
    config: MulticonditionConfig | None = None,
    *,
    reference_cond: str = "0p5",
    weight_mode: str = "none",
) -> tuple[pd.DataFrame, dict]:
    """Is the matched-core result a power effect or a composition effect?

    Draws random subsets of the same size from the published load and scores
    each one, which gives the distribution the core result should be read
    against.  Returns the per-draw table and a summary that locates the real
    core within it.
    """
    cfg = config or MulticonditionConfig()
    frame = conditions[reference_cond].reset_index(drop=True)
    core = common_core(conditions)
    rng = np.random.default_rng(cfg.subset_seed)

    draws = []
    for draw in range(cfg.n_subset_draws):
        index = rng.choice(len(frame), cfg.subset_size, replace=False)
        result = blocked_cv(
            frame.iloc[index].reset_index(drop=True),
            cfg,
            weight_mode=weight_mode,
            n_blocks=cfg.n_blocks_core,
            n_permutations=0,
        )
        draws.append(
            {
                "draw": draw,
                "spearman": result["spearman"],
                "mae": result["mae"],
                "baseline": result["baseline"],
                "dMAE_base": result["dMAE_base"],
            }
        )
    table = pd.DataFrame(draws)

    actual = blocked_cv(
        frame[frame["site_uniprot"].isin(core)].reset_index(drop=True),
        cfg,
        weight_mode=weight_mode,
        n_blocks=cfg.n_blocks_core,
        n_permutations=0,
    )
    return table, summarise_subset_reference(table, conditions, cfg, actual=actual)


def summarise_subset_reference(
    draws: pd.DataFrame,
    conditions: dict[str, pd.DataFrame],
    config: MulticonditionConfig | None = None,
    *,
    reference_cond: str = "0p5",
    weight_mode: str = "none",
    actual: dict | None = None,
) -> dict:
    """Locate the matched core within the random-subset distribution.

    Split out from :func:`random_subset_reference` so that a draws table
    read back from disk can be summarised without refitting the two hundred
    subsets behind it.
    """
    cfg = config or MulticonditionConfig()
    core = common_core(conditions)
    if actual is None:
        frame = conditions[reference_cond].reset_index(drop=True)
        actual = blocked_cv(
            frame[frame["site_uniprot"].isin(core)].reset_index(drop=True),
            cfg,
            weight_mode=weight_mode,
            n_blocks=cfg.n_blocks_core,
            n_permutations=0,
        )
    rho = draws["spearman"].to_numpy()
    improvement = draws["dMAE_base"].to_numpy()
    return {
        "reference_cond": reference_cond,
        "reference_substrate_ng": SUBSTRATE_NG[reference_cond],
        "n_draws": int(len(draws)),
        "subset_size": cfg.subset_size,
        "core_size": len(core),
        "spearman_mean": float(rho.mean()),
        "spearman_sd": float(rho.std()),
        "spearman_median": float(np.percentile(rho, 50)),
        "spearman_p5": float(np.percentile(rho, 5)),
        "spearman_p95": float(np.percentile(rho, 95)),
        "fraction_above_0.3": float((rho > 0.3).mean()),
        "fraction_above_0": float((rho > 0).mean()),
        "improvement_mean": float(improvement.mean()),
        "improvement_median": float(np.percentile(improvement, 50)),
        "core_spearman": actual["spearman"],
        "core_improvement": actual["dMAE_base"],
        "core_percentile": float((rho < actual["spearman"]).mean()),
    }


def response_dispersion(conditions: dict[str, pd.DataFrame]) -> pd.DataFrame:
    """Spread of the response within each load."""
    rows = []
    for cond in CONDITIONS:
        y = conditions[cond]["resp"].to_numpy(dtype=float)
        q25, q75 = np.percentile(y, [25, 75])
        rows.append(
            {
                "cond": cond,
                "substrate_ng": SUBSTRATE_NG[cond],
                "n": int(len(y)),
                "mean": float(y.mean()),
                "sd": float(y.std(ddof=1)),
                "min": float(y.min()),
                "max": float(y.max()),
                "iqr": float(q75 - q25),
                "at_or_below_0.10": int(np.sum(y <= 0.10)),
                "at_or_above_0.90": int(np.sum(y >= 0.90)),
            }
        )
    return pd.DataFrame(rows)


def measurability_confound(
    labels: pd.DataFrame,
    features: pd.DataFrame,
    conditions: dict[str, pd.DataFrame],
    config: MulticonditionConfig | None = None,
    *,
    reference_cond: str = "0p5",
) -> dict[str, pd.DataFrame]:
    """Test whether structure predicts measurability rather than response.

    If the core sites were simply the easy-to-measure ones, structure would
    track how often a site clears the confidence filter, and the core and
    non-core sites would differ structurally.  Three tables report that
    check; a fourth asks whether confidence tracks the response itself.
    """
    cfg = config or MulticonditionConfig()
    confident = labels[labels["conf"] > cfg.confidence_cutoff]
    per_condition_sites = {
        cond: set(confident[confident["cond"] == cond]["site_uniprot"])
        for cond in CONDITIONS
    }

    sites = features[["site_uniprot", *FEATURES]].copy()
    sites["n_cond_pass"] = [
        sum(1 for cond in CONDITIONS if site in per_condition_sites[cond])
        for site in sites["site_uniprot"]
    ]
    sites["mean_conf"] = sites["site_uniprot"].map(
        labels.groupby("site_uniprot")["conf"].mean()
    )

    structure_vs_measurability = pd.DataFrame(
        [
            {
                "feature": feature,
                "rho_n_cond_pass": spearmanr(sites[feature], sites["n_cond_pass"]).statistic,
                "p_n_cond_pass": spearmanr(sites[feature], sites["n_cond_pass"]).pvalue,
                "rho_mean_conf": spearmanr(sites[feature], sites["mean_conf"]).statistic,
                "p_mean_conf": spearmanr(sites[feature], sites["mean_conf"]).pvalue,
            }
            for feature in FEATURES
        ]
    )

    core = set(common_core(conditions))
    reference = conditions[reference_cond]
    non_core = set(reference["site_uniprot"]) - core

    core_vs_rest = pd.DataFrame(
        [
            {
                "feature": feature,
                "core_mean": float(sites[sites["site_uniprot"].isin(core)][feature].mean()),
                "rest_mean": float(
                    sites[sites["site_uniprot"].isin(non_core)][feature].mean()
                ),
                "mannwhitney_p": float(
                    mannwhitneyu(
                        sites[sites["site_uniprot"].isin(core)][feature],
                        sites[sites["site_uniprot"].isin(non_core)][feature],
                    ).pvalue
                ),
            }
            for feature in FEATURES
        ]
    )

    merged = confident[confident["cond"] == reference_cond].merge(
        features[["site_uniprot"]], on="site_uniprot"
    )
    conf_vs_response = spearmanr(merged["conf"], merged["resp"])

    core_response = reference[reference["site_uniprot"].isin(core)]["resp"]
    rest_response = reference[reference["site_uniprot"].isin(non_core)]["resp"]
    response_split = mannwhitneyu(core_response, rest_response)

    scalars = pd.DataFrame(
        [
            {
                "quantity": "spearman_confidence_vs_response",
                "value": float(conf_vs_response.statistic),
                "p": float(conf_vs_response.pvalue),
                "n": int(len(merged)),
            },
            {
                "quantity": "mannwhitney_core_vs_rest_response",
                "value": float(core_response.mean() - rest_response.mean()),
                "p": float(response_split.pvalue),
                "n": int(len(core_response) + len(rest_response)),
            },
            {
                "quantity": "bonferroni_threshold",
                "value": 0.05 / len(FEATURES),
                "p": float("nan"),
                "n": len(FEATURES),
            },
        ]
    )

    return {
        "site_measurability": sites,
        "structure_vs_measurability": structure_vs_measurability,
        "core_vs_rest": core_vs_rest,
        "scalars": scalars,
    }


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------


def _environment() -> dict[str, str]:
    """Record what produced these numbers.

    Permutation P values are counts of how many permuted fits beat the real
    one, so a value sitting a hair under the threshold can move across it
    between library versions even when the fitted model is identical to
    three decimals. Recording the environment lets a later run tell a real
    difference from a numerical one.
    """
    import platform

    import sklearn

    return {
        "python": platform.python_version(),
        "numpy": np.__version__,
        "pandas": pd.__version__,
        "scipy": __import__("scipy").__version__,
        "scikit-learn": sklearn.__version__,
        "platform": platform.platform(),
    }


@dataclass
class MulticonditionOutputs:
    """Everything the reanalysis writes, kept in memory for testing."""

    labels: pd.DataFrame
    coverage: pd.DataFrame
    diagnostics: pd.DataFrame
    replication: pd.DataFrame
    weighting: pd.DataFrame
    core: pd.DataFrame
    subset_draws: pd.DataFrame
    subset_summary: dict
    dispersion: pd.DataFrame
    confound: dict[str, pd.DataFrame] = field(default_factory=dict)


def run(
    archive_path: str | Path,
    features_path: str | Path,
    output_dir: str | Path,
    config: MulticonditionConfig | None = None,
    *,
    quiet: bool = False,
) -> MulticonditionOutputs:
    """Run the whole reanalysis and write every result table.

    Parameters mirror the manuscript's data availability statement: the
    parsed archive, the structural descriptors, and a directory to write to.
    """
    cfg = config or MulticonditionConfig()
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    raw = load_archive(archive_path)
    features = load_structure_features(features_path)

    labels = build_label_table(raw, cfg)
    conditions = assemble_conditions(labels, features, cfg)

    def staged(name: str, compute):
        """Compute a table, or reuse one already on disk.

        The permutation controls and the subset reference together fit
        several thousand forests, which on a small machine runs for well
        over an hour.  Each table is therefore written as soon as it is
        finished and reused on a later call, so an interrupted run resumes
        instead of starting over.  Delete a file to force it to be redone.
        """
        path = out / f"{name}.csv"
        if path.exists():
            if not quiet:
                print(f"  [reuse] {path.name}", flush=True)
            return pd.read_csv(path)
        if not quiet:
            print(f"  [run]   {name} ...", end="", flush=True)
            started = time.time()
        table = compute()
        table.to_csv(path, index=False)
        if not quiet:
            print(f" {time.time() - started:.0f}s", flush=True)
        return table

    labels.to_csv(out / "level_labels_with_weights.csv", index=False)
    coverage = staged("data_coverage", lambda: data_coverage(labels, features, cfg))
    diagnostics = staged("label_diagnostics", lambda: label_diagnostics(conditions, labels, cfg))
    dispersion = staged("response_distribution", lambda: response_dispersion(conditions))

    # The weighting sweep dominates the runtime, so each scheme is stored
    # separately and the combined table is assembled from the parts.
    per_mode = {
        mode: staged(
            f"weighting_{mode}",
            lambda mode=mode: replication_across_conditions(
                conditions, cfg, weight_modes=[mode]
            ),
        )
        for mode in cfg.weight_modes
    }
    weighting = pd.concat(
        [per_mode[mode] for mode in cfg.weight_modes], ignore_index=True
    )
    weighting.to_csv(out / "weighting_sensitivity.csv", index=False)
    replication = (
        weighting[weighting["weight_mode"].isin(["none", "inv_var"])]
        .reset_index(drop=True)
        .copy()
    )
    replication.to_csv(out / "condition_replication.csv", index=False)

    core = staged("matched_core", lambda: matched_core(conditions, cfg))
    subset_draws = staged(
        "random_subset_draws", lambda: random_subset_reference(conditions, cfg)[0]
    )
    subset_summary = summarise_subset_reference(subset_draws, conditions, cfg)

    confound = measurability_confound(labels, features, conditions, cfg)
    for name, table in confound.items():
        table.to_csv(out / f"confound_{name}.csv", index=False)

    manifest = {
        "environment": _environment(),
        "source": raw.get("_source"),
        "schema": raw.get("_schema"),
        "sha256_prefix": raw.get("_sha256_prefix"),
        "parse_issues": raw.get("_issues"),
        "config": {
            key: list(value) if isinstance(value, tuple) else value
            for key, value in cfg.__dict__.items()
        },
        "n_sites_with_structure": int(len(features)),
        "common_core": common_core(conditions),
        "subset_summary": subset_summary,
    }
    with open(out / "multicondition_manifest.json", "w", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2, sort_keys=True)

    return MulticonditionOutputs(
        labels=labels,
        coverage=coverage,
        diagnostics=diagnostics,
        replication=replication,
        weighting=weighting,
        core=core,
        subset_draws=subset_draws,
        subset_summary=subset_summary,
        dispersion=dispersion,
        confound=confound,
    )
