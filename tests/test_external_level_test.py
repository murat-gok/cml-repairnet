from __future__ import annotations

import math

import pandas as pd
import pytest

from cml_repairnet.revision.external_level_test import (
    PRESPECIFIED_ORDER,
    dose_response_summary,
    load_level_labels,
    model_ordering_agreement,
    model_ordering_exact_test,
    observed_level_orderings,
    observed_reproducibility_test,
    weighted_continuous_agreement,
)

SCAFFOLDS = ("casein", "hemoglobin", "collagen")
LEVELS = ("low", "mid", "high")
# Trabosh et al. 2026 Supplementary Table 4 values as stored in the repo.
VALUES = {
    ("casein", "low"): (1.34, 0.61),
    ("casein", "mid"): (0.97, 0.06),
    ("casein", "high"): (0.92, 0.01),
    ("hemoglobin", "low"): (0.89, 0.07),
    ("hemoglobin", "mid"): (0.87, 0.06),
    ("hemoglobin", "high"): (0.82, 0.02),
    ("collagen", "low"): (0.80, 0.06),
    ("collagen", "mid"): (0.65, 0.05),
    ("collagen", "high"): (0.60, 0.05),
}


@pytest.fixture()
def labels_path(tmp_path):
    rows = [
        {
            "scaffold": s,
            "condition_order": lv,
            "repair_fraction": VALUES[(s, lv)][0],
            "repair_sd": VALUES[(s, lv)][1],
            "replicates": 4,
        }
        for s in SCAFFOLDS
        for lv in LEVELS
    ]
    path = tmp_path / "bulk.csv"
    pd.DataFrame(rows).to_csv(path, index=False)
    return path


def aggregates(order: tuple[str, ...], *, drop: str | None = None) -> pd.DataFrame:
    """Build a one-configuration aggregate table predicting ``order``."""
    scores = {name: 0.9 - 0.1 * i for i, name in enumerate(order)}
    rows = []
    for scaffold, score in scores.items():
        if scaffold == drop:
            continue
        rows.append(
            {
                "scaffold": scaffold,
                "quality_mode": "strict",
                "domain_mode": "in_domain",
                "component_scheme": "stoichiometric",
                "model": "extra_trees",
                "aggregation": "sasa_weighted_mean",
                "predicted_repair": score,
            }
        )
    return pd.DataFrame(rows)


def test_loader_flags_impossible_repair_and_winsorizes(labels_path):
    frame = load_level_labels(labels_path)
    bad = frame[frame["implausible_gt_one"]]
    assert len(bad) == 1
    assert bad.iloc[0]["scaffold"] == "casein"
    assert bad.iloc[0]["repair_fraction_winsorized"] == 1.0
    # noisiest point must not dominate the inverse-variance weights
    assert bad.iloc[0]["weight"] < frame["weight"].median()


def test_observed_ordering_is_reproducible_across_all_levels(labels_path):
    frame = load_level_labels(labels_path)
    orderings = observed_level_orderings(frame)
    assert len(orderings) == 3
    assert orderings["matches_prespecified"].all()

    result = observed_reproducibility_test(frame)
    assert result["ordering_consistent_across_levels"] is True
    assert result["p_exact_prespecified_ordering"] == pytest.approx((1 / 6) ** 3)
    assert result["p_exact_any_consistent_ordering"] == pytest.approx(6 * (1 / 6) ** 3)
    assert "correlated" in result["correlation_caveat"]


def test_dose_response_is_monotone_but_underpowered(labels_path):
    frame = load_level_labels(labels_path)
    table, summary = dose_response_summary(frame)
    assert summary["n_scaffolds_monotone_decreasing"] == 3
    # a 3-scaffold sign test cannot go below 0.125
    assert summary["sign_test_one_sided_p"] == pytest.approx(0.125)
    assert summary["mean_weighted_slope_per_level"] < 0
    assert set(table["scaffold"]) == set(SCAFFOLDS)


def test_model_agreement_detects_match_and_mismatch(labels_path):
    _, matched = model_ordering_agreement(aggregates(PRESPECIFIED_ORDER))
    assert matched["n_configurations_complete"] == 1
    assert matched["fraction_matching"] == 1.0

    wrong = ("hemoglobin", "casein", "collagen")
    _, mismatched = model_ordering_agreement(aggregates(wrong))
    assert mismatched["n_configurations_matching_observed"] == 0
    assert mismatched["fraction_matching"] == 0.0


def test_exact_test_reports_one_sixth_floor(labels_path):
    result = model_ordering_exact_test(aggregates(PRESPECIFIED_ORDER))
    assert result["status"] == "evaluated"
    assert result["exact_match"] is True
    assert result["p_exact"] == pytest.approx(1 / 6)
    assert result["p_floor"] == pytest.approx(1 / 6)


def test_incomplete_coverage_is_not_reported_as_a_null_result(labels_path):
    result = model_ordering_exact_test(aggregates(PRESPECIFIED_ORDER, drop="collagen"))
    assert result["status"] == "incomplete_scaffold_coverage"
    assert result["missing_scaffolds"] == ["collagen"]
    assert result["p_exact"] is None
    assert "representation-validity failure" in result["interpretation"]


def test_continuous_agreement_respects_scaffold_clustering(labels_path):
    frame = load_level_labels(labels_path)
    result = weighted_continuous_agreement(aggregates(PRESPECIFIED_ORDER), frame)
    assert result["status"] == "evaluated"
    # nine observations, but only three independent clusters
    assert result["n_observations"] == 9
    assert result["n_independent_clusters"] == 3
    assert result["p_floor"] == pytest.approx(1 / math.factorial(3))
    assert result["p_exact_cluster_permutation"] >= result["p_floor"] - 1e-12
    assert result["weighted_pearson_r"] > 0
