import numpy as np
import pandas as pd

from cml_repairnet.models.stage3 import (
    _conformal_quantile,
    _ood_scores,
    assign_bsa_domains,
    generate_unique_spatial_partitions,
    partition_fingerprint,
    seeded_spatial_groups,
    summarize_uncertainty,
)


def _coordinate_frame() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "site_paper": [10, 100, 220, 300, 410, 560],
            "ca_x": [0.0, 0.1, 5.0, 5.1, 10.0, 10.1],
            "ca_y": [0.0, 0.1, 5.0, 5.1, 10.0, 10.1],
            "ca_z": [0.0, 0.1, 5.0, 5.1, 10.0, 10.1],
        }
    )


def test_assign_bsa_domains_uses_mature_numbering() -> None:
    domains = assign_bsa_domains(_coordinate_frame())
    assert domains.tolist() == [
        "domain_I", "domain_I", "domain_II", "domain_II", "domain_III", "domain_III"
    ]


def test_seeded_spatial_groups_are_reproducible() -> None:
    frame = _coordinate_frame()
    first = seeded_spatial_groups(frame, n_groups=3, seed=7)
    second = seeded_spatial_groups(frame, n_groups=3, seed=7)
    assert np.array_equal(first, second)
    assert np.unique(first).size == 3


def test_conformal_quantile_is_finite_and_conservative() -> None:
    residuals = np.array([0.01, 0.02, 0.05, 0.10, 0.30])
    q = _conformal_quantile(residuals, alpha=0.10)
    assert q >= np.quantile(residuals, 0.9)
    assert np.isfinite(q)


def test_ood_scores_flag_distant_observation() -> None:
    train = pd.DataFrame({"a": [0.0, 0.1, -0.1, 0.05, -0.05], "b": [0, 0.1, 0, -0.1, 0]})
    test = pd.DataFrame({"a": [0.0, 20.0], "b": [0.0, 20.0]})
    distances, percentiles, flags = _ood_scores(train, test, n_neighbors=2)
    assert distances[1] > distances[0]
    assert percentiles[1] >= percentiles[0]
    assert flags[1]


def test_uncertainty_summary_reports_coverage() -> None:
    frame = pd.DataFrame(
        {
            "covered": [True, False, True],
            "interval_width": [0.2, 0.3, 0.4],
            "ood_flag": [False, True, False],
            "absolute_error": [0.1, 0.5, 0.2],
            "ood_percentile": [0.1, 1.0, 0.2],
        }
    )
    summary = summarize_uncertainty(frame).iloc[0]
    assert summary["coverage"] == 2 / 3
    assert summary["ood_count"] == 1


def test_partition_fingerprint_is_label_invariant() -> None:
    first = np.array([0, 0, 1, 1, 2, 2])
    relabelled = np.array([8, 8, 4, 4, 1, 1])
    different = np.array([0, 1, 0, 1, 2, 2])
    assert partition_fingerprint(first) == partition_fingerprint(relabelled)
    assert partition_fingerprint(first) != partition_fingerprint(different)


def test_unique_partition_generator_returns_requested_count() -> None:
    rng = np.random.default_rng(42)
    frame = pd.DataFrame(
        {
            "ca_x": rng.normal(size=24),
            "ca_y": rng.normal(size=24),
            "ca_z": rng.normal(size=24),
        }
    )
    partitions = generate_unique_spatial_partitions(
        frame, n_groups=4, start_seed=10, target_count=4, max_candidates=200
    )
    assert len(partitions) == 4
    assert len({fingerprint for _, _, fingerprint in partitions}) == 4
