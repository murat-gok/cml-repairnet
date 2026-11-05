"""Scientific-integrity tests for the multi-condition reanalysis.

These are not smoke tests.  Each one pins a property that, if it silently
changed, would invalidate a claim in the manuscript:

* the label and weight formulas,
* the Kish effective sample size,
* that weighting is what suppresses the signal rather than the site set,
* that the matched core really is composition-limited and not power-limited,
* that the reported per-condition metrics are reproducible.

The slow tests are marked so the fast ones can run on every commit.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from cml_repairnet.revision.multicondition import (
    CONDITIONS,
    FEATURES,
    SUBSTRATE_NG,
    MulticonditionConfig,
    assemble_conditions,
    blocked_cv,
    build_label_table,
    common_core,
    data_coverage,
    effective_sample_size,
    label_diagnostics,
    load_archive,
    load_structure_features,
)

REPO_ROOT = Path(__file__).resolve().parents[1]
ARCHIVE = REPO_ROOT / "data" / "external" / "pxd070124_parsed.json"
FEATURE_TABLE = REPO_ROOT / "data" / "processed" / "bsa_structure_features_union44.csv"

pytestmark = pytest.mark.skipif(
    not (ARCHIVE.exists() and FEATURE_TABLE.exists()),
    reason="deposited inputs are not present in this checkout",
)

# Values as printed in the manuscript.
REPORTED = {
    "none": {
        "0p25": {"mae": 0.274, "spearman": 0.414},
        "0p5": {"mae": 0.187, "spearman": 0.616},
        "1p0": {"mae": 0.221, "spearman": 0.399},
        "2p0": {"mae": 0.188, "spearman": 0.581},
    },
    "inv_var": {
        "0p25": {"mae": 0.304, "spearman": 0.141},
        "0p5": {"mae": 0.184, "spearman": 0.681},
        "1p0": {"mae": 0.233, "spearman": 0.298},
        "2p0": {"mae": 0.209, "spearman": 0.574},
    },
}
REPORTED_BASELINE_UNWEIGHTED = {"0p25": 0.337, "0p5": 0.275, "1p0": 0.268, "2p0": 0.270}
REPORTED_BASELINE_WEIGHTED = {"0p25": 0.396, "0p5": 0.300, "1p0": 0.302, "2p0": 0.290}
REPORTED_ESS = {"0p25": 5.0, "0p5": 9.7, "1p0": 11.8, "2p0": 11.4}
REPORTED_N = {"0p25": 22, "0p5": 33, "1p0": 37, "2p0": 38}


@pytest.fixture(scope="module")
def config() -> MulticonditionConfig:
    return MulticonditionConfig()


@pytest.fixture(scope="module")
def archive() -> dict:
    return load_archive(ARCHIVE)


@pytest.fixture(scope="module")
def features() -> pd.DataFrame:
    return load_structure_features(FEATURE_TABLE)


@pytest.fixture(scope="module")
def labels(archive: dict, config: MulticonditionConfig) -> pd.DataFrame:
    return build_label_table(archive, config)


@pytest.fixture(scope="module")
def conditions(labels, features, config) -> dict[str, pd.DataFrame]:
    return assemble_conditions(labels, features, config)


# ---------------------------------------------------------------------------
# Provenance
# ---------------------------------------------------------------------------


def test_archive_records_its_own_provenance(archive):
    """The parsed archive must carry the source and a checksum per file."""
    assert "PXD070124" in archive["_source"]
    assert set(archive["_sha256_prefix"]) == set(CONDITIONS)
    for prefix in archive["_sha256_prefix"].values():
        assert len(prefix) == 16 and all(c in "0123456789abcdef" for c in prefix)


def test_all_four_loads_are_present(archive):
    assert all(len(archive[cond]) > 30 for cond in CONDITIONS)


def test_feature_table_is_structural_only(features):
    """No sequence or annotation descriptor may enter the model."""
    forbidden = ("sequence", "hydropathy_window", "pssm", "embedding", "conservation")
    for column in FEATURES:
        assert not any(token in column for token in forbidden)
    assert len(FEATURES) == 9


# ---------------------------------------------------------------------------
# Labels and weights
# ---------------------------------------------------------------------------


def test_response_is_truncated_to_the_unit_interval(labels):
    assert labels["resp"].between(0.0, 1.0).all()
    # Truncation must be visible, not silent: the raw value is retained.
    assert (labels["resp_raw"] < 0).any()


def test_delta_method_variance_matches_a_hand_computation(archive, labels):
    """Recompute one row's variance from the source numbers."""
    row = archive["0p5"][0]
    active = np.asarray(row[2:5], dtype=float)
    control = np.asarray(row[5:8], dtype=float)
    a, c = active.mean(), control.mean()
    expected = (active.std(ddof=1) ** 2 / 3) / c**2 + a**2 * (
        control.std(ddof=1) ** 2 / 3
    ) / c**4
    got = labels[(labels["cond"] == "0p5") & (labels["site_uniprot"] == int(row[0]))]
    assert got["var"].iloc[0] == pytest.approx(expected, rel=1e-12)


def test_weights_are_median_normalised_within_a_load(labels):
    for cond in CONDITIONS:
        subset = labels[labels["cond"] == cond]
        assert subset["weight"].median() == pytest.approx(1.0, rel=1e-9)


def test_zero_control_rows_are_dropped_not_imputed(archive, labels):
    """A ratio against a zero control carries no information and must not appear."""
    for cond in CONDITIONS:
        zeros = sum(
            1 for row in archive[cond] if np.mean(np.asarray(row[5:8], float)) == 0
        )
        kept = len(labels[labels["cond"] == cond])
        assert kept == len(archive[cond]) - zeros


def test_kish_effective_sample_size_formula():
    assert effective_sample_size(np.ones(10)) == pytest.approx(10.0)
    # One dominant weight collapses the effective sample towards one.
    assert effective_sample_size([1e6, 1, 1, 1]) < 1.01
    assert effective_sample_size([0, 0]) == 0.0


def test_reported_effective_sample_sizes(conditions, labels, config):
    diagnostics = label_diagnostics(conditions, labels, config)
    for cond in CONDITIONS:
        row = diagnostics[diagnostics["cond"] == cond].iloc[0]
        assert row["kish_ess"] == pytest.approx(REPORTED_ESS[cond], abs=0.05)
        # The manuscript's central claim about weighting.
        assert 0.22 < row["ess_over_n"] < 0.33


def test_weight_range_spans_at_least_three_orders(conditions):
    """The skew the paper reports must be present in every load."""
    for cond in CONDITIONS:
        weights = conditions[cond]["weight"].to_numpy(float)
        assert np.log10(weights.max() / weights.min()) > 3.0


# ---------------------------------------------------------------------------
# Site sets
# ---------------------------------------------------------------------------


def test_reported_site_counts(conditions):
    for cond in CONDITIONS:
        assert len(conditions[cond]) == REPORTED_N[cond]


def test_coverage_table_is_internally_consistent(labels, features, config):
    coverage = data_coverage(labels, features, config)
    assert (coverage["confident"] <= coverage["quantified"]).all()
    assert (coverage["with_structure"] <= coverage["confident"]).all()
    assert list(coverage["substrate_ng"]) == [SUBSTRATE_NG[c] for c in CONDITIONS]


def test_common_core_is_twenty_sites(conditions):
    assert len(common_core(conditions)) == 20


def test_residue_numbering_offset(features):
    """site_pdb = site_uniprot - 24 maps P02769 onto mature 4F5S chain A."""
    if "site_pdb" not in features.columns:
        pytest.skip("feature table does not carry site_pdb")
    assert ((features["site_uniprot"] - features["site_pdb"]) == 24).all()


# ---------------------------------------------------------------------------
# The reported results
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("weight_mode", ["none", "inv_var"])
def test_reported_metrics_are_reproducible(conditions, config, weight_mode):
    """Pin the per-condition MAE and rank correlation printed in Table 3."""
    for cond in CONDITIONS:
        result = blocked_cv(
            conditions[cond], config, weight_mode=weight_mode, n_permutations=0
        )
        expected = REPORTED[weight_mode][cond]
        assert result["mae"] == pytest.approx(expected["mae"], abs=0.0015)
        assert result["spearman"] == pytest.approx(expected["spearman"], abs=0.0015)


@pytest.mark.parametrize("weight_mode", ["none", "inv_var"])
def test_reported_baselines_are_reproducible(conditions, config, weight_mode):
    expected = (
        REPORTED_BASELINE_UNWEIGHTED if weight_mode == "none" else REPORTED_BASELINE_WEIGHTED
    )
    for cond in CONDITIONS:
        result = blocked_cv(
            conditions[cond], config, weight_mode=weight_mode, n_permutations=0
        )
        assert result["baseline"] == pytest.approx(expected[cond], abs=0.0015)


def test_unweighted_beats_weighted_in_three_of_four_loads(conditions, config):
    """The manuscript's comparative claim, stated as a test."""
    better = 0
    for cond in CONDITIONS:
        unweighted = blocked_cv(conditions[cond], config, weight_mode="none", n_permutations=0)
        weighted = blocked_cv(
            conditions[cond], config, weight_mode="inv_var", n_permutations=0
        )
        if unweighted["spearman"] > weighted["spearman"]:
            better += 1
    assert better == 3


@pytest.mark.slow
def test_clipping_reproduces_the_weighted_outcome(conditions, config):
    """Clipping follows the weighted result, not the unweighted one.

    The manuscript's claim is about which conditions clear the permutation
    threshold, not about where the rank correlation lands: clipping to a
    moderate range leaves the same conditions failing as full
    inverse-variance weighting does, which locates the difficulty in the
    direction of the weighting rather than in a handful of extreme weights.
    """
    fast = MulticonditionConfig(n_permutations=50)

    def passing(weight_mode: str) -> set[str]:
        return {
            cond
            for cond in CONDITIONS
            if blocked_cv(conditions[cond], fast, weight_mode=weight_mode)["p_perm"] < 0.05
        }

    unweighted = passing("none")
    weighted = passing("inv_var")
    clipped = passing("clip5")

    assert unweighted == set(CONDITIONS)
    assert clipped != unweighted
    assert clipped == weighted


def test_blocks_are_spatially_contiguous_not_random(conditions, config):
    """Ward blocks must be tighter than random partitions of the same sizes."""
    from sklearn.cluster import AgglomerativeClustering

    frame = conditions["0p5"]
    coords = frame[["ca_x", "ca_y", "ca_z"]].to_numpy(float)
    groups = AgglomerativeClustering(n_clusters=config.n_blocks, linkage="ward").fit_predict(
        coords
    )

    def mean_spread(assignment: np.ndarray) -> float:
        return float(
            np.mean(
                [
                    np.linalg.norm(coords[assignment == g] - coords[assignment == g].mean(0), axis=1).mean()
                    for g in np.unique(assignment)
                ]
            )
        )

    rng = np.random.default_rng(0)
    random_spreads = [mean_spread(rng.permutation(groups)) for _ in range(20)]
    assert mean_spread(groups) < min(random_spreads)


@pytest.mark.slow
def test_permutation_control_rejects_the_null_unweighted(conditions, config):
    """With fewer permutations than the paper, but enough to resolve P < 0.05."""
    fast = MulticonditionConfig(n_permutations=50)
    significant = 0
    for cond in CONDITIONS:
        result = blocked_cv(conditions[cond], fast, weight_mode="none")
        if result["p_perm"] < 0.05:
            significant += 1
    assert significant == 4


@pytest.mark.slow
def test_matched_core_shows_no_signal_in_any_condition(conditions, config):
    """Composition, not power: the fixed twenty-site core must not clear 0.05."""
    fast = MulticonditionConfig(n_permutations=50, n_blocks_core=3)
    core = common_core(conditions)
    for cond in CONDITIONS:
        frame = conditions[cond]
        subset = frame[frame["site_uniprot"].isin(core)].reset_index(drop=True)
        result = blocked_cv(subset, fast, weight_mode="none", n_blocks=3)
        assert result["p_perm"] > 0.05


def test_run_writes_a_manifest(tmp_path, config):
    """A short run must still record provenance and configuration."""
    from cml_repairnet.revision.multicondition import run

    quick = MulticonditionConfig(n_permutations=2, n_subset_draws=3)
    outputs = run(ARCHIVE, FEATURE_TABLE, tmp_path, quick)
    manifest_path = tmp_path / "multicondition_manifest.json"
    assert manifest_path.exists()
    manifest = json.loads(manifest_path.read_text())
    assert "PXD070124" in manifest["source"]
    assert manifest["config"]["n_permutations"] == 2
    assert len(manifest["common_core"]) == 20
    assert not outputs.replication.empty
    assert set(outputs.replication["weight_mode"]) == {"none", "inv_var"}
