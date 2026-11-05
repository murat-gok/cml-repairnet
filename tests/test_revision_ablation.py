from __future__ import annotations

import numpy as np
import pandas as pd

from cml_repairnet.revision.feature_ablation import (
    BFACTOR_COLUMNS,
    FeatureFamily,
    SpatialPartition,
    add_contact_flexibility_proxies,
    build_feature_families,
    evaluate_feature_family_partition,
)


def _synthetic_tables(n: int = 24) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    rng = np.random.default_rng(7)
    site = np.arange(1, n + 1)
    labels = pd.DataFrame(
        {
            "protein_id": ["P02769"] * n,
            "site_uniprot": site + 24,
            "site_paper": site,
            "spatial_group": np.repeat(np.arange(4), n // 4),
            "repair_published": np.clip(0.2 + 0.02 * site + rng.normal(0, 0.03, n), 0, 1),
            "repair_from_means_clipped": np.clip(0.2 + 0.02 * site, 0, 1),
            "label_weight": np.ones(n),
        }
    )
    sequence = pd.DataFrame(
        {
            "protein_id": ["P02769"] * n,
            "site_uniprot": site + 24,
            "distance_to_n_terminus": site,
            "window_mean_hydropathy": rng.normal(size=n),
            "window_cysteine_count": rng.integers(0, 3, size=n),
        }
    )
    structure = pd.DataFrame(
        {
            "protein_id": ["P02769"] * n,
            "site_uniprot": site + 24,
            "site_pdb": site,
            "structure_missing": np.zeros(n),
            "ca_x": np.repeat(np.arange(4), n // 4) * 20 + rng.normal(0, 0.2, n),
            "ca_y": rng.normal(0, 0.2, n),
            "ca_z": rng.normal(0, 0.2, n),
            "bfactor_ca": 10 + site,
            "neighbor_mean_bfactor_8a": 12 + site,
            "residue_sasa": 50 + site,
            "relative_sasa_lysine": np.clip(0.2 + site / 100, 0, 1.5),
            "distance_to_protein_centroid": np.abs(site - n / 2),
            "neighbor_hydrophobic_fraction_8a": rng.random(n),
            "neighbor_cysteine_count_8a": rng.integers(0, 2, size=n),
            "neighbor_aromatic_count_8a": rng.integers(0, 3, size=n),
            "ca_contact_count_8a": rng.integers(3, 12, size=n),
            "ca_contact_count_10a": rng.integers(8, 18, size=n),
            "ca_contact_count_12a": rng.integers(14, 25, size=n),
        }
    )
    return labels, sequence, structure


def test_transfer_safe_families_exclude_bfactors() -> None:
    labels, sequence, structure = _synthetic_tables()
    families = build_feature_families(labels, sequence, structure)
    assert set(BFACTOR_COLUMNS).isdisjoint(families["combined_transfer_safe"].feature_columns)
    assert set(BFACTOR_COLUMNS).isdisjoint(families["structure_transfer_safe"].feature_columns)
    assert set(BFACTOR_COLUMNS).issubset(families["combined_full"].feature_columns)


def test_contact_proxy_columns_are_finite() -> None:
    _, _, structure = _synthetic_tables()
    out = add_contact_flexibility_proxies(structure)
    proxy = out.filter(regex=r"^flex_proxy_")
    assert proxy.shape[1] == 6
    assert np.isfinite(proxy.to_numpy(dtype=float)).all()


def test_small_oof_ablation_runs() -> None:
    labels, sequence, structure = _synthetic_tables()
    family = build_feature_families(labels, sequence, structure)["combined_transfer_safe"]
    groups = labels["spatial_group"].to_numpy()
    partition = SpatialPartition(
        partition_id="synthetic_g4",
        groups=groups,
        method="synthetic",
        seed=1,
        fingerprint="synthetic",
        n_groups=4,
    )
    metrics, predictions = evaluate_feature_family_partition(
        family,
        partition,
        target="repair_published",
        model_names=("ridge",),
        n_jobs=1,
    )
    assert len(metrics) == 1
    assert len(predictions) == len(labels)
    assert np.isfinite(metrics.loc[0, "mae"])
    assert predictions["outer_fold"].nunique() == 4


def test_structure_only_proxy_family_is_isolated() -> None:
    labels, sequence, structure = _synthetic_tables()
    families = build_feature_families(labels, sequence, structure)
    family = families["structure_transfer_safe_plus_proxy"]
    assert set(BFACTOR_COLUMNS).isdisjoint(family.feature_columns)
    assert any(name.startswith("flex_proxy_") for name in family.feature_columns)
    assert "distance_to_n_terminus" not in family.feature_columns


def test_bfactor_and_proxy_diagnostic_families_exist() -> None:
    labels, sequence, structure = _synthetic_tables()
    families = build_feature_families(labels, sequence, structure)
    assert set(families["structure_bfactor_only"].feature_columns) == set(BFACTOR_COLUMNS)
    assert len(families["structure_proxy_only"].feature_columns) == 6
