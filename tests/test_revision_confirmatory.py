from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from cml_repairnet.revision.confirmatory import (
    run_multi_partition_permutation_control,
    run_nested_model_selection,
    run_primary_cluster_bootstrap,
)


def _tables(n: int = 24):
    rng = np.random.default_rng(123)
    site = np.arange(n)
    groups = np.repeat(np.arange(4), n // 4)
    signal = groups / 5 + rng.normal(0, 0.02, n)
    labels = pd.DataFrame({
        "protein_id": ["P02769"] * n,
        "site_uniprot": site + 25,
        "site_paper": site + 1,
        "spatial_group": groups,
        "repair_published": np.clip(0.15 + signal, 0, 1),
        "repair_from_means_clipped": np.clip(0.15 + signal, 0, 1),
        "label_weight": np.ones(n),
    })
    sequence = pd.DataFrame({
        "protein_id": ["P02769"] * n,
        "site_uniprot": site + 25,
        "distance_to_n_terminus": site,
        "window_mean_hydropathy": rng.normal(size=n),
        "window_cysteine_count": rng.integers(0, 3, n),
    })
    structure = pd.DataFrame({
        "protein_id": ["P02769"] * n,
        "site_uniprot": site + 25,
        "site_pdb": site + 1,
        "structure_missing": np.zeros(n),
        "ca_x": groups * 20 + rng.normal(0, 0.1, n),
        "ca_y": rng.normal(0, 0.1, n),
        "ca_z": rng.normal(0, 0.1, n),
        "bfactor_ca": 20 + signal * 10,
        "neighbor_mean_bfactor_8a": 21 + signal * 9,
        "residue_sasa": 30 + signal * 25,
        "relative_sasa_lysine": np.clip(0.2 + signal, 0, 1.2),
        "distance_to_protein_centroid": 10 - signal * 3,
        "neighbor_hydrophobic_fraction_8a": np.clip(0.3 + signal / 2, 0, 1),
        "neighbor_cysteine_count_8a": rng.integers(0, 2, n),
        "neighbor_aromatic_count_8a": rng.integers(0, 3, n),
        "ca_contact_count_8a": np.maximum(2, (12 - signal * 5).astype(int)),
        "ca_contact_count_10a": np.maximum(5, (18 - signal * 5).astype(int)),
        "ca_contact_count_12a": np.maximum(8, (24 - signal * 5).astype(int)),
    })
    return labels, sequence, structure


def test_multi_partition_permutation_runs(tmp_path: Path) -> None:
    labels, sequence, structure = _tables()
    out = run_multi_partition_permutation_control(
        labels, sequence, structure, tmp_path,
        feature_set="structure_transfer_safe",
        model_name="ridge",
        iterations=3,
        group_counts=(4,),
        random_partition_count=1,
        n_jobs=1,
    )
    assert len(out["partition_summary"]) == 2
    assert len(out["aggregate_draws"]) == 3
    assert np.isfinite(out["aggregate_summary"].loc[0, "mean_delta_mae"])


def test_nested_model_selection_runs(tmp_path: Path) -> None:
    labels, sequence, structure = _tables()
    out = run_nested_model_selection(
        labels, sequence, structure, tmp_path,
        feature_set="structure_transfer_safe",
        model_names=("ridge", "extra_trees"),
        fixed_model="extra_trees",
        group_counts=(4,),
        random_partition_count=0,
        n_jobs=1,
        bootstrap_iterations=50,
    )
    assert len(out["partition_metrics"]) == 1
    assert out["predictions"]["nested_y_pred"].notna().all()
    assert out["fold_selection"]["selected_model"].isin(["ridge", "extra_trees"]).all()


def test_primary_cluster_bootstrap_runs(tmp_path: Path) -> None:
    rows = []
    rng = np.random.default_rng(4)
    for i in range(24):
        y = i / 30
        group = str(i // 6)
        rows.append({
            "protein_id": "P02769", "site_uniprot": i + 25,
            "partition_id": "ward_g4", "outer_group": group, "outer_fold": i // 6 + 1,
            "feature_set": "structure_transfer_safe", "model": "extra_trees",
            "y_true": y, "y_pred": y + rng.normal(0, 0.03), "absolute_error": 0.0,
        })
        rows.append({
            "protein_id": "P02769", "site_uniprot": i + 25,
            "partition_id": "ward_g4", "outer_group": group, "outer_fold": i // 6 + 1,
            "feature_set": "null_mean", "model": "dummy_mean",
            "y_true": y, "y_pred": 0.4, "absolute_error": 0.0,
        })
    path = tmp_path / "pred.csv"
    pd.DataFrame(rows).to_csv(path, index=False)
    out = run_primary_cluster_bootstrap(
        path, tmp_path, iterations=20, sign_flip_iterations=100
    )
    assert out["observed"].loc[0, "delta_mae_null_minus_primary"] > 0
    assert set(out["summary"]["method"]) == {"site", "spatial_group_cluster"}
