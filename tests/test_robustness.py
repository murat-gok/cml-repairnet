import numpy as np
import pandas as pd

from cml_repairnet.evaluation.robustness import build_feature_view, paired_error_comparison
from cml_repairnet.models.baselines import _numeric_feature_columns


def test_site_pdb_is_excluded_from_features():
    frame = pd.DataFrame({
        "site_pdb": [1], "site_uniprot": [25], "ca_x": [0.0], "signal": [1.0]
    })
    assert _numeric_feature_columns(frame) == ["signal"]


def test_sequence_view_keeps_coordinates_only_for_grouping():
    labels = pd.DataFrame({
        "protein_id": ["P"], "site_uniprot": [10], "site_paper": [8],
        "repair_published": [0.5],
    })
    sequence = pd.DataFrame({
        "protein_id": ["P"], "site_uniprot": [10], "sequence_signal": [1.0]
    })
    structure = pd.DataFrame({
        "protein_id": ["P"], "site_uniprot": [10], "ca_x": [1.0], "ca_y": [2.0],
        "ca_z": [3.0], "structure_missing": [0], "residue_sasa": [50.0],
    })
    view = build_feature_view(labels, sequence, structure, "sequence")
    assert "sequence_signal" in view.columns
    assert "ca_x" in view.columns
    assert "residue_sasa" not in view.columns


def test_paired_error_comparison_detects_better_model():
    rows = []
    y = np.array([0.0, 0.25, 0.5, 0.75, 1.0])
    for idx, value in enumerate(y):
        common = {
            "protein_id": "P", "site_uniprot": idx + 1, "site_paper": idx + 1,
            "spatial_group": idx % 2, "y_true": value,
        }
        rows.append({**common, "model": "candidate", "y_pred": value + 0.01})
        rows.append({**common, "model": "dummy_mean", "y_pred": 0.5})
    result = paired_error_comparison(
        pd.DataFrame(rows), model="candidate", bootstrap_iterations=500,
        permutation_iterations=1000, seed=1,
    )
    assert result.absolute_improvement > 0
    assert result.relative_improvement_pct > 0
    assert result.n_sites == 5
