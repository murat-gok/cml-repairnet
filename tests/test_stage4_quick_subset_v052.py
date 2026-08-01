from __future__ import annotations

import pandas as pd

from cml_repairnet.atlas.sources import deterministic_quick_subset


def test_quick_subset_accepts_pre_annotation_accession_schema() -> None:
    frame = pd.DataFrame(
        {
            "candidate_id": [
                "human:P11111:K10",
                "human:P11111-2:K20",
                "mouse:Q22222:K30",
                "mouse:Q33333:K40",
            ],
            "species": ["human", "human", "mouse", "mouse"],
            "uniprot_accession_raw": ["P11111", "P11111-2", "Q22222", "Q33333"],
            "site_uniprot": [10, 20, 30, 40],
            "prm_validated": [True, False, True, False],
            "age_associated": [False, True, False, True],
            "evidence_score_raw": [3.0, 2.0, 2.5, 1.0],
        }
    )
    selected = deterministic_quick_subset(frame, max_proteins=2, max_sites=3)
    assert len(selected) <= 3
    assert "uniprot_accession_raw" in selected.columns
    assert "uniprot_accession" not in selected.columns
    assert "_quick_protein_key" not in selected.columns
    assert set(selected["species"]) == {"human", "mouse"}


def test_quick_subset_accepts_post_annotation_accession_schema() -> None:
    frame = pd.DataFrame(
        {
            "candidate_id": ["human:P11111:K10", "mouse:Q22222:K20"],
            "species": ["human", "mouse"],
            "uniprot_accession": ["P11111", "Q22222"],
            "site_uniprot": [10, 20],
            "prm_validated": [True, False],
            "age_associated": [False, True],
            "evidence_score_raw": [2.0, 1.0],
        }
    )
    selected = deterministic_quick_subset(frame, max_proteins=25, max_sites=60)
    pd.testing.assert_frame_equal(selected.reset_index(drop=True), frame.reset_index(drop=True))
