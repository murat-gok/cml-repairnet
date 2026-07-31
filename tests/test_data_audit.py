import pandas as pd
import pytest

from cml_repairnet.quality import DataAuditError, audit_stage1_tables


def _tables(n: int = 33) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    sites = list(range(25, 25 + n))
    labels = pd.DataFrame(
        {
            "protein_id": ["P02769"] * n,
            "site_uniprot": sites,
            "residue_uniprot": ["K"] * n,
            "repair_published": [0.5] * n,
            "repair_from_means_clipped": [0.45] * n,
        }
    )
    sequence = labels[["protein_id", "site_uniprot"]].copy()
    structure = labels[["protein_id", "site_uniprot"]].copy()
    structure["structure_missing"] = 0
    return labels, sequence, structure


def test_stage1_audit_passes_complete_tables() -> None:
    report = audit_stage1_tables(*_tables())
    assert report["status"] == "pass"
    assert report["n_sites"] == 33
    assert report["published_vs_reconstructed_mae"] == pytest.approx(0.05)


def test_stage1_audit_rejects_missing_structure() -> None:
    labels, sequence, structure = _tables()
    structure.loc[0, "structure_missing"] = 1
    with pytest.raises(DataAuditError, match="missing"):
        audit_stage1_tables(labels, sequence, structure)
