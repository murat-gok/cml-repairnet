from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd


class DataAuditError(ValueError):
    """Raised when a hard scientific data-integrity requirement is violated."""


def audit_stage1_tables(
    labels: pd.DataFrame,
    sequence_features: pd.DataFrame,
    structure_features: pd.DataFrame,
    *,
    expected_sites: int = 33,
) -> dict[str, Any]:
    """Audit the first CML-RepairNet benchmark tables.

    Hard failures protect the project from silent numbering errors, duplicate sites, invalid
    targets and incomplete structure mapping. Descriptive discrepancies are reported rather than
    hidden; in particular the published repair value and the ratio reconstructed from replicates
    are intentionally kept as separate quantities.
    """
    keys = ["protein_id", "site_uniprot"]
    for name, frame in {
        "labels": labels,
        "sequence_features": sequence_features,
        "structure_features": structure_features,
    }.items():
        missing = [column for column in keys if column not in frame.columns]
        if missing:
            raise DataAuditError(f"{name} is missing key columns: {missing}")
        if frame.duplicated(keys).any():
            duplicates = frame.loc[frame.duplicated(keys, keep=False), keys]
            raise DataAuditError(f"{name} contains duplicate sites:\n{duplicates}")

    if len(labels) != expected_sites:
        raise DataAuditError(f"Expected {expected_sites} labelled sites; found {len(labels)}")
    if set(labels["residue_uniprot"].dropna().astype(str)) != {"K"}:
        raise DataAuditError("All mapped UniProt benchmark residues must be lysine (K).")

    target = pd.to_numeric(labels["repair_published"], errors="coerce")
    if target.isna().any() or not target.between(0.0, 1.0).all():
        raise DataAuditError("repair_published must be finite and bounded in [0, 1].")

    labels_keys = set(map(tuple, labels[keys].to_numpy()))
    sequence_keys = set(map(tuple, sequence_features[keys].to_numpy()))
    structure_keys = set(map(tuple, structure_features[keys].to_numpy()))
    if labels_keys != sequence_keys:
        raise DataAuditError("Sequence-feature keys do not match label keys.")
    if labels_keys != structure_keys:
        raise DataAuditError("Structure-feature keys do not match label keys.")

    missing_structure = int(
        pd.to_numeric(structure_features["structure_missing"], errors="coerce").fillna(1).sum()
    )
    if missing_structure:
        raise DataAuditError(
            f"{missing_structure} benchmark sites are missing in the structure mapping."
        )

    reconstructed = pd.to_numeric(labels["repair_from_means_clipped"], errors="coerce")
    delta = (target - reconstructed).abs()
    report: dict[str, Any] = {
        "status": "pass",
        "n_sites": len(labels),
        "n_sequence_columns": int(sequence_features.shape[1]),
        "n_structure_columns": int(structure_features.shape[1]),
        "target_min": float(target.min()),
        "target_max": float(target.max()),
        "target_mean": float(target.mean()),
        "published_vs_reconstructed_mae": float(delta.mean()),
        "published_vs_reconstructed_max_abs_delta": float(delta.max()),
        "largest_delta_site_uniprot": int(labels.loc[delta.idxmax(), "site_uniprot"]),
        "structure_missing_sites": missing_structure,
        "pdb_numbering_column": "site_paper",
        "notes": [
            "Published Figure 3C repair remains the primary target.",
            "Replicate-derived repair is retained as a sensitivity target, not overwritten.",
            "PDB 4F5S uses mature BSA numbering; UniProt P02769 includes a 24-residue precursor.",
        ],
    }
    return report


def write_audit_report(report: dict[str, Any], json_path: Path, markdown_path: Path) -> None:
    json_path.parent.mkdir(parents=True, exist_ok=True)
    markdown_path.parent.mkdir(parents=True, exist_ok=True)
    json_path.write_text(json.dumps(report, indent=2), encoding="utf-8")

    lines = [
        "# CML-RepairNet Stage 1 Data Audit",
        "",
        f"**Status:** {report['status'].upper()}",
        "",
        f"- Labelled sites: {report['n_sites']}",
        f"- Sequence table columns: {report['n_sequence_columns']}",
        f"- Structure table columns: {report['n_structure_columns']}",
        f"- Published target range: {report['target_min']:.4f}–{report['target_max']:.4f}",
        f"- Published vs replicate-reconstructed MAE: "
        f"{report['published_vs_reconstructed_mae']:.4f}",
        f"- Largest absolute discrepancy: "
        f"{report['published_vs_reconstructed_max_abs_delta']:.4f} "
        f"(UniProt site {report['largest_delta_site_uniprot']})",
        f"- Missing structure sites: {report['structure_missing_sites']}",
        "",
        "## Interpretation safeguards",
        "",
    ]
    lines.extend(f"- {note}" for note in report["notes"])
    markdown_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
