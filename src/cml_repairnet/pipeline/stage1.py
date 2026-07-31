from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd

from cml_repairnet.data.download import download_public_bundle
from cml_repairnet.data.prepare_labels import build_bsa_labels
from cml_repairnet.features.sequence import build_sequence_features
from cml_repairnet.features.structure import build_structure_features
from cml_repairnet.quality import audit_stage1_tables, write_audit_report


def run_stage1(config: dict[str, Any], *, download: bool = True) -> dict[str, Any]:
    """Run the reproducible data-foundation stage of CML-RepairNet."""
    paths = config["paths"]
    labels_config = config["labels"]
    features_config = config["features"]

    article_source = Path(paths["article_source_xlsx"])
    fasta = Path(paths["bsa_fasta"])
    pdb = Path(paths["bsa_pdb"])
    labels_path = Path(paths["labels"])
    sequence_path = Path(paths["sequence_features"])
    structure_path = Path(paths["structure_features"])
    results_dir = Path(paths["results_dir"])

    required_raw = [article_source, fasta, pdb]
    if download and any(not path.exists() for path in required_raw):
        download_public_bundle(Path("data/raw"), include_pride_processed=False)

    missing = [str(path) for path in required_raw if not path.exists()]
    if missing:
        raise FileNotFoundError(
            "Stage 1 requires the following raw files. Run `cml-repairnet download` first: "
            + ", ".join(missing)
        )

    labels = build_bsa_labels(
        article_source,
        labels_path,
        fasta_path=fasta,
        precursor_offset=int(labels_config["precursor_offset"]),
        variance_floor=float(labels_config["variance_floor"]),
    )
    sequence = build_sequence_features(
        labels_path,
        fasta,
        sequence_path,
        radius=int(features_config["sequence_window_radius"]),
    )
    structure = build_structure_features(
        labels_path,
        pdb,
        structure_path,
        chain_id=str(features_config.get("pdb_chain", "A")),
        lysine_max_asa=float(features_config["lysine_max_asa_angstrom2"]),
        contact_radii=tuple(
            float(value) for value in features_config["contact_radii_angstrom"]
        ),
        pdb_site_column=str(features_config.get("pdb_site_column", "site_paper")),
    )

    report = audit_stage1_tables(labels, sequence, structure)
    write_audit_report(
        report,
        results_dir / "stage1_data_audit.json",
        results_dir / "stage1_data_audit.md",
    )
    return report
