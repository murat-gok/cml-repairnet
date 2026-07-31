from pathlib import Path

import pandas as pd

from cml_repairnet.features.sequence import build_sequence_features


def test_sequence_feature_center_is_lysine(tmp_path: Path) -> None:
    fasta = tmp_path / "test.fasta"
    fasta.write_text(">p\nAAAAKAAAA\n", encoding="utf-8")
    labels = tmp_path / "labels.csv"
    pd.DataFrame({"protein_id": ["p"], "site_uniprot": [5]}).to_csv(labels, index=False)
    output = tmp_path / "features.csv"
    frame = build_sequence_features(labels, fasta, output, radius=2)
    assert frame.loc[0, "sequence_window"] == "AAKAA"
    assert output.exists()
