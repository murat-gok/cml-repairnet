from pathlib import Path

import pandas as pd

from cml_repairnet.features.structure import build_structure_features


PDB_TEXT = """\
ATOM      1  N   LYS A  41       0.000   0.000   0.000  1.00 20.00           N  
ATOM      2  CA  LYS A  41       1.500   0.000   0.000  1.00 21.00           C  
ATOM      3  C   LYS A  41       2.200   1.300   0.000  1.00 22.00           C  
ATOM      4  O   LYS A  41       1.700   2.400   0.000  1.00 23.00           O  
ATOM      5  CB  LYS A  41       2.000  -0.900   1.100  1.00 21.00           C  
ATOM      6  N   ALA A  42       3.400   1.200   0.000  1.00 24.00           N  
ATOM      7  CA  ALA A  42       4.200   2.400   0.000  1.00 25.00           C  
ATOM      8  C   ALA A  42       5.700   2.100   0.000  1.00 26.00           C  
ATOM      9  O   ALA A  42       6.400   3.000   0.000  1.00 27.00           O  
ATOM     10  CB  ALA A  42       3.800   3.300   1.200  1.00 25.00           C  
TER      11      ALA A  42                                                      
END                                                                             
"""


def test_structure_features_use_mature_pdb_numbering(tmp_path: Path) -> None:
    labels = tmp_path / "labels.csv"
    pdb = tmp_path / "mini.pdb"
    output = tmp_path / "features.csv"
    # Paper/mature BSA position 41 corresponds to UniProt precursor position 65.
    pd.DataFrame(
        {"protein_id": ["P02769"], "site_paper": [41], "site_uniprot": [65]}
    ).to_csv(labels, index=False)
    pdb.write_text(PDB_TEXT, encoding="utf-8")

    frame = build_structure_features(labels, pdb, output)

    assert len(frame) == 1
    assert frame.loc[0, "site_uniprot"] == 65
    assert frame.loc[0, "site_pdb"] == 41
    assert frame.loc[0, "structure_missing"] == 0
    assert frame.loc[0, "bfactor_ca"] == 21.0
    assert frame.loc[0, "ca_contact_count_8a"] == 1
    assert output.exists()
