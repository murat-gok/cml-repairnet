from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from Bio.PDB import PDBParser, ShrakeRupley
from Bio.PDB.Polypeptide import protein_letters_3to1

HYDROPHOBIC = set("AILMFWVY")
AROMATIC = set("FWY")


def _residue_letter(residue_name: str) -> str:
    return protein_letters_3to1.get(residue_name.upper(), "X")


def _resolve_pdb_site(record: object, pdb_site_column: str) -> int:
    """Resolve the coordinate-system position used by the experimental structure.

    The CMLase paper reports mature BSA positions (for example K51), while UniProt P02769
    contains a 24-residue signal/propeptide precursor. PDB 4F5S is numbered in the mature
    583-residue coordinate system. Therefore the default must be ``site_paper`` rather than
    ``site_uniprot``. A fallback keeps unit tests and non-BSA inputs convenient.
    """
    if hasattr(record, pdb_site_column):
        return int(getattr(record, pdb_site_column))
    if hasattr(record, "site_pdb"):
        return int(record.site_pdb)
    if hasattr(record, "site_uniprot"):
        return int(record.site_uniprot)
    raise ValueError("No structure-coordinate site column was found in the labels table.")


def build_structure_features(
    labels_csv: Path,
    pdb_path: Path,
    output_csv: Path,
    chain_id: str = "A",
    lysine_max_asa: float = 205.0,
    contact_radii: tuple[float, ...] = (8.0, 10.0, 12.0),
    pdb_site_column: str = "site_paper",
) -> pd.DataFrame:
    """Build residue-level structural descriptors from a PDB coordinate file.

    Output rows retain the UniProt coordinate as the merge key and additionally report the
    actual PDB/mature-protein coordinate used for structure lookup.
    """
    labels = pd.read_csv(labels_csv)
    structure = PDBParser(QUIET=True).get_structure("protein", str(pdb_path))
    model = next(structure.get_models())
    if chain_id not in model:
        raise KeyError(f"Chain {chain_id!r} not found in {pdb_path}")
    chain = model[chain_id]

    ShrakeRupley().compute(model, level="R")

    residues = [residue for residue in chain if residue.id[0] == " " and "CA" in residue]
    if not residues:
        raise ValueError(f"No standard amino-acid residues with CA atoms found in {pdb_path}")
    residue_by_number = {int(residue.id[1]): residue for residue in residues}
    all_ca = np.vstack([residue["CA"].coord for residue in residues])
    centroid = all_ca.mean(axis=0)

    rows: list[dict[str, object]] = []
    for record in labels.itertuples(index=False):
        site_uniprot = int(record.site_uniprot)
        site_pdb = _resolve_pdb_site(record, pdb_site_column)
        residue = residue_by_number.get(site_pdb)
        if residue is None:
            rows.append(
                {
                    "protein_id": record.protein_id,
                    "site_uniprot": site_uniprot,
                    "site_pdb": site_pdb,
                    "structure_missing": 1,
                }
            )
            continue
        if _residue_letter(residue.resname) != "K":
            raise ValueError(
                f"PDB residue {site_pdb} is {residue.resname}, not lysine. "
                "Check precursor/mature-protein numbering."
            )

        ca = np.asarray(residue["CA"].coord, dtype=float)
        distances = np.linalg.norm(all_ca - ca, axis=1)
        neighbor_residues = [
            other
            for other in residues
            if other is not residue and np.linalg.norm(other["CA"].coord - ca) <= 8.0
        ]
        neighbor_letters = [_residue_letter(other.resname) for other in neighbor_residues]
        bfactor = float(residue["CA"].bfactor)
        neighbor_bfactors = [float(other["CA"].bfactor) for other in neighbor_residues]
        sasa = float(getattr(residue, "sasa", np.nan))

        row: dict[str, object] = {
            "protein_id": record.protein_id,
            "site_uniprot": site_uniprot,
            "site_pdb": site_pdb,
            "structure_missing": 0,
            "pdb_chain": chain_id,
            "ca_x": ca[0],
            "ca_y": ca[1],
            "ca_z": ca[2],
            "bfactor_ca": bfactor,
            "residue_sasa": sasa,
            "relative_sasa_lysine": min(sasa / lysine_max_asa, 1.5),
            "distance_to_protein_centroid": float(np.linalg.norm(ca - centroid)),
            "neighbor_mean_bfactor_8a": (
                float(np.mean(neighbor_bfactors)) if neighbor_bfactors else np.nan
            ),
            "neighbor_hydrophobic_fraction_8a": (
                float(np.mean([aa in HYDROPHOBIC for aa in neighbor_letters]))
                if neighbor_letters
                else 0.0
            ),
            "neighbor_cysteine_count_8a": neighbor_letters.count("C"),
            "neighbor_aromatic_count_8a": sum(aa in AROMATIC for aa in neighbor_letters),
        }
        for radius in contact_radii:
            row[f"ca_contact_count_{radius:g}a"] = int(
                np.sum((distances > 0) & (distances <= radius))
            )
        rows.append(row)

    features = pd.DataFrame(rows)
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    features.to_csv(output_csv, index=False)
    return features
