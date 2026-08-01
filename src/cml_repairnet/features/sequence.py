from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from Bio import SeqIO

AA20 = list("ACDEFGHIKLMNPQRSTVWY")
HYDROPATHY = {
    "A": 1.8,
    "C": 2.5,
    "D": -3.5,
    "E": -3.5,
    "F": 2.8,
    "G": -0.4,
    "H": -3.2,
    "I": 4.5,
    "K": -3.9,
    "L": 3.8,
    "M": 1.9,
    "N": -3.5,
    "P": -1.6,
    "Q": -3.5,
    "R": -4.5,
    "S": -0.8,
    "T": -0.7,
    "V": 4.2,
    "W": -0.9,
    "Y": -1.3,
}
AROMATIC = set("FWY")
BULKY_HYDROPHOBIC = set("FWLVIY")
ACIDIC = set("DE")
BASIC = set("KRH")
POLAR = set("STNQCY")


def read_single_fasta(path: Path) -> str:
    records = list(SeqIO.parse(str(path), "fasta"))
    if len(records) != 1:
        raise ValueError(f"Expected one FASTA record in {path}; found {len(records)}")
    return str(records[0].seq)


def _window(sequence: str, center_1based: int, radius: int) -> tuple[str, int]:
    center = center_1based - 1
    start = max(0, center - radius)
    end = min(len(sequence), center + radius + 1)
    left_pad = max(0, radius - center)
    right_pad = max(0, center + radius + 1 - len(sequence))
    return ("X" * left_pad) + sequence[start:end] + ("X" * right_pad), center


def _indicator(aa: str, group: set[str]) -> float:
    return float(aa in group)


def build_sequence_features(
    labels_csv: Path,
    fasta_path: Path,
    output_csv: Path,
    radius: int = 10,
) -> pd.DataFrame:
    labels = pd.read_csv(labels_csv)
    sequence = read_single_fasta(fasta_path)
    rows: list[dict[str, object]] = []

    for record in labels.itertuples(index=False):
        site = int(record.site_uniprot)
        window, center_zero = _window(sequence, site, radius)
        if sequence[center_zero] != "K":
            raise ValueError(f"Mapped site {site} is {sequence[center_zero]}, not lysine.")

        flank = window[:radius] + window[radius + 1 :]
        valid = [aa for aa in flank if aa in AA20]
        row: dict[str, object] = {
            "protein_id": record.protein_id,
            "site_uniprot": site,
            "sequence_window": window,
            "distance_to_n_terminus": site - 1,
            "distance_to_c_terminus": len(sequence) - site,
            "window_mean_hydropathy": float(np.mean([HYDROPATHY[a] for a in valid])),
            "window_aromatic_fraction": float(np.mean([a in AROMATIC for a in valid])),
            "window_bulky_hydrophobic_fraction": float(
                np.mean([a in BULKY_HYDROPHOBIC for a in valid])
            ),
            "window_acidic_fraction": float(np.mean([a in ACIDIC for a in valid])),
            "window_basic_fraction": float(np.mean([a in BASIC for a in valid])),
            "window_polar_fraction": float(np.mean([a in POLAR for a in valid])),
            "window_cysteine_count": flank.count("C"),
            "window_proline_count": flank.count("P"),
            "window_glycine_count": flank.count("G"),
        }
        for aa in AA20:
            row[f"composition_{aa}"] = flank.count(aa) / max(1, len(valid))

        for offset in range(-5, 6):
            if offset == 0:
                continue
            aa = (
                sequence[center_zero + offset] if 0 <= center_zero + offset < len(sequence) else "X"
            )
            prefix = f"pos_{offset:+d}"
            row[f"{prefix}_hydropathy"] = HYDROPATHY.get(aa, 0.0)
            row[f"{prefix}_aromatic"] = _indicator(aa, AROMATIC)
            row[f"{prefix}_bulky_hydrophobic"] = _indicator(aa, BULKY_HYDROPHOBIC)
            row[f"{prefix}_acidic"] = _indicator(aa, ACIDIC)
            row[f"{prefix}_basic"] = _indicator(aa, BASIC)
            row[f"{prefix}_cysteine"] = float(aa == "C")
            row[f"{prefix}_proline"] = float(aa == "P")
        rows.append(row)

    features = pd.DataFrame(rows)
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    features.to_csv(output_csv, index=False)
    return features
