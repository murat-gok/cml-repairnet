from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from Bio import SeqIO


def _read_figure_3b(path: Path) -> pd.DataFrame:
    raw = pd.read_excel(path, sheet_name="Figure 3B", header=None)
    data = raw.iloc[1:, :7].copy()
    data.columns = [
        "site_paper",
        "active_rep1",
        "active_rep2",
        "active_rep3",
        "control_rep1",
        "control_rep2",
        "control_rep3",
    ]
    return data.apply(pd.to_numeric, errors="coerce").dropna(subset=["site_paper"])


def _read_figure_3c(path: Path) -> pd.DataFrame:
    data = pd.read_excel(path, sheet_name="Figure 3C")
    data = data.iloc[:, :2].copy()
    data.columns = ["site_paper", "repair_published"]
    return data.apply(pd.to_numeric, errors="coerce").dropna(subset=["site_paper"])


def _read_bfactor(path: Path) -> pd.DataFrame:
    raw = pd.read_excel(path, sheet_name="Supplementary Figure 3-4")
    data = raw.iloc[:, :3].copy()
    data.columns = ["site_uniprot", "bfactor_source", "repair_percent_source"]
    for column in data.columns:
        data[column] = pd.to_numeric(data[column], errors="coerce")
    return data.dropna(subset=["site_uniprot", "bfactor_source"])


def _read_fasta_sequence(path: Path) -> str:
    records = list(SeqIO.parse(str(path), "fasta"))
    if len(records) != 1:
        raise ValueError(f"Expected one FASTA record in {path}; found {len(records)}")
    return str(records[0].seq)


def build_bsa_labels(
    source_xlsx: Path,
    output_csv: Path,
    fasta_path: Path | None = None,
    precursor_offset: int = 24,
    variance_floor: float = 1e-4,
) -> pd.DataFrame:
    """Reconstruct replicate-aware site-level repair labels from the article source data."""
    replicates = _read_figure_3b(source_xlsx)
    published = _read_figure_3c(source_xlsx)
    bfactor = _read_bfactor(source_xlsx)

    df = replicates.merge(published, on="site_paper", how="inner", validate="one_to_one")
    df["site_paper"] = df["site_paper"].astype(int)
    df["site_uniprot"] = df["site_paper"] + precursor_offset

    active_cols = [f"active_rep{i}" for i in range(1, 4)]
    control_cols = [f"control_rep{i}" for i in range(1, 4)]
    df["active_mean"] = df[active_cols].mean(axis=1)
    df["control_mean"] = df[control_cols].mean(axis=1)
    df["active_sd"] = df[active_cols].std(axis=1, ddof=1)
    df["control_sd"] = df[control_cols].std(axis=1, ddof=1)

    safe_control = df["control_mean"].replace(0, np.nan)
    df["repair_from_means_raw"] = 1.0 - (df["active_mean"] / safe_control)
    df["repair_from_means_clipped"] = df["repair_from_means_raw"].clip(0.0, 1.0)

    # Delta-method approximation for Var(1 - A/C), using replicate mean variances.
    var_a = (df["active_sd"] ** 2) / len(active_cols)
    var_c = (df["control_sd"] ** 2) / len(control_cols)
    c = safe_control
    a = df["active_mean"]
    df["repair_variance_approx"] = (var_a / (c**2)) + ((a**2) * var_c / (c**4))
    df["repair_variance_approx"] = df["repair_variance_approx"].fillna(1.0)
    df["label_weight_raw"] = 1.0 / (df["repair_variance_approx"] + variance_floor)
    median_weight = float(df["label_weight_raw"].median())
    df["label_weight"] = df["label_weight_raw"] / median_weight

    df = df.merge(
        bfactor[["site_uniprot", "bfactor_source"]],
        on="site_uniprot",
        how="left",
        validate="one_to_one",
    )

    if fasta_path is not None and fasta_path.exists():
        sequence = _read_fasta_sequence(fasta_path)
        df["residue_uniprot"] = [
            sequence[position - 1] if 1 <= position <= len(sequence) else "?"
            for position in df["site_uniprot"]
        ]
        invalid = df.loc[
            df["residue_uniprot"] != "K", ["site_paper", "site_uniprot", "residue_uniprot"]
        ]
        if not invalid.empty:
            raise ValueError(
                "Site-number mapping failed; expected lysine at every mapped UniProt position:\n"
                + invalid.to_string(index=False)
            )

    df["protein_id"] = "P02769"
    df["protein_name"] = "Bovine serum albumin"
    df["organism"] = "Bos taurus"
    df["assay_condition"] = "CMLase vs inactive control; overnight; n=3"
    df["label_source"] = "Trabosh et al. 2026, Figure 3B/3C source data"

    ordered = [
        "protein_id",
        "protein_name",
        "organism",
        "site_paper",
        "site_uniprot",
        "residue_uniprot",
        "repair_published",
        "repair_from_means_raw",
        "repair_from_means_clipped",
        "repair_variance_approx",
        "label_weight",
        "active_mean",
        "active_sd",
        "control_mean",
        "control_sd",
        *active_cols,
        *control_cols,
        "bfactor_source",
        "assay_condition",
        "label_source",
    ]
    ordered = [column for column in ordered if column in df.columns]
    df = df[ordered].sort_values("site_paper").reset_index(drop=True)

    output_csv.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_csv, index=False)
    return df
