from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

SPECIES_TO_ORGANISM = {"human": "hsapiens", "mouse": "mmusculus"}
SPECIES_TO_TAXON = {"human": 9606, "mouse": 10090}

REQUIRED_ATLAS_COLUMNS = {
    "candidate_id",
    "species",
    "organism_id",
    "uniprot_accession",
    "site_uniprot",
    "gene_names",
    "protein_name",
    "repairability_extra_trees",
    "repairability_lightgbm",
    "repairability_ridge",
    "consensus_rank_percentile",
    "evidence_score_raw",
    "transfer_reliability_score",
    "weighted_priority_score",
    "applicability_status",
    "prm_validated",
    "age_associated",
    "ortholog_overlap",
    "model_disagreement",
    "ood_percentile",
    "prediction_interval_width",
}


@dataclass(frozen=True)
class Stage4Inputs:
    atlas_all: Path
    atlas_in_domain: Path
    pareto: Path
    portfolio: Path
    watchlist: Path


def first_gene_symbol(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    text = str(value).strip()
    if not text:
        return ""
    # UniProt gene_names starts with the primary symbol and may then list synonyms.
    return text.split()[0]


def canonical_accession(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    return str(value).strip().split("-")[0]


def validate_atlas(frame: pd.DataFrame, *, name: str = "atlas") -> None:
    missing = sorted(REQUIRED_ATLAS_COLUMNS - set(frame.columns))
    if missing:
        raise ValueError(f"{name} is missing required columns: {missing}")
    if frame["candidate_id"].duplicated().any():
        duplicates = frame.loc[frame["candidate_id"].duplicated(), "candidate_id"].head().tolist()
        raise ValueError(f"{name} contains duplicate candidate_id values: {duplicates}")
    unknown_species = sorted(set(frame["species"].dropna()) - set(SPECIES_TO_ORGANISM))
    if unknown_species:
        raise ValueError(f"{name} contains unsupported species: {unknown_species}")


def load_atlas(path: Path, *, name: str | None = None) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(path)
    frame = pd.read_csv(path)
    validate_atlas(frame, name=name or path.name)
    frame = frame.copy()
    frame["gene_symbol"] = frame["gene_names"].map(first_gene_symbol)
    frame["canonical_accession"] = frame["uniprot_accession"].map(canonical_accession)
    frame["site_label"] = (
        frame["gene_symbol"].fillna("") + " K" + frame["site_uniprot"].astype(int).astype(str)
    )
    return frame


def unique_gene_symbols(frame: pd.DataFrame) -> list[str]:
    values = [str(value).strip() for value in frame.get("gene_symbol", pd.Series(dtype=str))]
    return sorted({value for value in values if value})


def unique_accessions(frame: pd.DataFrame) -> list[str]:
    series = frame.get("canonical_accession", frame.get("uniprot_accession", pd.Series(dtype=str)))
    return sorted({canonical_accession(value) for value in series if canonical_accession(value)})


def select_species(frame: pd.DataFrame, species: str) -> pd.DataFrame:
    if species not in SPECIES_TO_ORGANISM:
        raise ValueError(f"Unsupported species: {species}")
    return frame.loc[frame["species"] == species].copy()


def evidence_rich_mask(frame: pd.DataFrame) -> pd.Series:
    columns = [
        column
        for column in ("prm_validated", "age_associated", "ortholog_overlap")
        if column in frame
    ]
    if not columns:
        return pd.Series(False, index=frame.index)
    return frame[columns].fillna(False).astype(bool).any(axis=1)


def top_fraction_mask(frame: pd.DataFrame, column: str, fraction: float = 0.25) -> pd.Series:
    if not 0 < fraction <= 1:
        raise ValueError("fraction must be in (0, 1]")
    values = pd.to_numeric(frame[column], errors="coerce")
    valid = values.dropna()
    if valid.empty:
        return pd.Series(False, index=frame.index)
    threshold = valid.quantile(1 - fraction)
    return values >= threshold


def require_nonempty(values: Iterable[str], *, label: str) -> list[str]:
    result = [str(value).strip() for value in values if str(value).strip()]
    if not result:
        raise ValueError(f"No usable values for {label}")
    return result
