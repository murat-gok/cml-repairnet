from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pandas as pd
import yaml


@dataclass(frozen=True)
class Stage6Inputs:
    atlas_all: Path
    atlas_in_domain: Path
    pareto: Path
    portfolio: Path
    cross_species_pairs: Path
    enrichment: Path
    string_edges: Path
    string_nodes: Path
    network_manifest: Path
    stage5_modules: Path


@dataclass(frozen=True)
class Stage6Outputs:
    processed: Path
    results: Path
    figures: Path


def input_paths(root: Path) -> Stage6Inputs:
    return Stage6Inputs(
        atlas_all=root / "results/stage4/atlas/cml_repairability_atlas_all_sites.csv",
        atlas_in_domain=root / "results/stage4/atlas/cml_repairability_atlas_in_domain.csv",
        pareto=root / "results/stage4/atlas/cml_repairability_pareto_front.csv",
        portfolio=root / "results/stage4/atlas/cml_repairability_priority_portfolio.csv",
        cross_species_pairs=root / "data/processed/stage5/stage5_cross_species_pairs.csv",
        enrichment=root / "data/processed/stage5/stage5_gprofiler_enrichment_all.csv",
        string_edges=root / "data/processed/stage5/stage5_string_edges.csv",
        string_nodes=root / "data/processed/stage5/stage5_string_nodes.csv",
        network_manifest=root / "results/stage5/stage5_network_manifest.csv",
        stage5_modules=root / "data/processed/stage5/stage5_candidate_module_assignments.csv",
    )


def output_paths(root: Path) -> Stage6Outputs:
    processed = root / "data/processed/stage6"
    results = root / "results/stage6"
    figures = results / "figures"
    processed.mkdir(parents=True, exist_ok=True)
    results.mkdir(parents=True, exist_ok=True)
    figures.mkdir(parents=True, exist_ok=True)
    return Stage6Outputs(processed=processed, results=results, figures=figures)


def load_config(root: Path) -> dict[str, Any]:
    path = root / "configs/stage6.yaml"
    if not path.exists():
        raise FileNotFoundError(path)
    with path.open(encoding="utf-8") as handle:
        config = yaml.safe_load(handle)
    if not isinstance(config, dict):
        raise ValueError("Stage 6 config must be a mapping")
    return config


def read_csv(path: Path, *, required: tuple[str, ...] = ()) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(path)
    frame = pd.read_csv(path)
    missing = [column for column in required if column not in frame.columns]
    if missing:
        raise ValueError(f"{path.name} is missing required columns: {missing}")
    return frame


def first_gene(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    return str(value).strip().split()[0]


def add_labels(frame: pd.DataFrame) -> pd.DataFrame:
    result = frame.copy()
    if "gene_symbol" not in result:
        result["gene_symbol"] = result.get("gene_names", pd.Series(index=result.index, dtype=object)).map(first_gene)
    if "site_label" not in result and "site_uniprot" in result:
        site = pd.to_numeric(result["site_uniprot"], errors="coerce").astype("Int64").astype(str)
        result["site_label"] = result["gene_symbol"].fillna("") + " K" + site
    return result
