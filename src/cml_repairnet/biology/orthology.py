from __future__ import annotations

import re
from typing import Any

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from cml_repairnet.biology.clients import GProfilerClient, OfflineCacheMiss
from cml_repairnet.biology.schema import first_gene_symbol, select_species, unique_gene_symbols

ENSEMBL_ID = re.compile(r"^ENS[A-Z]*G\d+(?:\.\d+)?$", re.IGNORECASE)

ORTHOLOGY_COLUMNS = [
    "source_species",
    "target_species",
    "source_gene",
    "target_gene",
    "source_ensembl",
    "target_ensembl",
    "description",
    "namespaces",
]

PAIR_SUMMARY_COLUMNS = [
    "mapped_portfolio_pairs",
    "unique_human_genes",
    "unique_mouse_genes",
    "spearman_repairability",
    "spearman_priority",
    "median_repairability_abs_difference",
    "median_priority_abs_difference",
]


def _text(value: Any) -> str | None:
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return None
    text = str(value).strip()
    return text or None


def _target_symbol(item: dict[str, Any]) -> str | None:
    """Return the target-organism gene symbol from a g:Orth result.

    g:Orth uses ``converted`` for the input/source Ensembl identifier,
    ``ortholog_ensg`` for the target Ensembl identifier, and ``name`` for the
    target-organism gene symbol. Older or mocked responses may put a symbol in
    ``ortholog`` or ``converted``; those fallbacks are accepted only when they
    are not Ensembl-like identifiers.
    """
    for key in ("name", "ortholog_name", "ortholog"):
        value = _text(item.get(key))
        if value and not ENSEMBL_ID.match(value):
            return value
    converted = _text(item.get("converted"))
    if converted and not ENSEMBL_ID.match(converted):
        return converted
    return None


def flatten_ortholog_response(
    response: dict[str, Any],
    *,
    source_species: str,
    target_species: str,
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for item in response.get("result", []):
        source = _text(item.get("incoming") or item.get("query") or item.get("input"))
        target = _target_symbol(item)
        rows.append(
            {
                "source_species": source_species,
                "target_species": target_species,
                "source_gene": source,
                "target_gene": target,
                # Official g:Orth schema: converted=input/source Ensembl ID,
                # ortholog_ensg=target Ensembl ID, name=target gene symbol.
                "source_ensembl": _text(item.get("converted") or item.get("ensg")),
                "target_ensembl": _text(item.get("ortholog_ensg") or item.get("converted_ensg")),
                "description": _text(item.get("description")),
                "namespaces": _text(item.get("namespaces")),
            }
        )
    return pd.DataFrame(rows, columns=ORTHOLOGY_COLUMNS)


def _gene_level(frame: pd.DataFrame) -> pd.DataFrame:
    result = frame.copy()
    if "gene_symbol" not in result:
        result["gene_symbol"] = result["gene_names"].map(first_gene_symbol)
    result["gene_symbol"] = result["gene_symbol"].astype(str).str.strip()
    return (
        result.loc[result["gene_symbol"].ne("")]
        .sort_values("weighted_priority_score", ascending=False)
        .groupby("gene_symbol", as_index=False)
        .first()
    )


def _empty_summary() -> pd.DataFrame:
    return pd.DataFrame(columns=PAIR_SUMMARY_COLUMNS)


def cross_species_concordance(
    portfolio: pd.DataFrame,
    mapping: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    human = _gene_level(select_species(portfolio, "human"))
    mouse = _gene_level(select_species(portfolio, "mouse"))
    if mapping.empty:
        return pd.DataFrame(), _empty_summary()

    map_clean = mapping.dropna(subset=["source_gene", "target_gene"]).copy()
    map_clean["source_gene_key"] = map_clean["source_gene"].astype(str).str.upper()
    map_clean["target_gene_key"] = map_clean["target_gene"].astype(str).str.upper()
    human["gene_symbol_key"] = human["gene_symbol"].astype(str).str.upper()
    mouse["gene_symbol_key"] = mouse["gene_symbol"].astype(str).str.upper()

    merged = map_clean.merge(
        human,
        left_on="source_gene_key",
        right_on="gene_symbol_key",
        how="inner",
        suffixes=("", "_human"),
    ).merge(
        mouse,
        left_on="target_gene_key",
        right_on="gene_symbol_key",
        how="inner",
        suffixes=("_human", "_mouse"),
    )
    if merged.empty:
        return merged, _empty_summary()

    # A g:Orth query can return multiple target orthologs. Preserve them, but
    # remove exact repeated source-target pairs before concordance summaries.
    merged = merged.drop_duplicates(subset=["source_gene", "target_gene"]).copy()
    merged["repairability_abs_difference"] = (
        merged["repairability_extra_trees_human"] - merged["repairability_extra_trees_mouse"]
    ).abs()
    merged["priority_abs_difference"] = (
        merged["weighted_priority_score_human"] - merged["weighted_priority_score_mouse"]
    ).abs()
    merged["cross_species_score"] = (
        1.0 - merged["repairability_abs_difference"].clip(0, 1)
    ) * 0.6 + (1.0 - merged["priority_abs_difference"].clip(0, 1)) * 0.4
    merged = merged.sort_values("cross_species_score", ascending=False)

    rho_repair = (
        spearmanr(
            merged["repairability_extra_trees_human"],
            merged["repairability_extra_trees_mouse"],
        ).statistic
        if len(merged) >= 3
        else np.nan
    )
    rho_priority = (
        spearmanr(
            merged["weighted_priority_score_human"],
            merged["weighted_priority_score_mouse"],
        ).statistic
        if len(merged) >= 3
        else np.nan
    )
    summary = pd.DataFrame(
        [
            {
                "mapped_portfolio_pairs": len(merged),
                "unique_human_genes": merged["source_gene"].nunique(),
                "unique_mouse_genes": merged["target_gene"].nunique(),
                "spearman_repairability": rho_repair,
                "spearman_priority": rho_priority,
                "median_repairability_abs_difference": merged[
                    "repairability_abs_difference"
                ].median(),
                "median_priority_abs_difference": merged["priority_abs_difference"].median(),
            }
        ],
        columns=PAIR_SUMMARY_COLUMNS,
    )
    return merged, summary


def run_orthology(
    *,
    portfolio: pd.DataFrame,
    client: GProfilerClient,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    human_genes = unique_gene_symbols(select_species(portfolio, "human"))
    manifest: list[dict[str, Any]] = []
    try:
        response = client.orthologs(
            organism="hsapiens",
            target="mmusculus",
            query=human_genes,
            query_name="portfolio_human_to_mouse",
        )
        mapping = flatten_ortholog_response(
            response,
            source_species="human",
            target_species="mouse",
        )
        status, error = "completed", ""
    except OfflineCacheMiss as exc:
        mapping = pd.DataFrame(columns=ORTHOLOGY_COLUMNS)
        status, error = "offline_cache_missing", str(exc)
    except Exception as exc:  # noqa: BLE001
        mapping = pd.DataFrame(columns=ORTHOLOGY_COLUMNS)
        status, error = "failed", f"{type(exc).__name__}: {exc}"

    pairs, summary = cross_species_concordance(portfolio, mapping)
    valid_target_symbols = (
        int(mapping["target_gene"].notna().sum()) if "target_gene" in mapping else 0
    )
    manifest.append(
        {
            "direction": "human_to_mouse",
            "input_genes": len(human_genes),
            "mapping_rows": len(mapping),
            "valid_target_symbols": valid_target_symbols,
            "portfolio_pairs": len(pairs),
            "status": status,
            "error": error,
        }
    )
    manifest_frame = pd.DataFrame(manifest)
    if summary.empty:
        summary = _empty_summary()
    return mapping, pairs, manifest_frame, summary
