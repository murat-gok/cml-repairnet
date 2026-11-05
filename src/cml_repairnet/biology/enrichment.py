from __future__ import annotations

import json
from dataclasses import dataclass
from typing import Any

import numpy as np
import pandas as pd

from cml_repairnet.biology.clients import GProfilerClient, OfflineCacheMiss
from cml_repairnet.biology.schema import (
    SPECIES_TO_ORGANISM,
    evidence_rich_mask,
    select_species,
    top_fraction_mask,
    unique_gene_symbols,
)

DEFAULT_SOURCES = ["GO:BP", "GO:MF", "GO:CC", "REAC", "KEGG", "WP", "CORUM", "HPA"]


@dataclass(frozen=True)
class EnrichmentQuery:
    species: str
    name: str
    genes: list[str]
    background: list[str]
    description: str


def build_enrichment_queries(
    *,
    atlas_in_domain: pd.DataFrame,
    portfolio: pd.DataFrame,
    pareto: pd.DataFrame,
    top_fraction: float = 0.25,
) -> list[EnrichmentQuery]:
    queries: list[EnrichmentQuery] = []
    for species in SPECIES_TO_ORGANISM:
        background_frame = select_species(atlas_in_domain, species)
        background = unique_gene_symbols(background_frame)
        if not background:
            continue
        species_portfolio = select_species(portfolio, species)
        species_pareto = select_species(pareto, species)
        sets = [
            (
                "portfolio",
                unique_gene_symbols(species_portfolio),
                "Diversity-controlled Stage 4 priority portfolio",
            ),
            (
                "pareto",
                unique_gene_symbols(species_pareto),
                "In-domain Pareto-front proteins",
            ),
            (
                "evidence_rich",
                unique_gene_symbols(species_portfolio.loc[evidence_rich_mask(species_portfolio)]),
                "Portfolio proteins with PRM, age, or ortholog evidence",
            ),
            (
                "high_repairability",
                unique_gene_symbols(
                    background_frame.loc[
                        top_fraction_mask(
                            background_frame,
                            "repairability_extra_trees",
                            fraction=top_fraction,
                        )
                    ]
                ),
                f"Top {top_fraction:.0%} of in-domain Extra Trees predictions",
            ),
        ]
        for name, genes, description in sets:
            genes = sorted(set(genes) & set(background))
            if len(genes) < 2:
                continue
            queries.append(
                EnrichmentQuery(
                    species=species,
                    name=name,
                    genes=genes,
                    background=background,
                    description=description,
                )
            )
    return queries


def _extract_intersection_genes(item: dict[str, Any], response: dict[str, Any]) -> list[str]:
    intersections = item.get("intersections")
    if not intersections:
        return []
    # Current API returns a vector aligned to the mapped query genes. Preserve a
    # conservative fallback rather than guessing when metadata changes.
    metadata = response.get("meta", {}).get("genes_metadata", {})
    query_meta = metadata.get("query", {})
    query_name = str(item.get("query", "query_1"))
    info = query_meta.get(query_name, {}) if isinstance(query_meta, dict) else {}
    mapped = info.get("query") or info.get("ensgs") or []
    if len(mapped) != len(intersections):
        return []
    genes: list[str] = []
    for gene, hit in zip(mapped, intersections):
        if hit:
            genes.append(str(gene))
    return genes


def flatten_gprofiler_response(
    response: dict[str, Any],
    *,
    species: str,
    query_name: str,
    query_description: str,
    query_gene_count: int,
    background_gene_count: int,
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for item in response.get("result", []):
        term_size = float(item.get("term_size") or np.nan)
        effective_domain = float(item.get("effective_domain_size") or np.nan)
        precision = float(item.get("precision") or 0.0)
        background_fraction = (
            term_size / effective_domain
            if np.isfinite(term_size) and np.isfinite(effective_domain) and effective_domain > 0
            else np.nan
        )
        fold_enrichment = (
            precision / background_fraction
            if background_fraction and background_fraction > 0
            else np.nan
        )
        intersection_genes = _extract_intersection_genes(item, response)
        rows.append(
            {
                "species": species,
                "query_name": query_name,
                "query_description": query_description,
                "source": item.get("source"),
                "term_id": item.get("native"),
                "term_name": item.get("name"),
                "term_description": item.get("description"),
                "p_value_adjusted": item.get("p_value"),
                "significant": bool(item.get("significant", False)),
                "query_size_mapped": item.get("query_size"),
                "query_gene_count_input": query_gene_count,
                "background_gene_count_input": background_gene_count,
                "intersection_size": item.get("intersection_size"),
                "term_size": item.get("term_size"),
                "effective_domain_size": item.get("effective_domain_size"),
                "precision": item.get("precision"),
                "recall": item.get("recall"),
                "fold_enrichment": fold_enrichment,
                "intersection_genes": " | ".join(intersection_genes),
                "parents_json": json.dumps(item.get("parents", [])),
                "source_order": item.get("source_order"),
                "group_id": item.get("group_id"),
            }
        )
    return pd.DataFrame(rows)


def run_enrichment(
    *,
    atlas_in_domain: pd.DataFrame,
    portfolio: pd.DataFrame,
    pareto: pd.DataFrame,
    client: GProfilerClient,
    sources: list[str] | None = None,
    threshold: float = 0.05,
    no_iea: bool = False,
    top_fraction: float = 0.25,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    sources = sources or DEFAULT_SOURCES
    queries = build_enrichment_queries(
        atlas_in_domain=atlas_in_domain,
        portfolio=portfolio,
        pareto=pareto,
        top_fraction=top_fraction,
    )
    frames: list[pd.DataFrame] = []
    query_rows: list[dict[str, Any]] = []
    version_rows: list[dict[str, Any]] = []
    seen_versions: set[str] = set()
    for query in queries:
        organism = SPECIES_TO_ORGANISM[query.species]
        try:
            response = client.profile(
                organism=organism,
                query=query.genes,
                background=query.background,
                sources=sources,
                query_name=f"{query.species}_{query.name}",
                threshold=threshold,
                correction="fdr",
                no_iea=no_iea,
            )
            status = "completed"
            error = ""
            frame = flatten_gprofiler_response(
                response,
                species=query.species,
                query_name=query.name,
                query_description=query.description,
                query_gene_count=len(query.genes),
                background_gene_count=len(query.background),
            )
            if not frame.empty:
                frames.append(frame)
            if organism not in seen_versions:
                try:
                    versions = client.data_versions(organism)
                    version_rows.append(
                        {"organism": organism, "versions_json": json.dumps(versions)}
                    )
                except Exception as exc:  # noqa: BLE001  # version metadata must not fail the analysis
                    version_rows.append(
                        {"organism": organism, "versions_json": "", "error": str(exc)}
                    )
                seen_versions.add(organism)
        except OfflineCacheMiss as exc:
            status = "offline_cache_missing"
            error = str(exc)
        except Exception as exc:  # noqa: BLE001
            status = "failed"
            error = f"{type(exc).__name__}: {exc}"
        query_rows.append(
            {
                "species": query.species,
                "query_name": query.name,
                "query_description": query.description,
                "query_gene_count": len(query.genes),
                "background_gene_count": len(query.background),
                "status": status,
                "error": error,
                "genes": " | ".join(query.genes),
            }
        )
    enrichment = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
    query_manifest = pd.DataFrame(query_rows)
    versions = pd.DataFrame(version_rows)
    return enrichment, query_manifest, versions


def leading_enrichment_terms(
    enrichment: pd.DataFrame,
    *,
    top_per_source: int = 10,
) -> pd.DataFrame:
    if enrichment.empty:
        return enrichment.copy()
    significant = enrichment.loc[enrichment["significant"].fillna(False)].copy()
    if significant.empty:
        return significant
    significant["p_value_adjusted"] = pd.to_numeric(
        significant["p_value_adjusted"], errors="coerce"
    )
    significant["fold_enrichment"] = pd.to_numeric(significant["fold_enrichment"], errors="coerce")
    return (
        significant.sort_values(
            ["species", "query_name", "source", "p_value_adjusted", "fold_enrichment"],
            ascending=[True, True, True, True, False],
        )
        .groupby(["species", "query_name", "source"], as_index=False, group_keys=False)
        .head(top_per_source)
        .reset_index(drop=True)
    )
