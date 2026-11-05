from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import networkx as nx
import numpy as np
import pandas as pd

from cml_repairnet.biology.clients import OfflineCacheMiss, StringClient
from cml_repairnet.biology.schema import SPECIES_TO_TAXON, select_species, unique_gene_symbols


@dataclass(frozen=True)
class NetworkResult:
    edges: pd.DataFrame
    nodes: pd.DataFrame
    metrics: pd.DataFrame
    ppi_enrichment: pd.DataFrame
    manifest: pd.DataFrame


def normalize_string_edges(frame: pd.DataFrame, *, species: str, query_name: str) -> pd.DataFrame:
    if frame.empty:
        return pd.DataFrame(
            columns=[
                "species",
                "query_name",
                "node_a",
                "node_b",
                "string_id_a",
                "string_id_b",
                "combined_score",
                "nscore",
                "fscore",
                "pscore",
                "ascore",
                "escore",
                "dscore",
                "tscore",
            ]
        )
    rename = {
        "preferredName_A": "node_a",
        "preferredName_B": "node_b",
        "stringId_A": "string_id_a",
        "stringId_B": "string_id_b",
        "score": "combined_score",
    }
    result = frame.rename(columns=rename).copy()
    result.insert(0, "query_name", query_name)
    result.insert(0, "species", species)
    keep = [
        "species",
        "query_name",
        "node_a",
        "node_b",
        "string_id_a",
        "string_id_b",
        "combined_score",
        "nscore",
        "fscore",
        "pscore",
        "ascore",
        "escore",
        "dscore",
        "tscore",
    ]
    for column in keep:
        if column not in result:
            result[column] = np.nan
    result["edge_key"] = result.apply(
        lambda row: "--".join(sorted([str(row["node_a"]), str(row["node_b"])])), axis=1
    )
    return result[keep + ["edge_key"]].drop_duplicates("edge_key")


def graph_tables(
    edges: pd.DataFrame, *, species: str, query_name: str
) -> tuple[pd.DataFrame, pd.DataFrame]:
    graph = nx.Graph()
    for row in edges.itertuples(index=False):
        graph.add_edge(
            str(row.node_a),
            str(row.node_b),
            weight=float(row.combined_score) if pd.notna(row.combined_score) else 0.0,
        )
    if graph.number_of_nodes() == 0:
        empty_nodes = pd.DataFrame(
            columns=[
                "species",
                "query_name",
                "gene_symbol",
                "degree",
                "weighted_degree",
                "betweenness",
                "closeness",
                "eigenvector",
                "community",
            ]
        )
        metrics = pd.DataFrame(
            [
                {
                    "species": species,
                    "query_name": query_name,
                    "nodes": 0,
                    "edges": 0,
                    "density": 0.0,
                    "connected_components": 0,
                    "largest_component": 0,
                    "average_clustering": np.nan,
                    "transitivity": np.nan,
                    "assortativity": np.nan,
                }
            ]
        )
        return empty_nodes, metrics

    degree = dict(graph.degree())
    weighted_degree = dict(graph.degree(weight="weight"))
    betweenness = nx.betweenness_centrality(graph, weight=None, normalized=True)
    closeness = nx.closeness_centrality(graph)
    try:
        eigenvector = nx.eigenvector_centrality_numpy(graph, weight="weight")
    except Exception:  # noqa: BLE001
        eigenvector = {node: np.nan for node in graph.nodes}
    communities = (
        list(nx.community.greedy_modularity_communities(graph, weight="weight"))
        if graph.number_of_edges()
        else []
    )
    community_map: dict[str, int] = {}
    for index, community in enumerate(communities, start=1):
        for node in community:
            community_map[str(node)] = index
    nodes = pd.DataFrame(
        [
            {
                "species": species,
                "query_name": query_name,
                "gene_symbol": node,
                "degree": degree[node],
                "weighted_degree": weighted_degree[node],
                "betweenness": betweenness[node],
                "closeness": closeness[node],
                "eigenvector": eigenvector[node],
                "community": community_map.get(node, 0),
            }
            for node in graph.nodes
        ]
    ).sort_values(["degree", "betweenness"], ascending=False)
    components = [len(component) for component in nx.connected_components(graph)]
    try:
        assortativity = nx.degree_assortativity_coefficient(graph)
    except Exception:  # noqa: BLE001
        assortativity = np.nan
    metrics = pd.DataFrame(
        [
            {
                "species": species,
                "query_name": query_name,
                "nodes": graph.number_of_nodes(),
                "edges": graph.number_of_edges(),
                "density": nx.density(graph),
                "connected_components": len(components),
                "largest_component": max(components) if components else 0,
                "average_clustering": nx.average_clustering(graph, weight="weight")
                if graph.number_of_edges()
                else np.nan,
                "transitivity": nx.transitivity(graph) if graph.number_of_edges() else np.nan,
                "assortativity": assortativity,
            }
        ]
    )
    return nodes, metrics


def run_string_networks(
    *,
    portfolio: pd.DataFrame,
    pareto: pd.DataFrame,
    client: StringClient,
    required_score: int = 700,
) -> NetworkResult:
    edge_frames: list[pd.DataFrame] = []
    node_frames: list[pd.DataFrame] = []
    metric_frames: list[pd.DataFrame] = []
    ppi_frames: list[pd.DataFrame] = []
    manifest_rows: list[dict[str, Any]] = []
    for species in SPECIES_TO_TAXON:
        for query_name, frame in (("portfolio", portfolio), ("pareto", pareto)):
            species_frame = select_species(frame, species)
            genes = unique_gene_symbols(species_frame)
            if len(genes) < 2:
                continue
            try:
                raw_edges = client.network(
                    genes,
                    species=SPECIES_TO_TAXON[species],
                    query_name=f"{species}_{query_name}",
                    required_score=required_score,
                )
                edges = normalize_string_edges(raw_edges, species=species, query_name=query_name)
                nodes, metrics = graph_tables(edges, species=species, query_name=query_name)
                ppi = client.ppi_enrichment(
                    genes,
                    species=SPECIES_TO_TAXON[species],
                    query_name=f"{species}_{query_name}",
                )
                if not ppi.empty:
                    ppi.insert(0, "query_name", query_name)
                    ppi.insert(0, "species", species)
                edge_frames.append(edges)
                node_frames.append(nodes)
                metric_frames.append(metrics)
                ppi_frames.append(ppi)
                status = "completed"
                error = ""
            except OfflineCacheMiss as exc:
                status = "offline_cache_missing"
                error = str(exc)
            except Exception as exc:  # noqa: BLE001
                status = "failed"
                error = f"{type(exc).__name__}: {exc}"
            manifest_rows.append(
                {
                    "species": species,
                    "query_name": query_name,
                    "input_genes": len(genes),
                    "required_score": required_score,
                    "status": status,
                    "error": error,
                    "genes": " | ".join(genes),
                }
            )
    concat = lambda frames: (
        pd.concat([f for f in frames if not f.empty], ignore_index=True)
        if any(not f.empty for f in frames)
        else pd.DataFrame()
    )
    return NetworkResult(
        edges=concat(edge_frames),
        nodes=concat(node_frames),
        metrics=concat(metric_frames),
        ppi_enrichment=concat(ppi_frames),
        manifest=pd.DataFrame(manifest_rows),
    )
