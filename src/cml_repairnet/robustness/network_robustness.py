from __future__ import annotations

import numpy as np
import pandas as pd
import networkx as nx


def _graph(edges: pd.DataFrame, nodes: list[str], threshold: float) -> nx.Graph:
    graph = nx.Graph()
    graph.add_nodes_from(nodes)
    subset = edges.loc[pd.to_numeric(edges["combined_score"], errors="coerce") >= threshold]
    for row in subset.itertuples(index=False):
        graph.add_edge(str(row.node_a), str(row.node_b), weight=float(row.combined_score))
    return graph


def _metrics(graph: nx.Graph) -> dict[str, float | int]:
    components = list(nx.connected_components(graph)) if graph.number_of_nodes() else []
    largest = max((len(component) for component in components), default=0)
    return {
        "nodes": graph.number_of_nodes(),
        "edges": graph.number_of_edges(),
        "density": nx.density(graph) if graph.number_of_nodes() > 1 else 0.0,
        "components": len(components),
        "largest_component_nodes": largest,
        "largest_component_fraction": largest / graph.number_of_nodes() if graph.number_of_nodes() else 0.0,
        "isolated_nodes": sum(1 for node in graph if graph.degree(node) == 0),
        "mean_degree": float(np.mean([degree for _, degree in graph.degree()])) if graph.number_of_nodes() else 0.0,
        "transitivity": nx.transitivity(graph) if graph.number_of_edges() else 0.0,
    }


def run_network_robustness(
    edges: pd.DataFrame,
    nodes: pd.DataFrame,
    manifest: pd.DataFrame,
    config: dict,
    seed: int,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    required_edges = {"species", "query_name", "node_a", "node_b", "combined_score"}
    required_nodes = {"species", "query_name"}
    if required_edges - set(edges.columns):
        raise ValueError(f"STRING edge table missing: {sorted(required_edges - set(edges.columns))}")
    if required_nodes - set(nodes.columns):
        raise ValueError(f"STRING node table missing: {sorted(required_nodes - set(nodes.columns))}")
    node_column = "node" if "node" in nodes.columns else "gene_symbol" if "gene_symbol" in nodes.columns else None
    if node_column is None:
        raise ValueError("STRING node table must contain node or gene_symbol")
    required_manifest = {"species", "query_name", "genes"}
    if required_manifest - set(manifest.columns):
        raise ValueError(f"STRING manifest missing: {sorted(required_manifest - set(manifest.columns))}")
    thresholds = [float(value) for value in config.get("thresholds", [0.7, 0.8, 0.9])]
    threshold_rows = []
    hub_rows = []
    null_rows = []
    rng = np.random.default_rng(seed)

    for manifest_row in manifest.itertuples(index=False):
        species = str(manifest_row.species)
        query_name = str(manifest_row.query_name)
        manifest_nodes = {part.strip() for part in str(manifest_row.genes).split("|") if part.strip()}
        node_group = nodes.loc[(nodes["species"] == species) & (nodes["query_name"] == query_name)]
        metric_nodes = set(node_group[node_column].dropna().astype(str))
        node_names = sorted(manifest_nodes | metric_nodes)
        edge_group = edges.loc[(edges["species"] == species) & (edges["query_name"] == query_name)].copy()
        for threshold in thresholds:
            graph = _graph(edge_group, node_names, threshold)
            threshold_rows.append({"species": species, "query_name": query_name, "threshold": threshold, **_metrics(graph)})

        base = _graph(edge_group, node_names, thresholds[0])
        degree_order = [node for node, _ in sorted(base.degree(), key=lambda item: (-item[1], item[0]))]
        base_metrics = _metrics(base)
        for count in [int(value) for value in config.get("hub_removal_counts", [1, 3, 5])]:
            removed = degree_order[: min(count, len(degree_order))]
            graph = base.copy()
            graph.remove_nodes_from(removed)
            current = _metrics(graph)
            hub_rows.append({
                "species": species,
                "query_name": query_name,
                "removed_hub_count": len(removed),
                "removed_hubs": " | ".join(removed),
                **current,
                "edge_retention_fraction": current["edges"] / base_metrics["edges"] if base_metrics["edges"] else np.nan,
                "largest_component_retention_fraction": current["largest_component_nodes"] / base_metrics["largest_component_nodes"] if base_metrics["largest_component_nodes"] else np.nan,
            })

        iterations = int(config.get("rewiring_iterations", 250))
        swaps_per_edge = int(config.get("rewiring_swaps_per_edge", 5))
        observed = nx.transitivity(base) if base.number_of_edges() else 0.0
        null_values = []
        if base.number_of_edges() >= 4 and base.number_of_nodes() >= 4:
            for _ in range(iterations):
                rewired = base.copy()
                try:
                    nx.double_edge_swap(
                        rewired,
                        nswap=max(rewired.number_of_edges() * swaps_per_edge, 1),
                        max_tries=max(rewired.number_of_edges() * swaps_per_edge * 20, 100),
                        seed=int(rng.integers(0, 2**31 - 1)),
                    )
                    null_values.append(nx.transitivity(rewired))
                except nx.NetworkXAlgorithmError:
                    continue
        null_array = np.asarray(null_values, dtype=float)
        empirical_p = (
            (np.sum(null_array >= observed) + 1) / (len(null_array) + 1)
            if len(null_array) else np.nan
        )
        null_rows.append({
            "species": species,
            "query_name": query_name,
            "threshold": thresholds[0],
            "observed_transitivity": observed,
            "degree_preserving_null_iterations": len(null_array),
            "null_mean_transitivity": float(np.mean(null_array)) if len(null_array) else np.nan,
            "null_sd_transitivity": float(np.std(null_array, ddof=1)) if len(null_array) > 1 else np.nan,
            "empirical_upper_tail_p": empirical_p,
            "interpretation_scope": "degree-preserving topology test; not a random-gene-set test",
        })
    return pd.DataFrame(threshold_rows), pd.DataFrame(hub_rows), pd.DataFrame(null_rows)
