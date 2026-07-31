from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd
import yaml

from cml_repairnet.biology.annotations import run_uniprot_annotations
from cml_repairnet.biology.clients import GProfilerClient, StringClient, UniProtAnnotationClient
from cml_repairnet.biology.enrichment import leading_enrichment_terms, run_enrichment
from cml_repairnet.biology.figures import enrichment_bubble_plot, cross_species_scatter, module_heatmap, priority_landscape
from cml_repairnet.biology.modules import build_candidate_modules
from cml_repairnet.biology.networks import run_string_networks
from cml_repairnet.biology.orthology import run_orthology
from cml_repairnet.biology.report import write_stage5_report
from cml_repairnet.biology.schema import load_atlas


def stage5_paths(project_root: Path) -> dict[str, Path]:
    result_root = project_root / "results" / "stage5"
    processed = project_root / "data" / "processed" / "stage5"
    cache = project_root / "data" / "interim" / "stage5" / "api_cache"
    figures = result_root / "figures"
    stage4_atlas = project_root / "results" / "stage4" / "atlas"
    return {
        "portfolio_input": stage4_atlas / "cml_repairability_priority_portfolio.csv",
        "pareto_input": stage4_atlas / "cml_repairability_pareto_front.csv",
        "atlas_in_domain_input": stage4_atlas / "cml_repairability_atlas_in_domain.csv",
        "watchlist_input": stage4_atlas / "cml_repairability_domain_shift_watchlist.csv",
        "processed": processed,
        "results": result_root,
        "cache": cache,
        "figures": figures,
        "uniprot_annotations": processed / "stage5_uniprot_annotations.csv",
        "protein_families": processed / "stage5_protein_family_summary.csv",
        "annotation_manifest": result_root / "stage5_annotation_manifest.csv",
        "enrichment_all": processed / "stage5_gprofiler_enrichment_all.csv",
        "leading_enrichment": processed / "stage5_gprofiler_enrichment_leading.csv",
        "enrichment_manifest": result_root / "stage5_enrichment_manifest.csv",
        "enrichment_versions": result_root / "stage5_enrichment_data_versions.csv",
        "string_edges": processed / "stage5_string_edges.csv",
        "string_nodes": processed / "stage5_string_nodes.csv",
        "string_metrics": processed / "stage5_string_network_metrics.csv",
        "string_ppi": processed / "stage5_string_ppi_enrichment.csv",
        "network_manifest": result_root / "stage5_network_manifest.csv",
        "ortholog_mapping": processed / "stage5_ortholog_mapping.csv",
        "ortholog_pairs": processed / "stage5_cross_species_pairs.csv",
        "ortholog_summary": processed / "stage5_cross_species_summary.csv",
        "orthology_manifest": result_root / "stage5_orthology_manifest.csv",
        "module_assignments": processed / "stage5_candidate_module_assignments.csv",
        "module_profiles": processed / "stage5_candidate_module_profiles.csv",
        "module_selection": result_root / "stage5_module_k_selection.csv",
        "candidate_table": processed / "stage5_manuscript_candidate_table.csv",
        "audit_json": result_root / "stage5_audit.json",
        "summary_md": result_root / "stage5_biological_interpretation.md",
    }


def load_stage5_config(project_root: Path) -> dict[str, Any]:
    path = project_root / "configs" / "stage5.yaml"
    if not path.exists():
        return {}
    return yaml.safe_load(path.read_text(encoding="utf-8")) or {}




def _safe_read_csv(path: Path) -> pd.DataFrame:
    if not path.exists() or path.stat().st_size == 0:
        return pd.DataFrame()
    try:
        return pd.read_csv(path)
    except pd.errors.EmptyDataError:
        return pd.DataFrame()


def _ensure_dirs(paths: dict[str, Path]) -> None:
    for key in ("processed", "results", "cache", "figures"):
        paths[key].mkdir(parents=True, exist_ok=True)


def doctor_stage5(project_root: Path) -> dict[str, Any]:
    paths = stage5_paths(project_root)
    status = {
        "portfolio": {"path": str(paths["portfolio_input"]), "exists": paths["portfolio_input"].exists()},
        "pareto": {"path": str(paths["pareto_input"]), "exists": paths["pareto_input"].exists()},
        "atlas_in_domain": {"path": str(paths["atlas_in_domain_input"]), "exists": paths["atlas_in_domain_input"].exists()},
    }
    status["ready"] = all(item["exists"] for item in status.values())
    return status


def _inputs(paths: dict[str, Path]) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    return (
        load_atlas(paths["atlas_in_domain_input"], name="Stage 4 in-domain atlas"),
        load_atlas(paths["portfolio_input"], name="Stage 4 priority portfolio"),
        load_atlas(paths["pareto_input"], name="Stage 4 Pareto front"),
    )


def annotate_stage5(project_root: Path, *, offline: bool = False) -> dict[str, int]:
    paths = stage5_paths(project_root); _ensure_dirs(paths)
    _, portfolio, _ = _inputs(paths)
    client = UniProtAnnotationClient(paths["cache"] / "uniprot", offline=offline)
    annotations, families, manifest = run_uniprot_annotations(portfolio, client)
    annotations.to_csv(paths["uniprot_annotations"], index=False)
    families.to_csv(paths["protein_families"], index=False)
    manifest.to_csv(paths["annotation_manifest"], index=False)
    return {"annotations": len(annotations), "protein_families": len(families)}


def enrich_stage5(project_root: Path, *, offline: bool = False) -> dict[str, int]:
    paths = stage5_paths(project_root); _ensure_dirs(paths)
    config = load_stage5_config(project_root).get("enrichment", {})
    atlas, portfolio, pareto = _inputs(paths)
    client = GProfilerClient(paths["cache"] / "gprofiler", offline=offline)
    enrichment, manifest, versions = run_enrichment(
        atlas_in_domain=atlas,
        portfolio=portfolio,
        pareto=pareto,
        client=client,
        sources=config.get("sources"),
        threshold=float(config.get("threshold", 0.05)),
        no_iea=bool(config.get("no_iea", False)),
        top_fraction=float(config.get("high_repairability_fraction", 0.25)),
    )
    leading = leading_enrichment_terms(enrichment, top_per_source=int(config.get("top_per_source", 10)))
    enrichment.to_csv(paths["enrichment_all"], index=False)
    leading.to_csv(paths["leading_enrichment"], index=False)
    manifest.to_csv(paths["enrichment_manifest"], index=False)
    versions.to_csv(paths["enrichment_versions"], index=False)
    return {"terms": len(enrichment), "significant": int(enrichment.get("significant", pd.Series(dtype=bool)).fillna(False).astype(bool).sum()), "leading": len(leading)}


def network_stage5(project_root: Path, *, offline: bool = False) -> dict[str, int]:
    paths = stage5_paths(project_root); _ensure_dirs(paths)
    config = load_stage5_config(project_root).get("string", {})
    _, portfolio, pareto = _inputs(paths)
    client = StringClient(paths["cache"] / "string", offline=offline, min_delay_seconds=1.0)
    result = run_string_networks(
        portfolio=portfolio, pareto=pareto, client=client,
        required_score=int(config.get("required_score", 700)),
    )
    result.edges.to_csv(paths["string_edges"], index=False)
    result.nodes.to_csv(paths["string_nodes"], index=False)
    result.metrics.to_csv(paths["string_metrics"], index=False)
    result.ppi_enrichment.to_csv(paths["string_ppi"], index=False)
    result.manifest.to_csv(paths["network_manifest"], index=False)
    return {"edges": len(result.edges), "nodes": len(result.nodes)}


def orthology_stage5(project_root: Path, *, offline: bool = False) -> dict[str, int]:
    paths = stage5_paths(project_root); _ensure_dirs(paths)
    _, portfolio, _ = _inputs(paths)
    client = GProfilerClient(paths["cache"] / "gprofiler", offline=offline)
    mapping, pairs, manifest, summary = run_orthology(portfolio=portfolio, client=client)
    mapping.to_csv(paths["ortholog_mapping"], index=False)
    pairs.to_csv(paths["ortholog_pairs"], index=False)
    manifest.to_csv(paths["orthology_manifest"], index=False)
    summary.to_csv(paths["ortholog_summary"], index=False)
    return {
        "mapping_rows": len(mapping),
        "valid_target_symbols": int(mapping.get("target_gene", pd.Series(dtype=object)).notna().sum()),
        "portfolio_pairs": len(pairs),
    }


def modules_stage5(project_root: Path) -> dict[str, int]:
    paths = stage5_paths(project_root); _ensure_dirs(paths)
    config = load_stage5_config(project_root).get("modules", {})
    _, portfolio, _ = _inputs(paths)
    assignments, profiles, selection = build_candidate_modules(
        portfolio,
        random_state=int(config.get("random_state", 42)),
        max_k=int(config.get("max_k", 6)),
    )
    assignments.to_csv(paths["module_assignments"], index=False)
    profiles.to_csv(paths["module_profiles"], index=False)
    selection.to_csv(paths["module_selection"], index=False)
    # Manuscript-ready table carries module labels and optional UniProt annotations.
    candidate = assignments.copy()
    if paths["uniprot_annotations"].exists():
        annotations = _safe_read_csv(paths["uniprot_annotations"])
        if not annotations.empty and "uniprot_accession" in annotations:
            candidate = candidate.merge(annotations, on="uniprot_accession", how="left", suffixes=("", "_stage5"))
    candidate.sort_values("portfolio_rank").to_csv(paths["candidate_table"], index=False)
    return {"modules": len(profiles), "assignments": len(assignments)}


def figures_stage5(project_root: Path) -> dict[str, int]:
    paths = stage5_paths(project_root); _ensure_dirs(paths)
    made = 0
    if paths["leading_enrichment"].exists():
        enrichment = _safe_read_csv(paths["leading_enrichment"])
        if not enrichment.empty:
            enrichment_bubble_plot(enrichment, paths["figures"] / "stage5_enrichment_bubble"); made += 1
    if paths["ortholog_pairs"].exists():
        pairs = _safe_read_csv(paths["ortholog_pairs"])
        if not pairs.empty: cross_species_scatter(pairs, paths["figures"] / "stage5_cross_species_concordance"); made += 1
    if paths["module_profiles"].exists():
        profiles = _safe_read_csv(paths["module_profiles"])
        if not profiles.empty:
            module_heatmap(profiles, paths["figures"] / "stage5_candidate_modules"); made += 1
    _, portfolio, _ = _inputs(paths)
    priority_landscape(portfolio, paths["figures"] / "stage5_priority_landscape"); made += 1
    return {"figure_sets": made}


def report_stage5(project_root: Path) -> dict[str, Any]:
    paths = stage5_paths(project_root); _ensure_dirs(paths)
    if not paths["candidate_table"].exists():
        modules_stage5(project_root)
    return write_stage5_report(paths)


def run_stage5_all(project_root: Path, *, offline: bool = False) -> dict[str, Any]:
    result: dict[str, Any] = {"doctor": doctor_stage5(project_root)}
    if not result["doctor"]["ready"]:
        raise RuntimeError("Stage 5 prerequisites are missing. Run Stage 4 full atlas first.")
    result["annotate"] = annotate_stage5(project_root, offline=offline)
    result["enrich"] = enrich_stage5(project_root, offline=offline)
    result["network"] = network_stage5(project_root, offline=offline)
    result["orthology"] = orthology_stage5(project_root, offline=offline)
    result["modules"] = modules_stage5(project_root)
    result["figures"] = figures_stage5(project_root)
    result["report"] = report_stage5(project_root)
    return result
