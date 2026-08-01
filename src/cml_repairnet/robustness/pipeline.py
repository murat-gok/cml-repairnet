from __future__ import annotations

import json
from pathlib import Path

from .conservation import run_residue_conservation
from .enrichment_sensitivity import run_enrichment_sensitivity
from .figures import make_figures
from .network_robustness import run_network_robustness
from .ranking import run_rank_stability
from .report import build_report
from .schema import input_paths, load_config, output_paths, read_csv
from .tiers import run_candidate_tiers


def doctor_stage6(root: Path) -> dict:
    inputs = input_paths(root)
    status = {
        name: {"path": str(path), "exists": path.exists()} for name, path in inputs.__dict__.items()
    }
    status["config"] = {
        "path": str(root / "configs/stage6.yaml"),
        "exists": (root / "configs/stage6.yaml").exists(),
    }
    status["ready"] = all(item["exists"] for key, item in status.items() if key != "ready")
    return status


def rank_stage6(root: Path) -> dict:
    paths = input_paths(root)  # restored for rank_stage6
    outputs = output_paths(root)
    config = load_config(root)
    in_domain = read_csv(paths.atlas_in_domain, required=("candidate_id",))
    portfolio = read_csv(paths.portfolio, required=("candidate_id",))
    pareto = read_csv(paths.pareto, required=("candidate_id",))
    rank_config = dict(config["rank_stability"])
    rank_config["random_seed"] = config["random_seed"]
    stability, weights = run_rank_stability(in_domain, portfolio, pareto, rank_config)
    stability.to_csv(outputs.processed / "stage6_rank_stability.csv", index=False)
    weights.to_csv(outputs.processed / "stage6_priority_weight_sensitivity.csv", index=False)
    return {
        "candidates": len(stability),
        "iterations": rank_config["iterations"],
        "weight_scenarios": len(weights),
    }


def conservation_stage6(root: Path) -> dict:
    paths = input_paths(root)  # restored for conservation_stage6
    outputs = output_paths(root)
    config = load_config(root)
    pairs = read_csv(
        paths.cross_species_pairs, required=("source_gene", "sequence_human", "sequence_mouse")
    )
    conservation = run_residue_conservation(pairs, config["residue_conservation"])
    conservation.to_csv(outputs.processed / "stage6_residue_conservation.csv", index=False)
    return {
        "pairs": len(conservation),
        "conserved_lysine_pairs": int(conservation["residue_level_conserved_lysine"].sum()),
    }


def enrichment_stage6(root: Path) -> dict:
    paths = input_paths(root)  # restored for enrichment_stage6
    outputs = output_paths(root)
    config = load_config(root)
    enrichment = read_csv(paths.enrichment, required=("term_id", "intersection_genes"))
    sensitivity = run_enrichment_sensitivity(enrichment, config["enrichment_sensitivity"])
    sensitivity.to_csv(outputs.processed / "stage6_enrichment_driver_sensitivity.csv", index=False)
    fragile = (
        sensitivity["driver_sensitivity_class"]
        .isin(["single_driver", "leave_one_driver_fragile"])
        .sum()
        if not sensitivity.empty
        else 0
    )
    return {"terms": len(sensitivity), "fragile_terms": int(fragile)}


def network_stage6(root: Path) -> dict:
    paths = input_paths(root)  # restored for network_stage6
    outputs = output_paths(root)
    config = load_config(root)
    edges = read_csv(paths.string_edges, required=("node_a", "node_b", "combined_score"))
    nodes = read_csv(paths.string_nodes, required=("species", "query_name"))
    manifest = read_csv(paths.network_manifest, required=("species", "query_name", "genes"))
    thresholds, hubs, null = run_network_robustness(
        edges, nodes, manifest, config["network_robustness"], int(config["random_seed"])
    )
    thresholds.to_csv(outputs.processed / "stage6_network_threshold_sensitivity.csv", index=False)
    hubs.to_csv(outputs.processed / "stage6_network_hub_removal.csv", index=False)
    null.to_csv(outputs.processed / "stage6_network_degree_preserving_null.csv", index=False)
    return {"threshold_rows": len(thresholds), "hub_rows": len(hubs), "null_rows": len(null)}


def tiers_stage6(root: Path) -> dict:
    paths = input_paths(root)  # restored for tiers_stage6
    outputs = output_paths(root)
    config = load_config(root)
    atlas_all = read_csv(paths.atlas_all, required=("candidate_id", "applicability_status"))
    stability = read_csv(
        outputs.processed / "stage6_rank_stability.csv", required=("candidate_id",)
    )
    conservation = read_csv(outputs.processed / "stage6_residue_conservation.csv")
    modules = read_csv(paths.stage5_modules, required=("candidate_id",))
    tiers = run_candidate_tiers(atlas_all, stability, conservation, modules, config["tiers"])
    tiers.to_csv(outputs.processed / "stage6_candidate_tiers.csv", index=False)
    manuscript = tiers.loc[tiers["candidate_tier"] != "Watchlist"].copy()
    manuscript.to_csv(outputs.processed / "stage6_manuscript_candidate_table.csv", index=False)
    return {
        "all_sites": len(tiers),
        "tier_counts": tiers["candidate_tier"].value_counts().to_dict(),
    }


def figures_stage6(root: Path) -> dict:
    outputs = output_paths(root)
    stability = read_csv(outputs.processed / "stage6_rank_stability.csv")
    conservation = read_csv(outputs.processed / "stage6_residue_conservation.csv")
    enrichment = read_csv(outputs.processed / "stage6_enrichment_driver_sensitivity.csv")
    network = read_csv(outputs.processed / "stage6_network_threshold_sensitivity.csv")
    tiers = read_csv(outputs.processed / "stage6_candidate_tiers.csv")
    count = make_figures(stability, conservation, enrichment, network, tiers, outputs.figures)
    return {"figure_sets": count}


def report_stage6(root: Path) -> dict:
    outputs = output_paths(root)
    config = load_config(root)
    stability = read_csv(outputs.processed / "stage6_rank_stability.csv")
    weights = read_csv(outputs.processed / "stage6_priority_weight_sensitivity.csv")
    conservation = read_csv(outputs.processed / "stage6_residue_conservation.csv")
    enrichment = read_csv(outputs.processed / "stage6_enrichment_driver_sensitivity.csv")
    thresholds = read_csv(outputs.processed / "stage6_network_threshold_sensitivity.csv")
    hubs = read_csv(outputs.processed / "stage6_network_hub_removal.csv")
    null = read_csv(outputs.processed / "stage6_network_degree_preserving_null.csv")
    tiers = read_csv(outputs.processed / "stage6_candidate_tiers.csv")
    audit = {
        "stage6_version": "0.7.2",
        "stage6_completion_status": "complete",
        "in_domain_candidates": len(stability),
        "rank_iterations": int(config["rank_stability"]["iterations"]),
        "pareto_iterations": int(config["rank_stability"]["pareto_iterations"]),
        "residue_pairs": len(conservation),
        "residue_level_conserved_lysine_pairs": int(
            conservation["residue_level_conserved_lysine"].sum()
        ),
        "enrichment_terms_assessed": len(enrichment),
        "network_threshold_rows": len(thresholds),
        "candidate_tier_counts": tiers["candidate_tier"].value_counts().to_dict(),
        "interpretation": "Stage 6 is a prespecified computational robustness and manuscript-validation layer over frozen Stage 4-5 outputs; it is not experimental validation.",
        "rank_stability_scope": "Monte Carlo stress test under configured perturbations. Global top-k frequency is diagnostic; portfolio frequency recomputes complete species-specific nondominated fronts and reapplies Stage 4 species balance and the per-protein cap.",
        "portfolio_stability_scope": "Exact frozen Stage 4 selection order: species-specific nondominated front rank, weighted priority within front, species quota, and per-protein cap.",
        "baseline_portfolio_reproduction_count": int(weights["baseline_reproduction_count"].iloc[0])
        if len(weights)
        else 0,
        "baseline_portfolio_reproduction_fraction": float(
            weights["baseline_reproduction_fraction"].iloc[0]
        )
        if len(weights)
        else 0.0,
        "network_null_scope": "Degree-preserving topology null, not a degree-matched random-gene-set analysis.",
    }
    (outputs.results / "stage6_audit.json").write_text(
        json.dumps(audit, indent=2), encoding="utf-8"
    )
    build_report(
        stability,
        weights,
        conservation,
        enrichment,
        thresholds,
        hubs,
        null,
        tiers,
        audit,
        outputs.results / "stage6_robustness_report.md",
    )
    return audit


def run_stage6_all(root: Path) -> dict:
    doctor = doctor_stage6(root)
    if not doctor["ready"]:
        raise FileNotFoundError(f"Stage 6 inputs are incomplete: {doctor}")
    return {
        "doctor": doctor,
        "rank": rank_stage6(root),
        "conservation": conservation_stage6(root),
        "enrichment": enrichment_stage6(root),
        "network": network_stage6(root),
        "tiers": tiers_stage6(root),
        "figures": figures_stage6(root),
        "report": report_stage6(root),
    }
