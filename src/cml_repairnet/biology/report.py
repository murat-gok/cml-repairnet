from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd


def _count(path: Path) -> int:
    if not path.exists():
        return 0
    try:
        return len(pd.read_csv(path))
    except Exception:
        return 0


def build_stage5_audit(paths: dict[str, Path]) -> dict[str, Any]:
    audit = {
        "stage4_portfolio_sites": _count(paths["portfolio_input"]),
        "stage4_pareto_sites": _count(paths["pareto_input"]),
        "stage4_in_domain_sites": _count(paths["atlas_in_domain_input"]),
        "uniprot_annotations": _count(paths["uniprot_annotations"]),
        "protein_families": _count(paths["protein_families"]),
        "enrichment_terms_all": _count(paths["enrichment_all"]),
        "enrichment_terms_significant": 0,
        "string_edges": _count(paths["string_edges"]),
        "string_nodes": _count(paths["string_nodes"]),
        "ortholog_pairs": _count(paths["ortholog_pairs"]),
        "orthology_mapping_rows": _count(paths["ortholog_mapping"]),
        "candidate_modules": _count(paths["module_profiles"]),
        "module_assignments": _count(paths["module_assignments"]),
        "interpretation": (
            "Stage 5 performs biological annotation and hypothesis generation from the frozen "
            "Stage 4 computational atlas. Enrichment and network results are not independent "
            "experimental validation of CML repairability."
        ),
        "background_safeguard": (
            "Functional enrichment uses species-specific Stage 4 in-domain genes as the custom "
            "background, rather than the whole genome, to reduce candidate-universe bias."
        ),
    }
    if paths["enrichment_all"].exists():
        frame = _safe_read(paths["enrichment_all"])
        if "significant" in frame:
            audit["enrichment_terms_significant"] = int(
                frame["significant"].fillna(False).astype(bool).sum()
            )
    if paths["orthology_manifest"].exists():
        orth_manifest = _safe_read(paths["orthology_manifest"])
        if not orth_manifest.empty:
            audit["orthology_valid_target_symbols"] = int(
                pd.to_numeric(
                    orth_manifest.get("valid_target_symbols", pd.Series([0])), errors="coerce"
                )
                .fillna(0)
                .sum()
            )
            audit["orthology_status"] = str(orth_manifest.iloc[0].get("status", "unknown"))
    if paths["module_selection"].exists():
        selection = _safe_read(paths["module_selection"])
        if not selection.empty and "silhouette" in selection:
            scores = pd.to_numeric(selection["silhouette"], errors="coerce")
            if scores.notna().any():
                idx = scores.idxmax()
                audit["module_selected_k"] = int(selection.loc[idx, "k"])
                audit["module_best_silhouette"] = float(scores.loc[idx])
                audit["module_interpretation"] = (
                    "exploratory_low_separation"
                    if float(scores.loc[idx]) < 0.30
                    else "exploratory_moderate_separation"
                )

    manifests = {}
    for name in (
        "annotation_manifest",
        "enrichment_manifest",
        "network_manifest",
        "orthology_manifest",
    ):
        path = paths[name]
        if path.exists():
            manifests[name] = _safe_read(path).to_dict(orient="records")
    audit["external_service_manifests"] = manifests
    complete_core = all(
        paths[name].exists()
        for name in ("module_assignments", "module_profiles", "candidate_table")
    )
    audit["stage5_completion_status"] = "complete" if complete_core else "partial"
    return audit


def _safe_read(path: Path) -> pd.DataFrame:
    if not path.exists() or path.stat().st_size == 0:
        return pd.DataFrame()
    try:
        return pd.read_csv(path)
    except pd.errors.EmptyDataError:
        return pd.DataFrame()


def write_stage5_report(paths: dict[str, Path]) -> dict[str, Any]:
    audit = build_stage5_audit(paths)
    paths["audit_json"].parent.mkdir(parents=True, exist_ok=True)
    paths["audit_json"].write_text(json.dumps(audit, indent=2), encoding="utf-8")

    portfolio = _safe_read(paths["candidate_table"])
    enrichment = _safe_read(paths["leading_enrichment"])
    modules = _safe_read(paths["module_profiles"])
    orthology = _safe_read(paths["ortholog_summary"])

    lines = [
        "# CML-RepairNet Stage 5 Biological Interpretation Report",
        "",
        "## Scope",
        "",
        f"- Stage 4 priority sites analysed: {audit['stage4_portfolio_sites']}",
        f"- In-domain background sites: {audit['stage4_in_domain_sites']}",
        f"- Candidate modules: {audit['candidate_modules']}",
        f"- Significant enrichment terms: {audit['enrichment_terms_significant']}",
        f"- STRING network edges: {audit['string_edges']}",
        f"- Cross-species portfolio pairs: {audit['ortholog_pairs']}",
        "",
        "## Safeguards",
        "",
        f"- {audit['interpretation']}",
        f"- {audit['background_safeguard']}",
    ]
    if not portfolio.empty:
        lines += ["", "## Leading portfolio candidates", ""]
        columns = [
            c
            for c in [
                "portfolio_rank",
                "species",
                "gene_symbol",
                "site_uniprot",
                "repairability_extra_trees",
                "weighted_priority_score",
                "module_label",
            ]
            if c in portfolio
        ]
        lines.append(
            portfolio.sort_values("portfolio_rank").head(15)[columns].to_markdown(index=False)
        )
    if not enrichment.empty:
        lines += ["", "## Leading functional terms", ""]
        columns = [
            c
            for c in [
                "species",
                "query_name",
                "source",
                "term_id",
                "term_name",
                "p_value_adjusted",
                "fold_enrichment",
            ]
            if c in enrichment
        ]
        lines.append(enrichment.head(20)[columns].to_markdown(index=False))
    if not modules.empty:
        lines += ["", "## Candidate modules", ""]
        columns = [
            c
            for c in [
                "module_id",
                "module_label",
                "sites",
                "proteins",
                "human_sites",
                "mouse_sites",
                "repairability_extra_trees",
                "evidence_score_raw",
                "transfer_reliability_score",
            ]
            if c in modules
        ]
        lines.append(modules[columns].to_markdown(index=False))
    if not orthology.empty:
        lines += ["", "## Cross-species concordance", "", orthology.to_markdown(index=False)]
    paths["summary_md"].write_text("\n".join(lines) + "\n", encoding="utf-8")
    return audit
