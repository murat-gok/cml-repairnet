from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd

from cml_repairnet.atlas.modeling import score_candidate_features
from cml_repairnet.atlas.prioritize import prioritize_atlas
from cml_repairnet.atlas.sources import (
    deterministic_quick_subset,
    download_stage4_source_data,
    prepare_candidate_table,
    annotate_candidates_with_uniprot,
)
from cml_repairnet.atlas.structures import build_alphafold_features


def stage4_paths(project_root: Path) -> dict[str, Path]:
    data_root = project_root / "data"
    raw = data_root / "raw" / "stage4"
    interim = data_root / "interim" / "stage4"
    processed = data_root / "processed" / "stage4"
    results = project_root / "results" / "stage4"
    models = project_root / "models" / "stage4"
    return {
        "raw": raw,
        "interim": interim,
        "processed": processed,
        "results": results,
        "models": models,
        "source_manifest": raw / "stage4_source_manifest.json",
        "candidates": interim / "stage4_candidates.csv",
        "prm_evidence": interim / "stage4_prm_evidence.csv",
        "parse_diagnostics": results / "stage4_parse_diagnostics.csv",
        "annotated": interim / "stage4_candidates_uniprot.csv",
        "prm_mapping": interim / "stage4_prm_mapping.csv",
        "prm_mapping_excluded": results / "stage4_excluded_prm_mapping.csv",
        "sequence_excluded": results / "stage4_excluded_sequence.csv",
        "structure_features": processed / "stage4_candidate_structure_features.csv",
        "structure_excluded": results / "stage4_excluded_structure.csv",
        "scored": processed / "stage4_scored_candidates.csv",
        "atlas_dir": results / "atlas",
        "audit_json": results / "stage4_audit.json",
        "audit_md": results / "stage4_audit.md",
    }


def doctor_stage4(project_root: Path) -> dict[str, Any]:
    paths = stage4_paths(project_root)
    required_prior = {
        "bsa_labels": project_root / "data/processed/bsa_labels.csv",
        "bsa_structure_features": project_root / "data/processed/bsa_structure_features.csv",
    }
    status = {
        name: {"path": str(path), "exists": path.exists()}
        for name, path in required_prior.items()
    }
    status["stage3_model_freeze_recommended"] = {
        "path": str(project_root / "results/stage3"),
        "exists": (project_root / "results/stage3").exists(),
    }
    status["ready"] = all(item["exists"] for item in status.values() if isinstance(item, dict) and "exists" in item and not item.get("path", "").endswith("results/stage3"))
    paths["results"].mkdir(parents=True, exist_ok=True)
    (paths["results"] / "stage4_doctor.json").write_text(
        json.dumps(status, indent=2), encoding="utf-8"
    )
    return status


def download_stage4(project_root: Path, *, overwrite: bool = False) -> dict[str, Any]:
    paths = stage4_paths(project_root)
    return download_stage4_source_data(paths["raw"], overwrite=overwrite)


def prepare_stage4(
    project_root: Path,
    *,
    candidate_csv: Path | None = None,
    mode: str = "full",
) -> pd.DataFrame:
    paths = stage4_paths(project_root)
    workbooks = [
        paths["raw"] / "cellular_cml_sites.xlsx",
        paths["raw"] / "mouse_organ_cml_sites.xlsx",
    ]
    available = [path for path in workbooks if path.exists()]
    if not available and candidate_csv is None:
        raise FileNotFoundError(
            "No Stage 4 supplementary workbooks found. Run the download step or provide "
            "--candidate-csv."
        )
    candidates = prepare_candidate_table(
        available,
        paths["candidates"],
        paths["parse_diagnostics"],
        manual_candidate_csv=candidate_csv,
        prm_evidence_csv=paths["prm_evidence"],
    )
    if mode == "quick":
        candidates = deterministic_quick_subset(candidates)
        candidates.to_csv(paths["candidates"], index=False)
    return candidates


def annotate_stage4(project_root: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    paths = stage4_paths(project_root)
    return annotate_candidates_with_uniprot(
        paths["candidates"],
        paths["annotated"],
        paths["sequence_excluded"],
        paths["interim"] / "uniprot_cache",
        prm_evidence_csv=paths["prm_evidence"],
        prm_mapping_csv=paths["prm_mapping"],
        prm_mapping_excluded_csv=paths["prm_mapping_excluded"],
    )


def structures_stage4(
    project_root: Path,
    *,
    workers: int = 4,
    overwrite: bool = False,
    min_site_plddt: float = 70.0,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    paths = stage4_paths(project_root)
    return build_alphafold_features(
        paths["annotated"],
        paths["structure_features"],
        paths["structure_excluded"],
        paths["interim"] / "alphafold_cache",
        workers=workers,
        overwrite=overwrite,
        min_site_plddt=min_site_plddt,
    )


def score_stage4(project_root: Path, *, refit_model: bool = False) -> pd.DataFrame:
    paths = stage4_paths(project_root)
    return score_candidate_features(
        paths["structure_features"],
        paths["scored"],
        paths["models"],
        project_root / "data/processed/bsa_labels.csv",
        project_root / "data/processed/bsa_structure_features.csv",
        refit_model=refit_model,
    )


def prioritize_stage4(
    project_root: Path,
    *,
    top_n: int = 100,
    max_sites_per_protein: int = 2,
) -> dict[str, pd.DataFrame]:
    paths = stage4_paths(project_root)
    return prioritize_atlas(
        paths["scored"],
        paths["atlas_dir"],
        top_n=top_n,
        max_sites_per_protein=max_sites_per_protein,
    )


def _safe_count(path: Path) -> int:
    if not path.exists():
        return 0
    try:
        return len(pd.read_csv(path))
    except Exception:
        return 0


def audit_stage4(project_root: Path) -> dict[str, Any]:
    paths = stage4_paths(project_root)
    atlas_all = paths["atlas_dir"] / "cml_repairability_atlas_all_sites.csv"
    atlas_in = paths["atlas_dir"] / "cml_repairability_atlas_in_domain.csv"
    pareto = paths["atlas_dir"] / "cml_repairability_pareto_front.csv"
    portfolio = paths["atlas_dir"] / "cml_repairability_priority_portfolio.csv"
    watchlist = paths["atlas_dir"] / "cml_repairability_domain_shift_watchlist.csv"
    audit: dict[str, Any] = {
        "candidate_sites_parsed": _safe_count(paths["candidates"]),
        "sequence_valid_sites": _safe_count(paths["annotated"]),
        "sequence_excluded_sites": _safe_count(paths["sequence_excluded"]),
        "prm_peptides_parsed": _safe_count(paths["prm_evidence"]),
        "prm_peptides_mapped": _safe_count(paths["prm_mapping"]),
        "prm_peptides_excluded": _safe_count(paths["prm_mapping_excluded"]),
        "structure_eligible_sites": _safe_count(paths["structure_features"]),
        "structure_excluded_sites": _safe_count(paths["structure_excluded"]),
        "scored_sites": _safe_count(paths["scored"]),
        "atlas_all_sites": _safe_count(atlas_all),
        "atlas_in_domain_sites": _safe_count(atlas_in),
        "pareto_front_sites": _safe_count(pareto),
        "priority_portfolio_sites": _safe_count(portfolio),
        "domain_shift_watchlist_sites": _safe_count(watchlist),
        "interpretation": (
            "Stage 4 outputs are computational prioritization under domain shift. They are "
            "not experimental validation of repairability in human or mouse proteins."
        ),
        "transfer_safeguard": (
            "AlphaFold pLDDT is not used as crystallographic B-factor. B-factor-derived "
            "features are excluded from transfer scoring; pLDDT is used only for quality."
        ),
    }
    retryable_structure_failures = 0
    cache_incomplete_structure_rows = 0
    if paths["structure_excluded"].exists():
        structure_excluded = pd.read_csv(paths["structure_excluded"])
        reasons = structure_excluded.get(
            "structure_quality_reason", pd.Series("", index=structure_excluded.index)
        ).fillna("").astype(str)
        retryable_structure_failures = int(
            reasons.str.startswith("structure_processing_error:").sum()
        )
        cache_incomplete_structure_rows = int(
            reasons.eq("structure_not_eligible").sum()
        )
    audit["retryable_structure_failures"] = retryable_structure_failures
    audit["cache_incomplete_structure_rows"] = cache_incomplete_structure_rows
    audit["structure_completion_status"] = (
        "complete"
        if retryable_structure_failures == 0 and cache_incomplete_structure_rows == 0
        else "incomplete_retry_structures"
    )

    if atlas_all.exists():
        frame = pd.read_csv(atlas_all)
        audit["species_counts"] = frame["species"].value_counts(dropna=False).to_dict()
        audit["applicability_counts"] = frame["applicability_status"].value_counts().to_dict()
        audit["median_ood_percentile"] = float(frame["ood_percentile"].median())
        audit["median_interval_width"] = float(frame["prediction_interval_width"].median())
        audit["median_site_plddt"] = float(frame["site_plddt"].median())
        in_domain_count = int((frame["applicability_status"] == "in_domain").sum())
        pareto_count = int((frame["pareto_front"] == 1).sum())
        audit["pareto_fraction_of_in_domain"] = (
            pareto_count / in_domain_count if in_domain_count else None
        )
    paths["results"].mkdir(parents=True, exist_ok=True)
    paths["audit_json"].write_text(json.dumps(audit, indent=2, default=str), encoding="utf-8")
    lines = [
        "# CML-RepairNet Stage 4 Audit",
        "",
        f"- Candidate sites parsed: {audit['candidate_sites_parsed']}",
        f"- Sequence-valid sites: {audit['sequence_valid_sites']}",
        f"- PRM peptides parsed/mapped: {audit['prm_peptides_parsed']}/{audit['prm_peptides_mapped']}",
        f"- Structure-eligible sites: {audit['structure_eligible_sites']}",
        f"- Scored sites: {audit['scored_sites']}",
        f"- In-domain atlas sites: {audit['atlas_in_domain_sites']}",
        f"- Pareto-front sites: {audit['pareto_front_sites']}",
        f"- Priority portfolio sites: {audit['priority_portfolio_sites']}",
        f"- Domain-shift watchlist sites: {audit['domain_shift_watchlist_sites']}",
        f"- Retryable structure failures: {audit['retryable_structure_failures']}",
        f"- Cache-incomplete structure rows: {audit['cache_incomplete_structure_rows']}",
        f"- Structure completion status: {audit['structure_completion_status']}",
        "",
        "## Safeguards",
        "",
        f"- {audit['interpretation']}",
        f"- {audit['transfer_safeguard']}",
    ]
    paths["audit_md"].write_text("\n".join(lines) + "\n", encoding="utf-8")
    return audit


def run_stage4_all(
    project_root: Path,
    *,
    mode: str = "quick",
    candidate_csv: Path | None = None,
    workers: int = 2,
    overwrite_downloads: bool = False,
    overwrite_structures: bool = False,
    min_site_plddt: float = 70.0,
    top_n: int = 100,
) -> dict[str, Any]:
    doctor = doctor_stage4(project_root)
    if not doctor["ready"]:
        raise RuntimeError("Stage 4 prerequisites are missing. Run Stage 1–3 data preparation first.")
    paths = stage4_paths(project_root)
    if candidate_csv is None and not all(
        (paths["raw"] / name).exists()
        for name in ("cellular_cml_sites.xlsx", "mouse_organ_cml_sites.xlsx")
    ):
        download_stage4(project_root, overwrite=overwrite_downloads)
    prepare_stage4(project_root, candidate_csv=candidate_csv, mode=mode)
    annotate_stage4(project_root)
    structures_stage4(
        project_root,
        workers=workers,
        overwrite=overwrite_structures,
        min_site_plddt=min_site_plddt,
    )
    score_stage4(project_root)
    prioritize_stage4(project_root, top_n=top_n)
    return audit_stage4(project_root)
