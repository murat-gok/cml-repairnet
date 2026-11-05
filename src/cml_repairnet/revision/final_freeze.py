from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import pandas as pd

FINAL_VERSION = "0.12.0"


@dataclass(frozen=True)
class FinalFreezePaths:
    t5_atlas: Path
    t7_atlas: Path
    t8_enrichment: Path
    t8_network: Path
    stage4_portfolio: Path | None
    stage4_pareto: Path | None
    stage4_watchlist: Path | None
    confirmatory_report: Path | None
    external_report: Path | None


def first_existing(paths: Iterable[Path]) -> Path | None:
    for path in paths:
        if path.exists():
            return path
    return None


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def discover_paths(root: Path) -> FinalFreezePaths:
    root = root.resolve()
    phase2 = root / "results/revision_phase2"
    return FinalFreezePaths(
        t5_atlas=phase2 / "phase2_t5_atlas_interval_sensitivity.csv",
        t7_atlas=phase2 / "phase2_t7_compartment_layered_atlas.csv",
        t8_enrichment=phase2 / "phase2_t8_enrichment_background_robustness.csv",
        t8_network=phase2 / "phase2_t8_string_study_bias_summary.csv",
        stage4_portfolio=first_existing((
            root / "results/stage4/atlas/cml_repairability_priority_portfolio.csv",
            root / "supplementary_data/Table_S3_stage4_priority_portfolio.csv",
        )),
        stage4_pareto=first_existing((
            root / "results/stage4/atlas/cml_repairability_pareto_front.csv",
            root / "supplementary_data/Table_S4_stage4_pareto_front.csv",
        )),
        stage4_watchlist=first_existing((
            root / "results/stage4/atlas/cml_repairability_domain_shift_watchlist.csv",
            root / "results/stage4/atlas/cml_repairability_watchlist.csv",
            root / "supplementary_data/Table_S5_stage4_domain_shift_watchlist.csv",
        )),
        confirmatory_report=first_existing((
            root / "results/revision_confirmatory/revision_confirmatory_report.md",
            root / "results/revision_tests/revision_confirmatory_report.md",
        )),
        external_report=first_existing((
            root / "results/external_scaffold_validation/external_scaffold_report.md",
            root / "results/external_scaffold_validation/report.md",
        )),
    )


def _read(path: Path, required: tuple[str, ...] = ()) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(path)
    frame = pd.read_csv(path)
    missing = [column for column in required if column not in frame.columns]
    if missing:
        raise ValueError(f"{path} is missing required columns: {missing}")
    return frame


def _candidate_ids(frame: pd.DataFrame) -> set[str]:
    return set(frame["candidate_id"].dropna().astype(str))


def doctor_final_freeze(root: Path, output_dir: Path) -> dict:
    paths = discover_paths(root)
    required = {
        "t5_atlas": paths.t5_atlas,
        "t7_atlas": paths.t7_atlas,
        "t8_enrichment": paths.t8_enrichment,
        "t8_network": paths.t8_network,
    }
    optional = {
        "stage4_portfolio": paths.stage4_portfolio,
        "stage4_pareto": paths.stage4_pareto,
        "stage4_watchlist": paths.stage4_watchlist,
        "confirmatory_report": paths.confirmatory_report,
        "external_report": paths.external_report,
    }
    report: dict = {
        "version": FINAL_VERSION,
        "required": {
            name: {"path": str(path), "exists": path.exists()}
            for name, path in required.items()
        },
        "optional": {
            name: {"path": str(path) if path else None, "exists": bool(path and path.exists())}
            for name, path in optional.items()
        },
    }
    report["status"] = "PASS" if all(item["exists"] for item in report["required"].values()) else "FAIL"
    if report["status"] == "PASS":
        t5 = _read(paths.t5_atlas, ("candidate_id",))
        t7 = _read(paths.t7_atlas, ("candidate_id", "cmlase_response_score"))
        report["n_t5_candidates"] = len(t5)
        report["n_t7_candidates"] = len(t7)
        report["candidate_id_sets_equal"] = _candidate_ids(t5) == _candidate_ids(t7)
        report["candidate_id_duplicates_t5"] = int(t5["candidate_id"].duplicated().sum())
        report["candidate_id_duplicates_t7"] = int(t7["candidate_id"].duplicated().sum())
        if not report["candidate_id_sets_equal"] or report["candidate_id_duplicates_t5"] or report["candidate_id_duplicates_t7"]:
            report["status"] = "FAIL"
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "final_freeze_doctor.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    return report


def _merge_t5_t7(t5: pd.DataFrame, t7: pd.DataFrame) -> pd.DataFrame:
    t5_columns = [
        "candidate_id",
        "standard_residual_low",
        "standard_residual_high",
        "standard_residual_width",
        "density_ratio_weighted_low",
        "density_ratio_weighted_high",
        "density_ratio_weighted_width",
        "localized_distance_weighted_low",
        "localized_distance_weighted_high",
        "localized_distance_weighted_width",
        "descriptive_ood_inflated_low",
        "descriptive_ood_inflated_high",
        "descriptive_ood_inflated_width",
    ]
    missing = [column for column in t5_columns if column not in t5.columns]
    if missing:
        raise ValueError(f"T5 atlas missing final-freeze columns: {missing}")
    merged = t7.merge(t5[t5_columns], on="candidate_id", how="left", validate="one_to_one")
    if merged["standard_residual_width"].isna().any():
        raise ValueError("T5/T7 merge left candidates without sensitivity bands")

    merged["bsa_calibrated_sensitivity_lower"] = merged["standard_residual_low"]
    merged["bsa_calibrated_sensitivity_upper"] = merged["standard_residual_high"]
    merged["bsa_calibrated_sensitivity_width"] = merged["standard_residual_width"]
    merged["ood_inflated_descriptive_lower"] = merged["descriptive_ood_inflated_low"]
    merged["ood_inflated_descriptive_upper"] = merged["descriptive_ood_inflated_high"]
    merged["ood_inflated_descriptive_width"] = merged["descriptive_ood_inflated_width"]
    merged["response_score_interpretation"] = (
        "relative BSA-derived prioritization score; not an absolute repair fraction"
    )
    merged["uncertainty_interpretation"] = (
        "BSA-calibrated sensitivity band; no finite-sample proteome-transfer coverage guarantee"
    )
    merged["transfer_claim_boundary"] = (
        "hypothesis prioritization only; cross-protein site-level validity not established"
    )
    return merged


def _curated_columns(frame: pd.DataFrame) -> list[str]:
    preferred = [
        "candidate_id", "species", "uniprot_accession", "gene_symbol", "gene_names",
        "protein_name", "site_uniprot", "site_label", "evidence_types", "tissues_conditions",
        "prm_validated", "age_tested", "age_associated", "ortholog_overlap",
        "cmlase_response_score", "consensus_rank_percentile", "weighted_priority_score",
        "applicability_status", "ood_percentile", "structure_quality_status", "site_plddt",
        "bsa_calibrated_sensitivity_lower", "bsa_calibrated_sensitivity_upper",
        "bsa_calibrated_sensitivity_width", "ood_inflated_descriptive_lower",
        "ood_inflated_descriptive_upper", "ood_inflated_descriptive_width",
        "candidate_tier", "tier_reason", "pareto_front", "portfolio_selection_frequency",
        "pareto_membership_frequency", "accessibility_layer", "compartment_primary",
        "compartment_ambiguous", "biomedical_interpretation_class", "interpretation_note",
        "response_score_interpretation", "uncertainty_interpretation", "transfer_claim_boundary",
    ]
    return [column for column in preferred if column in frame.columns]


def _subset_by_ids(frame: pd.DataFrame, source: Path | None) -> pd.DataFrame:
    if source is None or not source.exists():
        return frame.iloc[0:0].copy()
    selection = _read(source, ("candidate_id",))
    ids = _candidate_ids(selection)
    return frame.loc[frame["candidate_id"].astype(str).isin(ids)].copy()


def build_final_freeze(root: Path, output_dir: Path) -> dict:
    doctor = doctor_final_freeze(root, output_dir)
    if doctor["status"] != "PASS":
        raise RuntimeError("Final-freeze doctor failed; inspect final_freeze_doctor.json")
    paths = discover_paths(root)
    t5 = _read(paths.t5_atlas, ("candidate_id",))
    t7 = _read(paths.t7_atlas, ("candidate_id", "cmlase_response_score", "candidate_tier"))
    final = _merge_t5_t7(t5, t7)
    curated = final[_curated_columns(final)].copy()

    output_dir.mkdir(parents=True, exist_ok=True)
    final_path = output_dir / "cml_repairnet_final_atlas_all_sites.csv"
    curated_path = output_dir / "cml_repairnet_final_atlas_curated.csv"
    final.to_csv(final_path, index=False)
    curated.to_csv(curated_path, index=False)

    manuscript = final.loc[final["candidate_tier"].astype(str) != "Watchlist"].copy()
    manuscript[_curated_columns(manuscript)].to_csv(
        output_dir / "cml_repairnet_final_manuscript_candidates.csv", index=False
    )
    _subset_by_ids(final, paths.stage4_portfolio)[_curated_columns(final)].to_csv(
        output_dir / "cml_repairnet_final_priority_portfolio.csv", index=False
    )
    _subset_by_ids(final, paths.stage4_pareto)[_curated_columns(final)].to_csv(
        output_dir / "cml_repairnet_final_pareto_candidates.csv", index=False
    )
    watch = _subset_by_ids(final, paths.stage4_watchlist)
    if watch.empty:
        watch = final.loc[final["applicability_status"].astype(str) == "out_of_domain"].copy()
    watch[_curated_columns(final)].to_csv(
        output_dir / "cml_repairnet_final_domain_shift_watchlist.csv", index=False
    )

    enrichment = _read(paths.t8_enrichment, ("background_robust_informative", "ontology_root_term"))
    informative = enrichment.loc[enrichment["background_robust_informative"].fillna(False).astype(bool)].copy()
    informative.to_csv(output_dir / "cml_repairnet_final_background_robust_enrichment.csv", index=False)
    network = _read(paths.t8_network)
    network.to_csv(output_dir / "cml_repairnet_final_string_sensitivity.csv", index=False)

    counts = {
        "all_sites": int(len(final)),
        "species": final["species"].value_counts(dropna=False).to_dict(),
        "applicability_status": final["applicability_status"].value_counts(dropna=False).to_dict(),
        "candidate_tier": final["candidate_tier"].value_counts(dropna=False).to_dict(),
        "biomedical_interpretation_class": final["biomedical_interpretation_class"].value_counts(dropna=False).to_dict(),
        "informative_background_robust_terms": int(len(informative)),
    }
    output_files = sorted(path for path in output_dir.iterdir() if path.is_file() and path.name != "final_freeze_manifest.json")
    manifest = {
        "version": FINAL_VERSION,
        "status": "FROZEN",
        "freeze_type": "annotation-and-claim freeze over unchanged Stage 4-6 scores, ranks, tiers and selections",
        "model_retrained": False,
        "scores_recalibrated": False,
        "portfolio_reselected": False,
        "primary_score_name": "cmlase_response_score",
        "primary_score_interpretation": "relative BSA-derived prioritization score; not an absolute repair fraction",
        "uncertainty_claim": "BSA-calibrated sensitivity bands only; no proteome-transfer coverage guarantee",
        "external_validation_claim": "not established; cross-scaffold analysis is an applicability stress test",
        "enrichment_claim": "limited mouse evidence-rich mitochondrial, myelin and broad metabolic context only",
        "counts": counts,
        "input_sha256": {
            name: sha256_file(path)
            for name, path in {
                "t5_atlas": paths.t5_atlas,
                "t7_atlas": paths.t7_atlas,
                "t8_enrichment": paths.t8_enrichment,
                "t8_network": paths.t8_network,
            }.items()
        },
        "output_sha256": {path.name: sha256_file(path) for path in output_files},
    }
    (output_dir / "final_freeze_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")

    report_lines = [
        "# CML-RepairNet final scientific freeze",
        "",
        f"Version: {FINAL_VERSION}",
        "",
        "## Frozen interpretation",
        "- The model is used for relative CMLase-response prioritization, not absolute repair-fraction prediction.",
        "- Reported sensitivity bands are calibrated within BSA and do not guarantee cross-protein coverage.",
        "- Candidate tiers describe computational priority; biomedical interpretation classes separately describe accessibility and delivery context.",
        "- Cross-scaffold results are an applicability stress test, not positive site-level external validation.",
        "- Functional enrichment is retained only when significant in at least two distinct background gene sets including the matched background, excluding ontology roots.",
        "",
        "## Frozen counts",
        f"- All atlas sites: {counts['all_sites']}",
        f"- Informative background-robust enrichment terms: {counts['informative_background_robust_terms']}",
        "",
        "## Manuscript title",
        "Structure-aware modeling of site-specific carboxymethyllysine repair: a framework for prioritizing age-related protein damage",
    ]
    (output_dir / "FINAL_SCIENTIFIC_FREEZE_REPORT.md").write_text("\n".join(report_lines) + "\n", encoding="utf-8")
    return manifest
