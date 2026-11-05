from __future__ import annotations

from pathlib import Path

import pandas as pd

from cml_repairnet.revision.final_freeze import build_final_freeze, doctor_final_freeze


def _write_fixture(root: Path) -> None:
    phase2 = root / "results/revision_phase2"
    atlas_dir = root / "results/stage4/atlas"
    phase2.mkdir(parents=True)
    atlas_dir.mkdir(parents=True)
    ids = ["human:P1:K1", "mouse:P2:K2"]
    t5 = pd.DataFrame({
        "candidate_id": ids,
        "standard_residual_low": [0.1, 0.2], "standard_residual_high": [0.8, 0.9], "standard_residual_width": [0.7, 0.7],
        "density_ratio_weighted_low": [0.0, 0.1], "density_ratio_weighted_high": [0.9, 1.0], "density_ratio_weighted_width": [0.9, 0.9],
        "localized_distance_weighted_low": [0.1, 0.2], "localized_distance_weighted_high": [0.85, 0.95], "localized_distance_weighted_width": [0.75, 0.75],
        "descriptive_ood_inflated_low": [0.0, 0.0], "descriptive_ood_inflated_high": [1.0, 1.0], "descriptive_ood_inflated_width": [1.0, 1.0],
    })
    t7 = pd.DataFrame({
        "candidate_id": ids, "species": ["human", "mouse"], "uniprot_accession": ["P1", "P2"],
        "gene_symbol": ["G1", "G2"], "site_uniprot": [1, 2], "site_label": ["G1 K1", "G2 K2"],
        "cmlase_response_score": [0.6, 0.7], "candidate_tier": ["Tier 1", "Watchlist"],
        "applicability_status": ["in_domain", "out_of_domain"], "biomedical_interpretation_class": ["secondary_hypothesis_only", "higher_context_plausibility"],
        "accessibility_layer": ["C_intracellular", "A_extracellular_or_surface"], "compartment_primary": ["cytosolic", "extracellular_secreted"],
    })
    t5.to_csv(phase2 / "phase2_t5_atlas_interval_sensitivity.csv", index=False)
    t7.to_csv(phase2 / "phase2_t7_compartment_layered_atlas.csv", index=False)
    pd.DataFrame({
        "species": ["mouse"], "query_name": ["evidence_rich"], "term_name": ["mitochondrion"],
        "background_robust_informative": [True], "ontology_root_term": [False]
    }).to_csv(phase2 / "phase2_t8_enrichment_background_robustness.csv", index=False)
    pd.DataFrame({"species": ["mouse"], "query_name": ["portfolio"], "variant": ["all_evidence"]}).to_csv(
        phase2 / "phase2_t8_string_study_bias_summary.csv", index=False
    )
    pd.DataFrame({"candidate_id": [ids[0]]}).to_csv(atlas_dir / "cml_repairability_priority_portfolio.csv", index=False)
    pd.DataFrame({"candidate_id": [ids[0]]}).to_csv(atlas_dir / "cml_repairability_pareto_front.csv", index=False)


def test_doctor_and_build(tmp_path: Path) -> None:
    _write_fixture(tmp_path)
    output = tmp_path / "results/final_scientific_freeze"
    doctor = doctor_final_freeze(tmp_path, output)
    assert doctor["status"] == "PASS"
    manifest = build_final_freeze(tmp_path, output)
    assert manifest["status"] == "FROZEN"
    assert manifest["model_retrained"] is False
    final = pd.read_csv(output / "cml_repairnet_final_atlas_all_sites.csv")
    assert "cmlase_response_score" in final
    assert "bsa_calibrated_sensitivity_lower" in final
    assert len(final) == 2
    informative = pd.read_csv(output / "cml_repairnet_final_background_robust_enrichment.csv")
    assert informative["term_name"].tolist() == ["mitochondrion"]


def test_doctor_fails_on_id_mismatch(tmp_path: Path) -> None:
    _write_fixture(tmp_path)
    path = tmp_path / "results/revision_phase2/phase2_t5_atlas_interval_sensitivity.csv"
    t5 = pd.read_csv(path)
    t5.loc[0, "candidate_id"] = "human:OTHER:K1"
    t5.to_csv(path, index=False)
    doctor = doctor_final_freeze(tmp_path, tmp_path / "out")
    assert doctor["status"] == "FAIL"
