from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from cml_repairnet.revision.phase2 import (
    TRANSFER_SAFE_FEATURES,
    classify_subcellular_location,
    run_t5_conformal_shift_sensitivity,
    run_t6_boundary_shrinkage,
    run_t7_compartment_layering,
    run_t8_background_and_network_sensitivity,
    stable_gene_set_fingerprint,
    is_ontology_root_term,
    weighted_quantile,
)


def test_weighted_quantile_respects_heavy_weight() -> None:
    values = np.array([0.1, 0.2, 0.9])
    assert weighted_quantile(values, np.array([1.0, 1.0, 20.0]), 0.5) == 0.9


def test_compartment_classification_is_conservative() -> None:
    secreted = classify_subcellular_location("Secreted, extracellular space")
    assert secreted["accessibility_layer"] == "A_extracellular_or_surface"
    mixed = classify_subcellular_location("Secreted. Cytoplasm. Nucleus.")
    assert mixed["accessibility_layer"] == "U_mixed_or_unknown"
    mitochondrial = classify_subcellular_location("Mitochondrial matrix")
    assert mitochondrial["accessibility_layer"] == "D_mitochondrial"


def _synthetic_inputs(tmp_path: Path) -> tuple[Path, Path, Path, Path]:
    rng = np.random.default_rng(7)
    n = 20
    labels = pd.DataFrame({"site_uniprot": np.arange(1, n + 1), "repair_published": np.linspace(0, 1, n)})
    structure = pd.DataFrame({"site_uniprot": np.arange(1, n + 1)})
    for index, column in enumerate(TRANSFER_SAFE_FEATURES):
        structure[column] = rng.normal(index, 1, n)
    predictions = pd.DataFrame(
        {
            "site_uniprot": np.arange(1, n + 1),
            "site_paper": np.arange(1, n + 1),
            "partition_id": "ward_g4",
            "feature_set": "structure_transfer_safe",
            "model": "extra_trees",
            "y_true": labels["repair_published"],
            "y_pred": 0.25 + 0.5 * labels["repair_published"] + rng.normal(0, 0.04, n),
        }
    )
    atlas = pd.DataFrame(
        {
            "candidate_id": [f"human:P{i}:K1" for i in range(40)],
            "species": "human",
            "gene_names": [f"G{i}" for i in range(40)],
            "uniprot_accession": [f"P{i}" for i in range(40)],
            "subcellular_location": ["Secreted" if i % 2 else "Cytoplasm" for i in range(40)],
            "repairability_extra_trees": rng.uniform(0.2, 0.9, 40),
            "ood_percentile": np.linspace(0, 1, 40),
            "applicability_status": ["in_domain"] * 20 + ["out_of_domain"] * 20,
            "sequence": ["MKKKAAAA"] * 40,
            "protein_length": 100,
            "evidence_score_norm": rng.uniform(0, 1, 40),
            "structure_quality_score": rng.uniform(0.5, 1, 40),
            "structure_status": "eligible",
        }
    )
    for index, column in enumerate(TRANSFER_SAFE_FEATURES):
        atlas[column] = rng.normal(index + 0.2, 1.3, 40)
    paths = []
    for name, frame in (
        ("labels.csv", labels),
        ("structure.csv", structure),
        ("predictions.csv", predictions),
        ("atlas.csv", atlas),
    ):
        path = tmp_path / name
        frame.to_csv(path, index=False)
        paths.append(path)
    return tuple(paths)  # type: ignore[return-value]


def test_t5_and_t6_end_to_end(tmp_path: Path) -> None:
    labels, structure, predictions, atlas = _synthetic_inputs(tmp_path)
    output = tmp_path / "out"
    t5 = run_t5_conformal_shift_sensitivity(labels, structure, predictions, atlas, output)
    assert set(t5["coverage"]["method"]) == {
        "standard_residual",
        "density_ratio_weighted",
        "localized_distance_weighted",
    }
    assert (output / "phase2_t5_atlas_interval_sensitivity.csv").exists()
    # Add repeated partitions to exercise aggregation in T6.
    pred = pd.read_csv(predictions)
    repeated = pd.concat(
        [pred.assign(partition_id=f"p{i}", y_pred=pred["y_pred"] + i * 0.002) for i in range(3)],
        ignore_index=True,
    )
    repeated.to_csv(predictions, index=False)
    t6 = run_t6_boundary_shrinkage(predictions, output, bootstrap_iterations=100)
    assert "raw_partition_mean" in set(t6["summary"]["prediction_variant"])
    assert len(t6["site"]) == 20


def test_t7_layering(tmp_path: Path) -> None:
    _, _, _, atlas = _synthetic_inputs(tmp_path)
    result = run_t7_compartment_layering(atlas, tmp_path / "out")
    assert "accessibility_layer" in result["layered"]
    assert not result["layered"]["hard_exclusion"].any()


def test_t8_online_term_robustness_first_pass(tmp_path: Path, monkeypatch) -> None:
    _, _, _, atlas = _synthetic_inputs(tmp_path)

    def fake_gprofiler_request(*args, **kwargs):
        return ([{
            "source": "GO:BP",
            "native": "GO:0000001",
            "name": "synthetic process",
            "p_value": 0.01,
            "significant": True,
            "query_size": 4,
            "intersection_size": 2,
            "term_size": 10,
            "effective_domain_size": 40,
            "precision": 0.5,
            "recall": 0.2,
        }], {"data_version": {"synthetic": "1"}})

    monkeypatch.setattr(
        "cml_repairnet.revision.phase2._gprofiler_request",
        fake_gprofiler_request,
    )
    out = tmp_path / "out_t8"
    result = run_t8_background_and_network_sensitivity(
        atlas,
        out,
        controls_per_query=2,
        run_online=True,
    )
    robustness = result["term_robustness"]
    assert "significant_background_list" in robustness.columns
    assert robustness["significant_background_list"].str.contains("matched").any()
    assert "distinct_backgrounds_tested" in robustness.columns
    assert "background_robust_informative" in robustness.columns
    assert (robustness["distinct_backgrounds_tested"] <= robustness["background_labels_tested"]).all()
    assert (out / "phase2_t8_enrichment_background_robustness.csv").exists()


def test_t8_background_fingerprint_and_root_filter() -> None:
    assert stable_gene_set_fingerprint(["B", "A", "A"]) == stable_gene_set_fingerprint(["A", "B"])
    assert is_ontology_root_term("KEGG", "KEGG:00000", "KEGG root term")
    assert is_ontology_root_term("WP", "WP:000000", "WIKIPATHWAYS")
    assert not is_ontology_root_term("KEGG", "KEGG:04146", "Peroxisome")
