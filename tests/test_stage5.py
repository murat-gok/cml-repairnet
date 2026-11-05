from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from cml_repairnet.biology.annotations import family_summary, normalize_uniprot_annotations
from cml_repairnet.biology.enrichment import build_enrichment_queries, flatten_gprofiler_response
from cml_repairnet.biology.modules import build_candidate_modules
from cml_repairnet.biology.networks import graph_tables, normalize_string_edges
from cml_repairnet.biology.orthology import cross_species_concordance, flatten_ortholog_response
from cml_repairnet.biology.pipeline import doctor_stage5, modules_stage5, report_stage5
from cml_repairnet.biology.schema import first_gene_symbol, load_atlas, validate_atlas


def atlas_frame(n: int = 12) -> pd.DataFrame:
    species = ["human"] * (n // 2) + ["mouse"] * (n - n // 2)
    genes = [f"H{i}" for i in range(n // 2)] + [f"M{i}" for i in range(n - n // 2)]
    values = np.linspace(0.25, 0.85, n)
    return pd.DataFrame(
        {
            "candidate_id": [f"{s}:P{i:05d}:K{i + 1}" for i, s in enumerate(species)],
            "species": species,
            "organism_id": [9606 if s == "human" else 10090 for s in species],
            "uniprot_accession": [f"P{i:05d}" for i in range(n)],
            "site_uniprot": np.arange(1, n + 1),
            "gene_names": [f"{g} SYN" for g in genes],
            "protein_name": [f"Protein {g}" for g in genes],
            "repairability_extra_trees": values,
            "repairability_lightgbm": np.clip(values + 0.03, 0, 1),
            "repairability_ridge": np.clip(values - 0.02, 0, 1),
            "consensus_rank_percentile": np.linspace(0.1, 1.0, n),
            "evidence_score_raw": np.tile([1.0, 2.0, 3.0], n // 3 + 1)[:n],
            "transfer_reliability_score": np.linspace(0.4, 0.95, n),
            "weighted_priority_score": np.linspace(0.3, 0.9, n),
            "applicability_status": ["in_domain"] * n,
            "prm_validated": [i % 4 == 0 for i in range(n)],
            "age_associated": [i % 3 == 0 for i in range(n)],
            "ortholog_overlap": [i % 2 == 0 for i in range(n)],
            "model_disagreement": np.linspace(0.02, 0.15, n),
            "ood_percentile": np.linspace(0.0, 0.6, n),
            "prediction_interval_width": np.linspace(0.5, 0.9, n),
            "site_plddt": np.linspace(75, 98, n),
            "relative_sasa_lysine": np.linspace(0.1, 0.8, n),
            "neighbor_cysteine_count_8a": np.arange(n) % 4,
            "neighbor_hydrophobic_fraction_8a": np.linspace(0.2, 0.7, n),
            "portfolio_rank": np.arange(1, n + 1),
            "evidence_score_norm": np.linspace(0.1, 1.0, n),
        }
    )


def test_gene_symbol_and_validation() -> None:
    frame = atlas_frame()
    validate_atlas(frame)
    assert first_gene_symbol("PMVK PMKI") == "PMVK"
    assert first_gene_symbol(np.nan) == ""


def test_load_atlas_adds_stage5_columns(tmp_path: Path) -> None:
    path = tmp_path / "atlas.csv"
    atlas_frame().to_csv(path, index=False)
    loaded = load_atlas(path)
    assert "gene_symbol" in loaded
    assert loaded.iloc[0]["site_label"].endswith("K1")


def test_enrichment_queries_use_custom_species_background() -> None:
    atlas = atlas_frame(16)
    atlas["gene_symbol"] = atlas["gene_names"].map(first_gene_symbol)
    portfolio = atlas.groupby("species", group_keys=False).head(4)
    pareto = atlas.groupby("species", group_keys=False).head(2)
    queries = build_enrichment_queries(atlas_in_domain=atlas, portfolio=portfolio, pareto=pareto)
    assert queries
    for query in queries:
        expected_background = set(atlas.loc[atlas["species"] == query.species, "gene_symbol"])
        assert set(query.background) == expected_background
        assert set(query.genes) <= expected_background


def test_flatten_gprofiler_response_computes_fold_enrichment() -> None:
    response = {
        "result": [
            {
                "source": "GO:BP",
                "native": "GO:1",
                "name": "Test",
                "description": "Test term",
                "p_value": 0.01,
                "significant": True,
                "query_size": 10,
                "intersection_size": 5,
                "term_size": 20,
                "effective_domain_size": 100,
                "precision": 0.5,
                "recall": 0.25,
                "parents": [],
                "source_order": 1,
                "group_id": 1,
            }
        ]
    }
    frame = flatten_gprofiler_response(
        response,
        species="human",
        query_name="portfolio",
        query_description="x",
        query_gene_count=10,
        background_gene_count=100,
    )
    assert frame.loc[0, "fold_enrichment"] == 2.5


def test_uniprot_annotation_family_summary() -> None:
    raw = pd.DataFrame(
        {
            "Entry": ["P1", "P2"],
            "Protein families": ["Alpha family; Beta family", "Alpha family"],
            "Gene Names": ["A", "B"],
        }
    )
    normalized = normalize_uniprot_annotations(raw)
    families = family_summary(normalized)
    alpha = families.loc[families["protein_family"] == "Alpha family"].iloc[0]
    assert alpha["proteins"] == 2


def test_string_edge_normalization_and_network_metrics() -> None:
    raw = pd.DataFrame(
        {
            "preferredName_A": ["A", "B", "A"],
            "preferredName_B": ["B", "C", "B"],
            "stringId_A": ["1.A", "1.B", "1.A"],
            "stringId_B": ["1.B", "1.C", "1.B"],
            "score": [0.9, 0.8, 0.9],
        }
    )
    edges = normalize_string_edges(raw, species="human", query_name="portfolio")
    nodes, metrics = graph_tables(edges, species="human", query_name="portfolio")
    assert len(edges) == 2
    assert len(nodes) == 3
    assert metrics.loc[0, "edges"] == 2


def test_ortholog_concordance() -> None:
    frame = atlas_frame(8)
    frame["gene_symbol"] = frame["gene_names"].map(first_gene_symbol)
    frame.loc[frame["species"] == "human", "gene_symbol"] = ["A", "B", "C", "D"]
    frame.loc[frame["species"] == "mouse", "gene_symbol"] = ["a", "b", "c", "d"]
    mapping = flatten_ortholog_response(
        {
            "result": [
                {"incoming": "A", "converted": "a"},
                {"incoming": "B", "converted": "b"},
                {"incoming": "C", "converted": "c"},
            ]
        },
        source_species="human",
        target_species="mouse",
    )
    pairs, summary = cross_species_concordance(frame, mapping)
    assert len(pairs) == 3
    assert summary.loc[0, "mapped_portfolio_pairs"] == 3


def test_candidate_modules_are_deterministic() -> None:
    frame = atlas_frame(24)
    a1, p1, s1 = build_candidate_modules(frame, random_state=7)
    a2, _p2, s2 = build_candidate_modules(frame, random_state=7)
    assert a1["module_id"].tolist() == a2["module_id"].tolist()
    assert len(p1) >= 2
    assert int(s1.loc[s1["silhouette"].idxmax(), "k"]) == int(
        s2.loc[s2["silhouette"].idxmax(), "k"]
    )


def test_stage5_offline_core_pipeline(tmp_path: Path) -> None:
    atlas_dir = tmp_path / "results/stage4/atlas"
    atlas_dir.mkdir(parents=True)
    frame = atlas_frame(20)
    frame.to_csv(atlas_dir / "cml_repairability_atlas_in_domain.csv", index=False)
    frame.head(10).to_csv(atlas_dir / "cml_repairability_priority_portfolio.csv", index=False)
    frame.head(5).to_csv(atlas_dir / "cml_repairability_pareto_front.csv", index=False)
    frame.iloc[0:0].to_csv(atlas_dir / "cml_repairability_domain_shift_watchlist.csv", index=False)
    assert doctor_stage5(tmp_path)["ready"]
    result = modules_stage5(tmp_path)
    assert result["assignments"] == 10
    audit = report_stage5(tmp_path)
    assert audit["stage5_completion_status"] == "complete"
    assert (tmp_path / "results/stage5/stage5_biological_interpretation.md").exists()


def test_ortholog_official_schema_uses_name_as_target_symbol() -> None:
    response = {
        "result": [
            {
                "incoming": "ACTB",
                "converted": "ENSG00000075624",
                "ortholog_ensg": "ENSMUSG00000029580",
                "name": "Actb",
                "description": "actin, beta",
                "namespaces": "HGNC,UNIPROT_GN",
            }
        ]
    }
    frame = flatten_ortholog_response(response, source_species="human", target_species="mouse")
    assert frame.loc[0, "source_gene"] == "ACTB"
    assert frame.loc[0, "target_gene"] == "Actb"
    assert frame.loc[0, "source_ensembl"] == "ENSG00000075624"
    assert frame.loc[0, "target_ensembl"] == "ENSMUSG00000029580"


def test_ortholog_concordance_is_case_insensitive() -> None:
    frame = atlas_frame(8)
    frame["gene_symbol"] = frame["gene_names"].map(first_gene_symbol)
    frame.loc[frame["species"] == "human", "gene_symbol"] = ["ACTB", "PMVK", "PFN1", "STAT1"]
    frame.loc[frame["species"] == "mouse", "gene_symbol"] = ["Actb", "Pmvk", "Pfn1", "Stat1"]
    response = {
        "result": [
            {"incoming": "ACTB", "converted": "ENSG1", "ortholog_ensg": "ENSMUSG1", "name": "Actb"},
            {"incoming": "PMVK", "converted": "ENSG2", "ortholog_ensg": "ENSMUSG2", "name": "Pmvk"},
            {"incoming": "PFN1", "converted": "ENSG3", "ortholog_ensg": "ENSMUSG3", "name": "Pfn1"},
        ]
    }
    mapping = flatten_ortholog_response(response, source_species="human", target_species="mouse")
    pairs, summary = cross_species_concordance(frame, mapping)
    assert len(pairs) == 3
    assert summary.loc[0, "mapped_portfolio_pairs"] == 3
