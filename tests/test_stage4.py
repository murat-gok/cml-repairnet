from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from cml_repairnet.atlas.prioritize import (
    add_pareto_ranking,
    add_priority_components,
    non_dominated_sort,
    select_diverse_portfolio,
)
from cml_repairnet.atlas.sources import (
    aggregate_candidate_evidence,
    extract_accessions,
    extract_sites,
    normalize_workbook,
)
from cml_repairnet.atlas.structures import (
    TRANSFER_FEATURES,
    AlphaFoldEntry,
    extract_site_features_from_structure,
)


def test_extract_accessions_and_sites() -> None:
    assert extract_accessions("sp|P02769|ALBU_BOVIN; Q9Y6R7-2") == ["P02769", "Q9Y6R7-2"]
    assert extract_sites("K58; Lys 92; 154") == [58, 92, 154]


def test_workbook_parser_and_aggregation(tmp_path: Path) -> None:
    workbook = tmp_path / "supp.xlsx"
    frame = pd.DataFrame(
        {
            "Majority protein IDs": ["sp|P12345|TEST_HUMAN", "sp|Q8N123|TEST_MOUSE"],
            "Positions within proteins": ["K10", "K20"],
            "Modification": ["CML", "CML"],
            "Sample": ["HUVEC", "MEF"],
            "PRM validated": ["yes", "no"],
        }
    )
    with pd.ExcelWriter(workbook) as writer:
        frame.to_excel(writer, sheet_name="CML sites", index=False)
    parsed, diagnostics = normalize_workbook(workbook, source_label="cellular")
    assert len(parsed) == 2
    assert set(parsed["species"]) == {"human", "mouse"}
    aggregated = aggregate_candidate_evidence(parsed)
    assert len(aggregated) == 2
    assert aggregated["candidate_id"].str.contains(":K").all()
    assert diagnostics[0].status == "parsed"


def _write_small_pdb(path: Path) -> None:
    residues = [
        (1, "CYS", 0.0, 0.0, 0.0, 91.0),
        (2, "LYS", 3.0, 0.0, 0.0, 88.0),
        (3, "PHE", 6.0, 0.0, 0.0, 85.0),
        (4, "ALA", 10.0, 0.0, 0.0, 80.0),
    ]
    lines = []
    for serial, (number, name, x, y, z, bfactor) in enumerate(residues, start=1):
        lines.append(
            f"ATOM  {serial:5d}  CA  {name:>3s} A{number:4d}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00{bfactor:6.2f}           C"
        )
    lines.extend(["TER", "END"])
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def test_structure_feature_extraction_is_transfer_safe(tmp_path: Path) -> None:
    pdb = tmp_path / "model.pdb"
    _write_small_pdb(pdb)
    candidates = pd.DataFrame(
        {
            "candidate_id": ["human:P12345:K2"],
            "species": ["human"],
            "uniprot_accession": ["P12345"],
            "site_uniprot": [2],
        }
    )
    entry = AlphaFoldEntry("P12345", "AF-P12345-F1", 1, 4, "file://unused", None, 4, 86.0)
    result = extract_site_features_from_structure(candidates, pdb, entry)
    assert result.loc[0, "structure_status"] == "eligible"
    assert result.loc[0, "neighbor_cysteine_count_8a"] == 1
    assert result.loc[0, "neighbor_aromatic_count_8a"] == 1
    assert set(TRANSFER_FEATURES).issubset(result.columns)
    assert "bfactor_ca" not in TRANSFER_FEATURES


def test_non_dominated_sort_and_portfolio() -> None:
    values = np.array([[1.0, 1.0], [0.8, 0.9], [0.9, 0.7], [0.2, 0.2]])
    ranks = non_dominated_sort(values)
    assert ranks[0] == 1
    assert ranks[-1] > 1

    frame = pd.DataFrame(
        {
            "candidate_id": ["a", "b", "c", "d"],
            "species": ["human", "human", "mouse", "mouse"],
            "uniprot_accession": ["P1", "P1", "Q1", "Q2"],
            "repairability_extra_trees": [0.9, 0.8, 0.7, 0.6],
            "repairability_lightgbm": [0.8, 0.7, 0.9, 0.5],
            "repairability_ridge": [0.85, 0.75, 0.8, 0.55],
            "consensus_rank_percentile": [1.0, 0.7, 1.0, 0.5],
            "evidence_score_raw": [3.0, 2.0, 3.0, 1.0],
            "structure_quality_score": [0.9, 0.8, 0.9, 0.75],
            "ood_percentile": [0.1, 0.2, 0.1, 0.3],
            "prediction_interval_width": [0.3, 0.4, 0.3, 0.5],
            "model_disagreement": [0.05, 0.04, 0.08, 0.03],
            "applicability_status": ["in_domain"] * 4,
        }
    )
    ranked = add_pareto_ranking(add_priority_components(frame))
    portfolio = select_diverse_portfolio(ranked, top_n=3, max_sites_per_protein=1)
    assert len(portfolio) == 3
    assert portfolio.groupby(["species", "uniprot_accession"]).size().max() == 1


def test_workbook_parser_detects_preamble_header(tmp_path: Path) -> None:
    workbook = tmp_path / "preamble.xlsx"
    raw = pd.DataFrame(
        [
            ["Supplementary Data 1", None, None, None, None],
            ["CML sites identified in HUVEC", None, None, None, None],
            [
                "Majority protein IDs",
                "Positions within proteins",
                "Modification",
                "Sample",
                "PRM validated",
            ],
            ["sp|P12345|TEST_HUMAN", "K10", "CML", "HUVEC", "yes"],
            ["sp|Q8N123|TEST_MOUSE", "K20", "CML", "MEF", "no"],
        ]
    )
    with pd.ExcelWriter(workbook) as writer:
        raw.to_excel(writer, sheet_name="CML sites", index=False, header=False)
    parsed, diagnostics = normalize_workbook(workbook, source_label="cellular")
    assert len(parsed) == 2
    assert diagnostics[0].header_row == 3
    assert diagnostics[0].parsed_sites == 2


def test_workbook_parser_detects_two_row_header(tmp_path: Path) -> None:
    workbook = tmp_path / "two_header.xlsx"
    raw = pd.DataFrame(
        [
            ["Notes", None, None, None],
            [
                "Protein information",
                "Modification information",
                "Modification information",
                "Evidence",
            ],
            ["UniProt accession", "Position in protein", "Type", "Sample"],
            ["P12345", 10, "CM(K)", "HUVEC"],
            ["Q8N123", 20, "CM(K)", "MEF"],
        ]
    )
    with pd.ExcelWriter(workbook) as writer:
        raw.to_excel(writer, sheet_name="results", index=False, header=False)
    parsed, diagnostics = normalize_workbook(workbook, source_label="cellular")
    assert len(parsed) == 2
    assert diagnostics[0].header_depth in {1, 2}
    assert diagnostics[0].parsed_sites == 2


def test_prepare_candidate_table_writes_diagnostics_on_failure(tmp_path: Path) -> None:
    from cml_repairnet.atlas.sources import prepare_candidate_table

    workbook = tmp_path / "invalid.xlsx"
    pd.DataFrame({"description": ["no sites here"]}).to_excel(workbook, index=False)
    diagnostics = tmp_path / "diagnostics.csv"
    try:
        prepare_candidate_table([workbook], tmp_path / "out.csv", diagnostics)
    except ValueError:
        pass
    else:
        raise AssertionError("Expected an empty-parser ValueError")
    assert diagnostics.exists()
    frame = pd.read_csv(diagnostics)
    assert "columns_preview" in frame.columns
