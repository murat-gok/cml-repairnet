from __future__ import annotations

from pathlib import Path

import pandas as pd

from cml_repairnet.atlas import sources
from cml_repairnet.atlas.sources import (
    annotate_candidates_with_uniprot,
    extract_sites,
    map_cml_peptide_to_sequence,
    normalize_workbook,
    prepare_candidate_table,
)


def test_site_parser_rejects_accession_digits() -> None:
    assert extract_sites("P01942") == []
    assert extract_sites("Q9JIF0") == []
    assert extract_sites("K58; Lys 92; 154") == [58, 92, 154]
    assert extract_sites(122) == [122]


def test_generic_cml_table_does_not_invent_peptide_or_prm(tmp_path: Path) -> None:
    workbook = tmp_path / "cellular.xlsx"
    raw = pd.DataFrame(
        [
            ["# discovery table", None, None, None],
            ["ID", "Genes", "Protein Descriptions", "Modified site"],
            ["Q9JIF0", "Prmt1", "Protein arginine N-methyltransferase 1", "K143; K145"],
        ]
    )
    with pd.ExcelWriter(workbook) as writer:
        raw.to_excel(writer, sheet_name="1_CM Sites_table_MEF", index=False, header=False)
    parsed, diagnostics = normalize_workbook(workbook, source_label="cellular")
    assert parsed["site_uniprot"].tolist() == [143, 145]
    assert not parsed["prm_validated"].any()
    assert parsed["peptide_sequence"].eq("").all()
    assert diagnostics[0].parsed_sites == 2


def test_age_sheet_uses_fdr_and_preserves_all_tested_sites(tmp_path: Path) -> None:
    workbook = tmp_path / "mouse.xlsx"
    raw = pd.DataFrame(
        [
            ["# Old vs Young heart", None, None, None, None],
            ["logFC", "P.Value", "adj.P.Val", "ID", "Modified site"],
            [1.2, 0.001, 0.01, "P12345", "K10, K20"],
            [-0.4, 0.2, 0.3, "Q8N123", "K30"],
            [0.5, 0.01, 0.02, "Q8N123", "ProtNterm"],
        ]
    )
    with pd.ExcelWriter(workbook) as writer:
        raw.to_excel(writer, sheet_name="1_Aging_Heart_O_vs_Y", index=False, header=False)
    parsed, diagnostics = normalize_workbook(
        workbook, source_label="mouse_organ", default_species="mouse"
    )
    assert len(parsed) == 3
    assert parsed["age_tested"].all()
    assert parsed.loc[parsed["uniprot_accession_raw"] == "P12345", "age_associated"].all()
    assert not parsed.loc[parsed["uniprot_accession_raw"] == "Q8N123", "age_associated"].any()
    assert set(parsed["tissue_or_condition"]) == {"Heart"}
    assert diagnostics[0].status == "parsed_age_sites"


def test_overlap_sheet_emits_human_and_mouse_rows(tmp_path: Path) -> None:
    workbook = tmp_path / "cellular.xlsx"
    frame = pd.DataFrame(
        {
            "HUVEC_ID": ["O14980"],
            "HUVEC_Residue": ["K"],
            "HUVEC_Modified site": ["K122"],
            "MEF_ID": ["Q6P5F9"],
            "MEF_Residue": ["K"],
            "MEF_Modified site": [122],
            "sequence_id": [98.0],
            "HUVEC_peptide": ["SDPTCVEKEKVYIGK"],
            "MEF_peptide": ["SDPTCVEKEKVYIGK"],
        }
    )
    with pd.ExcelWriter(workbook) as writer:
        frame.to_excel(writer, sheet_name="3_Overlap_after_alignment", index=False)
    parsed, diagnostics = normalize_workbook(workbook, source_label="cellular")
    assert len(parsed) == 2
    assert set(parsed["species"]) == {"human", "mouse"}
    assert parsed["ortholog_overlap"].all()
    assert set(parsed["site_uniprot"].astype(int)) == {122}
    assert diagnostics[0].status == "parsed_overlap_pairs"


def test_proteostasis_site_parser_uses_parenthetical_position(tmp_path: Path) -> None:
    workbook = tmp_path / "mouse.xlsx"
    frame = pd.DataFrame(
        {
            "ID": ["P68372", "Q8BGE6"],
            "Genes": ["Tubb4b", "Atg4b"],
            "Category": ["AUTOPHAGY", "AUTOPHAGY"],
            "Protein Descriptions": ["Tubulin", "ATG4B"],
            "Modified site": ["Heart, Kidney (Kidney324)", "Heart (ProtNterm)"],
        }
    )
    with pd.ExcelWriter(workbook) as writer:
        frame.to_excel(writer, sheet_name="6_CM_Sites_Proteostasis_Tissue", index=False)
    parsed, diagnostics = normalize_workbook(
        workbook, source_label="mouse", default_species="mouse"
    )
    assert len(parsed) == 1
    assert parsed.iloc[0]["site_uniprot"] == 324
    assert parsed.iloc[0]["tissue_or_condition"] == "Heart, Kidney"
    assert diagnostics[0].status == "parsed_proteostasis_sites"


def test_prm_peptide_is_mapped_during_annotation(tmp_path: Path, monkeypatch) -> None:
    discovery = tmp_path / "mouse_discovery.xlsx"
    prm_book = tmp_path / "prm.xlsx"
    pd.DataFrame({"ID": ["P01942"], "Modified site": ["K12"]}).to_excel(
        discovery, sheet_name="5_CM Sites_Tissues", index=False
    )
    prm_raw = pd.DataFrame(
        [
            ["# PRM panel", None, None, None],
            ["ID", "Genes", "Peptide sequence", "Charge state"],
            ["P01942", "Hba", "FLASVSTVLTSK[CML]YR", 2],
        ]
    )
    with pd.ExcelWriter(prm_book) as writer:
        prm_raw.to_excel(writer, sheet_name="7_PRM_list", index=False, header=False)

    candidates_csv = tmp_path / "stage4_candidates.csv"
    diagnostics_csv = tmp_path / "diagnostics.csv"
    prm_csv = tmp_path / "stage4_prm_evidence.csv"
    prepare_candidate_table(
        [discovery, prm_book],
        candidates_csv,
        diagnostics_csv,
        prm_evidence_csv=prm_csv,
    )
    prm = pd.read_csv(prm_csv)
    assert len(prm) == 1
    assert prm["site_uniprot"].isna().all()

    sequence = "MVLSGEDKSNIKAAWGKIGGHGAEYGAEALERMFASFPTTKTYFPHFDVSHGSAQVKGHGKKVADALASAAGHLDDLPGALSALSDLHAHKLRVDPVNFKLLSHCLLVTLASHHPADFTPAVHASLDKFLASVSTVLTSKYR"
    expected_site, status, _ = map_cml_peptide_to_sequence("FLASVSTVLTSK[CML]YR", sequence)
    assert status == "mapped"
    assert expected_site == 140

    def fake_fetch(accessions, cache_dir, **kwargs):
        return pd.DataFrame(
            [
                {
                    "uniprot_accession": "P01942",
                    "uniprot_entry_name": "HBA_MOUSE",
                    "protein_name": "Hemoglobin alpha",
                    "gene_names": "Hba",
                    "organism_id_uniprot": 10090,
                    "protein_length": len(sequence),
                    "sequence": sequence,
                    "reviewed_status": "reviewed",
                    "subcellular_location": "",
                }
            ]
        )

    monkeypatch.setattr(sources, "fetch_uniprot_records", fake_fetch)
    eligible, excluded = annotate_candidates_with_uniprot(
        candidates_csv,
        tmp_path / "eligible.csv",
        tmp_path / "excluded.csv",
        tmp_path / "cache",
        prm_evidence_csv=prm_csv,
        prm_mapping_csv=tmp_path / "prm_mapping.csv",
        prm_mapping_excluded_csv=tmp_path / "prm_mapping_excluded.csv",
    )
    mapped = eligible[eligible["site_uniprot"] == expected_site]
    assert len(mapped) == 1
    assert bool(mapped.iloc[0]["prm_validated"])
    assert "[CML]" in mapped.iloc[0]["peptide_examples"]
    assert not (excluded["site_uniprot"].astype(str) == "1942").any()
