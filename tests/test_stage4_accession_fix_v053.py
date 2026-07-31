from __future__ import annotations

from pathlib import Path

import pandas as pd

import cml_repairnet.atlas.sources as sources
from cml_repairnet.atlas.sources import (
    annotate_candidates_with_uniprot,
    extract_accessions,
    is_valid_uniprot_accession,
    normalize_workbook,
)


def test_strict_uniprot_accession_pattern_rejects_protein_words() -> None:
    text = (
        "ACTIVATING ASSOCIATED ATTACHMENT CONTAINING CYTOCHROME "
        "P68371 A0A0B4J2F0 Q9Y6R7-2"
    )
    assert extract_accessions(text) == ["P68371", "A0A0B4J2F0", "Q9Y6R7-2"]
    assert not is_valid_uniprot_accession("ACTIVATING")
    assert not is_valid_uniprot_accession("ASSOCIATED")
    assert is_valid_uniprot_accession("P68371")
    assert is_valid_uniprot_accession("A0A0B4J2F0")


def test_parser_prefers_accession_rich_column_over_protein_names(tmp_path: Path) -> None:
    workbook = tmp_path / "supp.xlsx"
    frame = pd.DataFrame(
        {
            "Protein": [
                "Elongation factor containing protein",
                "Cytochrome activating protein",
            ],
            "Database identifier": ["sp|P68371|TBA1_HUMAN", "A0A0B4J2F0"],
            "Position": ["K40", "K75"],
            "Modification": ["CML", "CML"],
            "Sample": ["HUVEC", "MEF"],
        }
    )
    with pd.ExcelWriter(workbook) as writer:
        frame.to_excel(writer, sheet_name="CML sites", index=False)

    parsed, diagnostics = normalize_workbook(workbook, source_label="cellular")
    assert parsed["uniprot_accession_raw"].tolist() == ["P68371", "A0A0B4J2F0"]
    assert diagnostics[0].accession_column == "Database identifier"


def test_annotate_filters_invalid_accessions_before_api(
    tmp_path: Path, monkeypatch
) -> None:
    candidates = pd.DataFrame(
        {
            "candidate_id": ["human:ACTIVATING:K2", "human:P68371:K2"],
            "species": ["human", "human"],
            "organism_id": [9606, 9606],
            "uniprot_accession_raw": ["ACTIVATING", "P68371"],
            "site_uniprot": [2, 2],
        }
    )
    candidates_csv = tmp_path / "candidates.csv"
    candidates.to_csv(candidates_csv, index=False)

    seen: list[str] = []

    def fake_fetch(accessions, cache_dir, **kwargs):
        seen.extend(accessions)
        return pd.DataFrame(
            [
                {
                    "uniprot_accession": "P68371",
                    "uniprot_entry_name": "TEST_HUMAN",
                    "protein_name": "Test protein",
                    "gene_names": "TEST",
                    "organism_id_uniprot": 9606,
                    "protein_length": 3,
                    "sequence": "AKA",
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
    )

    assert seen == ["P68371"]
    assert eligible["uniprot_accession"].tolist() == ["P68371"]
    invalid = excluded[excluded["uniprot_accession_raw"] == "ACTIVATING"]
    assert invalid["sequence_validation_reason"].iloc[0] == "invalid_uniprot_accession_format"
