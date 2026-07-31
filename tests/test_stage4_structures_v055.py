from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from cml_repairnet.atlas.structures import (
    AlphaFoldEntry,
    extract_site_features_from_structure,
    fetch_alphafold_metadata,
    select_entry_for_site,
)


def _entry(
    entry_id: str,
    start: int,
    end: int,
    sequence: str | None,
    confidence: float = 90.0,
) -> AlphaFoldEntry:
    return AlphaFoldEntry(
        accession="P04406",
        entry_id=entry_id,
        uniprot_start=start,
        uniprot_end=end,
        cif_url="https://example.org/model.cif",
        pdb_url=None,
        model_version=6,
        global_metric_value=confidence,
        sequence=sequence,
        model_accession=None,
    )


def test_select_entry_prefers_canonical_over_shorter_isoform() -> None:
    reference = "M" + "A" * 193 + "K" + "A" * 140
    canonical = _entry("AF-P04406-F1", 1, len(reference), reference, 95.0)
    isoform_sequence = reference[:293]
    isoform = _entry("AF-P04406-2-F1", 1, len(isoform_sequence), isoform_sequence, 99.0)

    selected = select_entry_for_site(
        [isoform, canonical],
        195,
        expected_accession="P04406",
        expected_sequence=reference,
    )

    assert selected is not None
    assert selected.entry_id == "AF-P04406-F1"


def test_select_entry_rejects_isoform_residue_conflict() -> None:
    reference = "M" + "A" * 9 + "K" + "A" * 20
    canonical = _entry("AF-P04406-F1", 1, len(reference), reference)
    isoform_sequence = reference[:10] + "R" + reference[11:]
    isoform = _entry("AF-P04406-2-F1", 1, len(reference), isoform_sequence, 99.0)

    selected = select_entry_for_site(
        [isoform, canonical],
        11,
        expected_accession="P04406",
        expected_sequence=reference,
    )

    assert selected is not None
    assert selected.entry_id == "AF-P04406-F1"


def test_fetch_metadata_supports_post_sunset_api_fields(tmp_path: Path) -> None:
    payload = [
        {
            "modelEntityId": "AF-P04406-F1",
            "sequenceStart": 5,
            "sequenceEnd": 14,
            "sequence": "AAAKAAAAAA",
            "cifUrl": "https://example.org/P04406.cif",
            "latestVersion": 6,
            "globalMetricValue": 91.5,
        }
    ]
    cache = tmp_path / "metadata"
    cache.mkdir()
    (cache / "P04406.json").write_text(json.dumps(payload), encoding="utf-8")

    entries = fetch_alphafold_metadata("P04406", cache)

    assert len(entries) == 1
    assert entries[0].entry_id == "AF-P04406-F1"
    assert entries[0].uniprot_start == 5
    assert entries[0].uniprot_end == 14
    assert entries[0].sequence == "AAAKAAAAAA"
    assert entries[0].model_accession == "P04406"


def _write_fragment_pdb(path: Path) -> None:
    lines = []
    # Coordinate numbering is local 1..300. Canonical site 250 in a fragment
    # beginning at 151 maps to local residue 100. Residue number 250 is an
    # intentional decoy that must not be selected.
    for number in range(1, 301):
        name = "LYS" if number == 100 else "ALA"
        x = float(number)
        lines.append(
            f"ATOM  {number:5d}  CA  {name:>3s} A{number:4d}    "
            f"{x:8.3f}{0.0:8.3f}{0.0:8.3f}  1.00{90.0:6.2f}           C"
        )
    lines.extend(["TER", "END"])
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def test_fragment_mapping_uses_entry_sequence_index_before_direct_number(tmp_path: Path) -> None:
    pdb = tmp_path / "fragment.pdb"
    _write_fragment_pdb(pdb)

    fragment_sequence = "A" * 99 + "K" + "A" * 200
    reference_sequence = "A" * 249 + "K" + "A" * 200
    candidates = pd.DataFrame(
        {
            "candidate_id": ["human:P12345:K250"],
            "species": ["human"],
            "uniprot_accession": ["P12345"],
            "site_uniprot": [250],
            "sequence": [reference_sequence],
        }
    )
    entry = AlphaFoldEntry(
        accession="P12345",
        entry_id="AF-P12345-F2",
        uniprot_start=151,
        uniprot_end=450,
        cif_url="file://unused",
        pdb_url=None,
        model_version=6,
        global_metric_value=90.0,
        sequence=fragment_sequence,
        model_accession="P12345",
    )

    result = extract_site_features_from_structure(candidates, pdb, entry)

    assert result.loc[0, "structure_status"] == "eligible"
    assert result.loc[0, "structure_residue_number"] == 100
    assert result.loc[0, "structure_numbering_mapping"] == "entry_sequence_index"
