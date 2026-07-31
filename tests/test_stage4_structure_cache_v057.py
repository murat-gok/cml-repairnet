from __future__ import annotations

import pandas as pd

from cml_repairnet.atlas.structures import (
    _is_retryable_structure_reason,
    _partition_cached_rows,
)


def _candidates() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "candidate_id": "human:P00001:K10",
                "species": "human",
                "uniprot_accession": "P00001",
                "site_uniprot": 10,
            },
            {
                "candidate_id": "human:P00001:K20",
                "species": "human",
                "uniprot_accession": "P00001",
                "site_uniprot": 20,
            },
        ]
    )


def test_cache_missing_full_run_site_is_processed() -> None:
    candidates = _candidates()
    quick_cache = pd.DataFrame(
        [
            {
                "candidate_id": "human:P00001:K10",
                "species": "human",
                "uniprot_accession": "P00001",
                "site_uniprot": 10,
                "structure_status": "eligible",
                "structure_reason": None,
            }
        ]
    )

    stable, pending = _partition_cached_rows(candidates, quick_cache)

    assert stable["candidate_id"].tolist() == ["human:P00001:K10"]
    assert pending["candidate_id"].tolist() == ["human:P00001:K20"]


def test_transient_failure_is_retried_but_permanent_row_is_cached() -> None:
    candidates = _candidates()
    cached = pd.DataFrame(
        [
            {
                "candidate_id": "human:P00001:K10",
                "species": "human",
                "uniprot_accession": "P00001",
                "site_uniprot": 10,
                "structure_status": "excluded",
                "structure_reason": "structure_processing_error:ConnectionError:dns",
            },
            {
                "candidate_id": "human:P00001:K20",
                "species": "human",
                "uniprot_accession": "P00001",
                "site_uniprot": 20,
                "structure_status": "excluded",
                "structure_reason": "alphafold_entry_not_found",
            },
        ]
    )

    stable, pending = _partition_cached_rows(candidates, cached)

    assert stable["candidate_id"].tolist() == ["human:P00001:K20"]
    assert pending["candidate_id"].tolist() == ["human:P00001:K10"]
    assert _is_retryable_structure_reason(
        "structure_processing_error:ConnectionError:dns"
    )
    assert not _is_retryable_structure_reason("alphafold_entry_not_found")
