from __future__ import annotations

import io
import json
import time
from collections.abc import Iterable
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, ClassVar

import pandas as pd
import requests
from requests import Response
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry


class OfflineCacheMiss(RuntimeError):
    pass


def _safe_key(text: str) -> str:
    return "".join(
        character if character.isalnum() or character in "-_." else "_" for character in text
    )


def build_retrying_session(*, total: int = 5, backoff_factor: float = 1.0) -> requests.Session:
    session = requests.Session()
    retry = Retry(
        total=total,
        connect=total,
        read=total,
        status=total,
        backoff_factor=backoff_factor,
        status_forcelist=(429, 500, 502, 503, 504),
        allowed_methods=frozenset({"GET", "POST"}),
        respect_retry_after_header=True,
    )
    adapter = HTTPAdapter(max_retries=retry)
    session.mount("https://", adapter)
    session.mount("http://", adapter)
    session.headers.update({"User-Agent": "CML-RepairNet/0.6.1 Stage5"})
    return session


@dataclass
class CachedHttpClient:
    cache_dir: Path
    offline: bool = False
    timeout: float = 90.0
    min_delay_seconds: float = 0.25

    def __post_init__(self) -> None:
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self.session = build_retrying_session()
        self._last_request = 0.0

    def _delay(self) -> None:
        elapsed = time.monotonic() - self._last_request
        remaining = self.min_delay_seconds - elapsed
        if remaining > 0:
            time.sleep(remaining)

    def _json_request(
        self,
        *,
        method: str,
        url: str,
        cache_name: str,
        params: dict[str, Any] | None = None,
        payload: dict[str, Any] | None = None,
    ) -> Any:
        cache_path = self.cache_dir / f"{_safe_key(cache_name)}.json"
        if cache_path.exists():
            return json.loads(cache_path.read_text(encoding="utf-8"))
        if self.offline:
            raise OfflineCacheMiss(f"Offline cache missing: {cache_path}")
        self._delay()
        response = self.session.request(
            method,
            url,
            params=params,
            json=payload,
            timeout=self.timeout,
        )
        self._last_request = time.monotonic()
        response.raise_for_status()
        data = response.json()
        cache_path.write_text(json.dumps(data, indent=2), encoding="utf-8")
        self._write_metadata(cache_path, response)
        return data

    def _tsv_request(
        self,
        *,
        method: str,
        url: str,
        cache_name: str,
        params: dict[str, Any] | None = None,
        data: dict[str, Any] | None = None,
    ) -> pd.DataFrame:
        cache_path = self.cache_dir / f"{_safe_key(cache_name)}.tsv"
        if cache_path.exists():
            text = cache_path.read_text(encoding="utf-8")
            return pd.read_csv(io.StringIO(text), sep="\t") if text.strip() else pd.DataFrame()
        if self.offline:
            raise OfflineCacheMiss(f"Offline cache missing: {cache_path}")
        self._delay()
        response = self.session.request(
            method,
            url,
            params=params,
            data=data,
            timeout=self.timeout,
        )
        self._last_request = time.monotonic()
        response.raise_for_status()
        cache_path.write_text(response.text, encoding="utf-8")
        self._write_metadata(cache_path, response)
        if not response.text.strip():
            return pd.DataFrame()
        return pd.read_csv(io.StringIO(response.text), sep="\t")

    def _write_metadata(self, cache_path: Path, response: Response) -> None:
        selected_headers = {
            key: value
            for key, value in response.headers.items()
            if key.lower()
            in {
                "x-uniprot-release",
                "x-uniprot-release-date",
                "etag",
                "last-modified",
                "content-type",
            }
        }
        metadata = {
            "retrieved_at_utc": datetime.now(timezone.utc).isoformat(),
            "url": response.url,
            "status_code": response.status_code,
            "headers": selected_headers,
        }
        cache_path.with_suffix(cache_path.suffix + ".meta.json").write_text(
            json.dumps(metadata, indent=2), encoding="utf-8"
        )

    def cached_metadata(self) -> list[dict[str, Any]]:
        rows: list[dict[str, Any]] = []
        for path in sorted(self.cache_dir.glob("*.meta.json")):
            try:
                row = json.loads(path.read_text(encoding="utf-8"))
                row["cache_file"] = path.name
                rows.append(row)
            except Exception:  # noqa: BLE001, S112
                continue
        return rows


class GProfilerClient(CachedHttpClient):
    base_url = "https://biit.cs.ut.ee/gprofiler/api"

    def profile(
        self,
        *,
        organism: str,
        query: list[str],
        background: list[str],
        sources: list[str],
        query_name: str,
        threshold: float = 0.05,
        correction: str = "fdr",
        no_iea: bool = False,
    ) -> dict[str, Any]:
        payload = {
            "organism": organism,
            "query": query,
            "sources": sources,
            "user_threshold": threshold,
            "all_results": True,
            "ordered": False,
            "combined": False,
            "measure_underrepresentation": False,
            "no_iea": no_iea,
            "domain_scope": "custom",
            "background": background,
            "significance_threshold_method": correction,
            "no_evidences": False,
        }
        cache_name = f"gost_{organism}_{query_name}_{len(query)}_{len(background)}"
        return self._json_request(
            method="POST",
            url=f"{self.base_url}/gost/profile/",
            cache_name=cache_name,
            payload=payload,
        )

    def orthologs(
        self,
        *,
        organism: str,
        target: str,
        query: list[str],
        query_name: str,
    ) -> dict[str, Any]:
        payload = {"organism": organism, "target": target, "query": query}
        return self._json_request(
            method="POST",
            url=f"{self.base_url}/orth/orth/",
            cache_name=f"gorth_{organism}_{target}_{query_name}_{len(query)}",
            payload=payload,
        )

    def data_versions(self, organism: str) -> dict[str, Any]:
        return self._json_request(
            method="GET",
            url=f"{self.base_url}/util/data_versions/",
            cache_name=f"gprofiler_versions_{organism}",
            params={"organism": organism},
        )


class StringClient(CachedHttpClient):
    base_url = "https://version-12-0.string-db.org/api"

    def map_identifiers(
        self,
        identifiers: list[str],
        *,
        species: int,
        query_name: str,
    ) -> pd.DataFrame:
        return self._tsv_request(
            method="POST",
            url=f"{self.base_url}/tsv/get_string_ids",
            cache_name=f"string_map_{species}_{query_name}_{len(identifiers)}",
            data={
                "identifiers": "\r".join(identifiers),
                "species": species,
                "echo_query": 1,
                "caller_identity": "CML-RepairNet",
            },
        )

    def network(
        self,
        identifiers: list[str],
        *,
        species: int,
        query_name: str,
        required_score: int = 700,
        network_type: str = "functional",
    ) -> pd.DataFrame:
        return self._tsv_request(
            method="POST",
            url=f"{self.base_url}/tsv/network",
            cache_name=f"string_network_{species}_{query_name}_{required_score}_{len(identifiers)}",
            data={
                "identifiers": "\r".join(identifiers),
                "species": species,
                "required_score": required_score,
                "network_type": network_type,
                "add_nodes": 0,
                "show_query_node_labels": 1,
                "caller_identity": "CML-RepairNet",
            },
        )

    def ppi_enrichment(
        self,
        identifiers: list[str],
        *,
        species: int,
        query_name: str,
    ) -> pd.DataFrame:
        return self._tsv_request(
            method="POST",
            url=f"{self.base_url}/tsv/ppi_enrichment",
            cache_name=f"string_ppi_{species}_{query_name}_{len(identifiers)}",
            data={
                "identifiers": "\r".join(identifiers),
                "species": species,
                "caller_identity": "CML-RepairNet",
            },
        )


class UniProtAnnotationClient(CachedHttpClient):
    base_url = "https://rest.uniprot.org/uniprotkb/search"
    fields: ClassVar[list[str]] = [
        "accession",
        "id",
        "gene_names",
        "protein_name",
        "protein_families",
        "keyword",
        "go_p",
        "go_c",
        "go_f",
        "xref_interpro",
        "xref_pfam",
        "cc_function",
        "cc_pathway",
        "cc_subcellular_location",
    ]

    def fetch(
        self, accessions: Iterable[str], *, query_name: str, chunk_size: int = 80
    ) -> pd.DataFrame:
        unique = sorted(
            {str(value).strip().split("-")[0] for value in accessions if str(value).strip()}
        )
        frames: list[pd.DataFrame] = []
        for index in range(0, len(unique), chunk_size):
            chunk = unique[index : index + chunk_size]
            query = " OR ".join(f"accession:{accession}" for accession in chunk)
            frame = self._tsv_request(
                method="GET",
                url=self.base_url,
                cache_name=f"uniprot_annotations_{query_name}_{index // chunk_size}_{len(chunk)}",
                params={
                    "query": f"({query})",
                    "format": "tsv",
                    "fields": ",".join(self.fields),
                    "size": 500,
                },
            )
            frames.append(frame)
        if not frames:
            return pd.DataFrame()
        return pd.concat(frames, ignore_index=True).drop_duplicates()
