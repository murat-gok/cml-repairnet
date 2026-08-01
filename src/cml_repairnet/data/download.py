from __future__ import annotations

import hashlib
import json
import re
from collections.abc import Iterable
from pathlib import Path
from urllib.parse import urljoin

import requests
from bs4 import BeautifulSoup
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

ARTICLE_SOURCE_URL = (
    "https://static-content.springer.com/esm/art%3A10.1038%2F"
    "s41467-026-75141-2/MediaObjects/41467_2026_75141_MOESM4_ESM.xlsx"
)
BSA_FASTA_URL = "https://rest.uniprot.org/uniprotkb/P02769.fasta"
BSA_PDB_URL = "https://files.rcsb.org/download/4F5S.pdb"
PRIDE_BASE_URL = "https://ftp.pride.ebi.ac.uk/pride/data/archive/2026/07/PXD070124/"
USER_AGENT = "CML-RepairNet/0.2 (+academic reproducibility; public-data downloader)"


def _session() -> requests.Session:
    retry = Retry(
        total=4,
        connect=4,
        read=4,
        backoff_factor=1.0,
        status_forcelist=(429, 500, 502, 503, 504),
        allowed_methods=frozenset({"GET"}),
    )
    session = requests.Session()
    session.headers.update({"User-Agent": USER_AGENT})
    session.mount("https://", HTTPAdapter(max_retries=retry))
    return session


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def download_file(url: str, destination: Path, timeout: int = 120) -> dict[str, object]:
    """Download atomically and return a provenance record."""
    destination.parent.mkdir(parents=True, exist_ok=True)
    partial = destination.with_suffix(destination.suffix + ".part")
    if partial.exists():
        partial.unlink()

    try:
        with _session().get(url, stream=True, timeout=timeout) as response:
            response.raise_for_status()
            with partial.open("wb") as handle:
                for chunk in response.iter_content(chunk_size=1024 * 1024):
                    if chunk:
                        handle.write(chunk)
        if not partial.exists() or partial.stat().st_size == 0:
            raise RuntimeError(f"Downloaded file is empty: {url}")
        partial.replace(destination)
    except Exception:
        partial.unlink(missing_ok=True)
        raise

    return {
        "url": url,
        "path": str(destination),
        "bytes": destination.stat().st_size,
        "sha256": sha256_file(destination),
    }


def list_pride_files(base_url: str = PRIDE_BASE_URL) -> list[str]:
    """List files from the public PRIDE HTTPS directory."""
    response = _session().get(base_url, timeout=60)
    response.raise_for_status()
    soup = BeautifulSoup(response.text, "html.parser")
    names: list[str] = []
    for anchor in soup.find_all("a", href=True):
        href = str(anchor["href"])
        name = href.rstrip("/").split("/")[-1]
        if name and not name.startswith("?") and name not in {"..", "."}:
            names.append(name)
    return sorted(set(names))


def select_small_processed_files(names: Iterable[str]) -> list[str]:
    patterns = [
        re.compile(r"processed_filt_quant\.csv$", re.IGNORECASE),
        re.compile(r"peptide_IDs_repository\.csv$", re.IGNORECASE),
        re.compile(r"Mapping_raw_files\.csv$", re.IGNORECASE),
        re.compile(r"\.sdrf\.tsv$", re.IGNORECASE),
    ]
    return sorted(name for name in names if any(pattern.search(name) for pattern in patterns))


def download_public_bundle(
    raw_dir: Path,
    include_pride_processed: bool = False,
    manifest_name: str = "download_manifest.json",
) -> list[dict[str, object]]:
    records = [
        download_file(ARTICLE_SOURCE_URL, raw_dir / "article_source_data.xlsx"),
        download_file(BSA_FASTA_URL, raw_dir / "P02769.fasta"),
        download_file(BSA_PDB_URL, raw_dir / "4F5S.pdb"),
    ]
    if include_pride_processed:
        try:
            pride_files = select_small_processed_files(list_pride_files())
        except requests.RequestException as exc:
            raise RuntimeError(
                "Could not list the PRIDE directory. Retry later or download processed files "
                "manually from PXD070124."
            ) from exc
        for filename in pride_files:
            records.append(download_file(urljoin(PRIDE_BASE_URL, filename), raw_dir / filename))

    manifest_path = raw_dir / manifest_name
    manifest_path.write_text(json.dumps(records, indent=2), encoding="utf-8")
    return records
