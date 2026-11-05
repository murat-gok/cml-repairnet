from __future__ import annotations

import hashlib
import json
import re
import time
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

NATURE_SUPPLEMENT_URLS: Mapping[str, str] = {
    "cellular_cml_sites": (
        "https://media.springernature.com/original/springer-static/esm/"
        "art%3A10.1038%2Fs41467-021-26982-6/MediaObjects/"
        "41467_2021_26982_MOESM4_ESM.xlsx"
    ),
    "mouse_organ_cml_sites": (
        "https://media.springernature.com/original/springer-static/esm/"
        "art%3A10.1038%2Fs41467-021-26982-6/MediaObjects/"
        "41467_2021_26982_MOESM5_ESM.xlsx"
    ),
}

UNIPROT_STREAM_URL = "https://rest.uniprot.org/uniprotkb/stream"
UNIPROT_FIELDS = (
    "accession,id,protein_name,gene_names,organism_id,length,sequence,reviewed,"
    "cc_subcellular_location"
)

# Official UniProtKB accession syntax. The previous generic ``[A-Z0-9]{10}``
# branch incorrectly treated ten-letter protein-name words such as ``ACTIVATING``
# and ``ASSOCIATED`` as accessions. UniProt's 10-character identifiers retain
# the same positional constraints as the 6-character format.
_ACCESSION_CORE_PATTERN = (
    r"(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|"
    r"[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})"
)
_ACCESSION_RE = re.compile(
    rf"(?<![A-Z0-9]){_ACCESSION_CORE_PATTERN}(?:-\d+)?(?![A-Z0-9])",
    re.IGNORECASE,
)
_ACCESSION_FULL_RE = re.compile(
    rf"^{_ACCESSION_CORE_PATTERN}(?:-\d+)?$",
    re.IGNORECASE,
)
_SITE_PREFIX_RE = re.compile(r"(?:K|LYS)\s*[-_:]?\s*(\d{1,5})", re.IGNORECASE)
_PURE_SITE_LIST_RE = re.compile(r"^[\s;,/|+KLYS0-9:_-]+$", re.IGNORECASE)
_PROTEOSTASIS_SITE_RE = re.compile(r"\((?:[A-Za-z _-]*?)(\d{1,5})\)")
_CML_PEPTIDE_MARKER_RE = re.compile(r"K\s*\[\s*CML\s*\]", re.IGNORECASE)

ACCESSION_ALIASES = (
    "id",
    "uniprot",
    "uniprot accessions",
    "uniprot accession",
    "uniprot id",
    "accession",
    "accessions",
    "protein id",
    "protein ids",
    "protein accessions",
    "majority protein ids",
    "leading razor protein",
    "protein",
    "proteins",
    "protein ids",
    "protein group ids",
)
SITE_ALIASES = (
    "site",
    "sites",
    "position",
    "positions",
    "positions within proteins",
    "position within protein",
    "positions in proteins",
    "position in protein",
    "modified site",
    "modified sites",
    "modified position",
    "modification position",
    "lysine position",
    "lys position",
    "amino acid position",
    "residue position",
    "cml site",
    "cm site",
)
PEPTIDE_ALIASES = (
    "sequence",
    "modified sequence",
    "peptide",
    "peptide sequence",
)
MODIFICATION_ALIASES = (
    "modification",
    "modifications",
    "site type",
    "ptm",
    "modification type",
)
TISSUE_ALIASES = (
    "tissue",
    "organ",
    "cell type",
    "sample",
    "condition",
    "experiment",
)
VALIDATION_ALIASES = (
    "prm",
    "validated",
    "validation",
    "heavy peptide",
    "targeted proteomics",
)


def requests_session() -> requests.Session:
    retry = Retry(
        total=5,
        connect=5,
        read=5,
        backoff_factor=0.8,
        status_forcelist=(429, 500, 502, 503, 504),
        allowed_methods=frozenset({"GET"}),
        respect_retry_after_header=True,
    )
    session = requests.Session()
    session.headers.update({"User-Agent": "CML-RepairNet/0.5.4 (research; public-data client)"})
    session.mount("https://", HTTPAdapter(max_retries=retry))
    return session


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def download_file(
    url: str,
    destination: Path,
    *,
    session: requests.Session | None = None,
    timeout: int = 120,
    overwrite: bool = False,
) -> dict[str, object]:
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists() and destination.stat().st_size > 0 and not overwrite:
        return {
            "url": url,
            "path": str(destination),
            "bytes": destination.stat().st_size,
            "sha256": sha256_file(destination),
            "cached": True,
        }
    session = session or requests_session()
    partial = destination.with_suffix(destination.suffix + ".part")
    response = session.get(url, timeout=timeout, stream=True)
    response.raise_for_status()
    with partial.open("wb") as handle:
        for chunk in response.iter_content(chunk_size=1024 * 1024):
            if chunk:
                handle.write(chunk)
    partial.replace(destination)
    return {
        "url": url,
        "path": str(destination),
        "bytes": destination.stat().st_size,
        "sha256": sha256_file(destination),
        "cached": False,
        "content_type": response.headers.get("content-type"),
        "etag": response.headers.get("etag"),
        "last_modified": response.headers.get("last-modified"),
    }


def download_stage4_source_data(
    raw_dir: Path,
    *,
    overwrite: bool = False,
) -> dict[str, object]:
    raw_dir.mkdir(parents=True, exist_ok=True)
    session = requests_session()
    entries: list[dict[str, object]] = []
    for name, url in NATURE_SUPPLEMENT_URLS.items():
        entries.append(
            {
                "name": name,
                **download_file(
                    url,
                    raw_dir / f"{name}.xlsx",
                    session=session,
                    overwrite=overwrite,
                ),
            }
        )
    manifest = {
        "dataset": "Di Sanzo et al. 2021 CML supplementary data",
        "doi": "10.1038/s41467-021-26982-6",
        "license": "CC BY 4.0 article supplementary material",
        "files": entries,
    }
    (raw_dir / "stage4_source_manifest.json").write_text(
        json.dumps(manifest, indent=2), encoding="utf-8"
    )
    return manifest


def _canonical_column(value: object) -> str:
    text = str(value).strip().lower()
    text = re.sub(r"[\n\r\t]+", " ", text)
    text = re.sub(r"[^a-z0-9]+", " ", text)
    return re.sub(r"\s+", " ", text).strip()


def _best_column(columns: Sequence[object], aliases: Sequence[str]) -> str | None:
    """Choose the most specific alias match without reverse-substring leakage.

    Earlier versions allowed a short column name such as ``ID`` or ``t`` to match
    longer aliases such as ``modified sequence`` or ``modification``. That caused
    accession columns to be interpreted as peptide/site-type columns. Matching is
    now limited to exact aliases, complete alias-token containment, or a sufficiently
    long alias occurring inside a composite header.
    """
    alias_set = {_canonical_column(alias) for alias in aliases if _canonical_column(alias)}
    best: tuple[float, str] | None = None
    for column in columns:
        original = str(column)
        canon = _canonical_column(column)
        canon_tokens = set(canon.split())
        for alias in alias_set:
            alias_tokens = set(alias.split())
            if canon == alias:
                score = 1000.0 + len(alias_tokens)
            elif alias_tokens and alias_tokens.issubset(canon_tokens):
                score = 100.0 + 10.0 * len(alias_tokens) + len(alias) / 1000.0
            elif len(alias) >= 4 and alias in canon:
                score = 50.0 + len(alias) / 1000.0
            else:
                continue
            candidate = (score, original)
            if best is None or candidate[0] > best[0]:
                best = candidate
    return best[1] if best else None


def is_valid_uniprot_accession(value: object) -> bool:
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return False
    return _ACCESSION_FULL_RE.fullmatch(str(value).strip().upper()) is not None


def normalize_uniprot_accession(value: object, *, keep_isoform: bool = True) -> str | None:
    if not is_valid_uniprot_accession(value):
        return None
    accession = str(value).strip().upper()
    return accession if keep_isoform else accession.split("-")[0]


def extract_accessions(value: object) -> list[str]:
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return []
    matches = [match.upper() for match in _ACCESSION_RE.findall(str(value).upper())]
    return list(dict.fromkeys(matches))


def extract_sites(value: object) -> list[int]:
    """Extract lysine positions while rejecting digits embedded in accessions.

    Numeric cells are accepted directly. Text cells require an explicit K/Lys
    prefix, except that bare numeric tokens are allowed inside an otherwise pure
    site list that already contains at least one explicit lysine position. Thus
    ``K58; Lys 92; 154`` is valid, while ``P01942`` is not interpreted as site 1942.
    """
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return []
    if isinstance(value, (int, np.integer)):
        return [int(value)] if int(value) > 0 else []
    if isinstance(value, float) and float(value).is_integer():
        return [int(value)] if int(value) > 0 else []
    text = str(value).strip()
    if not text:
        return []
    if re.fullmatch(r"\d{1,5}", text):
        site = int(text)
        return [site] if site > 0 else []
    explicit = [int(match) for match in _SITE_PREFIX_RE.findall(text)]
    sites = list(explicit)
    if explicit and _PURE_SITE_LIST_RE.fullmatch(text.upper()):
        for token in re.split(r"[;,/|+]", text):
            token = token.strip()
            if re.fullmatch(r"\d{1,5}", token):
                sites.append(int(token))
    return list(dict.fromkeys(site for site in sites if site > 0))


def parse_cml_peptide(value: object) -> tuple[str | None, int | None]:
    """Return the unmodified peptide and zero-based CML lysine offset."""
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return None, None
    raw = str(value).strip().upper()
    marker = _CML_PEPTIDE_MARKER_RE.search(raw)
    if marker is None:
        return None, None
    prefix_letters = re.sub(r"[^A-Z]", "", raw[: marker.start()])
    offset = len(prefix_letters)
    clean = _CML_PEPTIDE_MARKER_RE.sub("K", raw, count=1)
    clean = re.sub(r"\[[^\]]*\]", "", clean)
    clean = re.sub(r"[^A-Z]", "", clean)
    if not clean or offset >= len(clean) or clean[offset] != "K":
        return None, None
    return clean, offset


def map_cml_peptide_to_sequence(
    peptide: object, sequence: object
) -> tuple[int | None, str, str | None]:
    """Map a marked PRM peptide to a canonical sequence.

    Exact matching is preferred. If exact matching fails, an I/L-equivalent match
    is attempted because bottom-up proteomics cannot distinguish those residues.
    Multiple possible absolute lysine sites are reported as ambiguous.
    """
    clean, offset = parse_cml_peptide(peptide)
    if clean is None or offset is None:
        return None, "excluded", "missing_cml_marker"
    if not isinstance(sequence, str) or not sequence:
        return None, "excluded", "missing_uniprot_sequence"
    seq = re.sub(r"[^A-Z]", "", sequence.upper())

    def starts(haystack: str, needle: str) -> list[int]:
        found: list[int] = []
        position = haystack.find(needle)
        while position >= 0:
            found.append(position)
            position = haystack.find(needle, position + 1)
        return found

    positions = starts(seq, clean)
    method = "exact"
    if not positions:
        positions = starts(
            seq.translate(str.maketrans({"I": "J", "L": "J"})),
            clean.translate(str.maketrans({"I": "J", "L": "J"})),
        )
        method = "il_equivalent"
    sites = sorted({start + offset + 1 for start in positions})
    if not sites:
        return None, "excluded", "peptide_not_found_in_canonical_sequence"
    if len(sites) > 1:
        return None, "excluded", "ambiguous_multiple_peptide_matches"
    site = sites[0]
    if site < 1 or site > len(seq) or seq[site - 1] != "K":
        return None, "excluded", "mapped_residue_is_not_lysine"
    return site, "mapped", method


def _is_real_peptide(value: object) -> bool:
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return False
    text = str(value).strip().upper()
    if not text:
        return False
    clean = re.sub(r"\[[^\]]*\]", "", text)
    clean = re.sub(r"[^A-Z]", "", clean)
    return len(clean) >= 5 and set(clean).issubset(set("ACDEFGHIKLMNPQRSTVWY"))


def infer_species(text: str, default: str | None = None) -> tuple[str | None, int | None]:
    lowered = text.lower()
    human_tokens = ("huvec", "human", "homo sapiens", "9606")
    mouse_tokens = (
        "mef",
        "mouse",
        "murine",
        "mus musculus",
        "10090",
        "heart",
        "kidney",
        "liver",
    )
    human = any(token in lowered for token in human_tokens)
    mouse = any(token in lowered for token in mouse_tokens)
    if human and not mouse:
        return "human", 9606
    if mouse and not human:
        return "mouse", 10090
    if default == "human":
        return "human", 9606
    if default == "mouse":
        return "mouse", 10090
    return None, None


def _bool_from_text(value: object, tokens: Sequence[str]) -> bool:
    text = str(value).lower()
    return any(token in text for token in tokens)


def _pair_accessions_sites(accessions: list[str], sites: list[int]) -> list[tuple[str, int, str]]:
    if not accessions or not sites:
        return []
    if len(accessions) == len(sites):
        return [(accession, site, "paired") for accession, site in zip(accessions, sites)]
    if len(sites) == 1:
        return [(accession, sites[0], "shared_site") for accession in accessions]
    if len(accessions) == 1:
        return [(accessions[0], site, "shared_accession") for site in sites]
    # Ambiguous protein groups are intentionally narrowed to the leading accession.
    return [(accessions[0], site, "ambiguous_leading_accession") for site in sites]


@dataclass(frozen=True)
class ParseDiagnostics:
    workbook: str
    sheet: str
    rows: int
    parsed_sites: int
    accession_column: str | None
    site_column: str | None
    modification_column: str | None
    status: str
    header_row: int | None = None
    header_depth: int | None = None
    column_count: int | None = None
    columns_preview: str | None = None


def _unique_column_names(values: Sequence[object]) -> list[str]:
    """Return stable, non-empty and unique Excel column names."""
    result: list[str] = []
    seen: dict[str, int] = {}
    for index, value in enumerate(values):
        if value is None or (isinstance(value, float) and np.isnan(value)):
            base = f"unnamed_{index}"
        else:
            base = str(value).strip() or f"unnamed_{index}"
        count = seen.get(base, 0)
        seen[base] = count + 1
        result.append(base if count == 0 else f"{base}__{count + 1}")
    return result


def _combine_header_rows(first: Sequence[object], second: Sequence[object]) -> list[str]:
    combined: list[str] = []
    for index, (upper, lower) in enumerate(zip(first, second)):
        upper_text = "" if pd.isna(upper) else str(upper).strip()
        lower_text = "" if pd.isna(lower) else str(lower).strip()
        if lower_text and not lower_text.lower().startswith("unnamed"):
            if upper_text and _canonical_column(upper_text) != _canonical_column(lower_text):
                combined.append(f"{upper_text} {lower_text}".strip())
            else:
                combined.append(lower_text)
        elif upper_text:
            combined.append(upper_text)
        else:
            combined.append(f"unnamed_{index}")
    return _unique_column_names(combined)


def _header_alias_score(columns: Sequence[object]) -> tuple[int, bool, bool]:
    accession = _best_column(columns, ACCESSION_ALIASES) is not None
    site = _best_column(columns, SITE_ALIASES) is not None
    optional = sum(
        _best_column(columns, aliases) is not None
        for aliases in (PEPTIDE_ALIASES, MODIFICATION_ALIASES, TISSUE_ALIASES)
    )
    # Required columns dominate; optional columns break ties.
    score = 5 * int(accession) + 5 * int(site) + optional
    return score, accession, site


def _detect_header(raw: pd.DataFrame, *, max_scan_rows: int = 40) -> tuple[int, int, list[str]]:
    """Detect one- or two-row tabular headers beneath supplementary-file preambles."""
    if raw.empty:
        return 0, 1, []
    limit = min(max_scan_rows, len(raw))

    # Prefer a self-contained one-row header. This avoids combining a narrative
    # preamble row with the real header immediately below it.
    one_row_candidates: list[tuple[int, int, int, list[str]]] = []
    for row_index in range(limit):
        values = raw.iloc[row_index].tolist()
        columns = _unique_column_names(values)
        score, accession, site = _header_alias_score(columns)
        nonempty = sum(bool(pd.notna(value) and str(value).strip()) for value in values)
        if accession and site:
            one_row_candidates.append((score, int(nonempty), row_index, columns))
    if one_row_candidates:
        score, nonempty, row_index, columns = max(
            one_row_candidates, key=lambda item: (item[0], item[1], -item[2])
        )
        return row_index, 1, columns

    two_row_candidates: list[tuple[int, int, int, list[str]]] = []
    for row_index in range(limit - 1):
        upper = raw.iloc[row_index].tolist()
        lower = raw.iloc[row_index + 1].tolist()
        columns = _combine_header_rows(upper, lower)
        score, accession, site = _header_alias_score(columns)
        nonempty = sum(bool(pd.notna(value) and str(value).strip()) for value in lower)
        if accession and site:
            two_row_candidates.append((score, int(nonempty), row_index, columns))
    if two_row_candidates:
        score, nonempty, row_index, columns = max(
            two_row_candidates, key=lambda item: (item[0], item[1], -item[2])
        )
        return row_index, 2, columns

    # No alias-complete header was found. Keep the row with the strongest partial
    # evidence; content-based inference will inspect the columns below it.
    fallback: list[tuple[int, int, list[str]]] = []
    for row_index in range(limit):
        columns = _unique_column_names(raw.iloc[row_index].tolist())
        score, _, _ = _header_alias_score(columns)
        fallback.append((score, row_index, columns))
    score, row_index, columns = max(fallback, key=lambda item: (item[0], -item[1]))
    return row_index, 1, columns


def _content_column_score(series: pd.Series, extractor) -> tuple[int, float]:
    non_null = series.dropna()
    if non_null.empty:
        return 0, 0.0
    sample = non_null.head(300)
    parsed = sum(bool(extractor(value)) for value in sample)
    return parsed, parsed / len(sample)


def _infer_content_column(
    frame: pd.DataFrame,
    extractor,
    *,
    minimum_matches: int = 2,
    minimum_fraction: float = 0.10,
    excluded: Sequence[str] = (),
) -> str | None:
    best_name: str | None = None
    best_score = (0, 0.0)
    excluded_set = set(excluded)
    for column in frame.columns:
        name = str(column)
        if name in excluded_set:
            continue
        score = _content_column_score(frame[column], extractor)
        if score > best_score:
            best_name, best_score = name, score
    if best_score[0] >= minimum_matches and best_score[1] >= minimum_fraction:
        return best_name
    return None


def _accession_header_bonus(column: object) -> int:
    canon = _canonical_column(column)
    if not canon:
        return 0
    if "uniprot" in canon or "accession" in canon:
        return 50
    if any(
        token in canon
        for token in ("majority protein id", "protein ids", "protein id", "leading razor protein")
    ):
        return 30
    if canon in {"protein", "proteins", "protein name", "protein names"}:
        return 0
    return 5 if "protein" in canon else 0


def _select_accession_column(frame: pd.DataFrame) -> str | None:
    """Select an accession column using strict content evidence plus header hints.

    Header-only fuzzy matching can confuse a protein-name column with a UniProt
    identifier column. Valid accession count is therefore the primary criterion;
    explicit identifier-like headers only break ties.
    """
    candidates: list[tuple[int, float, int, str]] = []
    for column in frame.columns:
        name = str(column)
        matches, fraction = _content_column_score(frame[column], extract_accessions)
        if matches <= 0:
            continue
        candidates.append((matches, fraction, _accession_header_bonus(column), name))
    if not candidates:
        return None
    matches, fraction, _bonus, name = max(
        candidates,
        key=lambda item: (item[0], item[1], item[2], item[3]),
    )
    # A single valid identifier is acceptable for a small supplementary sheet,
    # while larger sheets naturally win through the match-count criterion.
    return name


def _read_supplementary_sheet(workbook: Path, sheet: str) -> tuple[pd.DataFrame, int, int]:
    raw = pd.read_excel(workbook, sheet_name=sheet, header=None, dtype=object)
    raw = raw.dropna(how="all").dropna(axis=1, how="all")
    if raw.empty:
        return pd.DataFrame(), 0, 1
    raw = raw.reset_index(drop=True)
    header_row, header_depth, columns = _detect_header(raw)
    data_start = header_row + header_depth
    frame = raw.iloc[data_start:].copy()
    frame.columns = columns
    frame = frame.dropna(how="all").reset_index(drop=True)
    return frame, header_row, header_depth


def _sheet_tissue(sheet: str) -> str:
    lowered = sheet.lower()
    for tissue in ("heart", "kidney", "liver"):
        if tissue in lowered:
            return tissue.title()
    if "huvec" in lowered:
        return "HUVEC"
    if "mef" in lowered:
        return "MEF"
    return sheet


def _proteostasis_site_and_tissues(value: object) -> tuple[list[int], str]:
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return [], ""
    text = str(value).strip()
    if "protnterm" in text.lower() or "prot nterm" in text.lower():
        return [], re.sub(r"\([^)]*\)", "", text).strip(" ,")
    sites = [int(match) for match in _PROTEOSTASIS_SITE_RE.findall(text)]
    tissues = re.sub(r"\([^)]*\)", "", text).strip(" ,")
    return list(dict.fromkeys(site for site in sites if site > 0)), tissues


def _numeric(record: pd.Series, column: str | None) -> float | None:
    if column is None:
        return None
    value = pd.to_numeric(record.get(column), errors="coerce")
    return float(value) if pd.notna(value) else None


def _base_evidence_row(**overrides: object) -> dict[str, object]:
    row: dict[str, object] = {
        "species": None,
        "organism_id": None,
        "uniprot_accession_raw": None,
        "site_uniprot": pd.NA,
        "source_label": "",
        "source_workbook": "",
        "source_sheet": "",
        "source_row": pd.NA,
        "pairing_status": "",
        "peptide_sequence": "",
        "tissue_or_condition": "",
        "prm_validated": False,
        "age_tested": False,
        "age_associated": False,
        "age_log2fc": np.nan,
        "age_p_value": np.nan,
        "age_adj_p_value": np.nan,
        "age_direction": "",
        "age_significant_tissue": "",
        "ortholog_overlap": False,
        "sequence_identity": np.nan,
        "evidence_type": "discovery_site",
    }
    row.update(overrides)
    return row


def normalize_workbook(
    workbook: Path,
    *,
    source_label: str,
    default_species: str | None = None,
) -> tuple[pd.DataFrame, list[ParseDiagnostics]]:
    """Normalize the Di Sanzo supplementary workbooks into evidence rows.

    Known sheet families are parsed explicitly because their semantics differ:
    discovery-site tables, human/mouse overlap tables, old-vs-young limma tables,
    proteostasis summaries, and PRM peptide panels. Generic parsing remains as a
    fallback for user-supplied tables.
    """
    excel = pd.ExcelFile(workbook)
    rows: list[dict[str, object]] = []
    diagnostics: list[ParseDiagnostics] = []

    for sheet in excel.sheet_names:
        try:
            frame, header_row, header_depth = _read_supplementary_sheet(workbook, sheet)
        except Exception as exc:  # noqa: BLE001
            diagnostics.append(
                ParseDiagnostics(
                    workbook.name,
                    sheet,
                    0,
                    0,
                    None,
                    None,
                    None,
                    f"read_error:{type(exc).__name__}",
                    None,
                    None,
                    None,
                    str(exc)[:300],
                )
            )
            continue

        preview = " | ".join(map(str, list(frame.columns)[:20]))
        lowered_sheet = sheet.lower()
        parsed_before = len(rows)
        accession_col = _select_accession_column(frame)
        site_col = _best_column(frame.columns, SITE_ALIASES)
        modification_col = _best_column(frame.columns, MODIFICATION_ALIASES)

        # PRM panels contain peptide-level validation but no absolute site. Preserve
        # them for sequence-based mapping during the UniProt annotation stage.
        if "prm" in lowered_sheet:
            peptide_col = _best_column(frame.columns, PEPTIDE_ALIASES)
            species, organism_id = infer_species(
                sheet, default="human" if "huvec" in lowered_sheet else "mouse"
            )
            if accession_col and peptide_col and species:
                for row_offset, (_, record) in enumerate(frame.iterrows(), start=1):
                    peptide = str(record.get(peptide_col, "")).strip()
                    clean, _ = parse_cml_peptide(peptide)
                    if clean is None:
                        continue  # unmodified PRM controls are not CML-site evidence
                    for accession in extract_accessions(record.get(accession_col)):
                        rows.append(
                            _base_evidence_row(
                                species=species,
                                organism_id=organism_id,
                                uniprot_accession_raw=accession,
                                site_uniprot=pd.NA,
                                source_label=source_label,
                                source_workbook=workbook.name,
                                source_sheet=sheet,
                                source_row=header_row + header_depth + row_offset + 1,
                                pairing_status="prm_requires_sequence_mapping",
                                peptide_sequence=peptide,
                                tissue_or_condition="HUVEC" if species == "human" else "PRM panel",
                                prm_validated=True,
                                evidence_type="prm_peptide",
                            )
                        )
            parsed_count = len(rows) - parsed_before
            diagnostics.append(
                ParseDiagnostics(
                    workbook.name,
                    sheet,
                    len(frame),
                    parsed_count,
                    accession_col,
                    None,
                    None,
                    "parsed_prm_peptides" if parsed_count else "parsed_zero_sites",
                    header_row + 1,
                    header_depth,
                    len(frame.columns),
                    preview,
                )
            )
            continue

        # The overlap table has two accession/site/peptide systems in each row.
        if "overlap" in lowered_sheet and {
            "HUVEC_ID",
            "HUVEC_Modified site",
            "MEF_ID",
            "MEF_Modified site",
        }.issubset(set(map(str, frame.columns))):
            for row_offset, (_, record) in enumerate(frame.iterrows(), start=1):
                identity = _numeric(record, "sequence_id")
                for species, organism_id, prefix, condition in (
                    ("human", 9606, "HUVEC", "HUVEC"),
                    ("mouse", 10090, "MEF", "MEF"),
                ):
                    accessions = extract_accessions(record.get(f"{prefix}_ID"))
                    sites = extract_sites(record.get(f"{prefix}_Modified site"))
                    peptide = str(record.get(f"{prefix}_peptide", "")).strip()
                    for accession, site, pairing in _pair_accessions_sites(accessions, sites):
                        rows.append(
                            _base_evidence_row(
                                species=species,
                                organism_id=organism_id,
                                uniprot_accession_raw=accession,
                                site_uniprot=site,
                                source_label=source_label,
                                source_workbook=workbook.name,
                                source_sheet=sheet,
                                source_row=header_row + header_depth + row_offset + 1,
                                pairing_status=f"ortholog_{pairing}",
                                peptide_sequence=peptide if _is_real_peptide(peptide) else "",
                                tissue_or_condition=condition,
                                ortholog_overlap=True,
                                sequence_identity=identity,
                                evidence_type="ortholog_overlap",
                            )
                        )
            parsed_count = len(rows) - parsed_before
            diagnostics.append(
                ParseDiagnostics(
                    workbook.name,
                    sheet,
                    len(frame),
                    parsed_count,
                    "HUVEC_ID | MEF_ID",
                    "HUVEC_Modified site | MEF_Modified site",
                    None,
                    "parsed_overlap_pairs" if parsed_count else "parsed_zero_sites",
                    header_row + 1,
                    header_depth,
                    len(frame.columns),
                    preview,
                )
            )
            continue

        # Old-vs-young limma tables are already CML-site tables. ``Modified site``
        # is the site column, not a modification-type field. FDR <= 0.05 defines
        # age association while every parsed K site remains marked as age-tested.
        if "aging_" in lowered_sheet and "_o_vs_y" in lowered_sheet:
            accession_col = _select_accession_column(frame)
            site_col = _best_column(frame.columns, SITE_ALIASES)
            tissue = _sheet_tissue(sheet)
            logfc_col = _best_column(frame.columns, ("logfc", "log2 fold change"))
            p_col = _best_column(frame.columns, ("p.value", "p value", "pvalue"))
            fdr_col = _best_column(frame.columns, ("adj.p.val", "adjusted p value", "fdr"))
            if accession_col and site_col:
                for row_offset, (_, record) in enumerate(frame.iterrows(), start=1):
                    accessions = extract_accessions(record.get(accession_col))
                    sites = extract_sites(record.get(site_col))
                    logfc = _numeric(record, logfc_col)
                    p_value = _numeric(record, p_col)
                    adj_p = _numeric(record, fdr_col)
                    significant = adj_p is not None and adj_p <= 0.05
                    direction = (
                        "old_higher"
                        if logfc is not None and logfc > 0
                        else "old_lower"
                        if logfc is not None and logfc < 0
                        else "no_change"
                    )
                    for accession, site, pairing in _pair_accessions_sites(accessions, sites):
                        rows.append(
                            _base_evidence_row(
                                species="mouse",
                                organism_id=10090,
                                uniprot_accession_raw=accession,
                                site_uniprot=site,
                                source_label=source_label,
                                source_workbook=workbook.name,
                                source_sheet=sheet,
                                source_row=header_row + header_depth + row_offset + 1,
                                pairing_status=pairing,
                                tissue_or_condition=tissue,
                                age_tested=True,
                                age_associated=significant,
                                age_log2fc=logfc,
                                age_p_value=p_value,
                                age_adj_p_value=adj_p,
                                age_direction=direction,
                                age_significant_tissue=tissue if significant else "",
                                evidence_type="age_test",
                            )
                        )
            parsed_count = len(rows) - parsed_before
            diagnostics.append(
                ParseDiagnostics(
                    workbook.name,
                    sheet,
                    len(frame),
                    parsed_count,
                    accession_col,
                    site_col,
                    None,
                    "parsed_age_sites" if parsed_count else "parsed_zero_sites",
                    header_row + 1,
                    header_depth,
                    len(frame.columns),
                    preview,
                )
            )
            continue

        # Proteostasis summary stores site as e.g. ``Heart, Kidney (Kidney324)``.
        if "proteostasis" in lowered_sheet:
            accession_col = _select_accession_column(frame)
            site_col = _best_column(frame.columns, SITE_ALIASES)
            if accession_col and site_col:
                for row_offset, (_, record) in enumerate(frame.iterrows(), start=1):
                    accessions = extract_accessions(record.get(accession_col))
                    sites, tissues = _proteostasis_site_and_tissues(record.get(site_col))
                    for accession, site, pairing in _pair_accessions_sites(accessions, sites):
                        rows.append(
                            _base_evidence_row(
                                species="mouse",
                                organism_id=10090,
                                uniprot_accession_raw=accession,
                                site_uniprot=site,
                                source_label=source_label,
                                source_workbook=workbook.name,
                                source_sheet=sheet,
                                source_row=header_row + header_depth + row_offset + 1,
                                pairing_status=pairing,
                                tissue_or_condition=tissues or sheet,
                                evidence_type="proteostasis_site",
                            )
                        )
            parsed_count = len(rows) - parsed_before
            diagnostics.append(
                ParseDiagnostics(
                    workbook.name,
                    sheet,
                    len(frame),
                    parsed_count,
                    accession_col,
                    site_col,
                    None,
                    "parsed_proteostasis_sites" if parsed_count else "parsed_zero_sites",
                    header_row + 1,
                    header_depth,
                    len(frame.columns),
                    preview,
                )
            )
            continue

        # Generic site-table fallback.
        accession_col = _select_accession_column(frame)
        site_col = _best_column(frame.columns, SITE_ALIASES)
        peptide_col = _best_column(frame.columns, PEPTIDE_ALIASES)
        tissue_col = _best_column(frame.columns, TISSUE_ALIASES)
        if accession_col is None:
            accession_col = _infer_content_column(frame, extract_accessions, minimum_matches=1)
        if site_col is None:
            site_col = _infer_content_column(
                frame, extract_sites, excluded=([accession_col] if accession_col else [])
            )
        validation_cols = [
            str(column)
            for column in frame.columns
            if any(alias in _canonical_column(column) for alias in VALIDATION_ALIASES)
        ]
        if accession_col is None or site_col is None:
            diagnostics.append(
                ParseDiagnostics(
                    workbook.name,
                    sheet,
                    len(frame),
                    0,
                    accession_col,
                    site_col,
                    modification_col,
                    "skipped_missing_required_columns",
                    header_row + 1,
                    header_depth,
                    len(frame.columns),
                    preview,
                )
            )
            continue

        for row_offset, (_, record) in enumerate(frame.iterrows(), start=1):
            row_text = " | ".join(str(value) for value in record.tolist() if pd.notna(value))
            species, organism_id = infer_species(
                f"{source_label} {sheet} {row_text}", default=default_species
            )
            if species is None:
                continue
            site_text = str(record.get(site_col, ""))
            if _bool_from_text(
                site_text, ("n-term", "nterm", "protnterm", "protein n termin")
            ) and not extract_sites(site_text):
                continue
            # Apply a modification-type restriction only to a distinct type column.
            if modification_col and modification_col != site_col:
                modification_text = str(record.get(modification_col, ""))
                if modification_text and not _bool_from_text(
                    modification_text.lower(),
                    ("cml", "carboxymethyl", "cm(k)", "cm-k", "cm k", "lys"),
                ):
                    continue
            accessions = extract_accessions(record.get(accession_col))
            sites = extract_sites(record.get(site_col))
            if not accessions or not sites:
                continue
            tissue_text = (
                str(record.get(tissue_col, "")).strip() if tissue_col else _sheet_tissue(sheet)
            )
            prm_validated = any(
                _bool_from_text(record.get(column, ""), ("yes", "true", "prm", "validated", "1"))
                for column in validation_cols
            )
            peptide = str(record.get(peptide_col, "")).strip() if peptide_col else ""
            if not _is_real_peptide(peptide):
                peptide = ""
            for accession, site, pairing_status in _pair_accessions_sites(accessions, sites):
                rows.append(
                    _base_evidence_row(
                        species=species,
                        organism_id=organism_id,
                        uniprot_accession_raw=accession,
                        site_uniprot=site,
                        source_label=source_label,
                        source_workbook=workbook.name,
                        source_sheet=sheet,
                        source_row=header_row + header_depth + row_offset + 1,
                        pairing_status=pairing_status,
                        peptide_sequence=peptide,
                        tissue_or_condition=tissue_text or sheet,
                        prm_validated=bool(prm_validated),
                        evidence_type="discovery_site",
                    )
                )
        parsed_count = len(rows) - parsed_before
        diagnostics.append(
            ParseDiagnostics(
                workbook.name,
                sheet,
                len(frame),
                parsed_count,
                accession_col,
                site_col,
                modification_col,
                "parsed" if parsed_count else "parsed_zero_sites",
                header_row + 1,
                header_depth,
                len(frame.columns),
                preview,
            )
        )
    return pd.DataFrame(rows), diagnostics


def aggregate_candidate_evidence(frame: pd.DataFrame) -> pd.DataFrame:
    if frame.empty:
        return frame.copy()
    frame = frame[pd.to_numeric(frame["site_uniprot"], errors="coerce").notna()].copy()
    if frame.empty:
        return frame
    frame["site_uniprot"] = pd.to_numeric(frame["site_uniprot"], errors="coerce").astype(int)
    key = ["species", "organism_id", "uniprot_accession_raw", "site_uniprot"]

    def join_unique(values: Iterable[object]) -> str:
        unique = sorted(
            {str(value).strip() for value in values if pd.notna(value) and str(value).strip()}
        )
        return " | ".join(unique)

    def join_peptides(values: Iterable[object]) -> str:
        return " | ".join(sorted({str(v).strip() for v in values if _is_real_peptide(v)}))

    aggregated = (
        frame.groupby(key, dropna=False)
        .agg(
            evidence_rows=("source_row", "size"),
            source_workbooks=("source_workbook", join_unique),
            source_sheets=("source_sheet", join_unique),
            evidence_types=("evidence_type", join_unique),
            tissues_conditions=("tissue_or_condition", join_unique),
            tissue_condition_count=(
                "tissue_or_condition",
                lambda x: len({str(v).strip() for v in x if pd.notna(v) and str(v).strip()}),
            ),
            prm_validated=("prm_validated", "max"),
            age_tested=("age_tested", "max"),
            age_associated=("age_associated", "max"),
            age_tested_tissues=("tissue_or_condition", join_unique),
            age_significant_tissues=("age_significant_tissue", join_unique),
            age_min_adj_p=("age_adj_p_value", "min"),
            age_max_abs_log2fc=(
                "age_log2fc",
                lambda x: pd.to_numeric(x, errors="coerce").abs().max(),
            ),
            age_directions=("age_direction", join_unique),
            ortholog_overlap=("ortholog_overlap", "max"),
            sequence_identity_max=("sequence_identity", "max"),
            peptide_examples=("peptide_sequence", join_peptides),
            ambiguous_mapping=(
                "pairing_status",
                lambda values: any(str(value).startswith("ambiguous") for value in values),
            ),
        )
        .reset_index()
    )
    # Only age-test evidence should populate the tested-tissue field.
    age_tissues = (
        frame[frame["age_tested"].fillna(False).astype(bool)]
        .groupby(key)["tissue_or_condition"]
        .agg(join_unique)
    )
    aggregated["age_tested_tissues"] = [
        age_tissues.get(tuple(row[k] for k in key), "") for _, row in aggregated.iterrows()
    ]
    aggregated["candidate_id"] = (
        aggregated["species"].astype(str)
        + ":"
        + aggregated["uniprot_accession_raw"].astype(str)
        + ":K"
        + aggregated["site_uniprot"].astype(int).astype(str)
    )
    aggregated["evidence_score_raw"] = (
        np.log1p(aggregated["evidence_rows"].astype(float))
        + 0.5 * np.log1p(aggregated["tissue_condition_count"].astype(float))
        + 1.0 * aggregated["prm_validated"].astype(float)
        + 0.5 * aggregated["age_associated"].astype(float)
        + 0.25 * aggregated["ortholog_overlap"].astype(float)
    )
    return aggregated.sort_values(["species", "uniprot_accession_raw", "site_uniprot"])


def prepare_candidate_table(
    workbooks: Sequence[Path],
    output_csv: Path,
    diagnostics_csv: Path,
    *,
    manual_candidate_csv: Path | None = None,
    prm_evidence_csv: Path | None = None,
) -> pd.DataFrame:
    raw_frames: list[pd.DataFrame] = []
    diagnostics: list[ParseDiagnostics] = []
    for workbook in workbooks:
        label = workbook.stem
        default = "mouse" if "mouse" in label.lower() or "organ" in label.lower() else None
        parsed, parsed_diagnostics = normalize_workbook(
            workbook, source_label=label, default_species=default
        )
        if not parsed.empty:
            raw_frames.append(parsed)
        diagnostics.extend(parsed_diagnostics)

    if manual_candidate_csv is not None:
        manual = pd.read_csv(manual_candidate_csv)
        required = {"species", "uniprot_accession_raw", "site_uniprot"}
        missing = required - set(manual.columns)
        if missing:
            raise ValueError(f"Manual candidate CSV is missing columns: {sorted(missing)}")
        if "organism_id" not in manual.columns:
            manual["organism_id"] = manual["species"].map({"human": 9606, "mouse": 10090})
        defaults = _base_evidence_row(
            source_row=np.arange(1, len(manual) + 1),
            source_workbook=manual_candidate_csv.name,
            source_sheet="manual",
            source_label="manual",
            pairing_status="manual",
            tissue_or_condition="manual",
            evidence_type="manual",
        )
        for column, value in defaults.items():
            if column not in manual.columns:
                manual[column] = value
        raw_frames.append(manual)

    diagnostics_csv.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame([diagnostic.__dict__ for diagnostic in diagnostics]).to_csv(
        diagnostics_csv, index=False
    )
    if not raw_frames:
        raise ValueError(
            "No CML site candidates were parsed. Diagnostics were written to "
            f"{diagnostics_csv}; alternatively provide --candidate-csv."
        )
    raw = pd.concat(raw_frames, ignore_index=True, sort=False)
    prm = raw[raw.get("evidence_type", "").eq("prm_peptide")].copy()
    if prm_evidence_csv is None:
        prm_evidence_csv = output_csv.with_name("stage4_prm_evidence.csv")
    prm_evidence_csv.parent.mkdir(parents=True, exist_ok=True)
    prm.to_csv(prm_evidence_csv, index=False)

    site_raw = raw[pd.to_numeric(raw["site_uniprot"], errors="coerce").notna()].copy()
    if site_raw.empty:
        raise ValueError(
            "Supplementary files contained no absolute CML lysine sites. "
            f"PRM-only evidence was written to {prm_evidence_csv}."
        )
    aggregated = aggregate_candidate_evidence(site_raw)
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    aggregated.to_csv(output_csv, index=False)
    return aggregated


def _chunked(values: Sequence[str], size: int) -> Iterable[list[str]]:
    for start in range(0, len(values), size):
        yield list(values[start : start + size])


def _empty_uniprot_records() -> pd.DataFrame:
    return pd.DataFrame(
        columns=[
            "uniprot_accession",
            "uniprot_entry_name",
            "protein_name",
            "gene_names",
            "organism_id_uniprot",
            "protein_length",
            "sequence",
            "reviewed_status",
            "subcellular_location",
        ]
    )


def fetch_uniprot_records(
    accessions: Sequence[str],
    cache_dir: Path,
    *,
    chunk_size: int = 80,
    pause_seconds: float = 0.1,
) -> pd.DataFrame:
    cache_dir.mkdir(parents=True, exist_ok=True)
    normalized = sorted(
        {
            accession
            for value in accessions
            if (accession := normalize_uniprot_accession(value, keep_isoform=False))
        }
    )
    if not normalized:
        return _empty_uniprot_records()

    missing: list[str] = []
    for accession in normalized:
        path = cache_dir / f"{accession}.json"
        if not path.exists():
            missing.append(accession)

    session = requests_session()
    for chunk in _chunked(missing, chunk_size):
        query = " OR ".join(f"accession:{accession}" for accession in chunk)
        response = session.get(
            UNIPROT_STREAM_URL,
            params={"query": f"({query})", "format": "tsv", "fields": UNIPROT_FIELDS},
            timeout=120,
        )
        try:
            response.raise_for_status()
        except requests.HTTPError as exc:
            detail = response.text.strip()[:1000]
            raise requests.HTTPError(
                f"UniProt request failed for validated accessions {chunk}: "
                f"HTTP {response.status_code}. {detail}",
                response=response,
            ) from exc
        if response.text.strip():
            from io import StringIO

            table = pd.read_csv(StringIO(response.text), sep="\t")
            for _, row in table.iterrows():
                accession = str(row.get("Entry", "")).upper()
                if not is_valid_uniprot_accession(accession):
                    continue
                record = {
                    "uniprot_accession": accession,
                    "uniprot_entry_name": row.get("Entry Name"),
                    "protein_name": row.get("Protein names"),
                    "gene_names": row.get("Gene Names"),
                    "organism_id_uniprot": row.get("Organism (ID)"),
                    "protein_length": row.get("Length"),
                    "sequence": row.get("Sequence"),
                    "reviewed_status": row.get("Reviewed"),
                    "subcellular_location": row.get("Subcellular location [CC]"),
                }
                (cache_dir / f"{accession}.json").write_text(
                    json.dumps(record, indent=2, default=str), encoding="utf-8"
                )
        time.sleep(pause_seconds)

    normalized_records: list[dict[str, object]] = []
    for accession in normalized:
        path = cache_dir / f"{accession}.json"
        if path.exists():
            normalized_records.append(json.loads(path.read_text(encoding="utf-8")))
    if not normalized_records:
        return _empty_uniprot_records()
    return pd.DataFrame(normalized_records)


def _append_pipe(existing: object, incoming: object) -> str:
    values: set[str] = set()
    for value in (existing, incoming):
        if value is None or (isinstance(value, float) and np.isnan(value)):
            continue
        values.update(part.strip() for part in str(value).split("|") if part.strip())
    return " | ".join(sorted(values))


def _recompute_evidence_score(frame: pd.DataFrame) -> pd.DataFrame:
    result = frame.copy()
    for column, default in (
        ("evidence_rows", 1),
        ("tissue_condition_count", 1),
        ("prm_validated", False),
        ("age_associated", False),
        ("ortholog_overlap", False),
    ):
        if column not in result.columns:
            result[column] = default
    result["evidence_score_raw"] = (
        np.log1p(pd.to_numeric(result["evidence_rows"], errors="coerce").fillna(1.0))
        + 0.5
        * np.log1p(pd.to_numeric(result["tissue_condition_count"], errors="coerce").fillna(1.0))
        + result["prm_validated"].fillna(False).astype(float)
        + 0.5 * result["age_associated"].fillna(False).astype(float)
        + 0.25 * result["ortholog_overlap"].fillna(False).astype(float)
    )
    return result


def _map_prm_evidence(
    prm: pd.DataFrame, records: pd.DataFrame
) -> tuple[pd.DataFrame, pd.DataFrame]:
    if prm.empty:
        empty = prm.copy()
        for column in (
            "uniprot_accession",
            "prm_mapping_status",
            "prm_mapping_reason",
            "mapped_site_uniprot",
        ):
            if column not in empty.columns:
                empty[column] = pd.Series(dtype="object")
        return empty.copy(), empty.copy()
    work = prm.copy()
    work["uniprot_accession"] = work["uniprot_accession_raw"].map(
        lambda value: normalize_uniprot_accession(value, keep_isoform=False)
    )
    work = work.merge(records, on="uniprot_accession", how="left")
    mapped_rows: list[dict[str, object]] = []
    excluded_rows: list[dict[str, object]] = []
    for _, row in work.iterrows():
        site, status, reason = map_cml_peptide_to_sequence(
            row.get("peptide_sequence"), row.get("sequence")
        )
        record = row.to_dict()
        record["prm_mapping_status"] = status
        record["prm_mapping_reason"] = reason
        record["mapped_site_uniprot"] = site
        expected = pd.to_numeric(row.get("organism_id"), errors="coerce")
        observed = pd.to_numeric(row.get("organism_id_uniprot"), errors="coerce")
        if (
            status == "mapped"
            and pd.notna(expected)
            and pd.notna(observed)
            and int(expected) != int(observed)
        ):
            record["prm_mapping_status"] = "excluded"
            record["prm_mapping_reason"] = "organism_mismatch"
            status = "excluded"
        (mapped_rows if status == "mapped" else excluded_rows).append(record)
    columns = list(work.columns) + [
        "prm_mapping_status",
        "prm_mapping_reason",
        "mapped_site_uniprot",
    ]
    mapped = pd.DataFrame(mapped_rows)
    excluded = pd.DataFrame(excluded_rows)
    if mapped.empty:
        mapped = pd.DataFrame(columns=columns)
    if excluded.empty:
        excluded = pd.DataFrame(columns=columns)
    return mapped, excluded


def annotate_candidates_with_uniprot(
    candidates_csv: Path,
    output_csv: Path,
    excluded_csv: Path,
    cache_dir: Path,
    *,
    prm_evidence_csv: Path | None = None,
    prm_mapping_csv: Path | None = None,
    prm_mapping_excluded_csv: Path | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    candidates = pd.read_csv(candidates_csv)
    if "uniprot_accession_raw" not in candidates.columns:
        raise ValueError(
            "Candidate table is missing 'uniprot_accession_raw'. Re-run Stage 4 prepare."
        )
    candidates["uniprot_accession"] = candidates["uniprot_accession_raw"].map(
        lambda value: normalize_uniprot_accession(value, keep_isoform=False)
    )
    invalid_format = candidates[candidates["uniprot_accession"].isna()].copy()
    if not invalid_format.empty:
        invalid_format["sequence_validation_status"] = "excluded"
        invalid_format["sequence_validation_reason"] = "invalid_uniprot_accession_format"
    valid_candidates = candidates[candidates["uniprot_accession"].notna()].copy()

    if prm_evidence_csv is None:
        prm_evidence_csv = candidates_csv.with_name("stage4_prm_evidence.csv")
    prm = pd.read_csv(prm_evidence_csv) if prm_evidence_csv.exists() else pd.DataFrame()
    if not prm.empty:
        prm["uniprot_accession"] = prm["uniprot_accession_raw"].map(
            lambda value: normalize_uniprot_accession(value, keep_isoform=False)
        )
        prm = prm[prm["uniprot_accession"].notna()].copy()

    all_accessions = valid_candidates["uniprot_accession"].tolist()
    if not prm.empty:
        all_accessions.extend(prm["uniprot_accession"].tolist())
    records = fetch_uniprot_records(all_accessions, cache_dir)

    mapped_prm, excluded_prm = (
        _map_prm_evidence(prm, records) if not prm.empty else (pd.DataFrame(), pd.DataFrame())
    )
    prm_mapping_csv = prm_mapping_csv or output_csv.with_name("stage4_prm_mapping.csv")
    prm_mapping_excluded_csv = prm_mapping_excluded_csv or excluded_csv.with_name(
        "stage4_excluded_prm_mapping.csv"
    )
    prm_mapping_csv.parent.mkdir(parents=True, exist_ok=True)
    mapped_prm.to_csv(prm_mapping_csv, index=False)
    excluded_prm.to_csv(prm_mapping_excluded_csv, index=False)

    # Merge mapped PRM evidence into existing candidates; create a new candidate
    # when a validated peptide identifies a site absent from discovery tables.
    if not mapped_prm.empty:
        for _, prm_row in mapped_prm.iterrows():
            accession = str(prm_row["uniprot_accession"])
            species = str(prm_row["species"])
            site = int(prm_row["mapped_site_uniprot"])
            mask = (
                (valid_candidates["species"].astype(str) == species)
                & (valid_candidates["uniprot_accession"].astype(str) == accession)
                & (pd.to_numeric(valid_candidates["site_uniprot"], errors="coerce") == site)
            )
            peptide = str(prm_row.get("peptide_sequence", ""))
            if mask.any():
                index = valid_candidates.index[mask][0]
                valid_candidates.at[index, "prm_validated"] = True
                current_evidence = pd.to_numeric(
                    valid_candidates.at[index, "evidence_rows"], errors="coerce"
                )
                valid_candidates.at[index, "evidence_rows"] = (
                    int(current_evidence) if pd.notna(current_evidence) else 0
                ) + 1
                valid_candidates.at[index, "source_workbooks"] = _append_pipe(
                    valid_candidates.at[index, "source_workbooks"], prm_row.get("source_workbook")
                )
                valid_candidates.at[index, "source_sheets"] = _append_pipe(
                    valid_candidates.at[index, "source_sheets"], prm_row.get("source_sheet")
                )
                valid_candidates.at[index, "evidence_types"] = _append_pipe(
                    valid_candidates.at[index, "evidence_types"]
                    if "evidence_types" in valid_candidates.columns
                    else "",
                    "prm_peptide",
                )
                valid_candidates.at[index, "tissues_conditions"] = _append_pipe(
                    valid_candidates.at[index, "tissues_conditions"],
                    prm_row.get("tissue_or_condition"),
                )
                valid_candidates.at[index, "tissue_condition_count"] = len(
                    [
                        part
                        for part in str(valid_candidates.at[index, "tissues_conditions"]).split("|")
                        if part.strip()
                    ]
                )
                valid_candidates.at[index, "peptide_examples"] = _append_pipe(
                    valid_candidates.at[index, "peptide_examples"], peptide
                )
            else:
                new_row = {column: np.nan for column in valid_candidates.columns}
                new_row.update(
                    {
                        "species": species,
                        "organism_id": int(prm_row["organism_id"]),
                        "uniprot_accession_raw": accession,
                        "uniprot_accession": accession,
                        "site_uniprot": site,
                        "evidence_rows": 1,
                        "source_workbooks": str(prm_row.get("source_workbook", "")),
                        "source_sheets": str(prm_row.get("source_sheet", "")),
                        "evidence_types": "prm_peptide",
                        "tissues_conditions": str(prm_row.get("tissue_or_condition", "PRM")),
                        "tissue_condition_count": 1,
                        "prm_validated": True,
                        "age_tested": False,
                        "age_associated": False,
                        "age_tested_tissues": "",
                        "age_significant_tissues": "",
                        "ortholog_overlap": False,
                        "peptide_examples": peptide,
                        "ambiguous_mapping": False,
                        "candidate_id": f"{species}:{accession}:K{site}",
                    }
                )
                valid_candidates = pd.concat(
                    [valid_candidates, pd.DataFrame([new_row])], ignore_index=True, sort=False
                )
        valid_candidates = _recompute_evidence_score(valid_candidates)

    if valid_candidates.empty:
        output_csv.parent.mkdir(parents=True, exist_ok=True)
        valid_candidates.to_csv(output_csv, index=False)
        invalid_format.to_csv(excluded_csv, index=False)
        raise ValueError("No valid UniProt accessions remain after strict format validation.")

    merged = valid_candidates.merge(records, on="uniprot_accession", how="left")
    merged["site_uniprot"] = pd.to_numeric(merged["site_uniprot"], errors="coerce").astype("Int64")

    def validate_row(row: pd.Series) -> tuple[str, str | None]:
        sequence = row.get("sequence")
        site = row.get("site_uniprot")
        if not isinstance(sequence, str) or not sequence:
            return "excluded", "missing_uniprot_sequence"
        if pd.isna(site) or int(site) < 1 or int(site) > len(sequence):
            return "excluded", "site_outside_canonical_sequence"
        if sequence[int(site) - 1].upper() != "K":
            return "excluded", f"canonical_residue_is_{sequence[int(site) - 1].upper()}_not_K"
        organism = pd.to_numeric(row.get("organism_id_uniprot"), errors="coerce")
        expected = pd.to_numeric(row.get("organism_id"), errors="coerce")
        if pd.notna(organism) and pd.notna(expected) and int(organism) != int(expected):
            return "excluded", "organism_mismatch"
        if "-" in str(row.get("uniprot_accession_raw", "")):
            return "warning", "isoform_mapped_to_canonical_sequence"
        return "eligible", None

    validation = merged.apply(validate_row, axis=1, result_type="expand")
    merged["sequence_validation_status"] = validation[0]
    merged["sequence_validation_reason"] = validation[1]
    eligible = merged[merged["sequence_validation_status"].isin(["eligible", "warning"])].copy()
    excluded_valid = merged[merged["sequence_validation_status"] == "excluded"].copy()
    excluded = pd.concat([invalid_format, excluded_valid], ignore_index=True, sort=False)
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    eligible.to_csv(output_csv, index=False)
    excluded.to_csv(excluded_csv, index=False)
    return eligible, excluded


def deterministic_quick_subset(
    frame: pd.DataFrame,
    *,
    max_proteins: int = 25,
    max_sites: int = 60,
) -> pd.DataFrame:
    """Return a deterministic, species-balanced quick subset.

    The prepare stage contains ``uniprot_accession_raw``; the canonical
    ``uniprot_accession`` column is only added during UniProt annotation.
    Accept both schemas so the helper is safe before and after annotation.
    Isoform suffixes are ignored only for protein-count and selection purposes;
    the original accession column is preserved in the returned frame.
    """
    if frame.empty:
        return frame.copy()

    accession_column = (
        "uniprot_accession"
        if "uniprot_accession" in frame.columns
        else "uniprot_accession_raw"
        if "uniprot_accession_raw" in frame.columns
        else None
    )
    if accession_column is None:
        raise ValueError(
            "Quick subset requires 'uniprot_accession_raw' before annotation "
            "or 'uniprot_accession' after annotation."
        )

    required = {
        "species",
        "prm_validated",
        "age_associated",
        "evidence_score_raw",
        "candidate_id",
    }
    missing = required - set(frame.columns)
    if missing:
        raise ValueError(f"Quick subset is missing columns: {sorted(missing)}")

    working = frame.copy()
    working["_quick_protein_key"] = (
        working[accession_column].astype(str).str.strip().str.upper().str.split("-").str[0]
    )

    if (
        len(working) <= max_sites
        and working["_quick_protein_key"].nunique(dropna=True) <= max_proteins
    ):
        return working.drop(columns="_quick_protein_key")

    ranked = working.sort_values(
        ["prm_validated", "age_associated", "evidence_score_raw", "candidate_id"],
        ascending=[False, False, False, True],
        kind="mergesort",
    )
    selected_proteins: set[str] = set()
    species_values = [value for value in ("human", "mouse") if value in set(ranked["species"])]
    other_species = sorted(set(ranked["species"].dropna()) - set(species_values))
    species_values.extend(other_species)

    per_species_limit = max(1, max_proteins // max(1, len(species_values)))
    for species in species_values:
        part = ranked[ranked["species"] == species]
        proteins: list[str] = []
        for accession in part["_quick_protein_key"].tolist():
            if not accession or accession == "NAN":
                continue
            if accession not in proteins:
                proteins.append(accession)
            if len(proteins) >= per_species_limit:
                break
        selected_proteins.update(proteins)

    # Fill unused capacity globally when one species has fewer proteins.
    if len(selected_proteins) < max_proteins:
        for accession in ranked["_quick_protein_key"].tolist():
            if not accession or accession == "NAN":
                continue
            selected_proteins.add(accession)
            if len(selected_proteins) >= max_proteins:
                break

    selected = ranked[ranked["_quick_protein_key"].isin(selected_proteins)].head(max_sites)
    return selected.drop(columns="_quick_protein_key").copy()
