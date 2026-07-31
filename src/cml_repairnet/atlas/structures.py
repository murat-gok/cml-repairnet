from __future__ import annotations

import json
import math
import re
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Sequence

import numpy as np
import pandas as pd
import requests
from Bio.PDB import MMCIFParser, PDBParser, ShrakeRupley
from Bio.PDB.Polypeptide import protein_letters_3to1

from cml_repairnet.atlas.sources import download_file, requests_session

ALPHAFOLD_API = "https://alphafold.ebi.ac.uk/api/prediction/{accession}"
HYDROPHOBIC = set("AILMFWVY")
AROMATIC = set("FWY")
TRANSFER_FEATURES = (
    "residue_sasa",
    "relative_sasa_lysine",
    "distance_to_protein_centroid",
    "neighbor_hydrophobic_fraction_8a",
    "neighbor_cysteine_count_8a",
    "neighbor_aromatic_count_8a",
    "ca_contact_count_8a",
    "ca_contact_count_10a",
    "ca_contact_count_12a",
)


_RETRYABLE_STRUCTURE_REASON_PREFIXES = (
    "structure_processing_error:",
)


def _is_retryable_structure_reason(value: object) -> bool:
    text = "" if value is None else str(value)
    return any(text.startswith(prefix) for prefix in _RETRYABLE_STRUCTURE_REASON_PREFIXES)


def _partition_cached_rows(
    candidates: pd.DataFrame,
    cached: pd.DataFrame | None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Split a per-protein cache into stable rows and rows that require work.

    Stage 4 quick and full runs may request different sites from the same protein.
    A cache file is therefore reusable only for candidate IDs that are present and
    non-retryable. Transient network/processing failures are deliberately retried
    on the next structures run.
    """
    required = candidates.copy()
    if cached is None or cached.empty or "candidate_id" not in cached.columns:
        return pd.DataFrame(), required

    current_ids = set(required["candidate_id"].astype(str))
    cached = cached[cached["candidate_id"].astype(str).isin(current_ids)].copy()
    reason = cached.get("structure_reason", pd.Series("", index=cached.index))
    retryable = reason.map(_is_retryable_structure_reason)
    stable = cached.loc[~retryable].drop_duplicates("candidate_id", keep="last")
    covered = set(stable["candidate_id"].astype(str))
    pending = required.loc[~required["candidate_id"].astype(str).isin(covered)].copy()
    return stable, pending


def _request_with_outer_retry(operation, *, attempts: int = 3, backoff_seconds: float = 2.0):
    """Retry request-layer failures beyond the HTTP adapter retry budget."""
    last_error = None
    for attempt in range(1, max(1, int(attempts)) + 1):
        try:
            return operation()
        except requests.RequestException as exc:
            last_error = exc
            if attempt >= attempts:
                raise
            time.sleep(float(backoff_seconds) * (2 ** (attempt - 1)))
    raise last_error  # pragma: no cover

_AF_ENTRY_ACCESSION = re.compile(
    r"^AF-(?P<accession>[A-Z0-9]{6,10}(?:-\d+)?)-F\d+(?:-v\d+)?$",
    re.IGNORECASE,
)


def _residue_letter(name: str) -> str:
    return protein_letters_3to1.get(str(name).upper(), "X")


def _clean_sequence(value: object) -> str | None:
    if value is None:
        return None
    sequence = re.sub(r"\s+", "", str(value)).upper()
    return sequence or None


def _accession_from_entry_id(entry_id: str) -> str | None:
    match = _AF_ENTRY_ACCESSION.match(str(entry_id).strip())
    return match.group("accession").upper() if match else None


@dataclass(frozen=True)
class AlphaFoldEntry:
    accession: str
    entry_id: str
    uniprot_start: int
    uniprot_end: int
    cif_url: str
    pdb_url: str | None
    model_version: int | None
    global_metric_value: float | None
    sequence: str | None = None
    model_accession: str | None = None
    source_database: str | None = None

    def contains(self, position: int) -> bool:
        return self.uniprot_start <= int(position) <= self.uniprot_end

    @property
    def span(self) -> int:
        return max(0, int(self.uniprot_end) - int(self.uniprot_start) + 1)

    def local_index(self, position: int) -> int:
        return int(position) - int(self.uniprot_start)


def fetch_alphafold_metadata(
    accession: str,
    cache_dir: Path,
    *,
    session: requests.Session | None = None,
    overwrite: bool = False,
) -> list[AlphaFoldEntry]:
    cache_dir.mkdir(parents=True, exist_ok=True)
    cache_path = cache_dir / f"{accession}.json"
    if cache_path.exists() and not overwrite:
        payload = json.loads(cache_path.read_text(encoding="utf-8"))
    else:
        session = session or requests_session()
        response = session.get(ALPHAFOLD_API.format(accession=accession), timeout=120)
        if response.status_code == 404:
            payload = []
        else:
            response.raise_for_status()
            payload = response.json()
        cache_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")

    entries: list[AlphaFoldEntry] = []
    for item in payload if isinstance(payload, list) else []:
        sequence = _clean_sequence(item.get("sequence") or item.get("uniprotSequence"))
        start = int(item.get("sequenceStart") or item.get("uniprotStart") or 1)
        end_value = item.get("sequenceEnd") or item.get("uniprotEnd")
        end = int(end_value or (start + len(sequence) - 1 if sequence else 10**9))
        cif_url = item.get("cifUrl") or item.get("bcifUrl")
        if not cif_url:
            continue
        metric = item.get("globalMetricValue")
        entry_id = str(item.get("modelEntityId") or item.get("entryId") or accession)
        model_accession = (
            _accession_from_entry_id(entry_id)
            or item.get("uniprotAccession")
            or item.get("uniprotId")
        )
        entries.append(
            AlphaFoldEntry(
                accession=str(accession).upper(),
                entry_id=entry_id,
                uniprot_start=start,
                uniprot_end=end,
                cif_url=str(cif_url),
                pdb_url=item.get("pdbUrl"),
                model_version=(int(item["latestVersion"]) if item.get("latestVersion") else None),
                global_metric_value=(float(metric) if metric is not None else None),
                sequence=sequence,
                model_accession=(str(model_accession).upper() if model_accession else None),
                source_database=(
                    str(item.get("sourceDatabase") or item.get("sourceDatabaseName"))
                    if item.get("sourceDatabase") or item.get("sourceDatabaseName")
                    else None
                ),
            )
        )
    return entries


def _entry_sequence_matches_reference(
    entry: AlphaFoldEntry,
    reference_sequence: str | None,
) -> bool | None:
    if not reference_sequence or not entry.sequence:
        return None
    reference = _clean_sequence(reference_sequence)
    if reference is None:
        return None
    start = int(entry.uniprot_start) - 1
    end = int(entry.uniprot_end)
    if start < 0 or end > len(reference):
        return False
    return reference[start:end] == entry.sequence


def _entry_site_residue_matches(
    entry: AlphaFoldEntry,
    position: int,
    reference_sequence: str | None,
) -> bool | None:
    if not entry.contains(position):
        return False
    if not reference_sequence or not entry.sequence:
        return None
    reference = _clean_sequence(reference_sequence)
    if reference is None or position < 1 or position > len(reference):
        return False
    local = entry.local_index(position)
    if local < 0 or local >= len(entry.sequence):
        return False
    return entry.sequence[local] == reference[position - 1]


def select_entry_for_site(
    entries: Sequence[AlphaFoldEntry],
    position: int,
    *,
    expected_accession: str | None = None,
    expected_sequence: str | None = None,
) -> AlphaFoldEntry | None:
    compatible = [entry for entry in entries if entry.contains(position)]
    if not compatible:
        return None

    expected_accession = str(expected_accession or "").upper() or None

    # If API-provided sequences are available, reject models that conflict with
    # the canonical UniProt residue at the candidate site. This prevents an
    # isoform model from being used with canonical numbering.
    site_compatible = [
        entry
        for entry in compatible
        if _entry_site_residue_matches(entry, position, expected_sequence) is not False
    ]
    if site_compatible:
        compatible = site_compatible

    def rank(entry: AlphaFoldEntry) -> tuple[object, ...]:
        model_accession = (
            str(entry.model_accession).upper()
            if entry.model_accession
            else _accession_from_entry_id(entry.entry_id)
        )
        exact_accession = bool(expected_accession and model_accession == expected_accession)
        full_sequence_match = _entry_sequence_matches_reference(entry, expected_sequence)
        site_match = _entry_site_residue_matches(entry, position, expected_sequence)
        edge_margin = min(
            int(position) - int(entry.uniprot_start),
            int(entry.uniprot_end) - int(position),
        )
        confidence = (
            float(entry.global_metric_value)
            if entry.global_metric_value is not None
            else -1.0
        )
        return (
            0 if exact_accession else 1,
            0 if full_sequence_match is True else 1 if full_sequence_match is None else 2,
            0 if site_match is True else 1 if site_match is None else 2,
            -edge_margin,
            -confidence,
            -entry.span,
            entry.entry_id,
        )

    return sorted(compatible, key=rank)[0]


def download_alphafold_structure(
    entry: AlphaFoldEntry,
    structure_dir: Path,
    *,
    session: requests.Session | None = None,
    overwrite: bool = False,
) -> Path:
    structure_dir.mkdir(parents=True, exist_ok=True)
    safe_entry = entry.entry_id.replace("/", "_")
    destination = structure_dir / f"{safe_entry}.cif"
    download_file(
        entry.cif_url,
        destination,
        session=session,
        overwrite=overwrite,
    )
    return destination


def _load_structure(path: Path):
    suffix = path.suffix.lower()
    if suffix in {".cif", ".mmcif"}:
        parser = MMCIFParser(QUIET=True)
    elif suffix in {".pdb", ".ent"}:
        parser = PDBParser(QUIET=True)
    else:
        raise ValueError(f"Unsupported structure format: {path}")
    structure = parser.get_structure(path.stem, str(path))
    return next(structure.get_models())


def _standard_residues(model) -> list[Any]:
    residues: list[Any] = []
    for chain in model:
        for residue in chain:
            if residue.id[0] == " " and "CA" in residue:
                residues.append(residue)
    return residues


def _resolve_residue_number(
    residues: Sequence[Any],
    site_uniprot: int,
    entry: AlphaFoldEntry,
) -> tuple[Any | None, int | None, str]:
    if not residues:
        return None, None, "no_residues"

    # The API sequence is the safest bridge between canonical UniProt numbering
    # and a coordinate file, including locally numbered long-protein fragments.
    if entry.sequence:
        local_index = entry.local_index(site_uniprot)
        if 0 <= local_index < len(entry.sequence) and local_index < len(residues):
            residue = residues[local_index]
            expected = entry.sequence[local_index]
            observed = _residue_letter(residue.resname)
            if observed == expected:
                return residue, int(residue.id[1]), "entry_sequence_index"

    by_number = {int(residue.id[1]): residue for residue in residues}
    min_number = min(by_number)

    if entry.uniprot_start > 1:
        local_number = min_number + int(site_uniprot) - int(entry.uniprot_start)
        if local_number in by_number:
            return by_number[local_number], local_number, "fragment_local_numbering"
        if site_uniprot in by_number:
            return by_number[site_uniprot], site_uniprot, "direct_uniprot_numbering"
    else:
        if site_uniprot in by_number:
            return by_number[site_uniprot], site_uniprot, "direct_uniprot_numbering"
        local_number = min_number + int(site_uniprot) - int(entry.uniprot_start)
        if local_number in by_number:
            return by_number[local_number], local_number, "fragment_local_numbering"

    return None, None, "site_not_in_structure"


def _neighbor_features(
    residue: Any,
    residues: Sequence[Any],
    centroid: np.ndarray,
    *,
    lysine_max_asa: float = 205.0,
) -> dict[str, object]:
    ca = np.asarray(residue["CA"].coord, dtype=float)
    all_ca = np.vstack([other["CA"].coord for other in residues])
    distances = np.linalg.norm(all_ca - ca, axis=1)
    neighbors = [
        other
        for other in residues
        if other is not residue and np.linalg.norm(other["CA"].coord - ca) <= 8.0
    ]
    letters = [_residue_letter(other.resname) for other in neighbors]
    confidences = [float(other["CA"].bfactor) for other in neighbors]
    site_plddt = float(residue["CA"].bfactor)
    sasa = float(getattr(residue, "sasa", np.nan))
    return {
        "ca_x": float(ca[0]),
        "ca_y": float(ca[1]),
        "ca_z": float(ca[2]),
        "site_plddt": site_plddt,
        "neighbor_mean_plddt_8a": float(np.mean(confidences)) if confidences else np.nan,
        "low_confidence_neighbor_fraction_8a": (
            float(np.mean(np.asarray(confidences) < 70.0)) if confidences else 0.0
        ),
        "residue_sasa": sasa,
        "relative_sasa_lysine": min(sasa / lysine_max_asa, 1.5) if np.isfinite(sasa) else np.nan,
        "distance_to_protein_centroid": float(np.linalg.norm(ca - centroid)),
        "neighbor_hydrophobic_fraction_8a": (
            float(np.mean([letter in HYDROPHOBIC for letter in letters])) if letters else 0.0
        ),
        "neighbor_cysteine_count_8a": int(letters.count("C")),
        "neighbor_aromatic_count_8a": int(sum(letter in AROMATIC for letter in letters)),
        "ca_contact_count_8a": int(np.sum((distances > 0) & (distances <= 8.0))),
        "ca_contact_count_10a": int(np.sum((distances > 0) & (distances <= 10.0))),
        "ca_contact_count_12a": int(np.sum((distances > 0) & (distances <= 12.0))),
    }


def extract_site_features_from_structure(
    candidates: pd.DataFrame,
    structure_path: Path,
    entry: AlphaFoldEntry,
) -> pd.DataFrame:
    model = _load_structure(structure_path)
    ShrakeRupley().compute(model, level="R")
    residues = _standard_residues(model)
    if not residues:
        raise ValueError(f"No standard CA-containing residues found in {structure_path}")
    centroid = np.vstack([residue["CA"].coord for residue in residues]).mean(axis=0)
    rows: list[dict[str, object]] = []
    for record in candidates.itertuples(index=False):
        site = int(record.site_uniprot)
        reference_sequence = _clean_sequence(getattr(record, "sequence", None))
        model_accession = (
            entry.model_accession or _accession_from_entry_id(entry.entry_id)
        )
        sequence_match = _entry_sequence_matches_reference(entry, reference_sequence)
        site_match = _entry_site_residue_matches(entry, site, reference_sequence)
        residue, structure_number, mapping = _resolve_residue_number(residues, site, entry)
        base = {
            "candidate_id": record.candidate_id,
            "species": record.species,
            "uniprot_accession": record.uniprot_accession,
            "site_uniprot": site,
            "alphafold_entry_id": entry.entry_id,
            "alphafold_model_accession": model_accession,
            "alphafold_uniprot_start": entry.uniprot_start,
            "alphafold_uniprot_end": entry.uniprot_end,
            "alphafold_model_version": entry.model_version,
            "alphafold_global_metric": entry.global_metric_value,
            "alphafold_sequence_match": sequence_match,
            "alphafold_site_residue_match": site_match,
            "structure_path": str(structure_path),
            "structure_source": "AlphaFoldDB",
            "confidence_kind": "pLDDT",
            "structure_residue_number": structure_number,
            "structure_numbering_mapping": mapping,
        }
        if sequence_match is False or site_match is False:
            rows.append(
                {
                    **base,
                    "structure_status": "excluded",
                    "structure_reason": "alphafold_sequence_conflicts_with_uniprot",
                }
            )
            continue
        if residue is None:
            rows.append({**base, "structure_status": "excluded", "structure_reason": mapping})
            continue
        letter = _residue_letter(residue.resname)
        if letter != "K":
            rows.append(
                {
                    **base,
                    "structure_status": "excluded",
                    "structure_reason": f"structure_residue_is_{letter}_not_K",
                }
            )
            continue
        rows.append(
            {
                **base,
                "structure_status": "eligible",
                "structure_reason": None,
                **_neighbor_features(residue, residues, centroid),
            }
        )
    return pd.DataFrame(rows)


def _process_accession(
    accession: str,
    candidates: pd.DataFrame,
    metadata_dir: Path,
    structure_dir: Path,
    per_protein_dir: Path,
    overwrite: bool,
) -> pd.DataFrame:
    output = per_protein_dir / f"{accession}.csv"
    cached = None
    if output.exists() and not overwrite:
        try:
            cached = pd.read_csv(output)
        except Exception:
            cached = None
    stable_cached, pending = (
        (pd.DataFrame(), candidates.copy())
        if overwrite
        else _partition_cached_rows(candidates, cached)
    )
    if pending.empty:
        return stable_cached.reset_index(drop=True)

    session = requests_session()
    try:
        entries = _request_with_outer_retry(
            lambda: fetch_alphafold_metadata(
                accession, metadata_dir, session=session, overwrite=overwrite
            )
        )
        if not entries:
            new_frame = pending[
                ["candidate_id", "species", "uniprot_accession", "site_uniprot"]
            ].copy()
            new_frame["structure_status"] = "excluded"
            new_frame["structure_reason"] = "alphafold_entry_not_found"
        else:
            reference_sequences = (
                pending.get("sequence", pd.Series(dtype=object))
                .dropna()
                .astype(str)
                .map(_clean_sequence)
                .dropna()
                .unique()
                .tolist()
            )
            reference_sequence = reference_sequences[0] if reference_sequences else None

            parts: list[pd.DataFrame] = []
            entry_groups: dict[str, tuple[AlphaFoldEntry, list[int]]] = {}
            missing_rows: list[dict[str, object]] = []
            for index, row in pending.iterrows():
                entry = select_entry_for_site(
                    entries,
                    int(row["site_uniprot"]),
                    expected_accession=accession,
                    expected_sequence=reference_sequence,
                )
                if entry is None:
                    missing_rows.append(
                        {
                            "candidate_id": row["candidate_id"],
                            "species": row["species"],
                            "uniprot_accession": accession,
                            "site_uniprot": int(row["site_uniprot"]),
                            "structure_status": "excluded",
                            "structure_reason": "no_sequence_compatible_alphafold_entry",
                        }
                    )
                    continue
                entry_groups.setdefault(entry.entry_id, (entry, []))[1].append(index)
            for entry, indices in entry_groups.values():
                structure_path = _request_with_outer_retry(
                    lambda entry=entry: download_alphafold_structure(
                        entry, structure_dir, session=session, overwrite=overwrite
                    )
                )
                parts.append(
                    extract_site_features_from_structure(
                        pending.loc[indices], structure_path, entry
                    )
                )
            if missing_rows:
                parts.append(pd.DataFrame(missing_rows))
            new_frame = (
                pd.concat(parts, ignore_index=True, sort=False)
                if parts
                else pd.DataFrame(missing_rows)
            )
    except Exception as exc:  # retain failure, but retry it automatically next run
        new_frame = pending[
            ["candidate_id", "species", "uniprot_accession", "site_uniprot"]
        ].copy()
        new_frame["structure_status"] = "excluded"
        new_frame["structure_reason"] = (
            f"structure_processing_error:{type(exc).__name__}:{exc}"
        )

    combined = pd.concat([stable_cached, new_frame], ignore_index=True, sort=False)
    current_ids = set(candidates["candidate_id"].astype(str))
    combined = combined[combined["candidate_id"].astype(str).isin(current_ids)]
    combined = combined.drop_duplicates("candidate_id", keep="last")
    output.parent.mkdir(parents=True, exist_ok=True)
    combined.to_csv(output, index=False)
    return combined.reset_index(drop=True)


def build_alphafold_features(
    candidates_csv: Path,
    output_csv: Path,
    excluded_csv: Path,
    cache_root: Path,
    *,
    workers: int = 4,
    overwrite: bool = False,
    min_site_plddt: float = 70.0,
    max_low_confidence_neighbor_fraction: float = 0.50,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    candidates = pd.read_csv(candidates_csv)
    metadata_dir = cache_root / "metadata"
    structure_dir = cache_root / "structures"
    per_protein_dir = cache_root / "features_by_protein"
    per_protein_dir.mkdir(parents=True, exist_ok=True)
    groups = list(candidates.groupby("uniprot_accession", sort=True))
    frames: list[pd.DataFrame] = []
    with ThreadPoolExecutor(max_workers=max(1, int(workers))) as executor:
        futures = {
            executor.submit(
                _process_accession,
                str(accession),
                group.copy(),
                metadata_dir,
                structure_dir,
                per_protein_dir,
                overwrite,
            ): str(accession)
            for accession, group in groups
        }
        for future in as_completed(futures):
            frames.append(future.result())
    features = (
        pd.concat(frames, ignore_index=True, sort=False) if frames else pd.DataFrame()
    )
    merged = candidates.merge(
        features,
        on=["candidate_id", "species", "uniprot_accession", "site_uniprot"],
        how="left",
        suffixes=("", "_structure"),
    )
    eligible_mask = merged["structure_status"].eq("eligible")
    quality_reason = np.where(
        pd.to_numeric(merged.get("site_plddt"), errors="coerce") < min_site_plddt,
        "site_plddt_below_threshold",
        np.where(
            pd.to_numeric(
                merged.get("low_confidence_neighbor_fraction_8a"), errors="coerce"
            )
            > max_low_confidence_neighbor_fraction,
            "low_confidence_neighborhood",
            "",
        ),
    )
    merged["structure_quality_status"] = np.where(
        eligible_mask & (quality_reason == ""), "eligible", "excluded"
    )
    merged["structure_quality_reason"] = np.where(
        ~eligible_mask,
        merged.get("structure_reason", "structure_not_eligible").fillna(
            "structure_not_eligible"
        ),
        quality_reason,
    )
    eligible = merged[merged["structure_quality_status"] == "eligible"].copy()
    excluded = merged[merged["structure_quality_status"] == "excluded"].copy()
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    eligible.to_csv(output_csv, index=False)
    excluded.to_csv(excluded_csv, index=False)
    return eligible, excluded
