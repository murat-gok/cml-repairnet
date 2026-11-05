from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import joblib
import numpy as np
import pandas as pd
import requests
from Bio.PDB import MMCIFParser, PDBParser, ShrakeRupley
from Bio.PDB.Polypeptide import protein_letters_3to1
from scipy.stats import kendalltau, spearmanr
from sklearn.base import clone
from sklearn.impute import SimpleImputer
from sklearn.metrics import mean_absolute_error
from sklearn.neighbors import NearestNeighbors
from sklearn.preprocessing import StandardScaler

from cml_repairnet.models.stage3 import stage3_frozen_estimators
from cml_repairnet.models.stage25 import build_sample_weights

HYDROPHOBIC = set("AILMFWVY")
AROMATIC = set("FWY")
DEFAULT_FEATURES = (
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
MODEL_WEIGHT_MODES: Mapping[str, str] = {
    "extra_trees": "none",
    "lightgbm": "clipped",
    "ridge": "none",
}


@dataclass(frozen=True)
class QualityThresholds:
    strict_site_plddt: float = 70.0
    strict_neighbor_plddt: float = 70.0
    relaxed_site_plddt: float = 50.0
    relaxed_neighbor_plddt: float = 50.0
    ood_percentile_cutoff: float = 0.95
    ood_neighbors: int = 5


def _residue_letter(name: str) -> str:
    return protein_letters_3to1.get(str(name).upper(), "X")


def _feature_position(feature: Mapping[str, Any], side: str) -> int | None:
    try:
        value = feature["location"][side].get("value")
        return int(value) if value is not None else None
    except (KeyError, TypeError, ValueError):
        return None


def select_mature_range(record: Mapping[str, Any]) -> tuple[int, int, str]:
    sequence = str(record.get("sequence", {}).get("value", ""))
    if not sequence:
        raise ValueError("UniProt record contains no canonical sequence")
    candidates: list[tuple[int, int, str]] = []
    for feature in record.get("features", []):
        if str(feature.get("type", "")).lower() != "chain":
            continue
        start = _feature_position(feature, "start")
        end = _feature_position(feature, "end")
        if start is None or end is None or end < start:
            continue
        description = str(feature.get("description", "canonical chain"))
        candidates.append((start, end, description))
    if candidates:
        start, end, description = max(candidates, key=lambda item: item[1] - item[0] + 1)
        return start, min(end, len(sequence)), description

    signal_end = 0
    for feature in record.get("features", []):
        if str(feature.get("type", "")).lower() == "signal":
            end = _feature_position(feature, "end")
            if end is not None:
                signal_end = max(signal_end, end)
    return signal_end + 1, len(sequence), "canonical sequence excluding annotated signal peptide"


def _request_json(url: str, timeout: float = 60.0) -> Any:
    response = requests.get(url, timeout=timeout, headers={"User-Agent": "CML-RepairNet/0.10"})
    response.raise_for_status()
    return response.json()


def _download(url: str, path: Path, timeout: float = 120.0) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists() and path.stat().st_size > 100:
        return path
    response = requests.get(url, timeout=timeout, headers={"User-Agent": "CML-RepairNet/0.10"})
    response.raise_for_status()
    path.write_bytes(response.content)
    return path


def fetch_uniprot_record(accession: str, cache_dir: Path) -> dict[str, Any]:
    path = cache_dir / "uniprot" / f"{accession}.json"
    if path.exists():
        return json.loads(path.read_text(encoding="utf-8"))
    record = _request_json(f"https://rest.uniprot.org/uniprotkb/{accession}.json")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(record, indent=2), encoding="utf-8")
    return record


def _choose_alphafold_entry(entries: Sequence[Mapping[str, Any]], accession: str) -> Mapping[str, Any]:
    exact = [
        entry for entry in entries
        if str(entry.get("uniprotAccession", entry.get("uniprotId", ""))).upper()
        == accession.upper()
    ]
    candidates = exact or list(entries)
    if not candidates:
        raise ValueError(f"No AlphaFold entry returned for {accession}")
    return max(
        candidates,
        key=lambda entry: int(entry.get("sequenceEnd", 0) or 0)
        - int(entry.get("sequenceStart", 1) or 1),
    )


def fetch_alphafold_structure(accession: str, cache_dir: Path) -> tuple[Path, dict[str, Any]]:
    metadata_path = cache_dir / "alphafold" / f"{accession}.json"
    if metadata_path.exists():
        entries = json.loads(metadata_path.read_text(encoding="utf-8"))
    else:
        entries = _request_json(f"https://alphafold.ebi.ac.uk/api/prediction/{accession}")
        metadata_path.parent.mkdir(parents=True, exist_ok=True)
        metadata_path.write_text(json.dumps(entries, indent=2), encoding="utf-8")
    if isinstance(entries, Mapping):
        entries = [entries]
    entry = dict(_choose_alphafold_entry(entries, accession))
    url = entry.get("pdbUrl") or entry.get("cifUrl") or entry.get("mmcifUrl")
    if not url:
        raise ValueError(f"AlphaFold entry for {accession} has no coordinate URL")
    suffix = ".cif" if "cif" in str(url).lower() else ".pdb"
    path = cache_dir / "alphafold" / f"{accession}{suffix}"
    return _download(str(url), path), entry


def _parse_structure(path: Path):
    if path.suffix.lower() in {".cif", ".mmcif"}:
        return MMCIFParser(QUIET=True).get_structure(path.stem, str(path))
    return PDBParser(QUIET=True).get_structure(path.stem, str(path))


def extract_transfer_safe_features(
    *,
    accession: str,
    scaffold: str,
    component: str,
    structure_path: Path,
    sequence: str,
    mature_start: int,
    mature_end: int,
    lysine_max_asa: float = 205.0,
) -> pd.DataFrame:
    structure = _parse_structure(structure_path)
    model = next(structure.get_models())
    chain = next(model.get_chains())
    ShrakeRupley().compute(model, level="R")
    residues = [res for res in chain if res.id[0] == " " and "CA" in res]
    by_number = {int(res.id[1]): res for res in residues}
    scope = [res for res in residues if mature_start <= int(res.id[1]) <= mature_end]
    if not scope:
        raise ValueError(f"No structural residues overlap mature range for {accession}")
    scope_ca = np.vstack([res["CA"].coord for res in scope]).astype(float)
    centroid = scope_ca.mean(axis=0)

    rows: list[dict[str, Any]] = []
    for position in range(max(1, mature_start), min(len(sequence), mature_end) + 1):
        if sequence[position - 1] != "K":
            continue
        residue = by_number.get(position)
        if residue is None or _residue_letter(residue.resname) != "K":
            rows.append(
                {
                    "scaffold": scaffold,
                    "component": component,
                    "uniprot_accession": accession,
                    "site_uniprot": position,
                    "structure_missing": True,
                }
            )
            continue
        ca = np.asarray(residue["CA"].coord, dtype=float)
        distances = np.linalg.norm(scope_ca - ca, axis=1)
        neighbors = [
            other for other in scope
            if other is not residue and np.linalg.norm(np.asarray(other["CA"].coord) - ca) <= 8.0
        ]
        neighbor_letters = [_residue_letter(other.resname) for other in neighbors]
        neighbor_plddt = [float(other["CA"].bfactor) for other in neighbors]
        sasa = float(getattr(residue, "sasa", np.nan))
        row: dict[str, Any] = {
            "scaffold": scaffold,
            "component": component,
            "uniprot_accession": accession,
            "site_uniprot": position,
            "structure_missing": False,
            "site_plddt": float(residue["CA"].bfactor),
            "neighbor_mean_plddt_8a": float(np.mean(neighbor_plddt)) if neighbor_plddt else np.nan,
            "residue_sasa": sasa,
            "relative_sasa_lysine": min(sasa / lysine_max_asa, 1.5),
            "distance_to_protein_centroid": float(np.linalg.norm(ca - centroid)),
            "neighbor_hydrophobic_fraction_8a": (
                float(np.mean([aa in HYDROPHOBIC for aa in neighbor_letters]))
                if neighbor_letters else 0.0
            ),
            "neighbor_cysteine_count_8a": int(neighbor_letters.count("C")),
            "neighbor_aromatic_count_8a": int(sum(aa in AROMATIC for aa in neighbor_letters)),
        }
        for radius in (8.0, 10.0, 12.0):
            row[f"ca_contact_count_{radius:g}a"] = int(
                np.sum((distances > 1e-8) & (distances <= radius))
            )
        rows.append(row)
    return pd.DataFrame(rows)


def prepare_external_sites(panel: pd.DataFrame, cache_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    all_rows: list[pd.DataFrame] = []
    manifest: list[dict[str, Any]] = []
    for record in panel.itertuples(index=False):
        accession = str(record.uniprot_accession)
        try:
            uniprot = fetch_uniprot_record(accession, cache_dir)
            sequence = str(uniprot["sequence"]["value"])
            mature_start, mature_end, chain_note = select_mature_range(uniprot)
            structure_path, af_entry = fetch_alphafold_structure(accession, cache_dir)
            frame = extract_transfer_safe_features(
                accession=accession,
                scaffold=str(record.scaffold),
                component=str(record.component),
                structure_path=structure_path,
                sequence=sequence,
                mature_start=mature_start,
                mature_end=mature_end,
            )
            frame["stoichiometric_weight"] = float(record.stoichiometric_weight)
            frame["equal_component_weight"] = float(record.equal_component_weight)
            all_rows.append(frame)
            manifest.append(
                {
                    "scaffold": record.scaffold,
                    "component": record.component,
                    "uniprot_accession": accession,
                    "status": "PASS",
                    "sequence_length": len(sequence),
                    "mature_start": mature_start,
                    "mature_end": mature_end,
                    "mature_range_note": chain_note,
                    "n_lysines": int((sequence[mature_start - 1:mature_end]).count("K")),
                    "n_structure_rows": len(frame),
                    "alphafold_entry_id": af_entry.get("entryId", af_entry.get("modelEntityId", "")),
                    "coordinate_file": str(structure_path),
                }
            )
        except Exception as exc:
            manifest.append(
                {
                    "scaffold": record.scaffold,
                    "component": record.component,
                    "uniprot_accession": accession,
                    "status": "FAIL",
                    "error": f"{type(exc).__name__}: {exc}",
                }
            )
    sites = pd.concat(all_rows, ignore_index=True) if all_rows else pd.DataFrame()
    return sites, pd.DataFrame(manifest)


def fit_ood_reference(
    train: pd.DataFrame,
    feature_columns: Sequence[str],
    n_neighbors: int = 5,
) -> dict[str, Any]:
    imputer = SimpleImputer(strategy="median")
    X = imputer.fit_transform(train[list(feature_columns)])
    scaler = StandardScaler().fit(X)
    Xz = scaler.transform(X)
    k = max(1, min(int(n_neighbors), len(Xz) - 1))
    nn = NearestNeighbors(n_neighbors=k + 1).fit(Xz)
    reference = nn.kneighbors(Xz, return_distance=True)[0][:, 1:].mean(axis=1)
    return {"imputer": imputer, "scaler": scaler, "reference": reference, "Xz": Xz, "k": k}


def score_ood(frame: pd.DataFrame, reference: Mapping[str, Any], feature_columns: Sequence[str]) -> pd.DataFrame:
    result = frame.copy()
    X = reference["imputer"].transform(result[list(feature_columns)])
    Xz = reference["scaler"].transform(X)
    nn = NearestNeighbors(n_neighbors=int(reference["k"])).fit(reference["Xz"])
    distance = nn.kneighbors(Xz, return_distance=True)[0].mean(axis=1)
    ref = np.asarray(reference["reference"], dtype=float)
    percentile = np.asarray([(ref <= value).mean() for value in distance], dtype=float)
    result["ood_distance"] = distance
    result["ood_percentile"] = percentile
    return result


def fit_and_score_models(
    *,
    labels: pd.DataFrame,
    structure_features: pd.DataFrame,
    external_sites: pd.DataFrame,
    target: str = "repair_published",
    feature_columns: Sequence[str] = DEFAULT_FEATURES,
    model_names: Sequence[str] = ("extra_trees", "lightgbm", "ridge"),
    seed: int = 20260727,
    model_output: Path | None = None,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    train = labels.merge(structure_features, on=["protein_id", "site_uniprot"], validate="one_to_one")
    missing = [column for column in feature_columns if column not in train or column not in external_sites]
    if missing:
        raise ValueError(f"Missing transfer-safe features: {missing}")
    y = train[target].to_numpy(dtype=float)
    X = train[list(feature_columns)].copy()
    X_external = external_sites[list(feature_columns)].copy()
    estimators = stage3_frozen_estimators(seed)
    fitted: dict[str, Any] = {}
    scored = external_sites.copy()
    for name in model_names:
        if name not in estimators:
            continue
        estimator = clone(estimators[name])
        weights = build_sample_weights(train, MODEL_WEIGHT_MODES.get(name, "none"))
        try:
            estimator.fit(X, y, model__sample_weight=weights)
        except (ImportError, ModuleNotFoundError) as exc:
            scored[f"pred_{name}"] = np.nan
            scored[f"model_error_{name}"] = str(exc)
            continue
        scored[f"pred_{name}"] = np.clip(estimator.predict(X_external), 0.0, 1.0)
        fitted[name] = estimator
    prediction_columns = [column for column in scored if column.startswith("pred_")]
    scored["pred_consensus"] = scored[prediction_columns].mean(axis=1, skipna=True)
    scored["pred_disagreement"] = scored[prediction_columns].std(axis=1, ddof=0, skipna=True)
    if model_output is not None:
        model_output.parent.mkdir(parents=True, exist_ok=True)
        joblib.dump({"models": fitted, "features": list(feature_columns), "seed": seed}, model_output)
    return scored, {"train": train, "fitted_models": fitted}


def add_quality_flags(frame: pd.DataFrame, thresholds: QualityThresholds) -> pd.DataFrame:
    result = frame.copy()
    result["quality_strict"] = (
        ~result["structure_missing"].fillna(True)
        & result["site_plddt"].ge(thresholds.strict_site_plddt)
        & result["neighbor_mean_plddt_8a"].ge(thresholds.strict_neighbor_plddt)
    )
    result["quality_relaxed"] = (
        ~result["structure_missing"].fillna(True)
        & result["site_plddt"].ge(thresholds.relaxed_site_plddt)
        & result["neighbor_mean_plddt_8a"].ge(thresholds.relaxed_neighbor_plddt)
    )
    result["in_domain"] = result["ood_percentile"].lt(thresholds.ood_percentile_cutoff)
    return result


def _component_site_weights(group: pd.DataFrame, component_scheme: str, site_scheme: str) -> np.ndarray:
    component_col = (
        "stoichiometric_weight" if component_scheme == "stoichiometric" else "equal_component_weight"
    )
    component_weight = group[component_col].astype(float).to_numpy()
    counts = group.groupby("component")["site_uniprot"].transform("count").astype(float).to_numpy()
    if site_scheme == "sasa":
        local = np.maximum(group["relative_sasa_lysine"].astype(float).to_numpy(), 0.01)
        denom = group.assign(_local=local).groupby("component")["_local"].transform("sum").to_numpy()
        return component_weight * local / np.maximum(denom, 1e-12)
    return component_weight / np.maximum(counts, 1.0)


def aggregate_predictions(
    frame: pd.DataFrame,
    *,
    models: Sequence[str] = ("extra_trees", "lightgbm", "ridge", "consensus"),
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    quality_modes = {
        "strict": frame["quality_strict"],
        "relaxed": frame["quality_relaxed"],
        "all_structured": ~frame["structure_missing"].fillna(True),
    }
    for quality_name, quality_mask in quality_modes.items():
        for domain_name, domain_mask in {
            "in_domain": frame["in_domain"],
            "all_domain": pd.Series(True, index=frame.index),
        }.items():
            eligible = frame[quality_mask & domain_mask].copy()
            for scaffold, group in eligible.groupby("scaffold"):
                total_sites = int((frame["scaffold"] == scaffold).sum())
                for component_scheme in ("stoichiometric", "equal_component"):
                    for model in models:
                        pred_col = f"pred_{model}"
                        if pred_col not in group or not group[pred_col].notna().any():
                            continue
                        for aggregation, site_scheme in (
                            ("mean", "uniform"),
                            ("sasa_weighted_mean", "sasa"),
                        ):
                            weights = _component_site_weights(group, component_scheme, site_scheme)
                            values = group[pred_col].to_numpy(dtype=float)
                            estimate = float(np.average(values, weights=weights))
                            rows.append(
                                {
                                    "scaffold": scaffold,
                                    "quality_mode": quality_name,
                                    "domain_mode": domain_name,
                                    "component_scheme": component_scheme,
                                    "model": model,
                                    "aggregation": aggregation,
                                    "predicted_repair": estimate,
                                    "n_eligible_sites": len(group),
                                    "n_total_sites": total_sites,
                                    "site_coverage": len(group) / total_sites if total_sites else np.nan,
                                    "n_components": group["component"].nunique(),
                                }
                            )
                        component_medians = group.groupby("component")[pred_col].median()
                        rows.append(
                            {
                                "scaffold": scaffold,
                                "quality_mode": quality_name,
                                "domain_mode": domain_name,
                                "component_scheme": component_scheme,
                                "model": model,
                                "aggregation": "median_of_component_medians",
                                "predicted_repair": float(component_medians.median()),
                                "n_eligible_sites": len(group),
                                "n_total_sites": total_sites,
                                "site_coverage": len(group) / total_sites if total_sites else np.nan,
                                "n_components": group["component"].nunique(),
                            }
                        )
                        accessible = group[
                            group["relative_sasa_lysine"]
                            >= group.groupby("component")["relative_sasa_lysine"].transform(
                                lambda values: values.quantile(0.75)
                            )
                        ]
                        if not accessible.empty:
                            weights = _component_site_weights(
                                accessible, component_scheme, "uniform"
                            )
                            rows.append(
                                {
                                    "scaffold": scaffold,
                                    "quality_mode": quality_name,
                                    "domain_mode": domain_name,
                                    "component_scheme": component_scheme,
                                    "model": model,
                                    "aggregation": "top_accessible_quartile_mean",
                                    "predicted_repair": float(
                                        np.average(accessible[pred_col], weights=weights)
                                    ),
                                    "n_eligible_sites": len(accessible),
                                    "n_total_sites": total_sites,
                                    "site_coverage": len(accessible) / total_sites if total_sites else np.nan,
                                    "n_components": accessible["component"].nunique(),
                                }
                            )
    return pd.DataFrame(rows)


def observed_summary(labels: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for scaffold, group in labels.groupby("scaffold"):
        raw = group["repair_fraction"].astype(float).to_numpy()
        rows.append(
            {
                "scaffold": scaffold,
                "observed_repair_median": float(np.median(raw)),
                "observed_repair_min": float(np.min(raw)),
                "observed_repair_max": float(np.max(raw)),
                "observed_repair_median_clipped": float(np.clip(np.median(raw), 0.0, 1.0)),
                "observed_rank": int(pd.Series(-np.median(raw)).rank(method="dense").iloc[0]),
                "n_conditions": len(group),
            }
        )
    result = pd.DataFrame(rows)
    result["observed_rank"] = result["observed_repair_median"].rank(
        ascending=False, method="dense"
    ).astype(int)
    return result.sort_values("observed_rank").reset_index(drop=True)


def compare_rankings(aggregates: pd.DataFrame, labels: pd.DataFrame) -> pd.DataFrame:
    observed = observed_summary(labels)
    key_columns = ["quality_mode", "domain_mode", "component_scheme", "model", "aggregation"]
    rows: list[dict[str, Any]] = []
    for keys, group in aggregates.groupby(key_columns, dropna=False):
        merged = (
            observed.merge(group, on="scaffold", how="inner")
            .sort_values("observed_rank")
            .reset_index(drop=True)
        )
        if len(merged) != len(observed):
            continue
        observed_values = merged["observed_repair_median_clipped"].to_numpy(dtype=float)
        predicted = merged["predicted_repair"].to_numpy(dtype=float)
        rho = float(spearmanr(observed_values, predicted).statistic)
        tau = float(kendalltau(observed_values, predicted).statistic)
        predicted_rank = pd.Series(predicted).rank(ascending=False, method="dense").astype(int)
        exact = bool(np.array_equal(predicted_rank.to_numpy(), merged["observed_rank"].to_numpy()))
        pairwise_total = 0
        pairwise_correct = 0
        for i in range(len(merged)):
            for j in range(i + 1, len(merged)):
                observed_sign = np.sign(observed_values[i] - observed_values[j])
                predicted_sign = np.sign(predicted[i] - predicted[j])
                pairwise_total += 1
                pairwise_correct += int(observed_sign == predicted_sign)
        row = dict(zip(key_columns, keys, strict=True))
        row.update(
            {
                "n_scaffolds": len(merged),
                "spearman": rho,
                "kendall_tau": tau,
                "exact_rank_match": exact,
                "pairwise_order_accuracy": pairwise_correct / pairwise_total,
                "calibration_mae_descriptive": float(mean_absolute_error(observed_values, predicted)),
                "minimum_site_coverage": float(merged["site_coverage"].min()),
                "interpretation": "protein-level external consistency; not site-level validation",
            }
        )
        rows.append(row)
    return pd.DataFrame(rows)


def bootstrap_rank_probability(
    sites: pd.DataFrame,
    labels: pd.DataFrame,
    *,
    model: str = "extra_trees",
    quality_mode: str = "strict",
    domain_mode: str = "in_domain",
    component_scheme: str = "stoichiometric",
    iterations: int = 2000,
    seed: int = 20260727,
) -> pd.DataFrame:
    quality_col = "quality_strict" if quality_mode == "strict" else "quality_relaxed"
    mask = sites[quality_col]
    if domain_mode == "in_domain":
        mask &= sites["in_domain"]
    work = sites[mask].copy()
    pred_col = f"pred_{model}"
    observed_order = observed_summary(labels)["scaffold"].tolist()
    rng = np.random.default_rng(seed)
    records: list[dict[str, Any]] = []
    required = set(observed_order)
    if set(work["scaffold"]) != required:
        return pd.DataFrame(
            [{
                "iteration": -1,
                "exact_rank_match": False,
                "status": "insufficient_scaffold_coverage",
                "available_scaffolds": ";".join(sorted(set(work["scaffold"]))),
            }]
        )
    for iteration in range(iterations):
        scaffold_scores: dict[str, float] = {}
        for scaffold, scaffold_group in work.groupby("scaffold"):
            component_scores: list[tuple[float, float]] = []
            weight_col = (
                "stoichiometric_weight"
                if component_scheme == "stoichiometric"
                else "equal_component_weight"
            )
            for _, component_group in scaffold_group.groupby("component"):
                sampled_idx = rng.choice(
                    component_group.index.to_numpy(), size=len(component_group), replace=True
                )
                sampled = component_group.loc[sampled_idx]
                sasa = np.maximum(sampled["relative_sasa_lysine"].to_numpy(dtype=float), 0.01)
                score = float(np.average(sampled[pred_col].to_numpy(dtype=float), weights=sasa))
                component_scores.append((score, float(component_group[weight_col].iloc[0])))
            values = np.asarray([item[0] for item in component_scores], dtype=float)
            weights = np.asarray([item[1] for item in component_scores], dtype=float)
            scaffold_scores[scaffold] = float(np.average(values, weights=weights))
        predicted_order = [
            name for name, _ in sorted(scaffold_scores.items(), key=lambda item: item[1], reverse=True)
        ]
        records.append(
            {
                "iteration": iteration,
                "exact_rank_match": predicted_order == observed_order,
                "predicted_order": ">".join(predicted_order),
                **{f"pred_{name}": scaffold_scores[name] for name in observed_order},
            }
        )
    result = pd.DataFrame(records)
    probability = float(result["exact_rank_match"].mean())
    result["exact_rank_probability"] = probability
    result["status"] = "PASS"
    return result


def data_fingerprint(paths: Iterable[Path]) -> str:
    digest = hashlib.sha256()
    for path in sorted(paths, key=lambda item: str(item)):
        digest.update(str(path).encode())
        digest.update(path.read_bytes())
    return digest.hexdigest()
