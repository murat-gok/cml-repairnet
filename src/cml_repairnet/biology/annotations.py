from __future__ import annotations

import re

import pandas as pd

from cml_repairnet.biology.clients import OfflineCacheMiss, UniProtAnnotationClient
from cml_repairnet.biology.schema import unique_accessions


def _find_column(frame: pd.DataFrame, candidates: tuple[str, ...]) -> str | None:
    normalized = {str(column).strip().lower(): column for column in frame.columns}
    for candidate in candidates:
        if candidate.lower() in normalized:
            return normalized[candidate.lower()]
    for key, original in normalized.items():
        if any(candidate.lower() in key for candidate in candidates):
            return original
    return None


def normalize_uniprot_annotations(frame: pd.DataFrame) -> pd.DataFrame:
    if frame.empty:
        return frame.copy()
    accession = _find_column(frame, ("Entry", "accession"))
    if accession is None:
        raise ValueError("UniProt annotation response does not contain accession/Entry")
    rename = {accession: "uniprot_accession"}
    mappings = {
        "Entry Name": "uniprot_entry_name_annotation",
        "Gene Names": "gene_names_annotation",
        "Protein names": "protein_name_annotation",
        "Protein families": "protein_families",
        "Keywords": "uniprot_keywords",
        "Gene Ontology (biological process)": "go_biological_process",
        "Gene Ontology (cellular component)": "go_cellular_component",
        "Gene Ontology (molecular function)": "go_molecular_function",
        "InterPro": "interpro_ids",
        "Pfam": "pfam_ids",
        "Function [CC]": "function_comment",
        "Pathway": "pathway_comment",
        "Subcellular location [CC]": "subcellular_location_annotation",
    }
    for source, target in mappings.items():
        column = _find_column(frame, (source,))
        if column:
            rename[column] = target
    result = frame.rename(columns=rename).copy()
    result["uniprot_accession"] = result["uniprot_accession"].astype(str).str.split("-").str[0]
    return result.drop_duplicates("uniprot_accession")


def family_summary(annotations: pd.DataFrame) -> pd.DataFrame:
    if annotations.empty or "protein_families" not in annotations:
        return pd.DataFrame(columns=["protein_family", "proteins", "accessions"])
    rows = []
    for row in annotations[["uniprot_accession", "protein_families"]].dropna().itertuples(index=False):
        text = str(row.protein_families)
        families = [part.strip() for part in re.split(r"[;|]", text) if part.strip()]
        for family in families:
            rows.append({"uniprot_accession": row.uniprot_accession, "protein_family": family})
    if not rows:
        return pd.DataFrame(columns=["protein_family", "proteins", "accessions"])
    frame = pd.DataFrame(rows).drop_duplicates()
    return (
        frame.groupby("protein_family")
        .agg(
            proteins=("uniprot_accession", "nunique"),
            accessions=("uniprot_accession", lambda x: " | ".join(sorted(set(x)))),
        )
        .reset_index()
        .sort_values(["proteins", "protein_family"], ascending=[False, True])
    )


def run_uniprot_annotations(
    portfolio: pd.DataFrame,
    client: UniProtAnnotationClient,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    accessions = unique_accessions(portfolio)
    try:
        raw = client.fetch(accessions, query_name="stage5_portfolio")
        annotations = normalize_uniprot_annotations(raw)
        families = family_summary(annotations)
        status, error = "completed", ""
    except OfflineCacheMiss as exc:
        annotations, families = pd.DataFrame(), pd.DataFrame()
        status, error = "offline_cache_missing", str(exc)
    except Exception as exc:
        annotations, families = pd.DataFrame(), pd.DataFrame()
        status, error = "failed", f"{type(exc).__name__}: {exc}"
    metadata = client.cached_metadata()
    manifest = pd.DataFrame([{
        "input_accessions": len(accessions),
        "annotated_accessions": len(annotations),
        "status": status,
        "error": error,
        "response_metadata_json": __import__("json").dumps(metadata),
    }])
    return annotations, families, manifest
