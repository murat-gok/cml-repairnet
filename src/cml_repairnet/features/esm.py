from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import torch
from Bio import SeqIO
from transformers import AutoModel, AutoTokenizer


def extract_residue_embeddings(
    labels_csv: Path,
    fasta_path: Path,
    output_npz: Path,
    model_name: str = "facebook/esm2_t12_35M_UR50D",
    device: str | None = None,
) -> None:
    """Extract frozen ESM-2 embeddings at labelled lysines.

    The first special token occupies index 0, so UniProt residue p maps to token index p.
    """
    labels = pd.read_csv(labels_csv)
    records = list(SeqIO.parse(str(fasta_path), "fasta"))
    if len(records) != 1:
        raise ValueError("Expected one FASTA record.")
    sequence = str(records[0].seq)

    selected_device = device or ("cuda" if torch.cuda.is_available() else "cpu")
    tokenizer = AutoTokenizer.from_pretrained(model_name)
    model = AutoModel.from_pretrained(model_name).to(selected_device).eval()
    tokens = tokenizer(sequence, return_tensors="pt", add_special_tokens=True)
    tokens = {key: value.to(selected_device) for key, value in tokens.items()}

    with torch.inference_mode():
        hidden = model(**tokens).last_hidden_state[0].cpu().numpy()

    sites = labels["site_uniprot"].astype(int).to_numpy()
    embeddings = np.vstack([hidden[site] for site in sites])
    output_npz.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        output_npz,
        protein_id=labels["protein_id"].to_numpy(),
        site_uniprot=sites,
        embeddings=embeddings,
        model_name=model_name,
    )
