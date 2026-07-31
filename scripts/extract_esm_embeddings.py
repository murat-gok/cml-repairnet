from pathlib import Path

from cml_repairnet.features.esm import extract_residue_embeddings
from cml_repairnet.utils.config import load_config

config = load_config()
paths = config["paths"]
extract_residue_embeddings(
    Path(paths["labels"]),
    Path(paths["bsa_fasta"]),
    Path("data/processed/bsa_esm_embeddings.npz"),
    model_name=config["features"]["esm_model"],
)
print("ESM embeddings saved.")
