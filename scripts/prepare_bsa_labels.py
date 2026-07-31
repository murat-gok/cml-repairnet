from pathlib import Path

from cml_repairnet.data.prepare_labels import build_bsa_labels
from cml_repairnet.utils.config import load_config

config = load_config()
paths = config["paths"]
labels = config["labels"]
frame = build_bsa_labels(
    Path(paths["article_source_xlsx"]),
    Path(paths["labels"]),
    fasta_path=Path(paths["bsa_fasta"]),
    precursor_offset=int(labels["precursor_offset"]),
    variance_floor=float(labels["variance_floor"]),
)
print(frame.head().to_string(index=False))
print(f"Prepared {len(frame)} labelled sites.")
