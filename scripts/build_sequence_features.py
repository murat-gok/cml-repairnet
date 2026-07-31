from pathlib import Path

from cml_repairnet.features.sequence import build_sequence_features
from cml_repairnet.utils.config import load_config

config = load_config()
paths = config["paths"]
radius = int(config["features"]["sequence_window_radius"])
frame = build_sequence_features(
    Path(paths["labels"]), Path(paths["bsa_fasta"]), Path(paths["sequence_features"]), radius=radius
)
print(f"Built {frame.shape[1] - 3} sequence-derived features for {len(frame)} sites.")
