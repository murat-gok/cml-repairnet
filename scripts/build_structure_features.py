from pathlib import Path

from cml_repairnet.features.structure import build_structure_features
from cml_repairnet.utils.config import load_config

config = load_config()
paths = config["paths"]
feature_config = config["features"]
frame = build_structure_features(
    Path(paths["labels"]),
    Path(paths["bsa_pdb"]),
    Path(paths["structure_features"]),
    lysine_max_asa=float(feature_config["lysine_max_asa_angstrom2"]),
    contact_radii=tuple(float(value) for value in feature_config["contact_radii_angstrom"]),
    chain_id=str(feature_config.get("pdb_chain", "A")),
    pdb_site_column=str(feature_config.get("pdb_site_column", "site_paper")),
)
print(f"Built {frame.shape[1] - 2} structure-derived features for {len(frame)} sites.")
