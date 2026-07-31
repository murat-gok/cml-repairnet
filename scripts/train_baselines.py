from pathlib import Path

from cml_repairnet.models.baselines import load_merged_features, run_nested_group_cv
from cml_repairnet.utils.config import load_config
from cml_repairnet.utils.reproducibility import set_global_seed

config = load_config()
seed = int(config["project"]["seed"])
set_global_seed(seed)
paths = config["paths"]
validation = config["validation"]
merged = load_merged_features(
    Path(paths["labels"]), Path(paths["sequence_features"]), Path(paths["structure_features"])
)
metrics, _ = run_nested_group_cv(
    merged,
    target=config["labels"]["target"],
    output_dir=Path(paths["results_dir"]),
    n_spatial_groups=int(validation["spatial_groups"]),
    bootstrap_iterations=int(validation["bootstrap_iterations"]),
    seed=seed,
)
print(metrics.to_string(index=False))
