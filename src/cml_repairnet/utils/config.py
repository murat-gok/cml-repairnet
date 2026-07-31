from __future__ import annotations

from pathlib import Path
from typing import Any

import yaml


def load_config(path: str | Path = "configs/base.yaml") -> dict[str, Any]:
    """Load a YAML configuration file."""
    config_path = Path(path)
    if not config_path.exists():
        raise FileNotFoundError(f"Configuration file not found: {config_path}")
    with config_path.open("r", encoding="utf-8") as handle:
        data = yaml.safe_load(handle)
    if not isinstance(data, dict):
        raise ValueError("Configuration root must be a mapping.")
    return data
