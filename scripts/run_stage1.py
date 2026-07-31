from __future__ import annotations

import argparse

from cml_repairnet.pipeline import run_stage1
from cml_repairnet.utils.config import load_config


def main() -> None:
    parser = argparse.ArgumentParser(description="Run the CML-RepairNet data-foundation stage.")
    parser.add_argument("--config", default="configs/base.yaml")
    parser.add_argument("--no-download", action="store_true")
    args = parser.parse_args()
    report = run_stage1(load_config(args.config), download=not args.no_download)
    print("Stage 1 PASS")
    for key, value in report.items():
        if key != "notes":
            print(f"{key}: {value}")


if __name__ == "__main__":
    main()
