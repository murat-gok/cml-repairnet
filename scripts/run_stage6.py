from __future__ import annotations

import argparse
import json
from pathlib import Path

from cml_repairnet.robustness.pipeline import (
    conservation_stage6,
    doctor_stage6,
    enrichment_stage6,
    figures_stage6,
    network_stage6,
    rank_stage6,
    report_stage6,
    run_stage6_all,
    tiers_stage6,
)


def project_root_from_script() -> Path:
    return Path(__file__).resolve().parents[1]


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="CML-RepairNet Stage 6 robustness analysis")
    parser.add_argument(
        "command",
        choices=(
            "doctor",
            "rank",
            "conservation",
            "enrichment",
            "network",
            "tiers",
            "figures",
            "report",
            "all",
        ),
    )
    parser.add_argument("--project-root", type=Path, default=project_root_from_script())
    return parser


def main() -> None:
    args = build_parser().parse_args()
    root = args.project_root.resolve()
    commands = {
        "doctor": doctor_stage6,
        "rank": rank_stage6,
        "conservation": conservation_stage6,
        "enrichment": enrichment_stage6,
        "network": network_stage6,
        "tiers": tiers_stage6,
        "figures": figures_stage6,
        "report": report_stage6,
        "all": run_stage6_all,
    }
    result = commands[args.command](root)
    print(json.dumps(result, indent=2, default=str))


if __name__ == "__main__":
    main()
