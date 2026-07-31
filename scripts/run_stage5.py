from __future__ import annotations

import argparse
import json
from pathlib import Path

from cml_repairnet.biology.pipeline import (
    annotate_stage5,
    doctor_stage5,
    enrich_stage5,
    figures_stage5,
    modules_stage5,
    network_stage5,
    orthology_stage5,
    report_stage5,
    run_stage5_all,
)


def project_root_from_script() -> Path:
    return Path(__file__).resolve().parents[1]


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="CML-RepairNet Stage 5 biological interpretation")
    parser.add_argument("command", choices=("doctor", "annotate", "enrich", "network", "orthology", "modules", "figures", "report", "all"))
    parser.add_argument("--project-root", type=Path, default=project_root_from_script())
    parser.add_argument("--offline", action="store_true", help="Use cached external-API responses only")
    return parser


def main() -> None:
    args = build_parser().parse_args()
    root = args.project_root.resolve()
    command = args.command
    if command == "doctor": result = doctor_stage5(root)
    elif command == "annotate": result = annotate_stage5(root, offline=args.offline)
    elif command == "enrich": result = enrich_stage5(root, offline=args.offline)
    elif command == "network": result = network_stage5(root, offline=args.offline)
    elif command == "orthology": result = orthology_stage5(root, offline=args.offline)
    elif command == "modules": result = modules_stage5(root)
    elif command == "figures": result = figures_stage5(root)
    elif command == "report": result = report_stage5(root)
    elif command == "all": result = run_stage5_all(root, offline=args.offline)
    else: raise AssertionError(command)
    print(json.dumps(result, indent=2, default=str))


if __name__ == "__main__":
    main()
