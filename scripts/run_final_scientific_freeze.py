from __future__ import annotations

import argparse
import json
from pathlib import Path

from cml_repairnet.revision.final_freeze import build_final_freeze, doctor_final_freeze


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Build the CML-RepairNet final scientific freeze outputs")
    parser.add_argument("command", choices=("doctor", "build"))
    parser.add_argument("--root", type=Path, default=Path("."))
    parser.add_argument("--output", type=Path, default=Path("results/final_scientific_freeze"))
    return parser


def main() -> None:
    args = build_parser().parse_args()
    root = args.root.resolve()
    output = args.output if args.output.is_absolute() else (root / args.output).resolve()
    if args.command == "doctor":
        result = doctor_final_freeze(root, output)
    else:
        result = build_final_freeze(root, output)
    print(json.dumps(result, indent=2, default=str))


if __name__ == "__main__":
    main()
