from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from cml_repairnet.virtual_screening import select_multifidelity_batch


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Select candidates for higher-fidelity in silico simulation."
    )
    parser.add_argument("--candidates", type=Path, required=True)
    parser.add_argument("--output", type=Path, default=Path("results/virtual_simulation_batch.csv"))
    parser.add_argument("--features", nargs="+", required=True)
    parser.add_argument("--batch-size", type=int, default=24)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    candidates = pd.read_csv(args.candidates)
    selected = select_multifidelity_batch(
        candidates,
        feature_columns=args.features,
        batch_size=args.batch_size,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    selected.to_csv(args.output, index=False)
    print(f"Wrote {len(selected)} candidates to {args.output}")


if __name__ == "__main__":
    main()
