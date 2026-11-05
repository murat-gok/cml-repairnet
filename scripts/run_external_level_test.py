"""Run the level-resolved external bulk-assay consistency test (revision A2).

Additive: reads the same inputs as the frozen external-scaffold pipeline and
writes to a separate results directory. Does not retrain or rescore anything.

Usage
-----
  python scripts/run_external_level_test.py
  python scripts/run_external_level_test.py --no-winsorize   # sensitivity run
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from cml_repairnet.revision.external_level_test import LevelTestPaths, run

DEFAULT_LABELS = Path("data/external/trabosh2026_bulk_protein_repair.csv")
DEFAULT_AGGREGATES = Path(
    "results/external_scaffold_validation/external_scaffold_aggregates.csv"
)
DEFAULT_OUT = Path("results/external_level_test")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--labels", type=Path, default=DEFAULT_LABELS)
    parser.add_argument("--aggregates", type=Path, default=DEFAULT_AGGREGATES)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT)
    parser.add_argument(
        "--no-winsorize",
        action="store_true",
        help="keep repair fractions above 1.0 as reported (sensitivity analysis)",
    )
    args = parser.parse_args()

    winsorize = not args.no_winsorize
    out_dir = args.out_dir if winsorize else args.out_dir.with_name(args.out_dir.name + "_raw")

    report = run(
        LevelTestPaths(labels=args.labels, aggregates=args.aggregates, out_dir=out_dir),
        winsorize=winsorize,
    )
    print(json.dumps(report, indent=2, ensure_ascii=False, default=str))
    print(f"\nOutputs written to: {out_dir}")


if __name__ == "__main__":
    main()
