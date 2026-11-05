#!/usr/bin/env python
"""Reproduce the multi-condition reanalysis reported in the JPR manuscript.

Usage
-----
    python scripts/run_multicondition_reanalysis.py

    python scripts/run_multicondition_reanalysis.py \
        --archive data/external/pxd070124_parsed.json \
        --features data/processed/bsa_structure_features_union44.csv \
        --output results/multicondition

Every table in the manuscript comes out of this one command.  The run takes
a few minutes: the per-load permutation controls and the random-subset
reference together fit roughly seven thousand forests.

Reported-value check
--------------------
Pass ``--check`` to compare the regenerated numbers against the values
printed in the manuscript.  The script exits non-zero if any of them moved,
which makes it usable as a regression test on the paper itself.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(REPO_ROOT / "src"))

from cml_repairnet.revision.multicondition import (  # noqa: E402
    CONDITIONS,
    SUBSTRATE_NG,
    MulticonditionConfig,
    run,
)

# Values as printed in the manuscript, keyed by condition.  Each is checked
# at the precision the manuscript prints it to.
REPORTED_UNWEIGHTED = {
    "0p25": {"mae": 0.274, "spearman": 0.414, "p_perm": 0.010},
    "0p5": {"mae": 0.187, "spearman": 0.616, "p_perm": 0.005},
    "1p0": {"mae": 0.221, "spearman": 0.399, "p_perm": 0.005},
    "2p0": {"mae": 0.188, "spearman": 0.581, "p_perm": 0.005},
}
REPORTED_WEIGHTED = {
    "0p25": {"mae": 0.304, "spearman": 0.141, "p_perm": 0.109},
    "0p5": {"mae": 0.184, "spearman": 0.681, "p_perm": 0.005},
    "1p0": {"mae": 0.233, "spearman": 0.298, "p_perm": 0.010},
    "2p0": {"mae": 0.209, "spearman": 0.574, "p_perm": 0.005},
}
REPORTED_ESS = {"0p25": 5.0, "0p5": 9.7, "1p0": 11.8, "2p0": 11.4}


def _print_table(title: str, frame, columns: list[str]) -> None:
    print()
    print("=" * 78)
    print(title)
    print("=" * 78)
    print(frame[columns].to_string(index=False))


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--archive",
        default="data/external/pxd070124_parsed.json",
        help="parsed PXD070124 quantification tables",
    )
    parser.add_argument(
        "--features",
        default="data/processed/bsa_structure_features_union44.csv",
        help="structural descriptors for the union of confident sites",
    )
    parser.add_argument(
        "--output", default="results/multicondition", help="directory for result tables"
    )
    parser.add_argument(
        "--permutations",
        type=int,
        default=None,
        help="override the number of site-structure permutations",
    )
    parser.add_argument(
        "--check",
        action="store_true",
        help="compare regenerated numbers against the manuscript and exit non-zero on drift",
    )
    args = parser.parse_args(argv)

    config = MulticonditionConfig()
    if args.permutations is not None:
        config = MulticonditionConfig(n_permutations=args.permutations)

    outputs = run(args.archive, args.features, args.output, config)

    _print_table(
        "Data coverage",
        outputs.coverage,
        ["substrate_ng", "quantified", "confident", "with_structure"],
    )
    _print_table(
        "Label diagnostics (manuscript Table 2)",
        outputs.diagnostics,
        ["substrate_ng", "n", "w_min", "w_max", "kish_ess", "ess_over_n", "truncated"],
    )
    _print_table(
        "Replication across loads (manuscript Table 3)",
        outputs.replication,
        ["substrate_ng", "weight_mode", "n", "mae", "baseline", "spearman", "dMAE", "p_perm"],
    )
    _print_table(
        "Weighting sensitivity",
        outputs.weighting,
        ["substrate_ng", "weight_mode", "mae", "spearman", "p_perm", "ess"],
    )
    _print_table(
        "Matched common core (manuscript Table 4)",
        outputs.core,
        ["substrate_ng", "n", "mae", "baseline", "spearman", "p_perm"],
    )
    _print_table(
        "Response dispersion",
        outputs.dispersion,
        ["substrate_ng", "n", "mean", "sd", "iqr", "at_or_below_0.10", "at_or_above_0.90"],
    )

    summary = outputs.subset_summary
    print()
    print("=" * 78)
    print("Random-subset reference")
    print("=" * 78)
    print(
        f"  {summary['n_draws']} draws of {summary['subset_size']} sites from the "
        f"{summary['reference_substrate_ng']} ng load"
    )
    print(
        f"  rho: mean {summary['spearman_mean']:.3f}  sd {summary['spearman_sd']:.3f}  "
        f"median {summary['spearman_median']:.3f}"
    )
    print(f"  fraction with rho > 0.3: {summary['fraction_above_0.3']:.0%}")
    print(
        f"  matched core: rho {summary['core_spearman']:.3f} "
        f"(percentile {summary['core_percentile']:.0%}), "
        f"improvement {summary['core_improvement']:.3f} "
        f"against a mean of {summary['improvement_mean']:.3f}"
    )

    confound = outputs.confound["structure_vs_measurability"]
    print()
    print("  measurability confound: smallest p across the nine descriptors is "
          f"{confound['p_n_cond_pass'].min():.3f} "
          f"(Bonferroni threshold {0.05 / len(confound):.4f})")

    print()
    print(f"Wrote {len(list(Path(args.output).glob('*')))} files to {args.output}")

    if not args.check:
        return 0

    print()
    print("=" * 78)
    print("Checking regenerated values against the manuscript")
    print("=" * 78)
    failures: list[str] = []

    def compare(label: str, got: float, want: float, decimals: int = 3) -> None:
        """Compare at the precision the manuscript actually reports.

        The effective sample sizes are printed to one decimal, so 9.7459 is
        exactly what "9.7" means; comparing it to four decimals would report
        a drift that is only rounding.
        """
        shown = round(got, decimals)
        drifted = abs(shown - want) > 10**-decimals / 2
        if drifted:
            failures.append(
                f"{label}: got {shown:.{decimals}f}, manuscript says {want:.{decimals}f}"
            )
        print(
            f"  {label:<34} got {shown:>8.{decimals}f}  paper {want:>8.{decimals}f}"
            f"  {'DRIFT' if drifted else 'ok'}"
        )

    for mode, reported in (("none", REPORTED_UNWEIGHTED), ("inv_var", REPORTED_WEIGHTED)):
        rows = outputs.replication[outputs.replication["weight_mode"] == mode]
        for cond in CONDITIONS:
            row = rows[rows["cond"] == cond].iloc[0]
            tag = f"{SUBSTRATE_NG[cond]}ng {mode}"
            for key in ("mae", "spearman", "p_perm"):
                compare(f"{tag} {key}", float(row[key]), reported[cond][key])

    for cond in CONDITIONS:
        row = outputs.diagnostics[outputs.diagnostics["cond"] == cond].iloc[0]
        compare(
            f"{SUBSTRATE_NG[cond]}ng Kish ESS",
            float(row["kish_ess"]),
            REPORTED_ESS[cond],
            decimals=1,
        )

    print()
    if failures:
        print(f"{len(failures)} value(s) no longer match the manuscript:")
        for failure in failures:
            print(f"  - {failure}")
        return 1
    print("All reported values reproduced.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
