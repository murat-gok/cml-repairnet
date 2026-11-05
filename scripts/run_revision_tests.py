from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import pandas as pd

from cml_repairnet.revision.feature_ablation import (
    BFACTOR_COLUMNS,
    build_feature_families,
    load_bsa_tables,
    run_feature_family_ablation,
    run_structure_permutation_control,
    write_revision_report,
)

DEFAULT_LABELS = Path("data/processed/bsa_labels.csv")
DEFAULT_SEQUENCE = Path("data/processed/bsa_sequence_features.csv")
DEFAULT_STRUCTURE = Path("data/processed/bsa_structure_features.csv")
DEFAULT_OUTPUT = Path("results/revision_tests")


def _add_common_paths(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--labels", type=Path, default=DEFAULT_LABELS)
    parser.add_argument("--sequence", type=Path, default=DEFAULT_SEQUENCE)
    parser.add_argument("--structure", type=Path, default=DEFAULT_STRUCTURE)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument(
        "--target",
        choices=("repair_published", "repair_from_means_clipped"),
        default="repair_published",
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run reviewer-driven CML-RepairNet revision tests"
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    doctor = subparsers.add_parser("doctor", help="Validate inputs and feature-family isolation")
    _add_common_paths(doctor)

    ablation = subparsers.add_parser("ablation", help="Run T1 feature-family ablation")
    _add_common_paths(ablation)
    ablation.add_argument("--mode", choices=("quick", "full"), default="quick")
    ablation.add_argument("--partitions", type=int, default=None)
    ablation.add_argument("--groups", type=int, nargs="+", default=[4])
    ablation.add_argument(
        "--models",
        nargs="+",
        choices=("extra_trees", "ridge", "lightgbm"),
        default=["extra_trees", "ridge", "lightgbm"],
    )
    ablation.add_argument("--n-jobs", type=int, default=1)
    ablation.add_argument("--bootstrap-iterations", type=int, default=10_000)

    permutation = subparsers.add_parser(
        "permutation", help="Run T2 permuted-structure negative control"
    )
    _add_common_paths(permutation)
    permutation.add_argument("--iterations", type=int, default=200)
    permutation.add_argument(
        "--feature-set",
        choices=("structure_transfer_safe", "combined_transfer_safe", "combined_transfer_safe_plus_proxy"),
        default="combined_transfer_safe",
    )
    permutation.add_argument(
        "--model", choices=("extra_trees", "ridge", "lightgbm"), default="extra_trees"
    )
    permutation.add_argument("--groups", type=int, default=4)
    permutation.add_argument(
        "--partition-method", choices=("ward", "kmeans_random"), default="ward"
    )
    permutation.add_argument("--partition-seed", type=int, default=20260729)
    permutation.add_argument("--n-jobs", type=int, default=1)

    all_cmd = subparsers.add_parser("all", help="Run doctor, T1 and T2")
    _add_common_paths(all_cmd)
    all_cmd.add_argument("--mode", choices=("quick", "full"), default="quick")
    all_cmd.add_argument("--partitions", type=int, default=None)
    all_cmd.add_argument("--groups", type=int, nargs="+", default=[4])
    all_cmd.add_argument(
        "--models",
        nargs="+",
        choices=("extra_trees", "ridge", "lightgbm"),
        default=["extra_trees", "ridge", "lightgbm"],
    )
    all_cmd.add_argument("--permutations", type=int, default=None)
    all_cmd.add_argument("--permutation-model", choices=("extra_trees", "ridge", "lightgbm"), default="extra_trees")
    all_cmd.add_argument("--n-jobs", type=int, default=1)
    all_cmd.add_argument("--bootstrap-iterations", type=int, default=10_000)
    return parser.parse_args()


def _load(args: argparse.Namespace) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    return load_bsa_tables(args.labels, args.sequence, args.structure)


def run_doctor(args: argparse.Namespace) -> dict[str, object]:
    labels, sequence, structure = _load(args)
    families = build_feature_families(labels, sequence, structure)
    if args.target not in labels.columns:
        raise ValueError(f"Target {args.target!r} is absent from labels")
    report = {
        "status": "PASS",
        "n_labels": len(labels),
        "n_sequence_rows": len(sequence),
        "n_structure_rows": len(structure),
        "target": args.target,
        "target_missing": int(labels[args.target].isna().sum()),
        "b_factor_columns_present": [column for column in BFACTOR_COLUMNS if column in structure.columns],
        "feature_sets": {
            name: {
                "n_features": len(family.feature_columns),
                "n_structure_features": len(family.structure_columns),
                "contains_bfactor": bool(set(BFACTOR_COLUMNS).intersection(family.feature_columns)),
            }
            for name, family in families.items()
        },
    }
    args.output.mkdir(parents=True, exist_ok=True)
    (args.output / "revision_doctor.json").write_text(
        json.dumps(report, indent=2, sort_keys=True), encoding="utf-8"
    )
    print(json.dumps(report, indent=2, sort_keys=True))
    if len(labels) != 33:
        print(
            f"WARNING: expected the frozen BSA benchmark to contain 33 sites, found {len(labels)}.",
            file=sys.stderr,
        )
    return report


def _partition_count(mode: str, explicit: int | None) -> int:
    if explicit is not None:
        if explicit < 1:
            raise ValueError("--partitions must be at least 1")
        return explicit
    return 5 if mode == "quick" else 50


def run_ablation(args: argparse.Namespace) -> dict[str, pd.DataFrame]:
    labels, sequence, structure = _load(args)
    outputs = run_feature_family_ablation(
        labels,
        sequence,
        structure,
        args.output,
        target=args.target,
        model_names=args.models,
        group_counts=args.groups,
        random_partition_count=_partition_count(args.mode, args.partitions),
        n_jobs=args.n_jobs,
        bootstrap_iterations=args.bootstrap_iterations,
    )
    print("\nT1 summary")
    print(outputs["summary"].to_string(index=False))
    print("\nT1 paired contrasts")
    print(outputs["pairwise_summary"].to_string(index=False))
    return outputs


def run_permutation(args: argparse.Namespace, *, iterations: int | None = None) -> dict[str, pd.DataFrame]:
    labels, sequence, structure = _load(args)
    n_iterations = iterations if iterations is not None else args.iterations
    outputs = run_structure_permutation_control(
        labels,
        sequence,
        structure,
        args.output,
        target=args.target,
        feature_set=getattr(args, "feature_set", "combined_transfer_safe"),
        model_name=getattr(args, "model", getattr(args, "permutation_model", "extra_trees")),
        iterations=n_iterations,
        n_groups=args.groups if isinstance(args.groups, int) else int(args.groups[0]),
        partition_method=getattr(args, "partition_method", "ward"),
        partition_seed=getattr(args, "partition_seed", 20260729),
        n_jobs=args.n_jobs,
    )
    print("\nT2 summary")
    print(outputs["summary"].to_string(index=False))
    return outputs


def main() -> None:
    args = parse_args()
    if args.command == "doctor":
        run_doctor(args)
        return
    if args.command == "ablation":
        outputs = run_ablation(args)
        write_revision_report(
            args.output,
            ablation_summary=outputs["summary"],
            pairwise_summary=outputs["pairwise_summary"],
        )
        return
    if args.command == "permutation":
        outputs = run_permutation(args)
        write_revision_report(args.output, permutation_summary=outputs["summary"])
        return
    if args.command == "all":
        run_doctor(args)
        ablation_outputs = run_ablation(args)
        permutations = args.permutations
        if permutations is None:
            permutations = 200 if args.mode == "quick" else 1000
        permutation_outputs = run_permutation(args, iterations=permutations)
        report = write_revision_report(
            args.output,
            ablation_summary=ablation_outputs["summary"],
            pairwise_summary=ablation_outputs["pairwise_summary"],
            permutation_summary=permutation_outputs["summary"],
        )
        print(f"\nReport: {report}")
        return
    raise RuntimeError(f"Unhandled command: {args.command}")


if __name__ == "__main__":
    main()
