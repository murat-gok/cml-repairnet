from __future__ import annotations

import argparse
from pathlib import Path

from cml_repairnet.revision.confirmatory import (
    run_multi_partition_permutation_control,
    run_nested_model_selection,
    run_primary_cluster_bootstrap,
    write_confirmatory_report,
)
from cml_repairnet.revision.feature_ablation import (
    load_bsa_tables,
    run_feature_family_ablation,
    write_revision_report,
)

DEFAULT_LABELS = Path("data/processed/bsa_labels.csv")
DEFAULT_SEQUENCE = Path("data/processed/bsa_sequence_features.csv")
DEFAULT_STRUCTURE = Path("data/processed/bsa_structure_features.csv")
DEFAULT_OUTPUT = Path("results/revision_confirmatory")
DEFAULT_PREDICTIONS = Path("results/revision_tests/revision_t1_ablation_oof_predictions.csv")

TARGETED_FAMILIES = [
    "structure_transfer_safe",
    "structure_transfer_safe_plus_proxy",
    "structure_full",
    "structure_full_plus_proxy",
    "structure_bfactor_only",
    "structure_proxy_only",
]


def add_data_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--labels", type=Path, default=DEFAULT_LABELS)
    parser.add_argument("--sequence", type=Path, default=DEFAULT_SEQUENCE)
    parser.add_argument("--structure", type=Path, default=DEFAULT_STRUCTURE)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument(
        "--target",
        choices=("repair_published", "repair_from_means_clipped"),
        default="repair_published",
    )


def add_partition_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--groups", type=int, nargs="+", default=[4])
    parser.add_argument("--random-partitions", type=int, default=5)
    parser.add_argument("--start-seed", type=int, default=20260729)
    parser.add_argument("--n-jobs", type=int, default=1)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run confirmatory reviewer tests T1b, T2-full, T3 and T4"
    )
    sub = parser.add_subparsers(dest="command", required=True)

    ablation = sub.add_parser("ablation", help="Run targeted T1b feature ablation")
    add_data_args(ablation)
    add_partition_args(ablation)
    ablation.add_argument(
        "--models", nargs="+", choices=("extra_trees", "lightgbm", "ridge"),
        default=["extra_trees", "lightgbm", "ridge"]
    )
    ablation.add_argument(
        "--families", nargs="+", default=TARGETED_FAMILIES,
        help="Feature families to evaluate; defaults to the targeted structure/B-factor/proxy panel",
    )
    ablation.add_argument("--bootstrap-iterations", type=int, default=10_000)

    permutation = sub.add_parser("permutation", help="Run T2 across all partitions")
    add_data_args(permutation)
    add_partition_args(permutation)
    permutation.add_argument("--iterations", type=int, default=2_000)
    permutation.add_argument(
        "--feature-set",
        choices=(
            "structure_transfer_safe",
            "structure_transfer_safe_plus_proxy",
            "combined_transfer_safe",
        ),
        default="structure_transfer_safe",
    )
    permutation.add_argument(
        "--model", choices=("extra_trees", "lightgbm", "ridge"), default="extra_trees"
    )
    permutation.add_argument("--checkpoint-every", type=int, default=25)
    permutation.add_argument("--no-resume", action="store_true")

    nested = sub.add_parser("nested", help="Run T3 nested spatial model selection")
    add_data_args(nested)
    add_partition_args(nested)
    nested.add_argument(
        "--feature-set",
        choices=("structure_transfer_safe", "structure_transfer_safe_plus_proxy", "structure_full"),
        default="structure_transfer_safe",
    )
    nested.add_argument(
        "--models", nargs="+", choices=("extra_trees", "lightgbm", "ridge"),
        default=["extra_trees", "lightgbm", "ridge"]
    )
    nested.add_argument(
        "--fixed-model", choices=("extra_trees", "lightgbm", "ridge"), default="extra_trees"
    )
    nested.add_argument("--bootstrap-iterations", type=int, default=10_000)

    bootstrap = sub.add_parser("bootstrap", help="Run T4 paired cluster bootstrap")
    bootstrap.add_argument("--predictions", type=Path, default=DEFAULT_PREDICTIONS)
    bootstrap.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    bootstrap.add_argument("--partition-id", default="ward_g4")
    bootstrap.add_argument("--feature-set", default="structure_transfer_safe")
    bootstrap.add_argument("--model", default="extra_trees")
    bootstrap.add_argument("--iterations", type=int, default=20_000)
    bootstrap.add_argument("--sign-flip-iterations", type=int, default=100_000)

    all_cmd = sub.add_parser("all", help="Run T1b, T2-full, T3 and T4")
    add_data_args(all_cmd)
    add_partition_args(all_cmd)
    all_cmd.add_argument(
        "--models", nargs="+", choices=("extra_trees", "lightgbm", "ridge"),
        default=["extra_trees", "lightgbm", "ridge"]
    )
    all_cmd.add_argument("--families", nargs="+", default=TARGETED_FAMILIES)
    all_cmd.add_argument("--permutations", type=int, default=2_000)
    all_cmd.add_argument("--checkpoint-every", type=int, default=25)
    all_cmd.add_argument("--no-resume", action="store_true")
    all_cmd.add_argument("--bootstrap-iterations", type=int, default=20_000)
    all_cmd.add_argument("--sign-flip-iterations", type=int, default=100_000)
    return parser.parse_args()


def load(args: argparse.Namespace):
    return load_bsa_tables(args.labels, args.sequence, args.structure)


def main() -> None:
    args = parse_args()
    if args.command == "bootstrap":
        out = run_primary_cluster_bootstrap(
            args.predictions,
            args.output,
            partition_id=args.partition_id,
            primary_feature_set=args.feature_set,
            primary_model=args.model,
            iterations=args.iterations,
            sign_flip_iterations=args.sign_flip_iterations,
        )
        report = write_confirmatory_report(
            args.output,
            t4_observed=out["observed"],
            t4_summary=out["summary"],
            t4_tests=out["tests"],
        )
        print(out["observed"].to_string(index=False))
        print(out["summary"].to_string(index=False))
        print(out["tests"].to_string(index=False))
        print(f"Report: {report}")
        return

    labels, sequence, structure = load(args)
    if args.command == "ablation":
        out = run_feature_family_ablation(
            labels, sequence, structure, args.output,
            target=args.target,
            model_names=args.models,
            family_names=args.families,
            group_counts=args.groups,
            random_partition_count=args.random_partitions,
            start_seed=args.start_seed,
            n_jobs=args.n_jobs,
            bootstrap_iterations=args.bootstrap_iterations,
        )
        report = write_revision_report(
            args.output,
            ablation_summary=out["summary"],
            pairwise_summary=out["pairwise_summary"],
        )
        print(out["summary"].to_string(index=False))
        print(out["pairwise_summary"].to_string(index=False))
        print(f"Report: {report}")
        return

    if args.command == "permutation":
        out = run_multi_partition_permutation_control(
            labels, sequence, structure, args.output,
            target=args.target,
            feature_set=args.feature_set,
            model_name=args.model,
            iterations=args.iterations,
            group_counts=args.groups,
            random_partition_count=args.random_partitions,
            start_seed=args.start_seed,
            n_jobs=args.n_jobs,
            resume=not args.no_resume,
            checkpoint_every=args.checkpoint_every,
        )
        report = write_confirmatory_report(args.output, t2_aggregate=out["aggregate_summary"])
        print(out["partition_summary"].to_string(index=False))
        print(out["aggregate_summary"].to_string(index=False))
        print(f"Report: {report}")
        return

    if args.command == "nested":
        out = run_nested_model_selection(
            labels, sequence, structure, args.output,
            target=args.target,
            feature_set=args.feature_set,
            model_names=args.models,
            fixed_model=args.fixed_model,
            group_counts=args.groups,
            random_partition_count=args.random_partitions,
            start_seed=args.start_seed,
            n_jobs=args.n_jobs,
            bootstrap_iterations=args.bootstrap_iterations,
        )
        report = write_confirmatory_report(
            args.output,
            t3_summary=out["summary"],
            t3_frequency=out["selection_frequency"],
        )
        print(out["partition_metrics"].to_string(index=False))
        print(out["summary"].to_string(index=False))
        print(out["selection_frequency"].to_string(index=False))
        print(f"Report: {report}")
        return

    if args.command == "all":
        t1 = run_feature_family_ablation(
            labels, sequence, structure, args.output,
            target=args.target,
            model_names=args.models,
            family_names=args.families,
            group_counts=args.groups,
            random_partition_count=args.random_partitions,
            start_seed=args.start_seed,
            n_jobs=args.n_jobs,
            bootstrap_iterations=args.bootstrap_iterations,
        )
        t2 = run_multi_partition_permutation_control(
            labels, sequence, structure, args.output,
            target=args.target,
            feature_set="structure_transfer_safe",
            model_name="extra_trees",
            iterations=args.permutations,
            group_counts=args.groups,
            random_partition_count=args.random_partitions,
            start_seed=args.start_seed,
            n_jobs=args.n_jobs,
            resume=not args.no_resume,
            checkpoint_every=args.checkpoint_every,
        )
        t3 = run_nested_model_selection(
            labels, sequence, structure, args.output,
            target=args.target,
            feature_set="structure_transfer_safe",
            model_names=args.models,
            fixed_model="extra_trees",
            group_counts=args.groups,
            random_partition_count=args.random_partitions,
            start_seed=args.start_seed,
            n_jobs=args.n_jobs,
            bootstrap_iterations=args.bootstrap_iterations,
        )
        predictions = args.output / "revision_t1_ablation_oof_predictions.csv"
        t4 = run_primary_cluster_bootstrap(
            predictions,
            args.output,
            iterations=args.bootstrap_iterations,
            sign_flip_iterations=args.sign_flip_iterations,
        )
        write_revision_report(
            args.output,
            ablation_summary=t1["summary"],
            pairwise_summary=t1["pairwise_summary"],
        )
        report = write_confirmatory_report(
            args.output,
            t2_aggregate=t2["aggregate_summary"],
            t3_summary=t3["summary"],
            t3_frequency=t3["selection_frequency"],
            t4_observed=t4["observed"],
            t4_summary=t4["summary"],
            t4_tests=t4["tests"],
        )
        print("\nT1b summary")
        print(t1["summary"].to_string(index=False))
        print("\nT2 aggregate")
        print(t2["aggregate_summary"].to_string(index=False))
        print("\nT3 summary")
        print(t3["summary"].to_string(index=False))
        print("\nT4 observed")
        print(t4["observed"].to_string(index=False))
        print(f"\nReport: {report}")
        return

    raise RuntimeError(args.command)


if __name__ == "__main__":
    main()
