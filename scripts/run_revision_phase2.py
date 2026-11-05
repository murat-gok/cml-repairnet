from __future__ import annotations

import argparse
import json
from pathlib import Path

from cml_repairnet.revision.phase2 import (
    discover_paths,
    run_phase2_doctor,
    run_t5_conformal_shift_sensitivity,
    run_t6_boundary_shrinkage,
    run_t7_compartment_layering,
    run_t8_background_and_network_sensitivity,
    write_phase2_report,
)

DEFAULT_OUTPUT = Path("results/revision_phase2")


def parser() -> argparse.ArgumentParser:
    main = argparse.ArgumentParser(
        description="Run CML-RepairNet Revision Phase II analyses T5–T8"
    )
    main.add_argument("--root", type=Path, default=Path("."))
    main.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    sub = main.add_subparsers(dest="command", required=True)

    sub.add_parser("doctor", help="Validate required and optional Phase II inputs")

    t5 = sub.add_parser("conformal", help="T5 covariate-shift interval sensitivity")
    t5.add_argument("--alpha", type=float, default=0.10)
    t5.add_argument("--partition-id", default="ward_g4")
    t5.add_argument("--density-ratio-clip", type=float, nargs=2, default=(0.05, 20.0))
    t5.add_argument("--ood-inflation-max", type=float, default=2.0)

    t6 = sub.add_parser("boundary", help="T6 boundary-shrinkage diagnostics")
    t6.add_argument("--low-threshold", type=float, default=0.10)
    t6.add_argument("--high-threshold", type=float, default=0.90)
    t6.add_argument("--bootstrap-iterations", type=int, default=20_000)

    sub.add_parser("compartments", help="T7 compartment/accessibility layering")

    t8 = sub.add_parser("enrichment", help="T8 multi-background enrichment and STRING sensitivity")
    t8.add_argument("--controls-per-query", type=int, default=5)
    t8.add_argument("--offline", action="store_true")
    t8.add_argument("--timeout", type=float, default=90.0)
    t8.add_argument("--retries", type=int, default=2)

    all_cmd = sub.add_parser("all", help="Run T5–T8 in sequence")
    all_cmd.add_argument("--alpha", type=float, default=0.10)
    all_cmd.add_argument("--bootstrap-iterations", type=int, default=20_000)
    all_cmd.add_argument("--controls-per-query", type=int, default=5)
    all_cmd.add_argument("--offline-enrichment", action="store_true")
    return main


def main() -> None:
    args = parser().parse_args()
    root = args.root.resolve()
    output = (root / args.output).resolve() if not args.output.is_absolute() else args.output
    paths = discover_paths(root)

    if args.command == "doctor":
        report = run_phase2_doctor(root, output)
        print(json.dumps(report, indent=2))
        return

    if args.command == "conformal":
        result = run_t5_conformal_shift_sensitivity(
            paths.labels,
            paths.structure,
            paths.predictions,
            paths.atlas,
            output,
            alpha=args.alpha,
            partition_id=args.partition_id,
            density_ratio_clip=tuple(args.density_ratio_clip),
            ood_inflation_max=args.ood_inflation_max,
        )
        report = write_phase2_report(output)
        print(result["coverage"].to_string(index=False))
        print(result["widths"].to_string(index=False))
        print(f"Report: {report}")
        return

    if args.command == "boundary":
        result = run_t6_boundary_shrinkage(
            paths.predictions,
            output,
            low_threshold=args.low_threshold,
            high_threshold=args.high_threshold,
            bootstrap_iterations=args.bootstrap_iterations,
        )
        report = write_phase2_report(output)
        print(result["summary"].to_string(index=False))
        print(result["strata"].to_string(index=False))
        print(f"Report: {report}")
        return

    if args.command == "compartments":
        result = run_t7_compartment_layering(
            paths.atlas,
            output,
            tiers_path=paths.tiers,
        )
        report = write_phase2_report(output)
        print(result["summary"].to_string(index=False))
        print(result["actionability"].to_string(index=False))
        print(f"Report: {report}")
        return

    if args.command == "enrichment":
        result = run_t8_background_and_network_sensitivity(
            paths.atlas,
            output,
            portfolio_path=paths.portfolio,
            pareto_path=paths.pareto,
            existing_enrichment_path=paths.enrichment,
            string_edges_path=paths.string_edges,
            controls_per_query=args.controls_per_query,
            run_online=not args.offline,
            timeout=args.timeout,
            retries=args.retries,
        )
        report = write_phase2_report(output)
        compact = (
            result["manifest"]
            .groupby(["species", "query_name", "background_name", "api_status"], as_index=False)
            .agg(query_gene_count=("query_gene_count", "first"), background_gene_count=("background_gene_count", "first"))
        )
        print(compact.to_string(index=False))
        print(result["term_robustness"].head(50).to_string(index=False))
        print(result["network_summary"].to_string(index=False))
        print(f"Report: {report}")
        return

    if args.command == "all":
        doctor = run_phase2_doctor(root, output)
        if doctor["status"] != "PASS":
            raise SystemExit("Phase II doctor failed; inspect results/revision_phase2/phase2_doctor.json")
        t5 = run_t5_conformal_shift_sensitivity(
            paths.labels,
            paths.structure,
            paths.predictions,
            paths.atlas,
            output,
            alpha=args.alpha,
        )
        t6 = run_t6_boundary_shrinkage(
            paths.predictions,
            output,
            bootstrap_iterations=args.bootstrap_iterations,
        )
        t7 = run_t7_compartment_layering(paths.atlas, output, tiers_path=paths.tiers)
        t8 = run_t8_background_and_network_sensitivity(
            paths.atlas,
            output,
            portfolio_path=paths.portfolio,
            pareto_path=paths.pareto,
            existing_enrichment_path=paths.enrichment,
            string_edges_path=paths.string_edges,
            controls_per_query=args.controls_per_query,
            run_online=not args.offline_enrichment,
        )
        report = write_phase2_report(output)
        print(t5["coverage"].to_string(index=False))
        print(t6["summary"].to_string(index=False))
        print(t7["actionability"].to_string(index=False))
        print(t8["term_robustness"].head(50).to_string(index=False))
        print(f"Report: {report}")


if __name__ == "__main__":
    main()
