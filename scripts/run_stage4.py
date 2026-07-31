from __future__ import annotations

import argparse
import json
from pathlib import Path

from cml_repairnet.atlas.pipeline import (
    annotate_stage4,
    audit_stage4,
    doctor_stage4,
    download_stage4,
    prepare_stage4,
    prioritize_stage4,
    run_stage4_all,
    score_stage4,
    structures_stage4,
)


def project_root_from_script() -> Path:
    return Path(__file__).resolve().parents[1]


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="CML-RepairNet Stage 4 human/mouse CML candidate atlas"
    )
    parser.add_argument(
        "command",
        choices=(
            "doctor",
            "download",
            "prepare",
            "annotate",
            "structures",
            "score",
            "prioritize",
            "audit",
            "all",
        ),
    )
    parser.add_argument("--project-root", type=Path, default=project_root_from_script())
    parser.add_argument("--mode", choices=("quick", "full"), default="quick")
    parser.add_argument("--candidate-csv", type=Path)
    parser.add_argument("--workers", type=int, default=2)
    parser.add_argument("--overwrite-downloads", action="store_true")
    parser.add_argument("--overwrite-structures", action="store_true")
    parser.add_argument("--refit-model", action="store_true")
    parser.add_argument("--min-site-plddt", type=float, default=70.0)
    parser.add_argument("--top-n", type=int, default=100)
    parser.add_argument("--max-sites-per-protein", type=int, default=2)
    return parser


def main() -> None:
    args = build_parser().parse_args()
    root = args.project_root.resolve()
    command = args.command
    if command == "doctor":
        result = doctor_stage4(root)
    elif command == "download":
        result = download_stage4(root, overwrite=args.overwrite_downloads)
    elif command == "prepare":
        result = prepare_stage4(root, candidate_csv=args.candidate_csv, mode=args.mode)
        print(result.head(10).to_string(index=False))
        print(f"\nPrepared candidates: {len(result)}")
        return
    elif command == "annotate":
        eligible, excluded = annotate_stage4(root)
        result = {"eligible": len(eligible), "excluded": len(excluded)}
    elif command == "structures":
        eligible, excluded = structures_stage4(
            root,
            workers=args.workers,
            overwrite=args.overwrite_structures,
            min_site_plddt=args.min_site_plddt,
        )
        result = {"eligible": len(eligible), "excluded": len(excluded)}
    elif command == "score":
        frame = score_stage4(root, refit_model=args.refit_model)
        result = {
            "scored": len(frame),
            "in_domain": int((frame["applicability_status"] == "in_domain").sum()),
        }
    elif command == "prioritize":
        frames = prioritize_stage4(
            root,
            top_n=args.top_n,
            max_sites_per_protein=args.max_sites_per_protein,
        )
        result = {name: len(frame) for name, frame in frames.items()}
    elif command == "audit":
        result = audit_stage4(root)
    elif command == "all":
        result = run_stage4_all(
            root,
            mode=args.mode,
            candidate_csv=args.candidate_csv,
            workers=args.workers,
            overwrite_downloads=args.overwrite_downloads,
            overwrite_structures=args.overwrite_structures,
            min_site_plddt=args.min_site_plddt,
            top_n=args.top_n,
        )
    else:  # pragma: no cover
        raise AssertionError(command)
    print(json.dumps(result, indent=2, default=str))


if __name__ == "__main__":
    main()
