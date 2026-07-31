from __future__ import annotations

import argparse
from pathlib import Path

from cml_repairnet.data.download import download_public_bundle
from cml_repairnet.data.prepare_labels import build_bsa_labels
from cml_repairnet.doctor import collect_doctor_report, format_doctor_report
from cml_repairnet.pipeline import run_stage1
from cml_repairnet.utils.config import load_config


def main() -> None:
    parser = argparse.ArgumentParser(prog="cml-repairnet")
    parser.add_argument("--config", default="configs/base.yaml")
    subparsers = parser.add_subparsers(dest="command", required=True)

    subparsers.add_parser("doctor", help="Check Python, dependencies and raw data.")

    download_parser = subparsers.add_parser("download", help="Download official public inputs.")
    download_parser.add_argument("--raw-dir", default="data/raw")
    download_parser.add_argument("--include-pride-processed", action="store_true")

    labels_parser = subparsers.add_parser("prepare-labels", help="Rebuild the 33-site labels.")
    labels_parser.add_argument("--source", default="data/raw/article_source_data.xlsx")
    labels_parser.add_argument("--fasta", default="data/raw/P02769.fasta")
    labels_parser.add_argument("--output", default="data/processed/bsa_labels.csv")

    stage1_parser = subparsers.add_parser(
        "stage1", help="Build and audit labels, sequence features and structure features."
    )
    stage1_parser.add_argument(
        "--no-download", action="store_true", help="Require raw files to exist locally."
    )

    args = parser.parse_args()
    if args.command == "doctor":
        print(format_doctor_report(collect_doctor_report()))
        return
    if args.command == "download":
        records = download_public_bundle(Path(args.raw_dir), args.include_pride_processed)
        for record in records:
            print(f"Downloaded {record['path']} ({record['bytes']} bytes)")
        return
    if args.command == "prepare-labels":
        frame = build_bsa_labels(Path(args.source), Path(args.output), Path(args.fasta))
        print(f"Prepared {len(frame)} labelled sites: {args.output}")
        return
    if args.command == "stage1":
        report = run_stage1(load_config(args.config), download=not args.no_download)
        print("Stage 1 completed and audited.")
        print(f"Sites: {report['n_sites']}")
        print(
            "Published/reconstructed target MAE: "
            f"{report['published_vs_reconstructed_mae']:.4f}"
        )
        print("Report: results/stage1_data_audit.md")
        return
    raise RuntimeError(f"Unhandled command: {args.command}")


if __name__ == "__main__":
    main()
