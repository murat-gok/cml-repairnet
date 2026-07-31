from __future__ import annotations

import argparse
from pathlib import Path

from cml_repairnet.data.download import download_public_bundle


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--raw-dir", default="data/raw")
    parser.add_argument("--include-pride-processed", action="store_true")
    args = parser.parse_args()
    records = download_public_bundle(Path(args.raw_dir), args.include_pride_processed)
    for record in records:
        print(f"Downloaded {record['path']} ({record['bytes']} bytes)")


if __name__ == "__main__":
    main()
