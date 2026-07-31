from __future__ import annotations

import importlib.util
import platform
import sys
from pathlib import Path
from typing import Any


def collect_doctor_report(project_root: Path = Path(".")) -> dict[str, Any]:
    required_packages = [
        "numpy",
        "pandas",
        "sklearn",
        "Bio",
        "openpyxl",
        "yaml",
        "requests",
    ]
    package_status = {
        package: importlib.util.find_spec(package) is not None for package in required_packages
    }
    raw_files = {
        "article_source_data": project_root / "data/raw/article_source_data.xlsx",
        "bsa_fasta": project_root / "data/raw/P02769.fasta",
        "bsa_pdb": project_root / "data/raw/4F5S.pdb",
    }
    python_ok = (3, 10) <= sys.version_info[:2] < (3, 13)
    return {
        "python": platform.python_version(),
        "python_supported": python_ok,
        "platform": platform.platform(),
        "packages": package_status,
        "raw_files": {name: path.exists() for name, path in raw_files.items()},
        "ready_for_stage1": python_ok
        and all(package_status.values())
        and all(path.exists() for path in raw_files.values()),
    }


def format_doctor_report(report: dict[str, Any]) -> str:
    lines = [
        "CML-RepairNet environment check",
        "=" * 32,
        f"Python: {report['python']} "
        f"({'OK' if report['python_supported'] else 'SUPPORTED RANGE: 3.10–3.12'})",
        f"Platform: {report['platform']}",
        "",
        "Packages:",
    ]
    lines.extend(
        f"  [{'OK' if installed else 'MISSING'}] {name}"
        for name, installed in report["packages"].items()
    )
    lines.append("")
    lines.append("Raw data:")
    lines.extend(
        f"  [{'OK' if exists else 'MISSING'}] {name}"
        for name, exists in report["raw_files"].items()
    )
    lines.extend(
        [
            "",
            "READY FOR STAGE 1" if report["ready_for_stage1"] else "NOT READY FOR STAGE 1",
        ]
    )
    if not all(report["raw_files"].values()):
        lines.append("Next command: cml-repairnet download")
    return "\n".join(lines)
