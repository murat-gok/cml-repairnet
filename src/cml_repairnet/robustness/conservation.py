from __future__ import annotations

import numpy as np
import pandas as pd
from Bio.Align import PairwiseAligner


def _map_position(alignment, position_zero_based: int, source_axis: int) -> int | None:
    source_blocks = alignment.aligned[source_axis]
    target_blocks = alignment.aligned[1 - source_axis]
    for (source_start, source_end), (target_start, target_end) in zip(source_blocks, target_blocks):
        if source_start <= position_zero_based < source_end:
            return int(target_start + (position_zero_based - source_start))
    return None


def _mapped_identity(alignment, sequence_a: str, sequence_b: str) -> tuple[int, int, float]:
    identical = 0
    aligned = 0
    for (a_start, a_end), (b_start, b_end) in zip(alignment.aligned[0], alignment.aligned[1]):
        length = int(a_end - a_start)
        aligned += length
        identical += sum(
            sequence_a[int(a_start) + offset] == sequence_b[int(b_start) + offset]
            for offset in range(length)
        )
    return identical, aligned, identical / aligned if aligned else np.nan


def run_residue_conservation(pairs: pd.DataFrame, config: dict) -> pd.DataFrame:
    required = {
        "source_gene",
        "target_gene",
        "uniprot_accession_human",
        "uniprot_accession_mouse",
        "site_uniprot_human",
        "site_uniprot_mouse",
        "sequence_human",
        "sequence_mouse",
        "repairability_extra_trees_human",
        "repairability_extra_trees_mouse",
    }
    missing = sorted(required - set(pairs.columns))
    if missing:
        raise ValueError(f"Cross-species table missing columns: {missing}")

    aligner = PairwiseAligner(mode="global")
    aligner.match_score = float(config.get("match_score", 2.0))
    aligner.mismatch_score = float(config.get("mismatch_score", -1.0))
    aligner.open_gap_score = float(config.get("open_gap_score", -10.0))
    aligner.extend_gap_score = float(config.get("extend_gap_score", -0.5))
    radius = int(config.get("local_window_radius", 7))
    rows = []

    for row in pairs.itertuples(index=False):
        human_sequence = str(row.sequence_human).strip().upper()
        mouse_sequence = str(row.sequence_mouse).strip().upper()
        human_site = int(row.site_uniprot_human)
        mouse_site = int(row.site_uniprot_mouse)
        if not (1 <= human_site <= len(human_sequence) and 1 <= mouse_site <= len(mouse_sequence)):
            raise ValueError(f"Site outside sequence for {row.source_gene}/{row.target_gene}")
        alignment = aligner.align(human_sequence, mouse_sequence)[0]
        mapped_mouse_zero = _map_position(alignment, human_site - 1, 0)
        mapped_human_zero = _map_position(alignment, mouse_site - 1, 1)
        mapped_mouse_site = None if mapped_mouse_zero is None else mapped_mouse_zero + 1
        mapped_human_site = None if mapped_human_zero is None else mapped_human_zero + 1
        reciprocal = mapped_mouse_site == mouse_site and mapped_human_site == human_site
        human_residue = human_sequence[human_site - 1]
        mouse_residue = mouse_sequence[mouse_site - 1]

        local_total = 0
        local_identical = 0
        start = max(0, human_site - 1 - radius)
        end = min(len(human_sequence), human_site + radius)
        for human_position in range(start, end):
            mouse_position = _map_position(alignment, human_position, 0)
            if mouse_position is None or not (0 <= mouse_position < len(mouse_sequence)):
                continue
            local_total += 1
            local_identical += human_sequence[human_position] == mouse_sequence[mouse_position]
        identical, aligned, global_identity = _mapped_identity(
            alignment, human_sequence, mouse_sequence
        )

        rows.append(
            {
                "human_gene": row.source_gene,
                "mouse_gene": row.target_gene,
                "human_accession": row.uniprot_accession_human,
                "mouse_accession": row.uniprot_accession_mouse,
                "human_site": human_site,
                "mouse_site": mouse_site,
                "human_residue": human_residue,
                "mouse_residue": mouse_residue,
                "mapped_mouse_site_from_human": mapped_mouse_site,
                "mapped_human_site_from_mouse": mapped_human_site,
                "reciprocal_site_alignment": reciprocal,
                "both_declared_sites_are_lysine": human_residue == "K" and mouse_residue == "K",
                "residue_level_conserved_lysine": reciprocal
                and human_residue == "K"
                and mouse_residue == "K",
                "same_uniprot_site_number": human_site == mouse_site,
                "site_number_difference": abs(human_site - mouse_site),
                "local_window_radius": radius,
                "local_aligned_positions": local_total,
                "local_sequence_identity": local_identical / local_total if local_total else np.nan,
                "global_aligned_positions": aligned,
                "global_sequence_identity": global_identity,
                "alignment_score": float(alignment.score),
                "repairability_human": float(row.repairability_extra_trees_human),
                "repairability_mouse": float(row.repairability_extra_trees_mouse),
                "repairability_abs_difference": abs(
                    float(row.repairability_extra_trees_human)
                    - float(row.repairability_extra_trees_mouse)
                ),
            }
        )
    return (
        pd.DataFrame(rows)
        .sort_values(
            [
                "residue_level_conserved_lysine",
                "local_sequence_identity",
                "repairability_abs_difference",
            ],
            ascending=[False, False, True],
        )
        .reset_index(drop=True)
    )
