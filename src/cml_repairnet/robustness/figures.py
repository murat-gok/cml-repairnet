from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


def _save(fig, path: Path) -> None:
    fig.tight_layout()
    fig.savefig(path.with_suffix(".png"), dpi=220, bbox_inches="tight")
    fig.savefig(path.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def make_figures(
    stability: pd.DataFrame,
    conservation: pd.DataFrame,
    enrichment: pd.DataFrame,
    network_thresholds: pd.DataFrame,
    tiers: pd.DataFrame,
    output_dir: Path,
) -> int:
    output_dir.mkdir(parents=True, exist_ok=True)
    count = 0

    fig, ax = plt.subplots(figsize=(8, 5))
    portfolio = stability.loc[stability["baseline_portfolio_member"] == 1]
    ax.scatter(portfolio["baseline_rank"], portfolio["portfolio_selection_frequency"], alpha=0.75)
    ax.set_xlabel("Stage 4 baseline rank (in-domain universe)")
    ax.set_ylabel("Diversity-constrained portfolio selection frequency")
    ax.set_title("Stage 6 exact Stage 4 portfolio-rule stability")
    _save(fig, output_dir / "stage6_rank_stability")
    count += 1

    if not conservation.empty:
        fig, ax = plt.subplots(figsize=(6, 6))
        ax.scatter(
            conservation["repairability_human"], conservation["repairability_mouse"], alpha=0.8
        )
        lower = min(
            conservation["repairability_human"].min(), conservation["repairability_mouse"].min()
        )
        upper = max(
            conservation["repairability_human"].max(), conservation["repairability_mouse"].max()
        )
        ax.plot([lower, upper], [lower, upper], linestyle="--")
        for row in conservation.itertuples(index=False):
            ax.annotate(
                row.human_gene, (row.repairability_human, row.repairability_mouse), fontsize=7
            )
        ax.set_xlabel("Human Extra Trees repairability")
        ax.set_ylabel("Mouse Extra Trees repairability")
        ax.set_title("Residue-level ortholog candidate concordance")
        _save(fig, output_dir / "stage6_residue_conservation")
        count += 1

    if not enrichment.empty:
        top = enrichment.sort_values("adjusted_p_value_stage5").head(25).copy()
        fig, ax = plt.subplots(figsize=(9, 7))
        labels = (top["species"] + " | " + top["source"] + " | " + top["term_name"]).str.slice(
            0, 80
        )
        ax.barh(
            range(len(top)),
            -top["adjusted_p_value_stage5"]
            .clip(lower=1e-300)
            .map(lambda value: __import__("math").log10(value)),
        )
        ax.set_yticks(range(len(top)), labels)
        ax.invert_yaxis()
        ax.set_xlabel("-log10(Stage 5 adjusted p-value)")
        ax.set_title("Significant terms assessed for driver fragility")
        _save(fig, output_dir / "stage6_enrichment_sensitivity")
        count += 1

    if not network_thresholds.empty:
        fig, ax = plt.subplots(figsize=(8, 5))
        for (species, query), group in network_thresholds.groupby(["species", "query_name"]):
            ax.plot(group["threshold"], group["edges"], marker="o", label=f"{species}:{query}")
        ax.set_xlabel("STRING combined-score threshold")
        ax.set_ylabel("Retained edges")
        ax.set_title("Network threshold robustness")
        ax.legend(fontsize=7)
        _save(fig, output_dir / "stage6_network_thresholds")
        count += 1

    fig, ax = plt.subplots(figsize=(7, 5))
    counts = (
        tiers["candidate_tier"]
        .value_counts()
        .reindex(["Tier 1", "Tier 2", "Tier 3", "Watchlist"], fill_value=0)
    )
    ax.bar(counts.index, counts.values)
    ax.set_ylabel("Candidate sites")
    ax.set_title("Stage 6 final candidate tiers")
    _save(fig, output_dir / "stage6_candidate_tiers")
    count += 1
    return count
