from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def _save(fig: plt.Figure, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(path.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(path.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def enrichment_bubble_plot(enrichment: pd.DataFrame, output: Path, *, top_n: int = 20) -> None:
    if enrichment.empty:
        return
    data = enrichment.loc[enrichment["significant"].fillna(False)].copy()
    if data.empty:
        return
    data["score"] = -np.log10(
        pd.to_numeric(data["p_value_adjusted"], errors="coerce").clip(lower=1e-300)
    )
    data = (
        data.sort_values(["score", "fold_enrichment"], ascending=False)
        .head(top_n)
        .sort_values("score")
    )
    fig, ax = plt.subplots(figsize=(9, max(5, 0.32 * len(data))))
    sizes = pd.to_numeric(data["intersection_size"], errors="coerce").fillna(1).clip(lower=1) * 18
    scatter = ax.scatter(data["score"], np.arange(len(data)), s=sizes, c=data["fold_enrichment"])
    context = (
        "["
        + data["species"].astype(str)
        + "/"
        + data["query_name"].astype(str)
        + "/"
        + data["source"].astype(str)
        + "] "
        + data["term_name"].astype(str)
    )
    ax.set_yticks(np.arange(len(data)), context.str.slice(0, 100))
    ax.set_xlabel("-log10 adjusted p-value")
    ax.set_title("Leading Stage 5 enrichment terms")
    fig.colorbar(scatter, ax=ax, label="Fold enrichment")
    _save(fig, output)


def cross_species_scatter(pairs: pd.DataFrame, output: Path) -> None:
    if pairs.empty:
        return
    fig, ax = plt.subplots(figsize=(7, 7))
    ax.scatter(
        pairs["repairability_extra_trees_human"], pairs["repairability_extra_trees_mouse"], s=50
    )
    low = min(ax.get_xlim()[0], ax.get_ylim()[0])
    high = max(ax.get_xlim()[1], ax.get_ylim()[1])
    ax.plot([low, high], [low, high], linestyle="--")
    ax.set_xlabel("Human Extra Trees repairability")
    ax.set_ylabel("Mouse Extra Trees repairability")
    ax.set_title("Cross-species ortholog concordance")
    for row in pairs.head(12).itertuples(index=False):
        ax.annotate(
            str(row.source_gene),
            (row.repairability_extra_trees_human, row.repairability_extra_trees_mouse),
            fontsize=7,
        )
    _save(fig, output)


def module_heatmap(profiles: pd.DataFrame, output: Path) -> None:
    if profiles.empty:
        return
    numeric = [column for column in MODULE_PLOT_COLUMNS if column in profiles]
    values = profiles[numeric].to_numpy(dtype=float)
    means = np.nanmean(values, axis=0)
    stds = np.nanstd(values, axis=0)
    stds[stds == 0] = 1
    z = (values - means) / stds
    fig, ax = plt.subplots(figsize=(10, max(3.5, 0.8 * len(profiles))))
    image = ax.imshow(z, aspect="auto")
    ax.set_xticks(np.arange(len(numeric)), numeric, rotation=45, ha="right")
    labels = profiles["module_id"].astype(str) + ": " + profiles["module_label"].astype(str)
    ax.set_yticks(np.arange(len(profiles)), labels)
    ax.set_title("Stage 5 candidate modules")
    fig.colorbar(image, ax=ax, label="Module z-score")
    _save(fig, output)


MODULE_PLOT_COLUMNS = [
    "repairability_extra_trees",
    "evidence_score_raw",
    "transfer_reliability_score",
    "model_disagreement",
    "site_plddt",
    "relative_sasa_lysine",
    "neighbor_cysteine_count_8a",
]


def priority_landscape(portfolio: pd.DataFrame, output: Path) -> None:
    if portfolio.empty:
        return
    fig, ax = plt.subplots(figsize=(8, 6))
    size = 25 + 90 * pd.to_numeric(portfolio["evidence_score_norm"], errors="coerce").fillna(0)
    scatter = ax.scatter(
        portfolio["repairability_extra_trees"],
        portfolio["transfer_reliability_score"],
        s=size,
        c=portfolio["weighted_priority_score"],
        alpha=0.8,
    )
    ax.set_xlabel("Extra Trees repairability")
    ax.set_ylabel("Transfer reliability")
    ax.set_title("Priority portfolio landscape")
    fig.colorbar(scatter, ax=ax, label="Weighted priority")
    for row in portfolio.nsmallest(12, "portfolio_rank").itertuples(index=False):
        gene = str(row.gene_names).split()[0]
        ax.annotate(
            f"{gene} K{int(row.site_uniprot)}",
            (row.repairability_extra_trees, row.transfer_reliability_score),
            fontsize=7,
        )
    _save(fig, output)
