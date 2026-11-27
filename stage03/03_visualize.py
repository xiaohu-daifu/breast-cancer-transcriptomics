"""
Visualise Stage 3 pathway enrichment results.

- Dot plot of enriched MSigDB Hallmark pathways across the four
  subtype signatures.
- Dot plot of enriched GO Biological Process terms for the Basal vs
  Luminal A comparison.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parent.parent

RESULTS_DIR = PROJECT_ROOT / "stage03" / "results"
FIGURES_DIR = PROJECT_ROOT / "stage03" / "figures"
ORA_PATH = RESULTS_DIR / "ora" / "ora_all_terms.csv"

FIGURES_DIR.mkdir(parents=True, exist_ok=True)

ADJ_P_THRESHOLD = 0.05
N_TOP_TERMS = 15


def display_label(contrast):
    label = contrast.replace("LuminalA", "Luminal A")
    label = label.replace("LuminalB", "Luminal B")
    return label.replace("_vs_", " vs ")


# Each subtype signature maps to a contrast + direction that captures
# genes up-regulated in that subtype.
SIGNATURES = [
    ("Basal", "Basal_vs_LuminalA", "up"),
    ("Luminal A", "Basal_vs_LuminalA", "down"),
    ("HER2", "HER2_vs_LuminalA", "up"),
    ("Luminal B", "LuminalB_vs_LuminalA", "up"),
]


def parse_overlap(overlap):
    """Parse '134/200' into (hits, size)."""
    hits, size = overlap.split("/")
    return int(hits), int(size)


def top_terms(ora, library, contrast, direction, n=N_TOP_TERMS):
    """Return the top n significant terms for a given query."""
    sub = ora[
        (ora["Library"] == library)
        & (ora["Contrast"] == contrast)
        & (ora["Direction"] == direction)
        & (ora["Adjusted P-value"] < ADJ_P_THRESHOLD)
    ].copy()

    sub[["hits", "size"]] = sub["Overlap"].apply(
        lambda x: pd.Series(parse_overlap(x))
    )
    sub["gene_ratio"] = sub["hits"] / sub["size"]
    sub["neg_log10_p"] = -np.log10(sub["Adjusted P-value"])

    return sub.sort_values("Adjusted P-value").head(n)


def draw_dotplot(ax, sub, title):
    """Draw a single dot plot on an axis."""
    sub = sub.sort_values("Adjusted P-value", ascending=False)

    scatter = ax.scatter(
        sub["gene_ratio"],
        sub["Term"],
        s=sub["hits"] * 3,
        c=sub["neg_log10_p"],
        cmap="viridis",
        alpha=0.85,
        edgecolors="black",
        linewidths=0.3,
    )

    ax.set_xlabel("Gene ratio")
    ax.set_title(title, fontsize=11)
    ax.tick_params(axis="y", labelsize=8)

    return scatter


ora = pd.read_csv(ORA_PATH)


# Hallmark dot plot (four subtype signatures)


fig, axes = plt.subplots(2, 2, figsize=(14, 12))

for ax, (subtype, contrast, direction) in zip(axes.flat, SIGNATURES):

    sub = top_terms(ora, "MSigDB_Hallmark_2020", contrast, direction)

    scatter = draw_dotplot(
        ax,
        sub,
        f"{subtype} signature\n({display_label(contrast)}, {direction})",
    )

fig.suptitle(
    "Enriched MSigDB Hallmark pathways across breast cancer subtypes",
    y=1.02,
    fontsize=13,
)

cbar = fig.colorbar(scatter, ax=axes, shrink=0.8, pad=0.02)
cbar.set_label("-log10(adjusted P-value)")

plt.tight_layout()
plt.savefig(
    FIGURES_DIR / "hallmark_dotplot.png",
    dpi=300,
    bbox_inches="tight",
)
plt.close()



# GO Biological Process dot plot (Basal vs Luminal A)


fig, axes = plt.subplots(1, 2, figsize=(14, 7))

for ax, direction, title in [
    (axes[0], "up", "Basal signature (Basal vs Luminal A, up)"),
    (axes[1], "down", "Luminal A signature (Basal vs Luminal A, down)"),
]:

    sub = top_terms(ora, "GO_Biological_Process_2021", "Basal_vs_LuminalA", direction)

    scatter = draw_dotplot(ax, sub, title)

fig.suptitle(
    "Enriched GO Biological Process terms (Basal vs Luminal A)",
    y=1.02,
    fontsize=13,
)

cbar = fig.colorbar(scatter, ax=axes, shrink=0.8, pad=0.02)
cbar.set_label("-log10(adjusted P-value)")

plt.tight_layout()
plt.savefig(
    FIGURES_DIR / "go_dotplot_basal_luminal.png",
    dpi=300,
    bbox_inches="tight",
)
plt.close()


# ---------------------------------------------------------------------------
# KEGG dot plot (four subtype signatures)
# ---------------------------------------------------------------------------

fig, axes = plt.subplots(2, 2, figsize=(14, 12))

for ax, (subtype, contrast, direction) in zip(axes.flat, SIGNATURES):

    sub = top_terms(ora, "KEGG_2021_Human", contrast, direction)

    scatter = draw_dotplot(
        ax,
        sub,
        f"{subtype} signature\n({display_label(contrast)}, {direction})",
    )

fig.suptitle(
    "Enriched KEGG pathways across breast cancer subtypes",
    y=1.02,
    fontsize=13,
)

cbar = fig.colorbar(scatter, ax=axes, shrink=0.8, pad=0.02)
cbar.set_label("-log10(adjusted P-value)")

plt.tight_layout()
plt.savefig(
    FIGURES_DIR / "kegg_dotplot.png",
    dpi=300,
    bbox_inches="tight",
)
plt.close()


print("Saved figures to", FIGURES_DIR)
