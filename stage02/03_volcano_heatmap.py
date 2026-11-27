"""
Generate publication figures for the Stage 2 differential-expression
results:

- Volcano plots for all six pairwise contrasts
- Heatmap of the top 50 differentially expressed genes

Inputs are the annotated DE tables produced by ``02_annotate_probes.R``
and the processed expression matrix from Stage 1.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


PROJECT_ROOT = Path(__file__).resolve().parent.parent

RESULTS_DIR = PROJECT_ROOT / "stage02" / "results"
FIGURES_DIR = PROJECT_ROOT / "stage02" / "figures"
EXPRESSION_PATH = PROJECT_ROOT / "data" / "GSE45827_expression.csv"
METADATA_PATH = PROJECT_ROOT / "data" / "GSE45827_metadata.csv"

FIGURES_DIR.mkdir(exist_ok=True)

FDR_THRESHOLD = 0.05
LOGFC_THRESHOLD = 1.0
N_TOP_GENES = 50
N_LABEL = 10

SUBTYPE_PALETTE = {
    "Basal": "#d62728",
    "HER2": "#ff7f0e",
    "Luminal A": "#1f77b4",
    "Luminal B": "#2ca02c",
}

CONTRASTS = [
    "Basal_vs_LuminalA",
    "Basal_vs_LuminalB",
    "Basal_vs_HER2",
    "LuminalB_vs_LuminalA",
    "HER2_vs_LuminalA",
    "HER2_vs_LuminalB",
]


def display_label(contrast):
    """Convert a contrast name to a human-readable label."""
    label = contrast.replace("LuminalA", "Luminal A")
    label = label.replace("LuminalB", "Luminal B")
    return label.replace("_vs_", " vs ")


def first_symbol(value):
    """Return the first gene symbol for a (possibly multi-mapped) probe."""
    if pd.isna(value):
        return None
    return str(value).split(";")[0].strip()


def save_figure(name):
    """Save and close the current figure."""
    plt.tight_layout()
    plt.savefig(FIGURES_DIR / name, dpi=300, bbox_inches="tight")
    plt.close()


# ---------------------------------------------------------------------------
# Volcano plots
# ---------------------------------------------------------------------------

for contrast in CONTRASTS:

    path = RESULTS_DIR / f"DE_{contrast}_annotated.csv"
    de = pd.read_csv(path)

    logfc = de["logFC"]
    fdr = de["adj.P.Val"]
    neg_log10_fdr = -np.log10(fdr.clip(lower=1e-300))

    is_sig = fdr < FDR_THRESHOLD
    is_up = is_sig & (logfc >= LOGFC_THRESHOLD)
    is_down = is_sig & (logfc <= -LOGFC_THRESHOLD)

    plt.figure(figsize=(8, 6))

    plt.scatter(
        logfc[~is_sig],
        neg_log10_fdr[~is_sig],
        s=8,
        c="lightgray",
        alpha=0.4,
        label="Not significant",
    )
    plt.scatter(
        logfc[is_up],
        neg_log10_fdr[is_up],
        s=8,
        c="#d62728",
        alpha=0.5,
        label="Upregulated",
    )
    plt.scatter(
        logfc[is_down],
        neg_log10_fdr[is_down],
        s=8,
        c="#1f77b4",
        alpha=0.5,
        label="Downregulated",
    )

    plt.axhline(
        -np.log10(FDR_THRESHOLD),
        color="gray",
        linestyle="--",
        linewidth=0.8,
    )
    plt.axvline(LOGFC_THRESHOLD, color="gray", linestyle="--", linewidth=0.8)
    plt.axvline(-LOGFC_THRESHOLD, color="gray", linestyle="--", linewidth=0.8)

    top = de[is_sig].reindex(
        de[is_sig]["logFC"].abs().sort_values(ascending=False).index
    ).head(N_LABEL)

    for _, row in top.iterrows():
        symbol = first_symbol(row["SYMBOL"])
        if symbol is not None:
            plt.annotate(
                symbol,
                (row["logFC"], neg_log10_fdr[row.name]),
                fontsize=8,
                alpha=0.9,
            )

    plt.xlabel("log2 fold change")
    plt.ylabel("-log10(FDR)")
    plt.title(f"{display_label(contrast)} (n = {de.shape[0]} probes)")
    plt.legend(loc="best", fontsize=8)

    save_figure(f"volcano_{contrast}.png")

# ---------------------------------------------------------------------------
# Heatmap of top 50 differentially expressed genes
# ---------------------------------------------------------------------------

expr = pd.read_csv(EXPRESSION_PATH, index_col=0)
metadata = pd.read_csv(METADATA_PATH).set_index("sample_id")

# Use the strongest contrast (Basal vs Luminal A) to select the top genes.
de = pd.read_csv(RESULTS_DIR / "DE_Basal_vs_LuminalA_annotated.csv")

sig = de[(de["adj.P.Val"] < FDR_THRESHOLD) & (de["logFC"].abs() >= LOGFC_THRESHOLD)]

top_genes = sig.reindex(
    sig["logFC"].abs().sort_values(ascending=False).index
).head(N_TOP_GENES)

top_probes = top_genes["Probe_ID"].tolist()
top_labels = [first_symbol(s) or p for p, s in zip(top_probes, top_genes["SYMBOL"])]

top_expr = expr.loc[top_probes]

# Z-score each gene across samples for visualisation.
top_expr_z = top_expr.sub(top_expr.mean(axis=1), axis=0).div(
    top_expr.std(axis=1), axis=0
)

# Order samples by subtype, then by sample id.
subtype_order = ["Luminal A", "Luminal B", "HER2", "Basal"]
sample_order = sorted(
    top_expr_z.columns,
    key=lambda s: (
        subtype_order.index(metadata.loc[s, "subtype"]),
        s,
    ),
)

top_expr_z = top_expr_z[sample_order]

row_colors = [SUBTYPE_PALETTE[metadata.loc[s, "subtype"]] for s in sample_order]

# Rename rows to gene symbols for readability.
top_expr_z.index = top_labels

grid = sns.clustermap(
    top_expr_z,
    col_cluster=False,
    row_cluster=True,
    cmap="vlag",
    center=0,
    figsize=(14, 10),
    col_colors=row_colors,
    yticklabels=True,
    xticklabels=False,
    cbar_kws={"label": "Z-scored expression"},
)

# Add a legend for the subtype colour bar.
for subtype, colour in SUBTYPE_PALETTE.items():
    grid.ax_col_dendrogram.bar(0, 0, color=colour, label=subtype)
grid.ax_col_dendrogram.legend(
    title="Subtype",
    loc="center left",
    bbox_to_anchor=(1.01, 0.5),
    fontsize=8,
)

grid.fig.suptitle(
    "Top 50 differentially expressed genes (Basal vs Luminal A)",
    y=1.03,
)

grid.savefig(FIGURES_DIR / "heatmap_top50_deg.png", dpi=300, bbox_inches="tight")
plt.close()
