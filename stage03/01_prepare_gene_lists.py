"""
Prepare gene lists for pathway enrichment from the Stage 2
differential-expression results.

Steps:
- Collapse Affymetrix probes to gene symbols (splitting multi-mappings).
- Define the background universe as all detected genes.
- Extract up- and down-regulated significant genes per contrast.
- Build a gene-level ranked list (by log2FC) per contrast for GSEA.

Outputs a single JSON manifest consumed by ``02_enrichment.py``.
"""

from pathlib import Path
import json

import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parent.parent

STAGE02_RESULTS = PROJECT_ROOT / "stage02" / "results"
OUT_DIR = PROJECT_ROOT / "stage03" / "results"
OUT_DIR.mkdir(parents=True, exist_ok=True)

FDR_THRESHOLD = 0.05
LOGFC_THRESHOLD = 1.0

CONTRASTS = [
    "Basal_vs_LuminalA",
    "Basal_vs_LuminalB",
    "Basal_vs_HER2",
    "LuminalB_vs_LuminalA",
    "HER2_vs_LuminalA",
    "HER2_vs_LuminalB",
]


def collapse_symbols(symbol_series):
    """Return the sorted set of gene symbols for a probe SYMBOL series."""
    symbols = set()
    for value in symbol_series:
        if pd.isna(value):
            continue
        for sym in str(value).split(";"):
            sym = sym.strip()
            if sym:
                symbols.add(sym)
    return sorted(symbols)


def gene_level_ranking(de):
    """Collapse probes to one log2FC per gene (max |log2FC|)."""
    rows = []
    for sym, logfc in zip(de["SYMBOL"], de["logFC"]):
        if pd.isna(sym):
            continue
        for s in str(sym).split(";"):
            s = s.strip()
            if s:
                rows.append((s, logfc))

    df = pd.DataFrame(rows, columns=["symbol", "logFC"])
    df["abs_logfc"] = df["logFC"].abs()
    df = df.sort_values("abs_logfc", ascending=False).drop_duplicates("symbol")

    return df.set_index("symbol")["logFC"]


# Background universe: all detected genes (all probes in a contrast table).
de_all = pd.read_csv(STAGE02_RESULTS / f"DE_{CONTRASTS[0]}_annotated.csv")
background = collapse_symbols(de_all["SYMBOL"])

gene_lists = {"background": background}

for contrast in CONTRASTS:

    de = pd.read_csv(STAGE02_RESULTS / f"DE_{contrast}_annotated.csv")

    significant = de[
        (de["adj.P.Val"] < FDR_THRESHOLD)
        & (de["logFC"].abs() >= LOGFC_THRESHOLD)
    ]

    up = collapse_symbols(significant[significant["logFC"] > 0]["SYMBOL"])
    down = collapse_symbols(significant[significant["logFC"] < 0]["SYMBOL"])

    gene_lists[f"{contrast}_up"] = up
    gene_lists[f"{contrast}_down"] = down

    ranked = gene_level_ranking(de).sort_values(ascending=False)
    gene_lists[f"{contrast}_ranked"] = [
        [gene, float(score)] for gene, score in ranked.items()
    ]

    print(f"{contrast}: up={len(up)}  down={len(down)}  ranked={len(ranked)}")

print(f"\nBackground universe: {len(background)} genes")

output_path = OUT_DIR / "gene_lists.json"

with open(output_path, "w") as f:
    json.dump(gene_lists, f, indent=2)

print(f"Saved to {output_path}")
