"""
Curate a compact, human-readable table of the top differentially
expressed genes per contrast for the Stage 2 deliverable.

The full DE tables are large and gitignored; this script distils them
into a small committed summary (top genes by |logFC| per contrast) that
is visible in the repository and feeds downstream pathway analysis.
"""

from pathlib import Path

import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parent.parent

RESULTS_DIR = PROJECT_ROOT / "stage02" / "results"
OUTPUT_PATH = RESULTS_DIR / "top_DEGs_per_contrast.csv"

FDR_THRESHOLD = 0.05
LOGFC_THRESHOLD = 1.0
N_TOP = 50

CONTRASTS = [
    "Basal_vs_LuminalA",
    "Basal_vs_LuminalB",
    "Basal_vs_HER2",
    "LuminalB_vs_LuminalA",
    "HER2_vs_LuminalA",
    "HER2_vs_LuminalB",
]


def display_label(contrast):
    label = contrast.replace("LuminalA", "Luminal A")
    label = label.replace("LuminalB", "Luminal B")
    return label.replace("_vs_", " vs ")


frames = []

for contrast in CONTRASTS:

    de = pd.read_csv(RESULTS_DIR / f"DE_{contrast}_annotated.csv")

    significant = de[
        (de["adj.P.Val"] < FDR_THRESHOLD)
        & (de["logFC"].abs() >= LOGFC_THRESHOLD)
    ].copy()

    significant["Direction"] = significant["logFC"].apply(
        lambda x: "Up" if x > 0 else "Down"
    )

    top = significant.reindex(
        significant["logFC"].abs().sort_values(ascending=False).index
    ).head(N_TOP)

    top.insert(0, "Contrast", display_label(contrast))

    frames.append(top)


result = pd.concat(frames, ignore_index=True)

columns = [
    "Contrast",
    "Direction",
    "Probe_ID",
    "SYMBOL",
    "ENTREZID",
    "GENENAME",
    "logFC",
    "adj.P.Val",
    "P.Value",
    "t",
    "B",
]

result = result[columns]

result.to_csv(OUTPUT_PATH, index=False)

print(f"Saved {len(result)} rows to {OUTPUT_PATH}")
print("\nRows per contrast:")
print(result["Contrast"].value_counts())
