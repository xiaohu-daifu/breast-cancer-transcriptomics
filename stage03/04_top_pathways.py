"""
Curate a compact, committed summary of the top enriched pathways per
subtype from the Stage 3 enrichment results.

The full ORA table is large and gitignored; this script distils it into
a small committed table (top 10 pathways per subtype per library) that is
visible in the repository.
"""

from pathlib import Path

import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parent.parent

ORA_PATH = PROJECT_ROOT / "stage03" / "results" / "ora" / "ora_all_terms.csv"
OUTPUT_PATH = (
    PROJECT_ROOT / "stage03" / "results" / "top_pathways_per_subtype.csv"
)

ADJ_P_THRESHOLD = 0.05
N_TOP = 10

SIGNATURES = [
    ("Basal", "Basal_vs_LuminalA", "up"),
    ("Luminal A", "Basal_vs_LuminalA", "down"),
    ("HER2", "HER2_vs_LuminalA", "up"),
    ("Luminal B", "LuminalB_vs_LuminalA", "up"),
]

LIBRARIES = [
    "MSigDB_Hallmark_2020",
    "KEGG_2021_Human",
    "GO_Biological_Process_2021",
]


def parse_overlap(overlap):
    hits, size = overlap.split("/")
    return int(hits), int(size)


ora = pd.read_csv(ORA_PATH)

rows = []

for subtype, contrast, direction in SIGNATURES:

    for library in LIBRARIES:

        sub = ora[
            (ora["Library"] == library)
            & (ora["Contrast"] == contrast)
            & (ora["Direction"] == direction)
            & (ora["Adjusted P-value"] < ADJ_P_THRESHOLD)
        ].copy()

        sub = sub.sort_values("Adjusted P-value").head(N_TOP)

        for _, r in sub.iterrows():

            hits, size = parse_overlap(r["Overlap"])

            rows.append({
                "Subtype": subtype,
                "Library": library,
                "Term": r["Term"],
                "Gene_count": hits,
                "Gene_set_size": size,
                "Gene_ratio": round(hits / size, 3),
                "Adjusted_P_value": r["Adjusted P-value"],
            })


result = pd.DataFrame(rows)

result.to_csv(OUTPUT_PATH, index=False)

print(f"Saved {len(result)} rows to {OUTPUT_PATH}")
print("\nRows per subtype:")
print(result["Subtype"].value_counts())
