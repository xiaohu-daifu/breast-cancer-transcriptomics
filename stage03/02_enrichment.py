"""
Run pathway enrichment on the Stage 3 gene lists using gseapy.

- Over-representation analysis (ORA) via Enrichr for GO, KEGG, and Hallmark.
- Gene set enrichment analysis (GSEA) via a pre-ranked gene list for Hallmark.

Both require network access to the Enrichr API. Requests are issued
per-library with retries to tolerate transient Enrichr failures.
"""

from pathlib import Path
import json
import time

import gseapy
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parent.parent

RESULTS_DIR = PROJECT_ROOT / "stage03" / "results"
ORA_DIR = RESULTS_DIR / "ora"
GSEA_DIR = RESULTS_DIR / "gsea"

ORA_DIR.mkdir(parents=True, exist_ok=True)
GSEA_DIR.mkdir(parents=True, exist_ok=True)

LIBRARIES = [
    "GO_Biological_Process_2021",
    "KEGG_2021_Human",
    "MSigDB_Hallmark_2020",
]

CONTRASTS = [
    "Basal_vs_LuminalA",
    "Basal_vs_LuminalB",
    "Basal_vs_HER2",
    "LuminalB_vs_LuminalA",
    "HER2_vs_LuminalA",
    "HER2_vs_LuminalB",
]

MAX_RETRIES = 3
RETRY_DELAY = 5


def enrichr_with_retry(genes, library):
    """Run a single-library Enrichr ORA query with retries."""
    last_error = None
    for attempt in range(1, MAX_RETRIES + 1):
        try:
            res = gseapy.enrichr(
                gene_list=genes,
                gene_sets=[library],
                organism="human",
                outdir=None,
            )
            return res.results
        except Exception as exc:
            last_error = exc
            print(f"  {library} attempt {attempt}/{MAX_RETRIES} failed: {exc}")
            if attempt < MAX_RETRIES:
                time.sleep(RETRY_DELAY)
    raise last_error


with open(RESULTS_DIR / "gene_lists.json") as f:
    gene_lists = json.load(f)



# Over-representation analysis (ORA)


ora_rows = []

for contrast in CONTRASTS:

    for direction in ["up", "down"]:

        genes = gene_lists[f"{contrast}_{direction}"]

        if not genes:
            continue

        for library in LIBRARIES:

            df = enrichr_with_retry(genes, library)

            df["Contrast"] = contrast
            df["Direction"] = direction
            df["Library"] = library

            ora_rows.append(df)

            print(f"ORA {contrast} {direction} {library}: {len(df)} terms")

            time.sleep(1)


ora = pd.concat(ora_rows, ignore_index=True)
ora.to_csv(ORA_DIR / "ora_all_terms.csv", index=False)

print(f"\nSaved ORA results ({len(ora)} rows) to {ORA_DIR / 'ora_all_terms.csv'}")



# Gene set enrichment analysis (GSEA) — Hallmark


gsea_rows = []

for contrast in CONTRASTS:

    ranked = pd.Series(dict(gene_lists[f"{contrast}_ranked"]))

    res = gseapy.prerank(
        rnk=ranked,
        gene_sets="MSigDB_Hallmark_2020",
        outdir=None,
        min_size=15,
        max_size=500,
    )

    df = res.res2d.copy()
    df["Contrast"] = contrast

    gsea_rows.append(df)

    print(f"GSEA {contrast}: {len(df)} terms")

    time.sleep(1)


gsea = pd.concat(gsea_rows, ignore_index=True)
gsea.to_csv(GSEA_DIR / "gsea_hallmark.csv", index=False)

print(f"\nSaved GSEA results ({len(gsea)} rows) to {GSEA_DIR / 'gsea_hallmark.csv'}")
