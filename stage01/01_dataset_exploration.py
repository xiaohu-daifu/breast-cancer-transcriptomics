"""
Explore the GEO dataset structure and extract, clean, and filter
sample metadata for the breast cancer analysis cohort.
"""

from pathlib import Path
import gzip
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parent.parent

data_path = PROJECT_ROOT / "data" / "GSE45827_series_matrix.txt.gz"
output_path = PROJECT_ROOT / "data" / "GSE45827_metadata.csv"


sample_ids = []
subtypes = []

with gzip.open(data_path, "rt", encoding="utf-8") as file:

    for line in file:

        if line.startswith("!Sample_geo_accession"):
            sample_ids = line.strip().split("\t")[1:]

        elif line.startswith("!Sample_characteristics_ch1"):
            values = line.strip().split("\t")[1:]

            if any("tumor subtype:" in value for value in values):
                subtypes = values


metadata = pd.DataFrame({
    "sample_id": sample_ids,
    "subtype": subtypes,
})

metadata["sample_id"] = metadata["sample_id"].str.strip('"')

metadata["subtype"] = (
    metadata["subtype"]
    .str.replace("tumor subtype: ", "", regex=False)
    .str.strip('"')
)

valid_subtypes = [
    "Basal",
    "Her2",
    "Luminal A",
    "Luminal B",
]

metadata = metadata[
    metadata["subtype"].isin(valid_subtypes)
].copy()


assert metadata["sample_id"].is_unique
assert metadata["subtype"].notna().all()
assert set(metadata["subtype"].unique()) == set(valid_subtypes)


print("Cohort size:", len(metadata))
print("\nSubtype counts:")
print(metadata["subtype"].value_counts())

metadata.to_csv(output_path, index=False)

print(f"\nSaved metadata to: {output_path}")