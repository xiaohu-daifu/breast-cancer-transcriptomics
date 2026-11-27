"""
Explore the GEO dataset structure and extract, clean, and filter
sample metadata for the breast cancer analysis cohort.

In addition to molecular subtype, this extracts the processing
``batch`` and ``hybridization`` technical variables, which are
needed for the batch-effect assessment performed during QC.
"""

import gzip
from pathlib import Path

import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parent.parent

data_path = PROJECT_ROOT / "data" / "GSE45827_series_matrix.txt.gz"
output_path = PROJECT_ROOT / "data" / "GSE45827_metadata.csv"


def clean_value(value):
    """Strip a GEO characteristic cell down to its bare value."""
    v = value.strip().strip('"')
    return v.split(":", 1)[1].strip() if ":" in v else v


sample_ids = []
subtype_values = []
batch_values = []
hybridization_values = []

with gzip.open(data_path, "rt", encoding="utf-8") as file:

    for line in file:

        if line.startswith("!Sample_geo_accession"):
            sample_ids = line.strip().split("\t")[1:]

        elif line.startswith("!Sample_characteristics_ch1"):
            values = line.strip().split("\t")[1:]

            if any("tumor subtype:" in value for value in values):
                subtype_values = values

            elif any("batch:" in value for value in values):
                batch_values = values

            elif any("hybridation:" in value for value in values):
                hybridization_values = values


metadata = pd.DataFrame({"sample_id": sample_ids})

metadata["sample_id"] = metadata["sample_id"].str.strip().str.strip('"')

metadata["subtype"] = [clean_value(value) for value in subtype_values]
metadata["batch"] = [clean_value(value) for value in batch_values]
metadata["hybridization"] = [
    clean_value(value) for value in hybridization_values
]

# Normalise the HER2 label to the conventional uppercase form.
metadata["subtype"] = metadata["subtype"].replace("Her2", "HER2")

valid_subtypes = [
    "Basal",
    "HER2",
    "Luminal A",
    "Luminal B",
]

metadata = metadata[
    metadata["subtype"].isin(valid_subtypes)
].copy()

metadata["batch"] = pd.to_numeric(metadata["batch"])
metadata["hybridization"] = pd.to_numeric(metadata["hybridization"])


assert metadata["sample_id"].is_unique
assert metadata["subtype"].notna().all()
assert set(metadata["subtype"].unique()) == set(valid_subtypes)
assert metadata[["batch", "hybridization"]].notna().all().all()


print("Cohort size:", len(metadata))
print("\nSubtype counts:")
print(metadata["subtype"].value_counts())

print("\nBatch by subtype:")
print(pd.crosstab(metadata["subtype"], metadata["batch"]))

print("\nHybridization by subtype:")
print(pd.crosstab(metadata["subtype"], metadata["hybridization"]))

metadata.to_csv(output_path, index=False)

print(f"\nSaved metadata to: {output_path}")
