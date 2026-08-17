from pathlib import Path
import gzip
import pandas as pd


data_path = Path("../data/GSE45827_series_matrix.txt.gz")


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


print(metadata.head())
print()
print(metadata["subtype"].value_counts(dropna=False))
valid_subtypes = [
    "Basal",
    "Her2",
    "Luminal A",
    "Luminal B",
]

metadata = metadata[metadata["subtype"].isin(valid_subtypes)].copy()

print(metadata.shape)
print(metadata["subtype"].value_counts())

# assert metadata["sample_id"].is_unique
# assert metadata["subtype"].notna().all()
# assert set(metadata["subtype"].unique()) == {
#     "Basal",
#     "Her2",
#     "Luminal A",
#     "Luminal B",
# }

print("Unique subtypes:")
print(metadata["subtype"].unique())

print("\nNumber of unique subtypes:")
print(metadata["subtype"].nunique())

print("\nCounts:")
print(metadata["subtype"].value_counts(dropna=False))

print("\nDuplicate sample IDs:")
print(metadata["sample_id"].duplicated().sum())

print("\nMissing values:")
print(metadata.isna().sum())

print(metadata.shape)
print(metadata["subtype"].value_counts())
print(metadata.isna().sum())
metadata.to_csv("../data/GSE45827_metadata.csv", index=False)