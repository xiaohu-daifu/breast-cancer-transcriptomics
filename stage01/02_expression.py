"""
Extract and validate the breast cancer expression matrix from
the GSE45827 GEO series matrix and select the analysis cohort.
"""

from pathlib import Path
from io import StringIO
import gzip
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parent.parent

data_path = PROJECT_ROOT / "data" / "GSE45827_series_matrix.txt.gz"
metadata_path = PROJECT_ROOT / "data" / "GSE45827_metadata.csv"
expression_path = PROJECT_ROOT / "data" / "GSE45827_expression.csv"


expression_lines = []

with gzip.open(data_path, "rt", encoding="utf-8") as file:

    in_matrix = False

    for line in file:

        if line.startswith("!series_matrix_table_begin"):
            in_matrix = True
            continue

        if line.startswith("!series_matrix_table_end"):
            break

        if in_matrix:
            expression_lines.append(line)


expression = pd.read_csv(
    StringIO("".join(expression_lines)),
    sep="\t"
)


metadata = pd.read_csv(metadata_path)


selected_samples = metadata["sample_id"].tolist()

assert metadata["sample_id"].isin(expression.columns).all()

expression = expression[
    ["ID_REF"] + selected_samples
]


assert list(expression.columns[1:]) == selected_samples
assert expression["ID_REF"].is_unique
assert not expression.iloc[:, 1:].isna().any().any()


print("Expression matrix shape:", expression.shape)


expression.to_csv(expression_path, index=False)

print(f"Saved expression matrix to: {expression_path}")