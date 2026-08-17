"""
Perform exploratory quality control of the cleaned GSE45827
breast cancer expression matrix.

The analysis includes:
- Expression value distribution
- Sample-wise expression distributions
- Low-variance gene filtering
- Principal component analysis (PCA)
- PCA visualization by molecular subtype
"""

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from sklearn.decomposition import PCA


PROJECT_ROOT = Path(__file__).resolve().parent.parent

expression_path = PROJECT_ROOT / "data" / "GSE45827_expression.csv"
metadata_path = PROJECT_ROOT / "data" / "GSE45827_metadata.csv"

FIGURES_DIR = PROJECT_ROOT / "figures"
FIGURES_DIR.mkdir(exist_ok=True)

expression = pd.read_csv(expression_path)
metadata = pd.read_csv(metadata_path).set_index("sample_id")

expression_values = expression.iloc[:, 1:]


values = expression_values.to_numpy().flatten()

plt.figure(figsize=(8, 5))
plt.hist(values, bins=100)

plt.xlabel("Expression value")
plt.ylabel("Frequency")
plt.title("GSE45827 expression value distribution")

plt.tight_layout()
plt.savefig(
    PROJECT_ROOT / "figures" / "expression_distribution.png",
    dpi=300,
    bbox_inches="tight",
)
plt.show()
plt.close



plt.figure(figsize=(12, 5))
plt.boxplot(
    expression_values.values,
    showfliers=False,
)

plt.xlabel("Samples")
plt.ylabel("Expression")
plt.title("Expression distribution across samples")

plt.tight_layout()
plt.savefig(
    PROJECT_ROOT / "figures" / "sample_expression_boxplot.png",
    dpi=300,
    bbox_inches="tight",
)
plt.show()
plt.close


X = expression_values.T

print("Expression matrix:", X.shape)

gene_variance = X.var(axis=0)

variance_threshold = 0.1

X = X.loc[:, gene_variance > variance_threshold]

print(
    f"After variance filtering "
    f"(threshold > {variance_threshold}): {X.shape}"
)



pca = PCA(n_components=2)

X_pca = pca.fit_transform(X)

pc1_variance = pca.explained_variance_ratio_[0]
pc2_variance = pca.explained_variance_ratio_[1]

print("\nVariance explained:")
print(f"PC1: {pc1_variance * 100:.1f}%")
print(f"PC2: {pc2_variance * 100:.1f}%")
print(f"PC1 + PC2: {(pc1_variance + pc2_variance) * 100:.1f}%")



sample_ids = expression.columns[1:]

pca_df = pd.DataFrame(
    X_pca,
    columns=["PC1", "PC2"],
    index=sample_ids,
)

pca_df["subtype"] = metadata.loc[sample_ids, "subtype"]



plt.figure(figsize=(8, 6))

for subtype in sorted(pca_df["subtype"].unique()):

    subset = pca_df[pca_df["subtype"] == subtype]

    plt.scatter(
        subset["PC1"],
        subset["PC2"],
        label=subtype,
    )

plt.xlabel(f"PC1 ({pc1_variance * 100:.1f}%)")
plt.ylabel(f"PC2 ({pc2_variance * 100:.1f}%)")
plt.title("PCA of GSE45827 breast cancer samples")
plt.legend(title="Subtype")

plt.tight_layout()
plt.savefig(
    PROJECT_ROOT / "figures" / "pca_subtypes.png",
    dpi=300,
    bbox_inches="tight",
)
plt.show()
plt.close