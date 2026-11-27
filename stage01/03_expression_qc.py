"""
Perform exploratory quality control of the cleaned GSE45827
breast cancer expression matrix.

The analysis includes:
- Expression value distribution
- Sample-wise expression distributions
- Low-variance gene filtering
- Principal component analysis (PCA) by subtype
- Heatmap of top variable genes
- Sample correlation heatmap
- Outlier detection
- Batch-variable assessment
"""

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from sklearn.decomposition import PCA


PROJECT_ROOT = Path(__file__).resolve().parent.parent

expression_path = PROJECT_ROOT / "data" / "GSE45827_expression.csv"
metadata_path = PROJECT_ROOT / "data" / "GSE45827_metadata.csv"

FIGURES_DIR = PROJECT_ROOT / "figures"
FIGURES_DIR.mkdir(exist_ok=True)


def save_figure(name):
    """Save and close the current figure into the figures directory."""
    plt.tight_layout()
    plt.savefig(
        FIGURES_DIR / name,
        dpi=300,
        bbox_inches="tight",
    )
    plt.close()


expression = pd.read_csv(expression_path)
metadata = pd.read_csv(metadata_path).set_index("sample_id")

expression_values = expression.iloc[:, 1:]

sample_ids = expression.columns[1:]

# ---------------------------------------------------------------------------
# 1. Expression value distribution
# ---------------------------------------------------------------------------

values = expression_values.to_numpy().flatten()

plt.figure(figsize=(8, 5))
plt.hist(values, bins=100)

plt.xlabel("Expression value")
plt.ylabel("Frequency")
plt.title("GSE45827 expression value distribution")"""
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
plt.title("GSE45827 expression value distribution")"""
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

<<<<<<< HEAD
save_figure("expression_distribution.png")
=======
plt.tight_layout()
plt.savefig(
    PROJECT_ROOT / "figures" / "expression_distribution.png",
    dpi=300,
    bbox_inches="tight",
)
plt.show()
plt.close()
<<<<<<< HEAD

>>>>>>> main

# ---------------------------------------------------------------------------
# 2. Sample-wise expression distributions
# ---------------------------------------------------------------------------

plt.figure(figsize=(12, 5))
plt.boxplot(
    expression_values.values,
    showfliers=False,
)

plt.xlabel("Samples")
plt.ylabel("Expression")
plt.title("Expression distribution across samples")

<<<<<<< HEAD
save_figure("sample_expression_boxplot.png")
=======
plt.tight_layout()
plt.savefig(
    PROJECT_ROOT / "figures" / "sample_expression_boxplot.png",
    dpi=300,
    bbox_inches="tight",
)
plt.show()
plt.close()
>>>>>>> main

# ---------------------------------------------------------------------------
# 3. Outlier detection (per-sample mean z-score)
# ---------------------------------------------------------------------------

sample_means = expression_values.mean(axis=0)

mean_z = (sample_means - sample_means.mean()) / sample_means.std()

outlier_flags = mean_z.abs() > 3

print("\nOutlier detection (per-sample mean expression, |z| > 3):")
print(f"Samples flagged: {int(outlier_flags.sum())}")

if outlier_flags.any():
    print(mean_z[outlier_flags])

# ---------------------------------------------------------------------------
# 4. Low-variance filtering
# ---------------------------------------------------------------------------

X = expression_values.T

print("\nExpression matrix:", X.shape)

gene_variance = X.var(axis=0)

variance_threshold = 0.1

X = X.loc[:, gene_variance > variance_threshold]

print(
    f"After variance filtering "
    f"(threshold > {variance_threshold}): {X.shape}"
)

# ---------------------------------------------------------------------------
# 5. Principal component analysis by subtype
# ---------------------------------------------------------------------------

pca = PCA(n_components=2)

X_pca = pca.fit_transform(X)

pc1_variance = pca.explained_variance_ratio_[0]
pc2_variance = pca.explained_variance_ratio_[1]

print("\nVariance explained:")
print(f"PC1: {pc1_variance * 100:.1f}%")
print(f"PC2: {pc2_variance * 100:.1f}%")
print(f"PC1 + PC2: {(pc1_variance + pc2_variance) * 100:.1f}%")

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

save_figure("pca_subtypes.png")

# ---------------------------------------------------------------------------
# 6. Batch-variable assessment
# ---------------------------------------------------------------------------

pca_df["batch"] = metadata.loc[sample_ids, "batch"].astype(str)
pca_df["hybridization"] = metadata.loc[
    sample_ids, "hybridization"
].astype(str)

plt.figure(figsize=(8, 6))

for batch in sorted(pca_df["batch"].unique()):

    subset = pca_df[pca_df["batch"] == batch]

    plt.scatter(
        subset["PC1"],
        subset["PC2"],
        label=f"batch {batch}",
    )

plt.xlabel(f"PC1 ({pc1_variance * 100:.1f}%)")
plt.ylabel(f"PC2 ({pc2_variance * 100:.1f}%)")
plt.title("PCA of GSE45827 breast cancer samples, coloured by batch")
plt.legend(title="Batch")

save_figure("pca_batch.png")

print("\nBatch variable confounding with subtype:")
print(
    pd.crosstab(
        metadata.loc[sample_ids, "subtype"],
        metadata.loc[sample_ids, "batch"],
    )
)
<<<<<<< HEAD

print("\nHybridization variable confounding with subtype:")
print(
    pd.crosstab(
        metadata.loc[sample_ids, "subtype"],
        metadata.loc[sample_ids, "hybridization"],
    )
)

# ---------------------------------------------------------------------------
# 7. Heatmap of top variable genes
# ---------------------------------------------------------------------------

N_TOP_GENES = 50

top_genes = gene_variance.nlargest(N_TOP_GENES).index

top_expression = expression_values.loc[top_genes]

# Z-score across samples for visualisation
top_expression_z = top_expression.sub(
    top_expression.mean(axis=1), axis=0
).div(
    top_expression.std(axis=1), axis=0
)

# Colour samples by subtype
subtype_palette = {
    "Basal": "#d62728",
    "HER2": "#ff7f0e",
    "Luminal A": "#1f77b4",
    "Luminal B": "#2ca02c",
}

sample_colors = metadata.loc[sample_ids, "subtype"].map(subtype_palette)

plot_df = top_expression_z.T

plot_df.columns = [
    str(col) for col in plot_df.columns
]

sns.clustermap(
    plot_df,
    col_cluster=True,
    row_cluster=True,
    row_colors=sample_colors,
    cmap="vlag",
    center=0,
    figsize=(12, 12),
    xticklabels=False,
    yticklabels=False,
    cbar_kws={"label": "Z-scored expression"},
)

plt.suptitle(
    "Top 50 variable genes across samples",
    y=1.02,
)

save_figure("heatmap_top50_variable_genes.png")

# ---------------------------------------------------------------------------
# 8. Sample correlation heatmap
# ---------------------------------------------------------------------------

sample_corr = X.T.corr()

sample_corr.index = sample_ids
sample_corr.columns = sample_ids

sns.clustermap(
    sample_corr,
    row_colors=sample_colors,
    col_colors=sample_colors,
    cmap="vlag",
    center=0,
    vmin=-1,
    vmax=1,
    figsize=(12, 12),
    xticklabels=False,
    yticklabels=False,
    cbar_kws={"label": "Pearson correlation"},
)

plt.suptitle(
    "Sample-to-sample expression correlation",
    y=1.02,
)

save_figure("sample_correlation.png")
=======
plt.show()
plt.close()

plt.tight_layout()
plt.savefig(
    PROJECT_ROOT / "figures" / "expression_distribution.png",
    dpi=300,
    bbox_inches="tight",
)
plt.show()
plt.close()
=======
>>>>>>> main



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
plt.close()


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
plt.close()
<<<<<<< HEAD
>>>>>>> main
=======

plt.tight_layout()
plt.savefig(
    PROJECT_ROOT / "figures" / "expression_distribution.png",
    dpi=300,
    bbox_inches="tight",
)
plt.show()
plt.close()



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
plt.close()


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
plt.close()
>>>>>>> main
