# Breast Cancer Transcriptomics

## Stage 1 — Data Acquisition, Cleaning & Exploration

### Goal

Acquire, clean, and explore breast cancer transcriptomic data from independent public cohorts as the foundation for downstream statistical analysis, biological interpretation, and machine learning.

---

## 1. Dataset Acquisition

### GEO: GSE45827

**Dataset:** `GSE45827_series_matrix.txt.gz`

**Source:** NCBI Gene Expression Omnibus (GEO)

**Accession:** GSE45827

**Title:** Expression data from Breast cancer subtypes

**Platform:** Affymetrix Human Genome U133 Plus 2.0 Array

**Experiment type:** Expression profiling by array

**Study design:** Primary invasive breast cancer samples collected at surgery before treatment.

The original cohort contains:

- 41 triple-negative/basal samples
- 30 HER2 samples
- 29 Luminal A samples
- 30 Luminal B samples
- 11 normal tissue samples
- 14 cell lines

### Original study design

> In a cohort study of primary invasive breast cancer (41 TN, 30 HER2, 29 Luminal A and 30 Luminal B) as well as 11 normal tissues samples and 14 cell lines, a tumor specimen at surgery before any patient treatment was obtained. Total RNA was extracted from all samples and the whole transcriptome was quantified with Affymetrix U133 Plus 2.0 Chips.

---

## 2. Raw Data Inspection

Inspected the GEO Series Matrix file to understand its structure and identify:

- Sample-level metadata
- GEO sample accession numbers
- Molecular subtype annotations
- Expression matrix boundaries
- Non-tumour samples and samples without molecular subtype annotations

---

## 3. Sample Metadata Extraction

Extracted GEO sample accession numbers and molecular subtype annotations from the GEO sample characteristics.

The molecular subtypes used for the analysis cohort are:

- Basal
- HER2
- Luminal A
- Luminal B

The metadata extraction and cohort filtering are implemented in:

`stage01/01_extract_metadata.py`

---

## 4. Metadata Cleaning

The extracted metadata were cleaned by:

- Removing GEO quotation marks from sample IDs and subtype labels
- Extracting the molecular subtype from the sample characteristics
- Excluding samples without one of the four molecular subtype labels
- Verifying that sample IDs are unique
- Verifying that all retained samples have valid subtype annotations

---

## 5. Analysis Cohort

After filtering, the analysis cohort contains **130 breast cancer samples**:

| Molecular subtype | Samples |
|---|---:|
| Basal | 41 |
| HER2 | 30 |
| Luminal A | 29 |
| Luminal B | 30 |
| **Total** | **130** |

The cleaned metadata are stored in:

`data/GSE45827_metadata.csv`

---

## 6. Expression Data

The expression matrix was extracted from the GEO Series Matrix file and restricted to the 130 samples in the analysis cohort.

The resulting matrix contains:

- **29,873 Affymetrix probe features**
- **130 breast cancer samples**
- **0 missing expression values**

The expression matrix contains Affymetrix probe IDs in `ID_REF` and sample-level expression values in the remaining columns.

The expression extraction and validation are implemented in:

`stage01/02_extract_expression.py`

The processed expression matrix is stored in:

`data/GSE45827_expression.csv`

### Expression/metadata validation

Verified that:

- All 130 metadata samples are present in the expression matrix
- Expression matrix columns match the metadata sample IDs exactly
- Probe IDs are unique
- No expression values are missing

---

## 7. Exploratory Data Analysis

Exploratory analysis was performed on the cleaned expression matrix.

### 7.1 Expression value distribution

The overall expression value distribution was examined using a histogram.

Expression values were predominantly distributed within an approximately **2–11** range, with a tail extending toward approximately 15.

Combined with the sample-level distributions, this is consistent with a processed microarray expression matrix rather than raw sequencing counts.

No additional log transformation was applied during this stage.

### 7.2 Sample-level expression distributions

Boxplots were used to compare expression distributions across the 130 samples.

The sample distributions were broadly consistent, with similar medians and interquartile ranges across samples. No sample showed an obvious global distributional shift requiring immediate exclusion.

### 7.3 Low-variance filtering

For exploratory PCA, genes/probes with very low variance across samples were removed.

A variance threshold of **0.1** was applied:

| | Features | Samples |
|---|---:|---:|
| Before filtering | 29,873 | 130 |
| After filtering | 28,086 | 130 |

This filtering was used specifically to reduce uninformative low-variance features for exploratory PCA.

### 7.4 Principal Component Analysis

PCA was performed on the variance-filtered expression matrix.

The first two principal components explained:

| Component | Variance explained |
|---|---:|
| PC1 | **17.4%** |
| PC2 | **6.0%** |
| PC1 + PC2 | **23.4%** |

The PCA revealed clear transcriptional structure corresponding to breast cancer molecular subtype.

In particular, Basal samples showed strong separation from the other subtypes along PC1, while Luminal A, Luminal B, and HER2 samples showed varying degrees of overlap.

This provides initial evidence that the expression matrix captures biologically meaningful molecular subtype structure.

The QC and PCA analysis are implemented in:

`stage01/03_expression_qc.py`

---

## 8. Stage 1 Interpretation

The initial exploratory analysis indicates that:

1. The GEO Series Matrix can be successfully parsed into a clean expression matrix and corresponding subtype metadata.
2. The 130-sample breast cancer cohort has complete expression data.
3. Sample-level expression distributions are broadly consistent across the cohort.
4. Low-variance filtering removes a relatively small proportion of features (1,787 of 29,873).
5. PCA reveals strong molecular structure associated with breast cancer subtype, particularly the separation of Basal samples.

These findings support proceeding to **differential expression analysis** to identify genes and transcriptional programs associated with molecular subtype.

---

## 9. Figures


### Expression value distribution

![Expression value distribution](figures/expression_distribution.png)

### Sample-level expression distributions

![Sample expression distributions](figures/sample_expression_boxplot.png)

### PCA by molecular subtype

![PCA of GSE45827 breast cancer samples](figures/pca_subtypes.png)

**Figure 3.** Principal component analysis of 130 GSE45827 breast cancer
samples using variance-filtered expression features. Samples are coloured
according to molecular subtype. PC1 explains 17.4% of the variance and
PC2 explains 6.0%.

---

## Stage 1 Status

- [x] Acquire GEO dataset
- [x] Inspect raw GEO structure
- [x] Extract sample IDs
- [x] Extract molecular subtype metadata
- [x] Clean metadata
- [x] Define analysis cohort
- [x] Save cleaned metadata
- [x] Extract expression matrix
- [x] Validate expression/metadata sample matching
- [x] Check missing expression values
- [x] Check probe ID uniqueness
- [x] Examine expression value distribution
- [x] Examine sample-level expression distributions
- [x] Perform low-variance filtering for exploratory PCA
- [x] Perform PCA
- [x] Examine PCA by molecular subtype


---

## Next Stage

### Stage 2 — Differential Expression

The next stage will identify genes differentially expressed between breast cancer molecular subtypes.

Planned comparisons include:

- Basal vs Luminal A
- Basal vs HER2
- Basal/TNBC vs non-Basal subtypes

The analysis will introduce:

- Log2 fold change
- Statistical significance testing
- Multiple-testing correction
- False discovery rate (FDR)
- Volcano plots
- Differential expression tables
- Heatmaps of significant genes