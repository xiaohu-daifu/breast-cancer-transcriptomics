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

**Study design:** Primary invasive breast cancer samples collected at surgery before treatment. The original cohort contains:

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

Extracted GEO sample accession numbers:

`GSM1116084` → `GSM1116238`

Extracted molecular subtype annotations from the GEO sample characteristics.

The relevant molecular subtypes are:

- Basal
- HER2
- Luminal A
- Luminal B

---

## 4. Metadata Cleaning

Cleaned the extracted metadata by:

- Removing GEO quotation marks from sample IDs and subtype labels
- Removing irrelevant sample annotations
- Excluding samples without one of the four molecular subtype labels
- Verifying sample IDs are unique
- Checking for missing subtype annotations

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

Extracted the expression matrix from the GEO Series Matrix file and restricted it to the 130 samples in the analysis cohort.

The resulting matrix contains:

- **29,873 unique probes**
- **130 breast cancer samples**
- **0 missing expression values**

The expression matrix contains Affymetrix probe IDs in `ID_REF` and sample-level expression values in the remaining columns.

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

### Planned analyses

- Expression distribution plots
- Sample-level quality control
- Principal component analysis (PCA)
- Sample clustering
- Subtype-level expression patterns
- Investigation of potential batch effects and other technical variation

### Key observations

*To be completed after exploratory analysis.*

### Figures

*Figures will be added here as the analysis progresses.*

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
- [ ] Perform expression QC
- [ ] Perform PCA
- [ ] Investigate potential batch effects
- [ ] Complete exploratory analysis