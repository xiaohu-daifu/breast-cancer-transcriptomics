# Breast Cancer Transcriptomic Subtyping & Immunotherapy Biomarker Discovery

A reproducible computational biology project investigating conserved transcriptomic
programs underlying breast cancer molecular subtypes and their relationship to
immune biology and potential immunotherapy response.

## Research Question

This project aims to answer the following biological question:

> **Which transcriptomic programs consistently distinguish breast cancer molecular
> subtypes across independent cohorts, and what do they reveal about disease
> biology, immune microenvironment, and potential immunotherapy response?**

The project progressively addresses four levels of biological questions:

1. Which transcriptomic programs distinguish Luminal A, Luminal B, HER2-enriched,
   and Basal-like/TNBC breast cancers?

2. What pathways and biological processes explain these transcriptomic differences?

3. Do the same transcriptomic programs appear across independent datasets, or are
   they specific to a particular cohort or technical platform?

4. Can validated transcriptomic programs help explain differences in the immune
   microenvironment and identify TNBCs with molecular features associated with
   immunotherapy response?


---

## Approach

The project is structured as a series of progressively more sophisticated
analyses:

```text
Breast cancer transcriptomic datasets
                │
                ▼
            _Stage 1_
       Data QC & exploration
                │
                ▼
            _Stage 2_
    Differential expression
                │
                ▼
            _Stage 3_
      Pathway interpretation
                │
                ▼
            _Stage 4_
Reproduce Thorsson 2018 (immune subtypes C1–C6)
                │
                ▼
            _Stage 5_
     Cross-dataset validation
                │
                ▼
            _Stage 6_
       Predictive modelling
                │
                ▼
            _Stage 7_
       Robustness analysis
                │
                ▼
            _Stage 8_
      Model interpretation
                │
                ▼
            _Stage 9_
 Immune & immunotherapy-related analysis

---

## Current status

- **Stages 1–3 complete** (data QC, differential expression, pathway
  enrichment) on the GSE45827 cohort — see `stage01/`, `stage02/`, `stage03/`.
- **Stage 4 in progress** — reproducing Thorsson et al. 2018 (immune subtypes
  C1–C6) with cross-cohort validation in METABRIC — see `stage04/`.

## Preprint goal

Stage 4 is being developed toward a preprint:

> **"Reproduction and cross-cohort validation of the TCGA immune subtypes
> (C1–C6) in breast cancer."**

The goal is to re-derive the C1–C6 immune subtypes from raw TCGA-BRCA
expression, quantify concordance with the published labels, cross-validate in
METABRIC (an independent, microarray cohort), and map the immune subtypes onto
PAM50 and prognosis.

### Note on the Stage 4 paper choice

The roadmap originally suggested reproducing TCGA (2012) or PAM50 (2009). We
instead selected **Thorsson et al. 2018** (*Immunity*), which classifies tumors
into six immune subtypes (C1–C6) from TCGA data. Candidate alternatives
(Bagaev et al. 2021; Sammut et al. 2022) were rejected after verification:
Bagaev is a commercial (BostonGene) product with paywalled, non-deposited data;
Sammut's primary data is EGA-controlled. Thorsson is fully open (free full text
and supplementary, TCGA data) and directly feeds the project's immune-focus.