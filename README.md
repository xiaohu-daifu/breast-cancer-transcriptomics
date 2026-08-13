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
    Published analysis reproduction
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