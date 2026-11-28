# Stage 4 — Reproduce Thorsson et al. 2018 (The Immune Landscape of Cancer)

## Goal

Reproduce the TCGA pan-cancer immune-subtype classification (C1–C6) from
Thorsson et al., *Immunity* 2018, and extend it to breast cancer with
cross-cohort validation in METABRIC : the foundation for a preprint.

## Preprint contribution

> **"Reproduction and cross-cohort validation of the TCGA immune subtypes (C1–C6) in breast cancer."**

- Re-derive C1–C6 from raw TCGA-BRCA expression and quantify concordance.
- Cross-validate C1–C6 in METABRIC (independent, microarray platform).
- Map C1–C6 onto PAM50 and prognosis in BRCA.

## Paper reference

Thorsson V, Gibbs DL, Brown SD, et al. "The Immune Landscape of Cancer."
*Immunity* 48(4):812–830 (2018). DOI 10.1016/j.immuni.2018.03.023.
Free full text + supplementary: Europe PMC **PMC5982584**.

## Data sources

| Data | Source | Status |
|---|---|---|
| TCGA-BRCA expression (RSEM) | cBioPortal `brca_tcga_pan_can_atlas_2018` | ✅ open |
| TCGA-BRCA clinical (PAM50, OS, TMB) | cBioPortal (same study) | ✅ open |
| METABRIC expression + clinical | cBioPortal `brca_metabric` | ✅ open |
| **C1–C6 gold-standard labels** | Thorsson supplementary (PMC5982584) | ⚠️ see note |

### C1–C6 labels (gold standard)

The per-sample immune-subtype assignments (C1–C6) are **not** in cBioPortal.
FireBrowse (which used to host them) is currently down. The reliable open
source is the **Thorsson 2018 supplementary tables** (free via Europe PMC
`PMC5982584`), which annotate the ~10,000 TCGA samples by immune subtype.
This needs to be extracted and mapped to the BRCA sample IDs.

## Method (to reproduce)

1. Compute **ssGSEA** scores for the ~160 immune signatures (from supplementary)
   across TCGA-BRCA.
2. Aggregate into the **5 immune determinants**.
3. Cluster into 6 subtypes; compare to published C1–C6 (adjusted Rand index).
4. Reproduce key associations (PAM50, survival, TMB/aneuploidy).
5. Extend: project onto METABRIC; map to PAM50; test survival.

## Scripts

```text
01_download_data.py   # cBioPortal clinical + expression (this stage)
```

## Reproducibility

```text
python stage04/01_download_data.py
```

(Subsequent scripts for signature scoring, clustering, and validation are
added in later phases.)

## Next steps

- Extract the ~160 immune signatures + C1–C6 labels from the Thorsson
  supplementary.
- Add signature-scoring (`02_signatures.py`) and clustering (`03_subtypes.py`).
