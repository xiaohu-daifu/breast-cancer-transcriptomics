"""
Download TCGA-BRCA (PanCancer Atlas) and METABRIC data from cBioPortal
for the Stage 4 reproduction of Thorsson et al. 2018 (immune subtypes).

Data pulled:
- Clinical (patient + sample) for both cohorts (PAM50 subtype, survival, TMB).
- mRNA expression for both cohorts.

The C1-C6 immune-subtype gold-standard labels are NOT on cBioPortal; they are
obtained separately from the Thorsson 2018 supplementary (Europe PMC PMC5982584).
"""

from pathlib import Path

import pandas as pd
import requests


CBIO = "https://www.cbioportal.org/api"

PROJECT_ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = PROJECT_ROOT / "data"
DATA_DIR.mkdir(exist_ok=True)

STUDIES = {
    "tcga_brca": {
        "study_id": "brca_tcga_pan_can_atlas_2018",
        "expression_profile": "brca_tcga_pan_can_atlas_2018_rna_seq_v2_mrna",
        "expression_sample_list": "brca_tcga_pan_can_atlas_2018_rna_seq_v2_mrna",
    },
    "metabric": {
        "study_id": "brca_metabric",
        "expression_profile": "brca_metabric_mrna",
        "expression_sample_list": "brca_metabric_mrna",
    },
}


def fetch_clinical(study_id, clinical_type):
    """Fetch a clinical table (PATIENT or SAMPLE) as a DataFrame."""
    url = f"{CBIO}/studies/{study_id}/clinical-data"
    params = {
        "clinicalDataType": clinical_type,
        "projection": "DETAILED",
    }
    res = requests.get(url, params=params, timeout=60)
    res.raise_for_status()
    rows = res.json()

    # Pivot wide (attributeId -> value) keyed by patientId/sampleId.
    key = "patientId" if clinical_type == "PATIENT" else "sampleId"
    df = pd.DataFrame(rows)
    if df.empty:
        return df
    wide = df.pivot_table(
        index=key,
        columns="clinicalAttributeId",
        values="value",
        aggfunc="first",
    )
    return wide.reset_index()


def fetch_expression(profile_id, sample_list_id):
    """
    Fetch the full mRNA expression matrix as a genes x samples DataFrame.

    NOTE: this is a large download (~20k genes x ~1000 samples for BRCA,
    ~2.5k samples for METABRIC). Run once; results are cached to disk.
    """
    url = f"{CBIO}/molecular-profiles/{profile_id}/molecular-data/fetch"
    params = {"projection": "SUMMARY"}
    body = {
        "sampleListId": sample_list_id,
    }
    res = requests.post(url, params=params, json=body, timeout=3600)
    res.raise_for_status()
    rows = res.json()

    df = pd.DataFrame(rows)
    matrix = df.pivot_table(
        index="entrezGeneId",
        columns="sampleId",
        values="value",
    )
    return matrix


def main():
    for name, cfg in STUDIES.items():

        out_dir = DATA_DIR / name
        out_dir.mkdir(exist_ok=True)

        study_id = cfg["study_id"]

        print(f"\n=== {name} ({study_id}) ===")

        for ctype in ["PATIENT", "SAMPLE"]:
            clinical = fetch_clinical(study_id, ctype)
            out_path = out_dir / f"clinical_{ctype.lower()}.csv"
            clinical.to_csv(out_path, index=False)
            print(f"  clinical {ctype}: {clinical.shape} -> {out_path}")

        expression = fetch_expression(
            cfg["expression_profile"],
            cfg["expression_sample_list"],
        )
        out_path = out_dir / "expression.csv.gz"
        expression.to_csv(out_path, compression="gzip")
        print(f"  expression: {expression.shape} -> {out_path}")


if __name__ == "__main__":
    main()
