# Data

## Source

All primary data in this project are derived from a single public GEO dataset:

- **Accession:** GSE45827
- **Title:** Expression data from Breast cancer subtypes
- **Source:** NCBI Gene Expression Omnibus (GEO)
- **Platform:** Affymetrix Human Genome U133 Plus 2.0 Array
- **URL:** https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE45827

## Download

The raw series matrix file is **not** committed to version control. To reproduce the analysis, download it from GEO:

```bash
wget -O data/GSE45827_series_matrix.txt.gz \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE45nnn/GSE45827/matrix/GSE45827_series_matrix.txt.gz"
```

Alternatively, download the raw `.CEL` archive for full pre-processing:

```bash
wget -O data/GSE45827_RAW.tar \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE45nnn/GSE45827/suppl/GSE45827_RAW.tar"
```

## Files

| File | Description | Status |
|---|---|---|
| `GSE45827_series_matrix.txt.gz` | Raw GEO series matrix (downloaded) | gitignored |
| `GSE45827_RAW.tar` | Raw Affymetrix `.CEL` archive (downloaded) | gitignored |
| `GSE45827_metadata.csv` | Cleaned sample metadata (subtype, batch, hybridization) | committed |
| `GSE45827_expression.csv` | Processed expression matrix (130 samples × 29,873 probes) | gitignored |

## Metadata columns

`GSE45827_metadata.csv` contains one row per tumour sample:

| Column | Description |
|---|---|
| `sample_id` | GEO sample accession (e.g. `GSM1116084`) |
| `subtype` | Molecular subtype (`Basal`, `HER2`, `Luminal A`, `Luminal B`) |
| `batch` | Processing batch (integer 1–10) |
| `hybridization` | Hybridization label (-1, 1, 2, 3) |

**Note:** the `batch` and `hybridization` variables are confounded with subtype (see `stage01/readme.md` §7.4), and this caveat must be carried through downstream analysis.
