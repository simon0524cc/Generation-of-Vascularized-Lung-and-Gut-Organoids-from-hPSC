# Data

## Source

Raw count matrices were generated in the study *Co-emergence of
mesoderm and endoderm lineages* (Miao et al., Cell, 2025) and are
deposited in GEO under accession **GSE250399**:

- Series: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE250399
- Direct download: https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE250399&format=file

## Download instructions

The analysis scripts expect one h5 count matrix per sample in
`data/raw/` (git-ignored; the files exceed GitHub's file-size limit).
Download the `filtered_feature_bc_matrix.h5` files listed below and
place them in `data/raw/` before running `scripts/01_preprocessing.R`.
The `raw_feature_bc_matrix.h5` files are not needed for the analyses
in this repository.

| Sample ID    | GEO sample | File to download |
|--------------|------------|------------------|
| Day3_B1      | GSM7978632 | GSM7978632_Day3_B1_filtered_feature_bc_matrix.h5 |
| Day3_B3      | GSM7978634 | GSM7978634_Day3_B3_filtered_feature_bc_matrix.h5 |
| Day7_vAFG_B1 | GSM7978636 | GSM7978636_Day7_vAFG_B1_filtered_feature_bc_matrix.h5 |
| Day7_vMHG_B3 | GSM7978639 | GSM7978639_Day7_vMHG_B3_filtered_feature_bc_matrix.h5 |
| Day21_vHLPO  | GSM7978640 | GSM7978640_Day21_vHLPO_filtered_feature_bc_matrix.h5 |
| Day21_vHIO   | GSM7978641 | GSM7978641_Day21_vHIO_filtered_feature_bc_matrix.h5 |
| Day21_vHCO   | GSM7978642 | GSM7978642_Day21_vHCO_filtered_feature_bc_matrix.h5 |

## Demo metadata

`metadata_demo.csv` is a lightweight per-sample table used to select
and label samples in the scripts. It contains no expression data.
