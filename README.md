# vHIO_vHLPO_scRNAseq Analysis

Single-cell RNA-seq analysis of vascularized human intestinal organoids
(vHIO / vHCO) and vascularized human lung organoids (vHLPO) at days 3, 7
and 21, accompanying the study *Co-development of mesoderm and endoderm
lineages* (Miao et al., Cell, 2025).

## Repository structure

    ├── README.md               # this file
    ├── LICENSE                 # MIT license
    ├── .gitignore
    ├── renv.lock               # R package version lockfile
    ├── vHIO_vHLPO_scRNAseq.Rproj
    ├── data/
    │   ├── README.md           # GEO download instructions (GSE250399)
    │   ├── metadata_demo.csv   # lightweight per-sample metadata
    │   └── raw/                # downloaded h5 matrices (git-ignored)
    ├── scripts/
    │   ├── 00_config.R         # shared parameters, annotation maps, helpers
    │   ├── 01_preprocessing.R  # QC, cell-cycle scoring, SCTransform
    │   ├── 02_clustering.R     # PCA, UMAP, Leiden clustering, marker export
    │   ├── 03_annotation.R     # cell-type annotation (incl. Day 21 refinements)
    │   └── 04_plotting.R       # paper-style figures and proportion tables
    ├── notebooks/              # optional interactive documentation
    └── outputs/                # analysis products (git-ignored)

Note: the template's `03_differential.R` step corresponds to the marker
analysis in this project; `FindAllMarkers` runs inside `02_clustering.R`,
so no separate differential-expression script is needed.

## System requirements

- macOS / Linux (untested on Windows)
- R 4.5.2
- ~20 GB free disk space and >= 32 GB RAM recommended (Day 21 samples are
  the largest)

Leiden clustering (Seurat's `algorithm = 4`) uses the `leidenbase` R
package, which is part of the locked environment; the scripts fall back
to Louvain clustering if it is unavailable.

Package versions are locked in `renv.lock` (generated with renv from a
clean R 4.5.2 session). To restore the exact environment:

```r
install.packages("renv")
renv::restore()   # run from the repository root
```

## Setup

1. Clone the repository and open `vHIO_vHLPO_scRNAseq.Rproj` in RStudio
   (or start R from the repository root).
2. Download the seven `filtered_feature_bc_matrix.h5` files listed in
   `data/README.md` from GEO (GSE250399) into `data/raw/`.
3. `renv::restore()` to install the locked package versions.

## Analysis workflow

Run the scripts in order; each script supports an optional comma-separated
sample list (`Rscript scripts/0X.R Day3_B1,Day3_B3`) and defaults to all
seven samples. The samples belong to two organoid projects:

| Project | Samples |
|---|---|
| vHLPO (vascularized human lung organoids) | `Day3_B1`, `Day7_vAFG_B1`, `Day21_vHLPO` |
| vHIO / vHCO (vascularized human intestinal organoids) | `Day3_B3`, `Day7_vMHG_B3`, `Day21_vHIO`, `Day21_vHCO` |

```bash
Rscript scripts/01_preprocessing.R    # QC, cell-cycle scoring, SCTransform
Rscript scripts/02_clustering.R       # PCA, UMAP, Leiden clustering, markers
Rscript scripts/03_annotation.R       # cell-type annotation
Rscript scripts/04_plotting.R         # paper figures and proportion tables
```

Each script writes a checkpoint rds under `outputs/<sample_id>/seurat_obj/`,
so the pipeline can be restarted from any module.

### Reproducing the manuscript figures

All manuscript figures are produced by `scripts/04_plotting.R`:

| Manuscript figure | `04_plotting.R` section | Output file |
|---|---|---|
| UMAPs, Figure 5a,b,d (vHLPO samples: `Day3_B1`, `Day7_vAFG_B1`, `Day21_vHLPO`) | 1. Annotated UMAP (`DimPlot`, `group.by = "Cell_Type_annotation"`) | `outputs/<sample_id>/figures/UMAP_CellType.pdf` |
| UMAPs, Figure 6a–d (vHIO/vHCO samples: `Day3_B3`, `Day7_vMHG_B3`, `Day21_vHIO`, `Day21_vHCO`) | 1. Annotated UMAP | `outputs/<sample_id>/figures/UMAP_CellType.pdf` |
| Cell-type proportions, Figure 5d and Figure 6e | 3. Cell-type proportions table | `outputs/<sample_id>/tables/Cell_Proportions_<sample_id>.csv` |
| Dot plots, Figure 5e–g and Figure 6F–I | 2. Grayscale dot plot (`DotPlot`) | `outputs/<sample_id>/figures/DotPlot_Grayscale_PaperStyle.pdf` |
| SMC / pericyte feature plots (Day 21 vHIO and vHCO) | 4. SMC / pericyte feature plot | `outputs/<sample_id>/figures/Featureplot_SMC_pericyte.pdf` |

Marker genes shown in the dot plots, the cell-type color maps and the
cell-type orders are defined per sample in `scripts/00_config.R`
(`paper_markers` / `dotplot_markers`, `my_colors`, `cell_order`,
`dotplot_order`). Positive cluster markers exported by
`scripts/02_clustering.R` are written to
`outputs/<sample_id>/tables/Markers_Leiden_Res<res>.csv`.

## Data availability

Raw count matrices are deposited in the Gene Expression Omnibus (GEO)
under accession **GSE250399** (MIAME-compliant):

https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE250399

Per-sample accessions are listed in `data/README.md` (GSM7978632,
GSM7978634, GSM7978636, GSM7978639, GSM7978640, GSM7978641, GSM7978642).
The h5 files themselves are git-ignored because they exceed GitHub's file
size limits.

## Code availability

This repository is archived on Zenodo with DOI
**10.5281/zenodo.21966164** (https://doi.org/10.5281/zenodo.21966164).
Please cite both the GitHub repository and the Zenodo record in the
references of the manuscript.

## License

This code is released under the MIT license (see `LICENSE`).

## Citation

If you use this code, please cite the associated paper (*Co-development of
mesoderm and endoderm lineages*, Miao et al., Cell, 2025) and the Zenodo
archive of this repository (https://doi.org/10.5281/zenodo.21966164).

