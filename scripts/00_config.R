# =============================================================================
# 00_config.R -- per-sample parameter table, annotation maps and shared helpers
# =============================================================================
# Source this file before running any pipeline module. It defines:
#   * sample_table : one row per sample with the parameters used in the paper
#   * cluster_annotations / cell_order / dotplot_order / my_colors /
#     paper_markers / dotplot_markers : per-sample annotation maps
#   * helper functions for reading, QC filtering and output paths
#
# All paths are resolved relative to the project root (the directory
# containing the .Rproj file) via the here package. Do not use setwd().

suppressPackageStartupMessages({
  library(here)
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(readr)
})

options(future.globals.maxSize = 2000 * 1024^2)
set.seed(1234)  # global seed, identical to the original per-sample scripts

# -----------------------------------------------------------------------------
# 1. Per-sample parameter table
# -----------------------------------------------------------------------------
# dims_umap / dims_neighbor / res_test / final_res are list columns
# (access with sample_table$dims_umap[[i]]).

sample_table <- data.frame(
  sample_id = c(
    "Day3_B1", "Day3_B3", "Day7_vAFG_B1", "Day7_vMHG_B3",
    "Day21_vHLPO", "Day21_vHIO", "Day21_vHCO"
  ),
  gsm = c(
    "GSM7978632", "GSM7978634", "GSM7978636", "GSM7978639",
    "GSM7978640", "GSM7978641", "GSM7978642"
  ),
  h5_file = c(
    "GSM7978632_Day3_B1_filtered_feature_bc_matrix.h5",
    "GSM7978634_Day3_B3_filtered_feature_bc_matrix.h5",
    "GSM7978636_Day7_vAFG_B1_filtered_feature_bc_matrix.h5",
    "GSM7978639_Day7_vMHG_B3_filtered_feature_bc_matrix.h5",
    "GSM7978640_Day21_vHLPO_filtered_feature_bc_matrix.h5",
    "GSM7978641_Day21_vHIO_filtered_feature_bc_matrix.h5",
    "GSM7978642_Day21_vHCO_filtered_feature_bc_matrix.h5"
  ),
  stage = c("Day 3", "Day 3", "Day 7", "Day 7", "Day 21", "Day 21", "Day 21"),
  max_features = c(10000, 10000, 11000, 10000, 10000, 10000, 10000),
  npcs = c(100, 100, 100, 100, 100, 200, 200),
  # width of the grayscale dot plot (inches); differs in the original scripts
  dotplot_width = c(12, 12, 12, 8, 12, 12, 12),
  # height of the resolution-comparison figure (inches)
  comparison_height = c(12, 6, 12, 12, 12, 12, 12),
  # whether the original script saved an unannotated UMAP.pdf after clustering
  save_umap_pdf = c(TRUE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE),
  stringsAsFactors = FALSE
)

sample_table$dims_umap <- I(list(
  1:30, 1:30, 1:30, 1:30, 1:30, 1:200, 1:200
))
sample_table$dims_neighbor <- I(list(
  1:50, 1:50, 1:50, 1:50, 1:100, 1:200, 1:200
))
sample_table$res_test <- I(list(
  c(0.6, 0.7, 0.8, 0.9, 1.0),
  c(0.6, 0.7, 0.8, 0.9, 1.0),
  c(0.6, 0.7, 0.8, 0.9, 1.0),
  c(0.6, 0.7, 0.8, 0.9, 1.0),
  c(0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0),
  c(0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0),
  c(0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0)
))
sample_table$final_res <- I(list(
  c(0.7), c(0.8), c(0.7, 0.8), c(0.7, 0.8), c(1.6, 2.0), c(0.7, 2.0), c(1.6)
))

get_sample <- function(sample_id) {
  sample_table[sample_table$sample_id == sample_id, , drop = FALSE]
}

# -----------------------------------------------------------------------------
# 2. Cluster -> cell-type annotation maps (copied verbatim from the originals)
# -----------------------------------------------------------------------------

cluster_annotations <- list(
  Day3_B1 = c(
    "1" = "Meso-endoderm progenitor",
    "2" = "Def. endoderm",
    "3" = "Def. endoderm",
    "4" = "Low-quality",
    "5" = "Lateral plate mesoderm",
    "6" = "Low-quality",
    "7" = "Def. endoderm",
    "8" = "Cardiac Mesoderm",
    "9" = "Meso-endoderm progenitor"
  ),
  Day3_B3 = c(
    "1" = "Meso-endoderm progenitor",
    "2" = "Def. endoderm",
    "3" = "Meso-endoderm progenitor",
    "4" = "Lateral plate mesoderm",
    "5" = "Lateral plate mesoderm",
    "6" = "Cardiac Mesoderm",
    "7" = "Def. endoderm",
    "8" = "Low-quality",
    "9" = "Cardiac Mesoderm",
    "10" = "Low-quality"
  ),
  Day7_vAFG_B1 = c(
    "1" = "Epithelium",
    "2" = "Epithelium",
    "3" = "Epithelium",
    "4" = "Endothelium",
    "5" = "Epithelium",
    "6" = "Mesenchymal cell",
    "7" = "non-differentiated cell",
    "8" = "Epithelium"
  ),
  Day7_vMHG_B3 = c(
    "1" = "Mesenchymal cell",
    "2" = "Mesenchymal cell",
    "3" = "Epithelium",
    "4" = "Epithelium",
    "5" = "Endothelium",
    "6" = "Mesenchymal cell",
    "7" = "Endothelium",
    "8" = "Mesenchymal cell",
    "9" = "Epithelium",
    "10" = "Epithelium",
    "11" = "Mesenchymal cell"
  ),
  Day21_vHLPO = c(
    "1" = "Pro. lung fibroblast",
    "2" = "Lung fibroblast",
    "3" = "Neural crest",
    "4" = "Lung fibroblast",
    "5" = "Intestinal-like epithelium",
    "6" = "Low-quality",
    "7" = "Lung fibroblast",
    "8" = "Epithelial progenitor-2",
    "9" = "Endothelium",
    "10" = "Lung fibroblast",
    "11" = "Lung fibroblast",
    "12" = "Pro. lung fibroblast",
    "13" = "Low-quality",
    "14" = "Lung fibroblast",
    "15" = "Epithelial progenitor-1",
    "16" = "Pro. lung fibroblast",
    "17" = "Low-quality",
    "18" = "Placode-lineage",
    "19" = "Pro. lung fibroblast",
    "20" = "Pericyte-like",
    "21" = "Endothelium",
    "22" = "Epithelial progenitor-2"
  ),
  Day21_vHIO = c(
    "1" = "Fibroblast",
    "2" = "Fibroblast",
    "3" = "Fibroblast",
    "4" = "Fibroblast",
    "5" = "Fibroblast",
    "6" = "Fibroblast",
    "7" = "Low-quality",
    "8" = "Pericyte",
    "9" = "Pro. fibroblast",
    "10" = "Epithelium",
    "11" = "Fibroblast",
    "12" = "Fibroblast",
    "13" = "Lymphatic EC",
    "14" = "Fibroblast",
    "15" = "Pro. fibroblast",
    "16" = "Endothelium",
    "17" = "Pericyte",
    "18" = "Epithelium",
    "19" = "Fibroblast",
    "20" = "Pro. fibroblast"
  ),
  Day21_vHCO = c(
    "1" = "Fibroblast",
    "2" = "Fibroblast",
    "3" = "Fibroblast",
    "4" = "Fibroblast",
    "5" = "Pericyte",
    "6" = "Fibroblast",
    "7" = "Pro. fibroblast",
    "8" = "Endothelium",
    "9" = "Epithelium",
    "10" = "Fibroblast",
    "11" = "Fibroblast",
    "12" = "Fibroblast",
    "13" = "Low-quality",
    "14" = "Epithelium",
    "15" = "Lymphatic EC",
    "16" = "Fibroblast",
    "17" = "Fibroblast",
    "18" = "Pro. fibroblast"
  )
)

# -----------------------------------------------------------------------------
# 3. Cell-type factor orders (UMAP legend) and dot-plot row orders
# -----------------------------------------------------------------------------

cell_order <- list(
  Day3_B1 = c("Cardiac Mesoderm", "Meso-endoderm progenitor",
              "Lateral plate mesoderm", "Def. endoderm", "Low-quality"),
  Day3_B3 = c("Cardiac Mesoderm", "Meso-endoderm progenitor",
              "Lateral plate mesoderm", "Def. endoderm", "Low-quality"),
  Day7_vAFG_B1 = c("Epithelium", "Endothelium", "Mesenchymal cell",
                   "non-differentiated cell"),
  Day7_vMHG_B3 = c("Epithelium", "Endothelium", "Mesenchymal cell"),
  Day21_vHLPO = c("Lung fibroblast", "Pro. lung fibroblast", "Endothelium",
                  "Epithelial progenitor-1", "Epithelial progenitor-2",
                  "Placode-lineage", "Pericyte-like", "Neural crest",
                  "Intestinal-like epithelium", "Low-quality"),
  Day21_vHIO = c("Fibroblast", "Pro. fibroblast", "Endothelium",
                 "Lymphatic EC", "Epithelium", "Smooth muscle cell",
                 "Pericyte", "Low-quality"),
  Day21_vHCO = c("Fibroblast", "Pro. fibroblast", "Endothelium",
                 "Lymphatic EC", "Epithelium", "Smooth muscle cell",
                 "Pericyte", "Immune cell", "Low-quality")
)

dotplot_order <- list(
  Day3_B1 = c("Cardiac Mesoderm", "Meso-endoderm progenitor",
              "Lateral plate mesoderm", "Def. endoderm", "Low-quality"),
  Day3_B3 = c("Cardiac Mesoderm", "Meso-endoderm progenitor",
              "Lateral plate mesoderm", "Def. endoderm", "Low-quality"),
  Day7_vAFG_B1 = c("non-differentiated cell", "Epithelium", "Endothelium",
                   "Mesenchymal cell"),
  Day7_vMHG_B3 = c("Epithelium", "Endothelium", "Mesenchymal cell"),
  Day21_vHLPO = c("Epithelial progenitor-1", "Epithelial progenitor-2",
                  "Endothelium", "Lung fibroblast", "Pro. lung fibroblast",
                  "Pericyte-like", "Neural crest", "Placode-lineage",
                  "Intestinal-like epithelium", "Low-quality"),
  Day21_vHIO = c("Epithelium", "Endothelium", "Lymphatic EC", "Fibroblast",
                 "Pro. fibroblast", "Pericyte", "Smooth muscle cell",
                 "Low-quality"),
  Day21_vHCO = c("Epithelium", "Endothelium", "Lymphatic EC", "Fibroblast",
                 "Pro. fibroblast", "Pericyte", "Smooth muscle cell",
                 "Immune cell", "Low-quality")
)

# -----------------------------------------------------------------------------
# 4. Per-sample color maps (identical to the original scripts)
# -----------------------------------------------------------------------------

my_colors <- list(
  Day3_B1 = c(
    "Meso-endoderm progenitor"   = "#AAC7AD",
    "Lateral plate mesoderm"     = "#FE891E",
    "Cardiac Mesoderm"           = "#DB3494",
    "Def. endoderm"              = "#78A1CD",
    "Low-quality"                = "#BFBFC0"
  ),
  Day3_B3 = c(
    "Meso-endoderm progenitor"   = "#AAC7AD",
    "Lateral plate mesoderm"     = "#FE891E",
    "Cardiac Mesoderm"           = "#DB3494",
    "Def. endoderm"              = "#78A1CD",
    "Low-quality"                = "#BFBFC0"
  ),
  Day7_vAFG_B1 = c(
    "non-differentiated cell"    = "#FCE0D2",
    "Mesenchymal cell"           = "#EF4D24",
    "Epithelium"                 = "#3C7AB6",
    "Endothelium"                = "#36704D"
  ),
  Day7_vMHG_B3 = c(
    "Mesenchymal cell"           = "#EF4D24",
    "Epithelium"                 = "#3C7AB6",
    "Endothelium"                = "#36704D"
  ),
  Day21_vHLPO = c(
    "Endothelium"                = "#36704D",
    "Pro. lung fibroblast"       = "#E89FC6",
    "Lung fibroblast"            = "#C7212A",
    "Epithelial progenitor-1"    = "#C9EBFD",
    "Epithelial progenitor-2"    = "#0B51AE",
    "Pericyte-like"              = "#D1817E",
    "Neural crest"               = "#532C8B",
    "Placode-lineage"            = "#F1B379",
    "Intestinal-like epithelium" = "#1A1A1C",
    "Low-quality"                = "#BFBFC0"
  ),
  Day21_vHIO = c(
    "Epithelium"                 = "#3C7AB6",
    "Endothelium"                = "#36704D",
    "Pro. fibroblast"            = "#E89FC6",
    "Fibroblast"                 = "#C7212A",
    "Smooth muscle cell"         = "#F7E89E",
    "Pericyte"                   = "#D1817E",
    "Lymphatic EC"               = "#65A73E",
    "Low-quality"                = "#BFBFC0"
  ),
  Day21_vHCO = c(
    "Epithelium"                 = "#3C7AB6",
    "Endothelium"                = "#36704D",
    "Pro. fibroblast"            = "#E89FC6",
    "Fibroblast"                 = "#C7212A",
    "Smooth muscle cell"         = "#F7E89E",
    "Pericyte"                   = "#D1817E",
    "Lymphatic EC"               = "#65A73E",
    "Immune cell"                = "#109992",
    "Low-quality"                = "#BFBFC0"
  )
)

# -----------------------------------------------------------------------------
# 5. Paper marker genes (Fig. S3B order, filtered to genes present in data at
#    run time). paper_markers feeds the Annotation_Check_DotPlot; the grayscale
#    dot plot uses dotplot_markers[[sid]] when defined, otherwise paper_markers.
# -----------------------------------------------------------------------------

paper_markers <- list(
  Day3_B1 = c(
    "NANOG", "POU5F1", "TBXT", "MESP1", "NPPB", "MYL7", "TNNC1", "MYL9",
    "NODAL", "MIXL1", "EOMES", "GSC", "HAND1", "BMP4", "FOXF1", "ISL1",
    "FOXA1", "FOXA2", "SOX17", "HHEX", "TOP2A", "MKI67"
  ),
  Day3_B3 = c(
    "NANOG", "POU5F1", "TBXT", "MESP1", "NPPB", "MYL7", "TNNC1", "MYL9",
    "NODAL", "MIXL1", "EOMES", "GSC", "HAND1", "BMP4", "FOXF1", "ISL1",
    "FOXA1", "FOXA2", "SOX17", "HHEX", "TOP2A", "MKI67"
  ),
  Day7_vAFG_B1 = c(
    "NANOG", "POU5F1", "SOX2", "IRX1", "OTX2", "ISL1", "SIX1", "PDX1",
    "ONECUT1", "ALB", "TTR", "PECAM1", "CDH5", "PDGFRA", "COL1A1"
  ),
  Day7_vMHG_B3 = c(
    "NANOG", "POU5F1", "SOX2", "CDX2", "PECAM1", "CDH5", "PDGFRA", "COL1A1"
  ),
  Day21_vHLPO = c(
    "NANOG", "POU5F1", "EPCAM", "KRT19", "NKX2-1", "SOX9", "TPPP3", "LEF1",
    "PECAM1", "CDH5", "HPGD", "KIT", "PDGFRA", "TBX4", "FOXF1", "MKI67",
    "TOP2A", "PDGFRB", "LAMC3", "FOXD3", "SOX10", "PAX3", "SIX1",
    "CDH17", "CDX2"
  ),
  Day21_vHIO = c(
    "NANOG", "POU5F1", "EPCAM", "CDX2", "APOA1", "LCN15", "PECAM1", "CDH5",
    "LYVE1", "PROX1", "STAB2", "PDGFRA", "COL1A1", "MKI67", "TOP2A",
    "PDGFRB", "CSPG4", "TAGLN", "CNN1"
  ),
  Day21_vHCO = c(
    "NANOG", "POU5F1", "EPCAM", "CDX2", "APOA1", "LCN15", "PECAM1", "CDH5",
    "LYVE1", "PROX1", "STAB2", "PDGFRA", "COL1A1", "MKI67", "TOP2A",
    "PDGFRB", "CSPG4", "TAGLN", "CNN1"
  )
)

# Extra markers appended to the grayscale dot plot only (Day21_vHCO adds
# immune markers after the immune re-annotation; the annotation-check dot plot
# in that script still uses paper_markers alone).
dotplot_markers <- list(
  Day21_vHCO = c(
    "NANOG", "POU5F1", "EPCAM", "CDX2", "APOA1", "LCN15", "PECAM1", "CDH5",
    "LYVE1", "PROX1", "STAB2", "PDGFRA", "COL1A1", "MKI67", "TOP2A",
    "PDGFRB", "CSPG4", "TAGLN", "CNN1", "CD84", "SPI1"
  )
)

# -----------------------------------------------------------------------------
# 6. Fixed thresholds and gene sets shared across samples
# -----------------------------------------------------------------------------

qc_defaults <- list(min_features = 500, max_umi = 100000, max_mt = 25)

immune_genes <- c("CD14", "CD163", "CD84", "SPI1")
immune_threshold <- 0.01  # low threshold on purpose (immune cells are rare)
smc_genes <- c("PDGFRB", "CSPG4", "TAGLN", "CNN1")

# -----------------------------------------------------------------------------
# 7. Helper functions (extracted from the original scripts)
# -----------------------------------------------------------------------------

# Read one 10x h5 matrix into a Seurat object with QC metadata columns.
read_10x_single <- function(h5_path, sample_id) {
  mtx <- Read10X_h5(h5_path)
  # Handle the possible list structure of multi-modal h5 files
  if (is.list(mtx)) {
    mtx <- if ("Gene Expression" %in% names(mtx)) mtx[["Gene Expression"]] else mtx[[1]]
  }
  obj <- CreateSeuratObject(counts = mtx, project = sample_id, assay = "RNA",
                            min.cells = 3, min.features = 200)
  obj$sample <- sample_id
  obj[["percent.mt"]]   <- PercentageFeatureSet(obj, pattern = "^MT-")
  obj[["percent.ribo"]] <- PercentageFeatureSet(obj, pattern = "^RP[SL]")
  # Prefix cell barcodes with the sample id
  RenameCells(obj, new.names = paste0(sample_id, "_", colnames(obj)))
}

# QC filter: 500 <= nFeature <= max_features, nCount <= 100000, MT <= 25%.
qc_filter <- function(x, min_features = 500, max_features = 10000,
                      max_umi = 100000, max_mt = 25) {
  subset(x, subset = nFeature_RNA >= min_features & nFeature_RNA <= max_features &
           nCount_RNA   <= max_umi       & percent.mt <= max_mt)
}

# -----------------------------------------------------------------------------
# 8. Output paths (project-root relative via here)
# -----------------------------------------------------------------------------

rds_path <- function(sample_id, fname) {
  here::here("outputs", sample_id, "seurat_obj", fname)
}

fig_path <- function(sample_id, fname) {
  here::here("outputs", sample_id, "figures", fname)
}

tab_path <- function(sample_id, fname) {
  here::here("outputs", sample_id, "tables", fname)
}

# Parse the command line: "Rscript scripts/0X.R Day3_B1,Day3_B3" runs a subset;
# no argument runs all samples in sample_table.
parse_sample_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) > 0) strsplit(args[1], ",")[[1]] else sample_table$sample_id
}
