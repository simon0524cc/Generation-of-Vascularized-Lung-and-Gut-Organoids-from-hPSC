#!/usr/bin/env Rscript
# =============================================================================
# 01_preprocessing.R -- read, QC, cell-cycle scoring and SCTransform
# =============================================================================
# Usage:
#   Rscript scripts/01_preprocessing.R                 # all seven samples
#   Rscript scripts/01_preprocessing.R Day3_B1         # one sample
#   Rscript scripts/01_preprocessing.R Day3_B1,Day3_B3 # several samples
#
# Produces: outputs/<sample_id>/seurat_obj/01_preprocessed.rds
# Depends on: 00_config.R and data/raw/<h5_file> (see data/README.md).

here::i_am("scripts/01_preprocessing.R")
source(here::here("scripts", "00_config.R"))

for (sid in parse_sample_args()) {
  row <- get_sample(sid)
  message("===> Processing sample: ", sid)

  # 1) Read 10x data from the local copy of the GEO h5 file
  h5_path <- here::here("data", "raw", row$h5_file)
  stopifnot(file.exists(h5_path))
  obj <- read_10x_single(h5_path, sid)

  # 2) QC filtering (thresholds identical to the original scripts)
  obj <- qc_filter(obj, max_features = row$max_features)
  summary(obj$nFeature_RNA)
  summary(obj$nCount_RNA)
  summary(obj$percent.mt)
  message("QC done, cells remaining: ", ncol(obj))

  # 3) Cell-cycle scoring on log-normalized RNA.
  # The paper regresses cell-cycle scores and MT percentage in SCTransform.
  data("cc.genes.updated.2019", package = "Seurat")
  DefaultAssay(obj) <- "RNA"
  obj <- NormalizeData(obj, verbose = FALSE)
  obj <- CellCycleScoring(obj,
                          s.features   = cc.genes.updated.2019$s.genes,
                          g2m.features = cc.genes.updated.2019$g2m.genes,
                          set.ident = FALSE)

  # 4) SCTransform with regression of cell cycle and mitochondrial content
  obj <- SCTransform(obj,
                     vars.to.regress = c("S.Score", "G2M.Score", "percent.mt"),
                     method = "glmGamPoi",
                     variable.features.n = 3000,
                     conserve.memory = TRUE,
                     verbose = TRUE)

  out <- rds_path(sid, "01_preprocessed.rds")
  dir.create(dirname(out), showWarnings = FALSE, recursive = TRUE)
  saveRDS(obj, out)
  message("===> Saved: ", out)
}

message("===> 01_preprocessing.R finished")
