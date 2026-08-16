#!/usr/bin/env Rscript
# =============================================================================
# 03_annotation.R -- cell-type annotation with sample-specific refinements
# =============================================================================
# Usage: Rscript scripts/03_annotation.R [sample_id1,sample_id2,...]
#
# Produces per sample:
#   outputs/<sample_id>/seurat_obj/03_annotated.rds
#   outputs/<sample_id>/figures/Annotation_Check_DotPlot.pdf
# Sample-specific steps (order matches the original scripts):
#   Day21_vHLPO : CDH17/CDX2 split of Intestinal-like epithelium
#   Day21_vHCO  : immune module scoring -> "Immune cell",
#                 per-gene immune FeaturePlots, SMC subclustering
#   Day21_vHIO  : pericyte subclustering -> "Smooth muscle cell"
# Depends on: 02_clustering.R.

here::i_am("scripts/03_annotation.R")
source(here::here("scripts", "00_config.R"))

for (sid in parse_sample_args()) {
  row <- get_sample(sid)
  message("===> Processing sample: ", sid)

  in_rds <- rds_path(sid, "02_clustered.rds")
  stopifnot(file.exists(in_rds))
  # Create the per-sample output directories (idempotent)
  for (d in c(dirname(rds_path(sid, "checkpoint.rds")),
              dirname(fig_path(sid, "figure.pdf")),
              dirname(tab_path(sid, "table.csv")))) {
    dir.create(d, showWarnings = FALSE, recursive = TRUE)
  }
  obj <- readRDS(in_rds)

  # 1) Marker check dot plot (RdBu palette, before annotation)
  markers <- unique(paper_markers[[sid]])
  markers <- markers[markers %in% rownames(obj)]
  DefaultAssay(obj) <- "RNA"
  obj <- NormalizeData(obj, verbose = FALSE)
  p_check <- DotPlot(obj, features = markers, group.by = "seurat_clusters") +
    RotatedAxis() +
    scale_color_distiller(palette = "RdBu") +
    ggtitle("Check Clusters against Paper Markers")
  print(p_check)
  ggsave(fig_path(sid, "Annotation_Check_DotPlot.pdf"), p_check,
         width = 12, height = 8)

  # 2) Base annotation from the cluster map
  ann <- cluster_annotations[[sid]]
  obj$Cell_Type_annotation <- unname(ann[as.character(obj$seurat_clusters)])

  # 3) Sample-specific refinements
  if (sid == "Day21_vHLPO") {
    # Split Intestinal-like epithelium: cells expressing CDH17 or CDX2 stay;
    # the remainder is demoted to Epithelial progenitor-2.
    expr_threshold <- 0
    gene_exprs <- FetchData(obj, vars = c("CDH17", "CDX2"))
    current_intestinal <- rownames(obj@meta.data)[
      obj$Cell_Type_annotation == "Intestinal-like epithelium"]
    true_intestinal <- rownames(gene_exprs)[
      rownames(gene_exprs) %in% current_intestinal &
        (gene_exprs$CDH17 > expr_threshold | gene_exprs$CDX2 > expr_threshold)]
    progenitor2_cells <- setdiff(current_intestinal, true_intestinal)
    obj$Cell_Type_annotation[progenitor2_cells] <- "Epithelial progenitor-2"
    message("\n--- Intestinal epithelium split report ---")
    message("Cells originally annotated as Intestinal-like epithelium: ",
            length(current_intestinal))
    message("Kept as Intestinal-like epithelium: ", length(true_intestinal))
    message("Demoted to Epithelial progenitor-2: ", length(progenitor2_cells), "\n")
  }

  if (sid == "Day21_vHCO") {
    # (a) Per-gene immune marker FeaturePlots (executed in the original script)
    for (gene in c("CD163", "CD84", "SPI1")) {
      p <- FeaturePlot(obj, features = gene, order = TRUE)
      ggsave(fig_path(sid, paste0(gene, "_feature_plot.pdf")), p,
             width = 6, height = 5, dpi = 300)
    }

    # (b) Immune module scoring: cells above the threshold become "Immune cell"
    immune_genes_use <- immune_genes[immune_genes %in% rownames(obj)]
    message("Genes used for the immune score: ",
            paste(immune_genes_use, collapse = ", "))
    if ("Immune_Score1" %in% colnames(obj@meta.data)) {
      obj$Immune_Score1 <- NULL  # remove a stale score column before re-scoring
    }
    obj <- AddModuleScore(obj, features = list(immune_genes_use),
                          name = "Immune_Score")
    immune_barcodes <- rownames(obj@meta.data)[
      obj$Immune_Score1 > immune_threshold]
    message("Immune cells identified at threshold > ", immune_threshold,
            ": ", length(immune_barcodes))
    obj$Cell_Type_annotation[
      rownames(obj@meta.data) %in% immune_barcodes] <- "Immune cell"

    # (c) Verification UMAP pair (saving kept commented out, as in the original)
    Idents(obj) <- "Cell_Type_annotation"
    p_score_umap <- FeaturePlot(obj, features = "Immune_Score1",
                                reduction = "umap", pt.size = 0.8) +
      scale_color_viridis_c(option = "plasma") +
      ggtitle("A. Immune Module Score Distribution") +
      theme(plot.title = element_text(hjust = 0.5))
    p_highlight_umap <- DimPlot(obj, group.by = "Cell_Type_annotation",
                                cells.highlight = immune_barcodes,
                                cols.highlight = "#E31A1C",
                                cols = "grey90", pt.size = 1) +
      ggtitle(paste0("B. Highlighted Immune Cells (n = ",
                     length(immune_barcodes), ")")) +
      theme(plot.title = element_text(hjust = 0.5)) +
      NoLegend()
    combined_plot <- p_score_umap | p_highlight_umap
    print(combined_plot)
    # ggsave(fig_path(sid, "Immune_Cells_UMAP_Verification.pdf"),
    #        combined_plot, width = 12, height = 5.5)
  }

  if (sid %in% c("Day21_vHIO", "Day21_vHCO")) {
    # Pericyte subclustering -> Smooth muscle cell.
    # Resolution and cluster picks differ between the two samples.
    smc_res <- if (sid == "Day21_vHIO") 0.4 else 0.6
    smc_clusters <- if (sid == "Day21_vHIO") c("2", "3") else c("0", "3", "4", "5")
    Idents(obj) <- "Cell_Type_annotation"
    peri_obj <- subset(obj, idents = "Pericyte")
    peri_obj <- FindVariableFeatures(peri_obj, selection.method = "vst",
                                     nfeatures = 2000)
    peri_obj <- ScaleData(peri_obj)
    peri_obj <- RunPCA(peri_obj, features = VariableFeatures(peri_obj))
    peri_obj <- FindNeighbors(peri_obj, dims = 1:15)
    peri_obj <- FindClusters(peri_obj, resolution = smc_res)
    FeaturePlot(peri_obj, features = c("TAGLN", "CNN1"))
    VlnPlot(peri_obj, features = c("TAGLN", "CNN1"), pt.size = 0)
    p_umap <- DimPlot(peri_obj, reduction = "umap", label = FALSE,
                      pt.size = 2, repel = TRUE)
    print(p_umap)
    smc_barcodes_subcluster <- rownames(peri_obj@meta.data)[
      peri_obj$seurat_clusters %in% smc_clusters]
    obj$Cell_Type_annotation <- as.character(obj$Cell_Type_annotation)
    obj$Cell_Type_annotation[smc_barcodes_subcluster] <- "Smooth muscle cell"
  }

  # 4) Factor levels for plotting and final Idents
  obj$Cell_Type_annotation <- factor(obj$Cell_Type_annotation,
                                     levels = cell_order[[sid]])
  Idents(obj) <- "Cell_Type_annotation"
  table(Idents(obj))

  out <- rds_path(sid, "03_annotated.rds")
  dir.create(dirname(out), showWarnings = FALSE, recursive = TRUE)
  saveRDS(obj, out)
  message("===> Saved: ", out)
}

message("===> 03_annotation.R finished")
