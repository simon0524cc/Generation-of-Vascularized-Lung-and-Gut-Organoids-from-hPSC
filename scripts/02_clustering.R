#!/usr/bin/env Rscript
# =============================================================================
# 02_clustering.R -- PCA, UMAP, Leiden clustering scan and marker export
# =============================================================================
# Usage: Rscript scripts/02_clustering.R [sample_id1,sample_id2,...]
#
# Produces per sample:
#   outputs/<sample_id>/seurat_obj/02_clustered.rds
#   outputs/<sample_id>/figures/UMAP_Leiden_Res_<res>.pdf (resolution scan)
#   outputs/<sample_id>/figures/Comparison_Resolutions.pdf
#   outputs/<sample_id>/tables/Markers_Leiden_Res<res>.csv (final resolutions)
#   outputs/<sample_id>/figures/UMAP.pdf (only where save_umap_pdf is TRUE)
# Depends on: 01_preprocessing.R.

here::i_am("scripts/02_clustering.R")
source(here::here("scripts", "00_config.R"))

for (sid in parse_sample_args()) {
  row <- get_sample(sid)
  message("===> Processing sample: ", sid)

  in_rds <- rds_path(sid, "01_preprocessed.rds")
  stopifnot(file.exists(in_rds))
  # Create the per-sample output directories (idempotent)
  for (d in c(dirname(rds_path(sid, "checkpoint.rds")),
              dirname(fig_path(sid, "figure.pdf")),
              dirname(tab_path(sid, "table.csv")))) {
    dir.create(d, showWarnings = FALSE, recursive = TRUE)
  }
  obj <- readRDS(in_rds)

  # 1) PCA on the SCT assay, then UMAP and SNN graph
  obj <- RunPCA(obj, npcs = row$npcs, seed.use = 42, assay = "SCT",
                verbose = FALSE)
  obj <- RunUMAP(obj, dims = row$dims_umap[[1]], reduction = "pca",
                 verbose = FALSE)
  obj <- FindNeighbors(obj, dims = row$dims_neighbor[[1]], reduction = "pca",
                       verbose = FALSE)
  snn_graph <- grep("snn$", names(obj@graphs), value = TRUE)[1]

  # 2) Resolution scan: Leiden (algorithm = 4) with Louvain fallback
  #    (the paper requires Leiden; without python leidenalg we fall back).
  #    random.seed = 1 reproduces the original analysis: the original
  #    scripts did not pass random.seed, and Seurat resets its default
  #    of 0 to 1 before running Leiden.
  plot_list <- list()
  for (res in row$res_test[[1]]) {
    obj <- tryCatch({
      FindClusters(obj, graph.name = snn_graph, resolution = res,
                   algorithm = 4, random.seed = 1, verbose = FALSE)
    }, error = function(e) {
      message("Leiden failed (leidenalg not installed?), falling back to Louvain.")
      FindClusters(obj, graph.name = snn_graph, resolution = res,
                   algorithm = 1, verbose = FALSE)
    })
    num_clust <- length(unique(Idents(obj)))
    p <- DimPlot(obj, reduction = "umap", label = TRUE, label.size = 5) +
      ggtitle(paste0("Leiden | Res: ", res, " | Clusters: ", num_clust)) +
      theme(plot.title = element_text(hjust = 0.5))
    ggsave(fig_path(sid, paste0("UMAP_Leiden_Res_", res, ".pdf")),
           p, width = 7, height = 6)
    plot_list[[as.character(res)]] <- p
  }
  final_plot <- wrap_plots(plot_list, ncol = 3)
  ggsave(fig_path(sid, "Comparison_Resolutions.pdf"), final_plot,
         width = 18, height = row$comparison_height)

  # 3) Final resolutions: set clusters and export positive markers
  #    (random.seed = 1, matching the original analysis -- see step 2)
  for (res in row$final_res[[1]]) {
    obj <- FindClusters(obj, graph.name = snn_graph, resolution = res,
                        algorithm = 4, random.seed = 1, verbose = FALSE)
    all_markers <- FindAllMarkers(obj, only.pos = TRUE, min.pct = 0.25,
                                  logfc.threshold = 0.25)
    write_csv(all_markers, tab_path(sid, paste0("Markers_Leiden_Res", res, ".csv")))
  }

  # 4) Optional unannotated UMAP (only where the original script saved one)
  #    NOTE: ggsave() must receive a ggplot object; the original scripts
  #    passed the Seurat object directly, which fails in recent ggplot2.
  if (row$save_umap_pdf) {
    p_umap <- DimPlot(obj, reduction = "umap", label = TRUE, label.size = 5)
    ggsave(fig_path(sid, "UMAP.pdf"), p_umap, width = 10, height = 7)
  }

  # 5) Checkpoint (keeps the clustering state of the last final resolution,
  #    matching the original scripts)
  out <- rds_path(sid, "02_clustered.rds")
  dir.create(dirname(out), showWarnings = FALSE, recursive = TRUE)
  saveRDS(obj, out)
  message("===> Saved: ", out)
}

message("===> 02_clustering.R finished")
