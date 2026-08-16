#!/usr/bin/env Rscript
# =============================================================================
# 04_plotting.R -- paper-style figures and cell-proportion tables
# =============================================================================
# Usage: Rscript scripts/04_plotting.R [sample_id1,sample_id2,...]
#
# Produces per sample:
#   outputs/<sample_id>/figures/UMAP_CellType.pdf
#   outputs/<sample_id>/figures/DotPlot_Grayscale_PaperStyle.pdf
#   outputs/<sample_id>/tables/Cell_Proportions_<sample_id>.csv
#   Day21_vHIO / Day21_vHCO additionally:
#     outputs/<sample_id>/figures/Featureplot_SMC_pericyte.pdf
# Depends on: 03_annotation.R.

here::i_am("scripts/04_plotting.R")
source(here::here("scripts", "00_config.R"))

for (sid in parse_sample_args()) {
  row <- get_sample(sid)
  message("===> Processing sample: ", sid)

  in_rds <- rds_path(sid, "03_annotated.rds")
  stopifnot(file.exists(in_rds))
  # Create the per-sample output directories (idempotent)
  for (d in c(dirname(rds_path(sid, "checkpoint.rds")),
              dirname(fig_path(sid, "figure.pdf")),
              dirname(tab_path(sid, "table.csv")))) {
    dir.create(d, showWarnings = FALSE, recursive = TRUE)
  }
  obj <- readRDS(in_rds)

  # 1) Final annotated UMAP with the per-sample color map
  p1 <- DimPlot(obj, reduction = "umap", group.by = "Cell_Type_annotation",
                cols = my_colors[[sid]], label = FALSE, pt.size = 2,
                repel = TRUE) +
    NoAxes() +
    theme(plot.title = element_blank()) +
    ggtitle(NULL) +
    guides(color = guide_legend(ncol = 1, override.aes = list(size = 15),
                                keyheight = 5))
  print(p1)
  ggsave(fig_path(sid, "UMAP_CellType.pdf"), p1, width = 14, height = 12)

  # 2) Grayscale dot plot (paper style: gray90 -> black)
  dot_markers <- if (!is.null(dotplot_markers[[sid]])) {
    dotplot_markers[[sid]]
  } else {
    paper_markers[[sid]]
  }
  dot_markers <- unique(dot_markers)
  dot_markers <- dot_markers[dot_markers %in% rownames(obj)]
  # Seurat 5.4's CellsByIdentities() duplicates cells whose identity is NA
  # (x[NA == level] yields NA placeholders), which crashes DotPlot() with
  # "duplicate 'row.names'". Label unannotated clusters explicitly for display.
  dot_ids <- obj$Cell_Type_annotation
  dot_ids <- factor(dot_ids, levels = c(dotplot_order[[sid]], "Unannotated"))
  dot_ids[is.na(dot_ids)] <- "Unannotated"
  obj$dotplot_clusters <- droplevels(factor(
    dot_ids, levels = rev(c(dotplot_order[[sid]], "Unannotated"))))
  Idents(obj) <- obj$dotplot_clusters  # DotPlot reads Idents(); keep them NA-free
  p2 <- DotPlot(obj, features = dot_markers, group.by = "dotplot_clusters",
                dot.scale = 8, dot.min = 0.01) +
    scale_color_gradient(low = "gray90", high = "black") +
    theme_bw(base_size = 12, base_family = "sans") +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", linewidth = 0.8, fill = NA),
      panel.background = element_rect(fill = "white"),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,
                                 color = "black", face = "italic", size = 11),
      axis.text.y = element_text(color = "black", face = "bold", size = 11),
      axis.title = element_blank(),
      legend.position = "right",
      legend.title = element_text(size = 10, face = "bold"),
      legend.text = element_text(size = 9)
    )
  print(p2)
  ggsave(fig_path(sid, "DotPlot_Grayscale_PaperStyle.pdf"), p2,
         width = row$dotplot_width, height = 4, dpi = 300)

  # 3) Cell-type proportions table
  cell_counts <- table(obj$Cell_Type_annotation)
  prop_df <- data.frame(
    Stage = row$stage,
    Cell_Type = names(cell_counts),
    Count = as.numeric(cell_counts),
    Percentage = as.numeric(prop.table(cell_counts) * 100)
  )
  prop_df <- prop_df[order(prop_df$Percentage, decreasing = TRUE), ]
  prop_df$Percentage <- round(prop_df$Percentage, 2)
  print(prop_df)
  write.csv(prop_df, tab_path(sid, paste0("Cell_Proportions_", sid, ".csv")),
            row.names = FALSE)

  # 4) SMC / pericyte feature plot (Day21_vHIO and Day21_vHCO only)
  if (sid %in% c("Day21_vHIO", "Day21_vHCO")) {
    p_all <- FeaturePlot(obj, features = smc_genes, ncol = 2, pt.size = 0.5,
                         order = TRUE, min.cutoff = "q10")
    print(p_all)
    ggsave(fig_path(sid, "Featureplot_SMC_pericyte.pdf"), p_all,
           width = 8, height = 8)
  }
}

message("===> 04_plotting.R finished")
