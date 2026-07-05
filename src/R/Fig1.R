# Fig1.R — assemble Figure 1: global transcriptome remodelling landscape
# Run from repo root under UTF-8 locale. Reads results/ produced by 01-07 + S1-S8.
source("src/R/00_config.R")
suppressPackageStartupMessages({
  library(ggplot2); library(patchwork); library(dplyr); library(scales)
  library(ragg); library(ggridges)
})

save_full <- function(p, path_stem, w, h) {
  grDevices::cairo_pdf(paste0(path_stem, ".pdf"),
                       width = w / 25.4, height = h / 25.4, family = "Helvetica")
  print(p); dev.off()
  svglite::svglite(paste0(path_stem, ".svg"), width = w / 25.4, height = h / 25.4)
  print(p); dev.off()
  ragg::agg_tiff(paste0(path_stem, ".tiff"),
                 width = w / 25.4, height = h / 25.4, units = "in", res = 600)
  print(p); dev.off()
  invisible(NULL)
}

{
  dir_f1 <- file.path(DIR_RESULTS, "figures", "F1_global_landscape")
  load(file.path(DIR_RESULTS, "normalized", "normalized_data.RData"))

  tp_levels    <- c("0h", "3h", "6h", "9h", "12h")
  stage_labels <- c("GV", "GVBD", "MI-6", "MI-9", "MII")

  mrna_mat <- matrix(NA_real_, nrow = length(gene_ids_filt), ncol = 5,
                     dimnames = list(gene_ids_filt, tp_levels))
  te_mat   <- matrix(NA_real_, nrow = length(gene_ids_filt), ncol = 5,
                     dimnames = list(gene_ids_filt, tp_levels))
  for (i in seq_along(tp_levels)) {
    tp      <- tp_levels[i]
    pm      <- pair_map[pair_map$time_point == tp, ]
    m_cols  <- pm$transcriptome_id
    te_cols <- grep(paste0("^", tp, "_rep"), colnames(te_A), value = TRUE)
    mrna_mat[, i] <- rowMeans(log1p(norm_trans_A[gene_ids_filt, m_cols, drop=FALSE]), na.rm=TRUE)
    te_mat[, i]   <- rowMeans(te_A[gene_ids_filt, te_cols, drop=FALSE], na.rm=TRUE)
  }
  mrna_delta <- sweep(mrna_mat, 1, mrna_mat[, "0h"], "-")
  te_delta   <- sweep(te_mat,   1, te_mat[,   "0h"], "-")

  # F1A (already saved; re-build from data)
  stage_df2 <- data.frame(
    x = 1:5,
    stage = c("GV", "GVBD", "MI (6h)", "MI (9h)", "MII"),
    time  = c("0 h", "3 h", "6 h", "9 h", "12 h"),
    n_rep = c("n=5", "n=5", "n=5", "n=3", "n=5"),
    color = unname(col_timepoint), y_circ = 1.0
  )
  arrow_df2 <- data.frame(x1 = 1:4 + 0.18, x2 = 2:5 - 0.18, y1 = 1.0, y2 = 1.0)
  omics_df2 <- data.frame(
    x = rep(1:5, 2),
    y = c(rep(0.45, 5), rep(0.10, 5)),
    omics = rep(c("T","TE"), each = 5)
  )
  p1a2 <- ggplot() +
    annotate("segment", x = 0.5, xend = 0.5, y = 1.35, yend = 1.15,
             arrow = arrow(length = unit(0.1, "cm"), type = "closed"),
             colour = "grey40", linewidth = 0.7) +
    annotate("text", x = 0.5, y = 1.46, label = "HCG\n(t=0h)",
             size = 2.2, hjust = 0.5, colour = "grey40", family = "Helvetica") +
    geom_segment(data = arrow_df2,
                 aes(x = x1, xend = x2, y = y1, yend = y2),
                 arrow = arrow(length = unit(0.09, "cm"), type = "open"),
                 colour = "grey60", linewidth = 0.5) +
    geom_point(data = stage_df2, aes(x = x, y = y_circ, colour = stage),
               size = 9, shape = 16, alpha = 0.85) +
    geom_text(data = stage_df2, aes(x = x, y = y_circ, label = stage),
              size = 2.2, colour = "white", fontface = "bold", family = "Helvetica") +
    geom_text(data = stage_df2, aes(x = x, y = 0.73, label = time),
              size = 2.4, colour = "grey30", family = "Helvetica") +
    geom_text(data = stage_df2, aes(x = x, y = 0.58, label = n_rep),
              size = 2.1, colour = "grey50", family = "Helvetica") +
    annotate("segment", x = 0.7, xend = 5.3, y = 0.45, yend = 0.45,
             colour = col_omics[["Transcriptome"]], linewidth = 0.9) +
    annotate("segment", x = 0.7, xend = 5.3, y = 0.10, yend = 0.10,
             colour = col_omics[["Translatome"]], linewidth = 0.9) +
    geom_point(data = omics_df2, aes(x = x, y = y,
               colour = omics), size = 2.6, shape = 21, fill = "white", stroke = 1.0) +
    annotate("text", x = 0.35, y = 0.45, label = "Transcriptome",
             size = 2.1, hjust = 1, colour = col_omics[["Transcriptome"]], fontface = "bold") +
    annotate("text", x = 0.35, y = 0.10, label = "Translatome",
             size = 2.1, hjust = 1, colour = col_omics[["Translatome"]], fontface = "bold") +
    scale_colour_manual(
      values = c(setNames(stage_df2$color, stage_df2$stage),
                 "T" = col_omics[["Transcriptome"]],
                 "TE" = col_omics[["Translatome"]])
    ) +
    coord_cartesian(xlim = c(-0.5, 6.0), ylim = c(-0.1, 1.6)) +
    labs(title = "Experimental design") +
    theme_void(base_family = "Helvetica") +
    theme(legend.position = "none",
          plot.title = element_text(size = 8, face = "bold", hjust = 0))

  # F1B
  f1b_df <- data.frame(
    tp = factor(stage_labels, levels = stage_labels),
    median = apply(mrna_mat, 2, median, na.rm = TRUE),
    q25 = apply(mrna_mat, 2, quantile, 0.25, na.rm = TRUE),
    q75 = apply(mrna_mat, 2, quantile, 0.75, na.rm = TRUE)
  )
  p1b2 <- ggplot(f1b_df, aes(x = tp, y = median, group = 1)) +
    geom_ribbon(aes(ymin = q25, ymax = q75),
                fill = col_omics[["Transcriptome"]], alpha = 0.18) +
    geom_line(colour = col_omics[["Transcriptome"]], linewidth = 0.9) +
    geom_point(colour = col_omics[["Transcriptome"]], size = 2.2,
               shape = 21, fill = "white", stroke = 1.1) +
    labs(x = "Stage", y = "mRNA (log1p)",
         title = "Global mRNA trajectory",
         subtitle = sprintf("Median ± IQR, n=%d genes", nrow(mrna_mat))) +
    theme_oocyte(base_size = 8)

  # F1C
  f1c_mrna <- data.frame(tp = factor(stage_labels, levels = stage_labels),
    median = apply(mrna_delta, 2, median, na.rm = TRUE),
    q25    = apply(mrna_delta, 2, quantile, 0.25, na.rm = TRUE),
    q75    = apply(mrna_delta, 2, quantile, 0.75, na.rm = TRUE),
    omics  = "Transcriptome")
  f1c_te <- data.frame(tp = factor(stage_labels, levels = stage_labels),
    median = apply(te_delta, 2, median, na.rm = TRUE),
    q25    = apply(te_delta, 2, quantile, 0.25, na.rm = TRUE),
    q75    = apply(te_delta, 2, quantile, 0.75, na.rm = TRUE),
    omics  = "Translatome")
  f1c_df <- rbind(f1c_mrna, f1c_te)
  f1c_df$omics <- factor(f1c_df$omics, levels = c("Transcriptome", "Translatome"))

  p1c2 <- ggplot(f1c_df, aes(x = tp, y = median, group = omics)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60", linewidth = 0.4) +
    geom_ribbon(aes(ymin = q25, ymax = q75, fill = omics), alpha = 0.15, colour = NA) +
    geom_line(aes(colour = omics), linewidth = 0.9) +
    geom_point(aes(colour = omics), size = 2.2, shape = 21, fill = "white", stroke = 1.1) +
    scale_colour_manual(values = col_omics,
                        labels = c(Transcriptome = "mRNA", Translatome = "TE"),
                        name = NULL) +
    scale_fill_manual(values = col_omics, guide = "none") +
    labs(x = "Stage", y = expression(Delta ~ "from GV"),
         title = "mRNA vs TE trajectories",
         subtitle = "Median ± IQR delta from GV") +
    theme_oocyte(base_size = 8) +
    theme(legend.position = "inside",
          legend.position.inside = c(0.03, 0.08),
          legend.justification = c(0, 0),
          legend.key.size = unit(3.5, "mm"))

  # F1D
  set.seed(params$seed)
  samp_g <- sample(rownames(mrna_mat), min(8000L, nrow(mrna_mat)))
  library(ggridges)
  f1d_df <- do.call(rbind, lapply(seq_along(tp_levels), function(i) {
    data.frame(stage = factor(stage_labels[i], levels = rev(stage_labels)),
               mrna  = mrna_mat[samp_g, i])
  }))
  p1d2 <- ggplot(f1d_df, aes(x = mrna, y = stage, fill = stage, colour = stage)) +
    geom_density_ridges(alpha = 0.70, scale = 1.4, linewidth = 0.5,
                        quantile_lines = TRUE, quantiles = 2) +
    scale_fill_manual(values = rev(col_timepoint), guide = "none") +
    scale_colour_manual(values = rev(col_timepoint), guide = "none") +
    labs(x = "mRNA (log1p)", y = NULL,
         title = "mRNA distribution shift",
         subtitle = "Progressive leftward shift") +
    theme_oocyte(base_size = 8)

  # 2-row layout: A (design, full width, short) / B + C + D (equal thirds)
  p1_full <- (p1a2 / (p1b2 | p1c2 | p1d2)) +
    plot_layout(heights = c(0.5, 1)) +
    plot_annotation(
      tag_levels = list(c("A", "B", "C", "D")),
      theme = theme(plot.tag = element_text(size = 10, face = "bold",
                                             family = "Helvetica"))
    ) &
    theme(plot.title = element_text(size = 8, face = "bold"),
          plot.subtitle = element_text(size = 6.2, colour = "grey40"),
          plot.margin = margin(3, 5, 3, 4))

  save_full(p1_full,
            file.path(dir_f1, "F1_main_assembly"),
            w = 183, h = 118)
  cat("F1 main assembly (A+B+C+D) saved\n")
}
