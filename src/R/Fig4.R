# Fig4.R — assemble Figure 4: translational-engagement compensation
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
  dir_f4  <- file.path(DIR_RESULTS, "figures", "F4_tbr_landscape")
  dir_sup <- file.path(DIR_RESULTS, "supplementary")

  # Re-load p4a (re-run from script data)
  s1  <- read.csv(file.path(dir_sup, "S1_tbr_cs_matrix.csv"), stringsAsFactors = FALSE)
  s6  <- read.csv(file.path(dir_sup, "S6_threshold_sensitivity.csv"), stringsAsFactors = FALSE)
  s4f <- read.csv(file.path(dir_sup, "S4_functional_cpe_features.csv"), stringsAsFactors = FALSE)

  col_quad <- c(
    "Stable / Minor Change" = "#A7A9AC",
    "TE-Only Activation"    = unname(col_omics["Translatome"]),
    "Compensatory"          = "#00857C",
    "Coordinated Clearance" = unname(col_omics["Transcriptome"]),
    "Partial Buffering"     = "#E07B39"
  )

  # F4A scatter (simplified for assembly)
  sc2 <- s1[is.finite(s1$tbr) & is.finite(s1$cs) &
              abs(s1$tbr) < 10 & abs(s1$cs) < 10, ]
  sc2$quadrant <- factor(sc2$quadrant, levels = names(col_quad))

  p4a2 <- ggplot(sc2, aes(x = cs, y = tbr, colour = quadrant)) +
    geom_point(alpha = 0.08, size = 0.3, shape = 16) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey55", linewidth = 0.3) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey55", linewidth = 0.3) +
    scale_colour_manual(values = col_quad, name = NULL,
                         labels = function(x) sub("Partial Buffering", "Partial Compensation", x),
                         guide = guide_legend(override.aes = list(alpha = 0.9, size = 2),
                                               ncol = 2)) +
    labs(x = "CS (Compensation Score)", y = "TBR (log₂)",
         title = "Compensation landscape") +
    theme_oocyte(base_size = 8) +
    theme(legend.position = "bottom",
          legend.text = element_text(size = 5.6),
          legend.key.size = unit(2.6, "mm"),
          legend.margin = margin(0, 0, 0, 0))

  # F4B — all genes with finite TBR (n=19,452); headline TBR>0 = 57.7%
  tbr_all <- s1$tbr[is.finite(s1$tbr)]
  tbr_cap <- pmax(-8, pmin(8, tbr_all))
  pct_pos <- round(100 * mean(tbr_all > 0), 1)
  p4b2 <- ggplot(data.frame(tbr = tbr_cap), aes(x = tbr)) +
    geom_histogram(binwidth = 0.5, fill = col_omics[["Translatome"]], alpha = 0.8,
                   colour = "white", linewidth = 0.2) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40", linewidth = 0.5) +
    annotate("text", x = 6, y = Inf, vjust = 2,
             label = paste0("TBR > 0\n", pct_pos, "%"),
             size = 2.6, colour = col_omics[["Translatome"]], fontface = "bold") +
    labs(x = "TBR (log₂)", y = "Number of genes",
         title = "TBR distribution") +
    theme_oocyte(base_size = 8)

  # F4C
  sc3 <- s1[is.finite(s1$tbr) & is.finite(s1$delta_mrna), ]
  sc3$tbr_cap   <- pmax(-6, pmin(6, sc3$tbr))
  sc3$quadrant  <- factor(sc3$quadrant, levels = names(col_quad))
  m2   <- lm(tbr_cap ~ delta_mrna, data = sc3)
  beta2 <- round(coef(m2)[2], 3)
  r2_2  <- round(summary(m2)$r.squared, 3)

  p4c2 <- ggplot(sc3, aes(x = delta_mrna, y = tbr_cap, colour = quadrant)) +
    geom_point(alpha = 0.10, size = 0.3, shape = 16) +
    geom_smooth(data = sc3, aes(x = delta_mrna, y = tbr_cap),
                method = "lm", formula = y ~ x,
                colour = "grey25", linewidth = 0.9, se = TRUE,
                fill = "grey25", alpha = 0.15, inherit.aes = FALSE) +
    annotate("text", x = 5.8, y = -5.6,
             label = paste0("β=", beta2, " · R²=", r2_2),
             hjust = 1, size = 2.5, colour = "grey25") +
    scale_colour_manual(values = col_quad, guide = "none") +
    scale_x_continuous(limits = c(-6.5, 6.5)) +
    scale_y_continuous(limits = c(-6.5, 6.5)) +
    labs(x = "ΔmRNA (log₂)", y = "TBR (log₂)",
         title = "ΔmRNA weakly tracks\ntranslational compensation") +
    theme_oocyte(base_size = 8)

  # F4D
  quad_n2 <- s1 %>%
    filter(!is.na(quadrant)) %>% count(quadrant) %>%
    mutate(pct = n / sum(n) * 100)
  quad_n2$quadrant <- factor(quad_n2$quadrant,
                              levels = names(col_quad)[names(col_quad) %in% quad_n2$quadrant])
  quad_n2$lab <- sprintf("%s  %.1f%%",
                         sub("Partial Buffering", "Partial Compensation", quad_n2$quadrant),
                         quad_n2$pct)
  p4d2 <- ggplot(quad_n2, aes(x = 1, y = pct, fill = quadrant)) +
    geom_col(width = 0.75, alpha = 0.88) +
    geom_text(aes(label = lab),
              position = position_stack(vjust = 0.5),
              size = 1.95, lineheight = 0.85, colour = "white", fontface = "bold") +
    scale_fill_manual(values = col_quad, guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0.03, 0.03))) +
    coord_cartesian(clip = "off") +
    labs(x = NULL, y = "% genes", title = "Quadrant\nproportions") +
    theme_oocyte(base_size = 8) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

  # F4E (threshold-sensitivity heatmap; S6 columns: cs_threshold / tbr_threshold / fisher_or)
  if ("tbr_threshold" %in% names(s6) && "fisher_or" %in% names(s6)) {
    s6$tbr_threshold <- factor(s6$tbr_threshold, levels = sort(unique(s6$tbr_threshold)))
    s6$cs_threshold  <- factor(s6$cs_threshold,  levels = sort(unique(s6$cs_threshold)))
    p4e2 <- ggplot(s6, aes(x = cs_threshold, y = tbr_threshold, fill = fisher_or)) +
      geom_tile(colour = "white", linewidth = 0.4) +
      geom_text(aes(label = sprintf("%.2f", fisher_or)), size = 2.3, colour = "grey20") +
      scale_fill_gradient2(low = "white", mid = "#C8D8F0",
                            high = col_omics[["Translatome"]], midpoint = 1.5,
                            name = "Fisher\nOR") +
      labs(x = "CS threshold", y = "TBR threshold",
           title = "CPE enrichment\nthreshold robustness") +
      theme_oocyte(base_size = 8) +
      theme(legend.key.size = unit(3.2, "mm"),
            legend.title = element_text(size = 6))
  } else {
    p4e2 <- ggplot() +
      annotate("text", x = 0.5, y = 0.5, label = "F4E: see F4E_threshold_heatmap.pdf",
               size = 3) + theme_void()
  }

  # F4F
  s1$gene_id_clean  <- sub("\\.\\d+$", "", s1$gene_id)
  df_f2 <- merge(s1, s4f, by = "gene_id_clean", all.x = TRUE)
  df_f2 <- df_f2[!is.na(df_f2$n_functional_cpe) & !is.na(df_f2$tbr) &
                    is.finite(df_f2$tbr), ]
  df_f2$cpe_g <- cut(df_f2$n_functional_cpe, breaks = c(-Inf, 0, 1, 2, Inf),
                      labels = c("0", "1", "2", "≥3"))
  dose_df2 <- df_f2 %>%
    group_by(cpe_g) %>%
    summarise(mean_tbr = mean(tbr, na.rm = TRUE),
              se_tbr   = sd(tbr, na.rm = TRUE) / sqrt(n()),
              n        = n(), .groups = "drop")
  p4f2 <- ggplot(dose_df2, aes(x = cpe_g, y = mean_tbr, group = 1)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey55", linewidth = 0.4) +
    geom_ribbon(aes(ymin = mean_tbr - se_tbr, ymax = mean_tbr + se_tbr),
                fill = col_omics[["Translatome"]], alpha = 0.18, colour = NA) +
    geom_line(colour = col_omics[["Translatome"]], linewidth = 1.0) +
    geom_point(size = 2.5, colour = col_omics[["Translatome"]],
               shape = 21, fill = "white", stroke = 1.2) +
    labs(x = "Functional CPEs", y = "Mean TBR (± SE)",
         title = "CPE–TBR\ndose response") +
    theme_oocyte(base_size = 8)

  # F4 shows the compensation phenomenon only (A-D); all CPE analysis is in F6,
  # and the CS/TBR threshold sensitivity is Supplementary Fig. S6.
  p4_full <- (p4a2 | p4b2) / (p4c2 | p4d2) +
    plot_layout(heights = c(1, 1)) +
    plot_annotation(
      tag_levels = list(c("A", "B", "C", "D")),
      theme = theme(plot.tag = element_text(size = 10, face = "bold",
                                             family = "Helvetica"))
    ) &
    theme(plot.title = element_text(size = 8, face = "bold"),
          plot.subtitle = element_text(size = 6.2, colour = "grey40"),
          plot.margin = margin(3, 5, 3, 4))

  save_full(p4_full, file.path(dir_f4, "F4_main_assembly"), w = 183, h = 140)
  cat("F4 main assembly (A–F) saved\n")
}
