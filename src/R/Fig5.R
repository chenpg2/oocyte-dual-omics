# Fig5_full_assembly.R — Batch 5
#
# 新增 panel：
#   F5D 补偿波基因规模气泡图（Wave × 区间 × 基因数）
#   F5E 各 Wave vTE 均值时序线图（分面）
# 全图拼接：F5A + F5B / F5D + F5E / F5C（热图独立）
#
# 输入：results/trd/compensation_waves.csv
#        results/trd/rate_matrix.csv
#        results/trd/interval_rate_stats.csv
# 输出：results/figures/F5_trd_waves/F5[D/E]_*.pdf
#        results/figures/F5_trd_waves/F5_main_assembly.pdf/.tiff

source("src/R/00_config.R")
library(ggplot2)
library(patchwork)
library(dplyr)
suppressPackageStartupMessages({
  library(ragg); library(ComplexHeatmap); library(circlize)
  library(grid); library(ggplotify)
})

dir_fig <- file.path(DIR_RESULTS, "figures", "F5_trd_waves")
dir.create(dir_fig, recursive = TRUE, showWarnings = FALSE)

save_pub <- function(p, name, w = 88, h = 80, dir = dir_fig) {
  grDevices::cairo_pdf(file.path(dir, paste0(name, ".pdf")),
                       width = w/25.4, height = h/25.4, family = "Helvetica")
  print(p); dev.off()
  svglite::svglite(file.path(dir, paste0(name, ".svg")),
                   width = w/25.4, height = h/25.4); print(p); dev.off()
  agg_tiff(file.path(dir, paste0(name, ".tiff")),
           width = w/25.4, height = h/25.4, units = "in", res = 600)
  print(p); dev.off()
}

# ============================================================
# 1. 数据加载
# ============================================================

waves    <- read.csv(file.path(DIR_RESULTS, "trd", "compensation_waves.csv"),
                     stringsAsFactors = FALSE)
rate_mat <- read.csv(file.path(DIR_RESULTS, "trd", "rate_matrix.csv"),
                     check.names = FALSE, stringsAsFactors = FALSE)
rate_stats <- read.csv(file.path(DIR_RESULTS, "trd", "interval_rate_stats.csv"),
                        stringsAsFactors = FALSE)

interval_order  <- c("GV→GVBD", "GVBD→MI6", "MI6→MI9", "MI9→MII")
interval_labels <- c("GV→GVBD", "GVBD→MI-6", "MI-6→MI-9", "MI-9→MII")
interval_sign   <- c("Positive", "Negative", "Positive", "Negative")
vte_cols        <- paste0("vTE_", interval_order)

wave_order  <- c("Wave 1 (GV→GVBD)", "Wave 2 (GVBD→MI6)",
                 "Wave 3 (MI6→MI9)",  "Wave 4 (MI9→MII)")
wave_labels <- c("W1: GV→GVBD", "W2: GVBD→MI-6",
                 "W3: MI-6→MI-9", "W4: MI-9→MII")
col_wave <- c(
  "W1: GV→GVBD"   = "#B42F37",
  "W2: GVBD→MI-6" = "#E07B39",
  "W3: MI-6→MI-9" = "#2D4B8E",
  "W4: MI-9→MII"  = "#00857C"
)
col_pos <- col_omics[["Translatome"]]
col_neg <- col_omics[["Transcriptome"]]
col_sign <- c("Positive" = col_pos, "Negative" = col_neg)

# ============================================================
# 2. F5D — 气泡时序图（每 Wave 在每区间的活跃基因数）
# ============================================================

# 每个基因属于某 Wave（max_vTE_interval 即其峰值区间）
waves_filt <- waves[waves$wave %in% wave_order, ]
waves_filt$wave_label <- wave_labels[match(waves_filt$wave, wave_order)]

# 交叉统计：每 wave_label × max_vTE_interval 的基因数
# 注意 max_vTE_interval 格式：如 "GV→GVBD"，需与 interval_order 对齐
# 先重新计算：对每个 wave，统计各区间的 mean vTE（from rate_matrix）
rownames(rate_mat) <- rate_mat$gene_id
common_genes <- intersect(waves_filt$gene_id, rownames(rate_mat))
waves_match  <- waves_filt[waves_filt$gene_id %in% common_genes, ]
vte_sub      <- as.matrix(rate_mat[waves_match$gene_id, vte_cols])

# 构建气泡数据：wave × interval 的平均 vTE + 基因数
bubble_df <- do.call(rbind, lapply(seq_along(wave_order), function(wi) {
  wv    <- wave_order[wi]
  wl    <- wave_labels[wi]
  idx   <- which(waves_match$wave == wv)
  n_g   <- length(idx)
  do.call(rbind, lapply(seq_along(interval_order), function(ii) {
    data.frame(
      wave      = wl,
      interval  = interval_labels[ii],
      mean_vTE  = mean(vte_sub[idx, ii], na.rm = TRUE),
      n_genes   = n_g,
      sign      = interval_sign[ii],
      stringsAsFactors = FALSE
    )
  }))
}))
bubble_df$wave     <- factor(bubble_df$wave,     levels = wave_labels)
bubble_df$interval <- factor(bubble_df$interval, levels = interval_labels)

# 气泡大小编码：|mean_vTE| × n_genes / 1000
bubble_df$bubble_size <- abs(bubble_df$mean_vTE) * sqrt(bubble_df$n_genes / 100)

p5d <- ggplot(bubble_df, aes(x = interval, y = wave,
                               size = bubble_size, color = sign)) +
  geom_point(alpha = 0.80) +
  geom_text(aes(label = sprintf("%.3f", mean_vTE)),
            size = 1.8, color = "gray20", vjust = -1.3, fontface = "plain") +
  scale_color_manual(values = col_sign, name = "Direction") +
  scale_size_continuous(range = c(2, 10), guide = "none") +
  labs(x = "Maturation interval",
       y = "Compensation wave",
       title = "Wave activity across intervals",
       subtitle = "Bubble size ∝ |mean vTE| × √n · colour = direction") +
  theme_oocyte(base_size = 8) +
  theme(legend.key.size = unit(3.5, "mm"),
        axis.text.x = element_text(angle = 25, hjust = 1, size = 6.5))

save_pub(p5d, "F5D_wave_bubble", w = 100, h = 72)
cat("F5D saved\n")

# ============================================================
# 3. F5E — 各 Wave vTE 均值时序线图（分面）
# ============================================================

wave_vte_df <- do.call(rbind, lapply(seq_along(wave_order), function(wi) {
  wv  <- wave_order[wi]
  wl  <- wave_labels[wi]
  idx <- which(waves_match$wave == wv)
  do.call(rbind, lapply(seq_along(interval_order), function(ii) {
    vals <- vte_sub[idx, ii]
    data.frame(
      wave      = wl,
      interval  = interval_labels[ii],
      int_num   = ii,
      mean_vTE  = mean(vals, na.rm = TRUE),
      se_vTE    = sd(vals, na.rm = TRUE) / sqrt(length(idx)),
      sign      = interval_sign[ii],
      stringsAsFactors = FALSE
    )
  }))
}))
wave_vte_df$wave     <- factor(wave_vte_df$wave,     levels = wave_labels)
wave_vte_df$interval <- factor(wave_vte_df$interval, levels = interval_labels)

p5e <- ggplot(wave_vte_df, aes(x = interval, y = mean_vTE,
                                group = wave, color = wave, fill = wave)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray55", linewidth = 0.4) +
  geom_ribbon(aes(ymin = mean_vTE - se_vTE, ymax = mean_vTE + se_vTE),
              color = NA, alpha = 0.18) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.2, shape = 21, fill = "white", stroke = 1.0) +
  scale_color_manual(values = col_wave, name = NULL) +
  scale_fill_manual(values  = col_wave, name = NULL, guide = "none") +
  facet_wrap(~ wave, nrow = 1, scales = "fixed",
             labeller = labeller(wave = function(x) sub(":.*", "", x))) +
  scale_x_discrete(labels = c("GV→GVBD" = "→GVBD", "GVBD→MI-6" = "→MI-6",
                              "MI-6→MI-9" = "→MI-9", "MI-9→MII" = "→MII")) +
  labs(x = "Interval (→ end stage)", y = expression("Mean vTE (log"[2]*"/3h)"),
       title = "Per-wave vTE kinetics",
       subtitle = "Each wave peaks at its assigned interval (mean ± SE)") +
  theme_oocyte(base_size = 7.5) +
  theme(legend.position = "none",
        strip.text = element_text(size = 6.5, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 5.6),
        panel.spacing = unit(1.6, "mm"))

save_pub(p5e, "F5E_wave_timeseries", w = 140, h = 72)
cat("F5E saved\n")

# ============================================================
# 4. 重建 F5A / F5B 用于拼接
# ============================================================

rate_stats$interval <- factor(rate_stats$interval, levels = interval_order)
rate_stats$sign     <- interval_sign[match(rate_stats$interval, interval_order)]
n_genes_rm          <- nrow(rate_mat)
rate_stats$se       <- rate_stats$sd_vTE / sqrt(n_genes_rm)
f5a_df              <- rate_stats[rate_stats$interval %in% interval_order, ]
f5a_df$x_num        <- as.integer(f5a_df$interval)
f5a_df$ymin_se      <- f5a_df$mean_vTE - f5a_df$se
f5a_df$ymax_se      <- f5a_df$mean_vTE + f5a_df$se

p5a <- ggplot(f5a_df, aes(x = interval, y = mean_vTE, group = 1)) +
  geom_hline(yintercept = 0, color = "gray55", linewidth = 0.5, linetype = "dashed") +
  annotate("rect", xmin = 0.5, xmax = 1.5, ymin = -Inf, ymax = Inf, fill = col_pos, alpha = 0.05) +
  annotate("rect", xmin = 1.5, xmax = 2.5, ymin = -Inf, ymax = Inf, fill = col_neg, alpha = 0.05) +
  annotate("rect", xmin = 2.5, xmax = 3.5, ymin = -Inf, ymax = Inf, fill = col_pos, alpha = 0.05) +
  annotate("rect", xmin = 3.5, xmax = 4.5, ymin = -Inf, ymax = Inf, fill = col_neg, alpha = 0.05) +
  geom_ribbon(aes(ymin = ymin_se, ymax = ymax_se), fill = "gray65", alpha = 0.30, color = NA) +
  geom_line(color = "gray30", linewidth = 0.7) +
  geom_point(aes(color = sign), size = 3.5, shape = 21, fill = "white", stroke = 1.5) +
  geom_text(aes(label = sprintf("%+.3f", mean_vTE)),
            vjust = -1.0, size = 2.5, color = "gray20") +
  scale_color_manual(values = col_sign, guide = "none") +
  scale_x_discrete(labels = interval_labels) +
  scale_y_continuous(expand = expansion(mult = c(0.12, 0.20))) +
  labs(x = "Interval", y = expression("Mean vTE (log"[2]*")"),
       title = "Biphasic vTE oscillation",
       subtitle = sprintf("n=%d genes | +/−/+/− pattern", n_genes_rm)) +
  theme_oocyte(base_size = 8) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1, size = 6.3))

# F5B (violin)
set.seed(params$seed)
samp_idx <- sample(nrow(rate_mat), min(8000L, nrow(rate_mat)))
f5b_long <- do.call(rbind, lapply(seq_along(interval_order), function(i) {
  intv  <- interval_order[i]
  col_v <- paste0("vTE_", intv)
  if (!col_v %in% colnames(rate_mat)) return(NULL)
  data.frame(interval = factor(intv, levels = interval_order),
             x_label  = interval_labels[i],
             vTE      = rate_mat[samp_idx, col_v],
             sign     = interval_sign[i],
             stringsAsFactors = FALSE)
}))
f5b_long <- f5b_long[is.finite(f5b_long$vTE), ]
f5b_long$x_label <- factor(f5b_long$x_label, levels = interval_labels)
med_df <- f5b_long %>% group_by(x_label, sign) %>%
  summarise(med = median(vTE, na.rm = TRUE), .groups = "drop")

p5b <- ggplot(f5b_long, aes(x = x_label, y = vTE)) +
  geom_hline(yintercept = 0, color = "gray55", linewidth = 0.4, linetype = "dashed") +
  geom_violin(aes(fill = sign, color = sign), trim = TRUE, alpha = 0.35,
              linewidth = 0.4, scale = "width") +
  geom_boxplot(width = 0.10, outlier.shape = NA, linewidth = 0.5,
               color = "gray25", fill = "white", alpha = 0.8) +
  geom_point(data = med_df, aes(x = x_label, y = med, color = sign),
             size = 2.5, shape = 18, inherit.aes = FALSE) +
  scale_fill_manual(values = col_sign, guide = "none") +
  scale_color_manual(values = col_sign, guide = "none") +
  coord_cartesian(ylim = c(-2.5, 2.5)) +
  labs(x = "Interval", y = expression("vTE (log"[2]*"/3h)"),
       title = "Population vTE distributions",
       subtitle = "Median flips sign each interval") +
  theme_oocyte(base_size = 8) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1, size = 6.3))

# ============================================================
# 4c. F5C — per-gene wave heatmap (ComplexHeatmap -> ggplot via ggplotify)
# ============================================================

hm_wave_order <- c("Wave 1 (GV→GVBD)", "Wave 2 (GVBD→MI6)",
                   "Wave 3 (MI6→MI9)",  "Wave 4 (MI9→MII)")
hm_col_lab    <- c("GV→\nGVBD", "GVBD→\nMI-6", "MI-6→\nMI-9", "MI-9→\nMII")
hm_vte_cols   <- paste0("vTE_", interval_order)

rownames(rate_mat) <- rate_mat$gene_id
hm_waves <- waves[waves$wave %in% hm_wave_order, ]
set.seed(params$seed)
hm_samp <- do.call(rbind, lapply(hm_wave_order, function(w) {
  s <- hm_waves[hm_waves$wave == w, ]
  s[sample(nrow(s), min(1200L, nrow(s))), ]
}))
hm_samp <- hm_samp[hm_samp$gene_id %in% rownames(rate_mat), ]
hm_mat  <- as.matrix(rate_mat[hm_samp$gene_id, hm_vte_cols])
colnames(hm_mat) <- hm_col_lab
hm_wf   <- factor(hm_samp$wave, levels = hm_wave_order, labels = c("W1", "W2", "W3", "W4"))

hm_max    <- as.numeric(quantile(abs(hm_mat), 0.98, na.rm = TRUE))
hm_colfun <- colorRamp2(c(-hm_max, 0, hm_max),
                        c(col_omics[["Transcriptome"]], "white", col_omics[["Translatome"]]))
hm_wavecol <- c("W1" = "#B42F37", "W2" = "#E07B39", "W3" = "#2D4B8E", "W4" = "#00857C")

ht_wave <- Heatmap(
  hm_mat, col = hm_colfun, name = "vTE",
  row_split = hm_wf, cluster_rows = FALSE, cluster_columns = FALSE,
  show_row_names = FALSE, show_column_names = TRUE,
  column_names_gp = gpar(fontsize = 6), column_names_rot = 0, column_names_centered = TRUE,
  left_annotation = rowAnnotation(Wave = hm_wf, col = list(Wave = hm_wavecol),
                                  show_legend = FALSE, annotation_name_gp = gpar(fontsize = 6),
                                  simple_anno_size = unit(2.5, "mm")),
  row_title_gp = gpar(fontsize = 6.5, fontface = "bold"),
  row_gap = unit(1, "mm"), border = FALSE, use_raster = TRUE, raster_quality = 4,
  heatmap_legend_param = list(title = "vTE", title_gp = gpar(fontsize = 6.5, fontface = "bold"),
                              labels_gp = gpar(fontsize = 6), grid_height = unit(2.5, "mm"),
                              grid_width = unit(2.5, "mm"), direction = "horizontal",
                              legend_width = unit(20, "mm")),
  column_title = "vTE by wave",
  column_title_gp = gpar(fontsize = 7.5, fontface = "bold"))

p5c <- ggplotify::as.ggplot(
  grid::grid.grabExpr(draw(ht_wave, heatmap_legend_side = "bottom",
                           padding = unit(c(1, 1, 1, 1), "mm"))))

# ============================================================
# 5. 全图拼接 F5_main_assembly (A B C top, D E bottom)
# ============================================================

design_f5 <- "AABBCC\nAABBCC\nDDDEEE"

p5_full <- (p5a + p5b + p5c + p5d + p5e) +
  plot_layout(design = design_f5, heights = c(1, 1, 1.25)) +
  plot_annotation(
    tag_levels = list(c("A", "B", "C", "D", "E")),
    theme = theme(plot.tag = element_text(size = 10, face = "bold"))
  ) &
  theme(plot.subtitle = element_text(size = 6.3, colour = "grey40"),
        plot.margin = margin(2, 4, 2, 3),
        plot.tag = element_text(size = 10, face = "bold"))

save_pub(p5_full, "F5_main_assembly", w = 183, h = 190)
cat("F5_main_assembly saved\n")

cat("\n=== Batch 5 / F5 完成 ===\n")
cat("输出目录:", dir_fig, "\n")
