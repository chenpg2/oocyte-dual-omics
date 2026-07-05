# Fig3_full_assembly.R — Batch 5
#
# 新增 panel：
#   F3D 环形比例图（五轨迹基因占比）
#   F3E TBR × 轨迹小提琴图
#   F3F CS  × 轨迹箱线图
# 最终拼接：
#   Row1: F3A（相位图）| F3D（环形）
#   Row2: F3E（TBR violin）| F3F（CS boxplot）
#   Row3: F3B（时序线图，宽）
#   Note: F3C（双层热图）作为独立 PDF，不嵌入 patchwork
#
# 输入：results/normalized/normalized_data.RData
#        results/tta/cluster_assignments.csv
#        results/supplementary/S1_tbr_cs_matrix.csv
# 输出：results/figures/F3_trajectories/F3[D/E/F]_*.pdf
#        results/figures/F3_trajectories/F3_main_assembly.pdf/.tiff

source("src/R/00_config.R")
library(ggplot2)
library(patchwork)
library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(grid)
suppressPackageStartupMessages(library(ragg))

dir_fig <- file.path(DIR_RESULTS, "figures", "F3_trajectories")
dir.create(dir_fig, recursive = TRUE, showWarnings = FALSE)

# ---- 通用导出函数（PDF + SVG + TIFF）----
save_pub <- function(p, name, w = 88, h = 80, dir = dir_fig) {
  path_pdf  <- file.path(dir, paste0(name, ".pdf"))
  path_svg  <- file.path(dir, paste0(name, ".svg"))
  path_tiff <- file.path(dir, paste0(name, ".tiff"))
  grDevices::cairo_pdf(path_pdf, width = w/25.4, height = h/25.4, family = "Helvetica")
  print(p); dev.off()
  svglite::svglite(path_svg, width = w/25.4, height = h/25.4); print(p); dev.off()
  agg_tiff(path_tiff, width = w/25.4, height = h/25.4, units = "in", res = 600)
  print(p); dev.off()
}

# ============================================================
# 1. 数据加载
# ============================================================

load(file.path(DIR_RESULTS, "normalized", "normalized_data.RData"))
tta <- read.csv(file.path(DIR_RESULTS, "tta", "cluster_assignments.csv"),
                stringsAsFactors = FALSE)
s1  <- read.csv(file.path(DIR_RESULTS, "supplementary", "S1_tbr_cs_matrix.csv"),
                stringsAsFactors = FALSE)

tp_levels    <- c("0h", "3h", "6h", "9h", "12h")
stage_labels <- c("GV", "GVBD", "MI-6", "MI-9", "MII")

col_traj <- c(
  "TE-Only Activation"          = "#00857C",
  "Late Compensatory Buffering" = "#2D4B8E",
  "Mild Coordinated Clearance"  = "#A7A9AC",
  "Coordinated Clearance"       = "#6E6E6E",
  "Deep Coordinated Clearance"  = "#B42F37"
)
traj_order <- names(col_traj)
traj_short <- c(
  "TE-Only Activation"          = "C1\nTE-Only",
  "Late Compensatory Buffering" = "C2\nCompen.",
  "Mild Coordinated Clearance"  = "C5\nMild Clear.",
  "Coordinated Clearance"       = "C3\nClearance",
  "Deep Coordinated Clearance"  = "C4\nDeep Clear."
)

# S1 已含 trajectory_type 列，直接筛选
df_s1 <- s1[!is.na(s1$trajectory_type) & s1$trajectory_type %in% traj_order, ]
df_s1$trajectory_type <- factor(df_s1$trajectory_type, levels = traj_order)

# ============================================================
# 2. F3D — 环形比例图（coord_polar + geom_col）
# ============================================================

traj_n <- tta %>%
  filter(trajectory_type %in% traj_order) %>%
  count(trajectory_type) %>%
  mutate(
    trajectory_type = factor(trajectory_type, levels = traj_order),
    pct    = n / sum(n) * 100,
    label  = sprintf("%s\n%.1f%%\n(n=%d)", traj_short[as.character(trajectory_type)], pct, n),
    ymax   = cumsum(pct),
    ymin   = lag(ymax, default = 0),
    ymid   = (ymax + ymin) / 2
  )

p3d <- ggplot(traj_n, aes(ymax = ymax, ymin = ymin,
                            xmax = 4, xmin = 2.5, fill = trajectory_type)) +
  geom_rect() +
  geom_text(aes(x = 5.1, y = ymid, label = label, color = trajectory_type),
            size = 1.9, lineheight = 0.95, fontface = "bold", inherit.aes = FALSE) +
  scale_fill_manual(values  = col_traj, guide = "none") +
  scale_color_manual(values = col_traj, guide = "none") +
  coord_polar(theta = "y", start = 0) +
  xlim(c(0.5, 6.5)) +
  labs(title = "Trajectory composition",
       subtitle = sprintf("n = %d genes total", sum(traj_n$n))) +
  theme_void(base_size = 8) +
  theme(
    plot.title    = element_text(size = 8, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 7,  hjust = 0.5, color = "gray40"),
    plot.margin   = margin(2, 2, 2, 2, "mm")
  )

save_pub(p3d, "F3D_trajectory_donut", w = 80, h = 80)
cat("F3D saved\n")

# ============================================================
# 3. F3E — TBR × 轨迹小提琴图
# ============================================================

med_e <- df_s1 %>%
  group_by(trajectory_type) %>%
  summarise(med = median(tbr, na.rm = TRUE), .groups = "drop")

p3e <- ggplot(df_s1, aes(x = trajectory_type, y = tbr,
                           fill = trajectory_type, color = trajectory_type)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray55", linewidth = 0.4) +
  geom_hline(yintercept = 0.3, linetype = "dotted", color = "gray45", linewidth = 0.35) +
  geom_violin(trim = TRUE, alpha = 0.35, linewidth = 0.35, scale = "width") +
  geom_boxplot(width = 0.10, outlier.shape = NA, linewidth = 0.45,
               color = "gray25", fill = "white", alpha = 0.85) +
  geom_point(data = med_e, aes(x = trajectory_type, y = med, color = trajectory_type),
             size = 2.5, shape = 18, inherit.aes = FALSE) +
  scale_fill_manual(values  = col_traj, guide = "none") +
  scale_color_manual(values = col_traj, guide = "none") +
  scale_x_discrete(labels = traj_short) +
  coord_cartesian(ylim = c(-3, 6)) +
  labs(x = NULL,
       y = "TBR (log₂)",
       title = "TBR distribution per trajectory",
       subtitle = "C2 highest TBR; C3/C4 negative (coordinated decline)") +
  theme_oocyte(base_size = 7.5) +
  theme(axis.text.x = element_text(size = 6, lineheight = 1.05))

save_pub(p3e, "F3E_tbr_by_trajectory", w = 100, h = 72)
cat("F3E saved\n")

# ============================================================
# 4. F3F — CS × 轨迹箱线图
# ============================================================

med_f <- df_s1 %>%
  group_by(trajectory_type) %>%
  summarise(med = median(cs, na.rm = TRUE), .groups = "drop")

p3f <- ggplot(df_s1, aes(x = trajectory_type, y = cs,
                           fill = trajectory_type)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray55", linewidth = 0.4) +
  geom_hline(yintercept = 1, linetype = "dotted", color = "gray45", linewidth = 0.35) +
  geom_boxplot(outlier.shape = NA, linewidth = 0.45,
               alpha = 0.70, width = 0.55) +
  geom_point(data = med_f, aes(x = trajectory_type, y = med),
             size = 2.5, shape = 18, color = "gray20", inherit.aes = FALSE) +
  scale_fill_manual(values = col_traj, guide = "none") +
  scale_x_discrete(labels = traj_short) +
  coord_cartesian(ylim = c(-2, 5)) +
  labs(x = NULL,
       y = "Clearance Severity (CS = −ΔmRNA)",
       title = "mRNA clearance per trajectory",
       subtitle = "C4 (Deep Clearance) highest CS; C1 (TE-Only) minimal") +
  theme_oocyte(base_size = 7.5) +
  theme(axis.text.x = element_text(size = 6, lineheight = 1.05))

save_pub(p3f, "F3F_cs_by_trajectory", w = 100, h = 72)
cat("F3F saved\n")

# ============================================================
# 5. 重建 F3A（相位图）和 F3B（时序线图）用于拼接
# ============================================================

## --- F3A 重建 ---
mrna_by_gene <- matrix(NA_real_, nrow = length(gene_ids_filt), ncol = 5,
                        dimnames = list(gene_ids_filt, tp_levels))
te_by_gene   <- matrix(NA_real_, nrow = length(gene_ids_filt), ncol = 5,
                        dimnames = list(gene_ids_filt, tp_levels))
for (i in seq_along(tp_levels)) {
  tp      <- tp_levels[i]
  pm      <- pair_map[pair_map$time_point == tp, ]
  m_cols  <- pm$transcriptome_id
  te_cols <- grep(paste0("^", tp, "_rep"), colnames(te_A), value = TRUE)
  mrna_by_gene[, i] <- rowMeans(log1p(norm_trans_A[gene_ids_filt, m_cols, drop = FALSE]), na.rm = TRUE)
  te_by_gene[, i]   <- rowMeans(te_A[gene_ids_filt, te_cols, drop = FALSE], na.rm = TRUE)
}
mrna_delta <- sweep(mrna_by_gene, 1, mrna_by_gene[, "0h"], "-")
te_delta   <- sweep(te_by_gene,   1, te_by_gene[,   "0h"], "-")
tta_map    <- setNames(tta$trajectory_type, tta$gene_id)

traj_short_pp <- c(
  "TE-Only Activation"          = "C1",
  "Late Compensatory Buffering" = "C2",
  "Mild Coordinated Clearance"  = "C5",
  "Coordinated Clearance"       = "C3",
  "Deep Coordinated Clearance"  = "C4"
)
pp_df <- do.call(rbind, lapply(traj_order, function(traj) {
  genes <- intersect(names(tta_map)[tta_map == traj], rownames(mrna_delta))
  data.frame(trajectory = traj, stage = stage_labels, stage_num = seq_along(stage_labels),
             mrna_d = colMeans(mrna_delta[genes, , drop = FALSE], na.rm = TRUE),
             te_d   = colMeans(te_delta[genes,   , drop = FALSE], na.rm = TRUE),
             stringsAsFactors = FALSE)
}))
pp_df$trajectory <- factor(pp_df$trajectory, levels = traj_order)

arr_df <- do.call(rbind, lapply(traj_order, function(traj) {
  sub <- pp_df[pp_df$trajectory == traj, ][order(pp_df$stage_num[pp_df$trajectory == traj]), ]
  n   <- nrow(sub)
  data.frame(trajectory = traj,
             x = sub$mrna_d[-n], y = sub$te_d[-n],
             xend = sub$mrna_d[-1], yend = sub$te_d[-1],
             stringsAsFactors = FALSE)
}))
arr_df$trajectory <- factor(arr_df$trajectory, levels = traj_order)
label_end <- pp_df[pp_df$stage == "MII", ]
label_gv  <- pp_df[pp_df$stage == "GV",  ]

p3a <- ggplot() +
  geom_hline(yintercept = 0, linetype = "dotted", color = "gray65", linewidth = 0.35) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "gray65", linewidth = 0.35) +
  geom_segment(data = arr_df,
               aes(x = x, y = y, xend = xend, yend = yend, color = trajectory),
               linewidth = 0.85,
               arrow = arrow(length = unit(2.0, "mm"), type = "closed", ends = "last"),
               lineend = "round") +
  geom_point(data = pp_df,
             aes(x = mrna_d, y = te_d, color = trajectory),
             size = 1.8, shape = 21, fill = "white", stroke = 0.9) +
  geom_point(data = label_gv,
             aes(x = mrna_d, y = te_d, color = trajectory),
             size = 2.5, shape = 19) +
  ggrepel::geom_text_repel(
    data = label_end,
    aes(x = mrna_d, y = te_d, label = traj_short_pp[as.character(trajectory)],
        color = trajectory),
    size = 2.2, fontface = "bold", linewidth = 0.3,
    box.padding = 0.3, point.padding = 0.2, max.overlaps = 20, show.legend = FALSE) +
  scale_color_manual(values = col_traj, name = NULL,
                     labels = setNames(names(traj_short_pp), traj_order)) +
  labs(x = expression(Delta~"mRNA from GV"),
       y = expression(Delta~"TE from GV"),
       title = "Phase portrait: mRNA–TE trajectory") +
  theme_oocyte(base_size = 7.5) +
  theme(legend.position = "none")

## --- F3B 重建（简化版，用于拼接） ---
omics_colors <- c("mRNA delta" = col_omics[["Transcriptome"]],
                   "TE delta"   = col_omics[["Translatome"]])
omics_lty    <- c("mRNA delta" = "solid", "TE delta" = "dashed")
omics_shape  <- c("mRNA delta" = 21L,     "TE delta" = 24L)

rows <- list()
for (traj in traj_order) {
  genes <- intersect(names(tta_map)[tta_map == traj], rownames(mrna_delta))
  n_g   <- length(genes)
  for (i in seq_along(stage_labels)) {
    m_v <- mrna_delta[genes, i]; t_v <- te_delta[genes, i]
    rows <- c(rows, list(
      data.frame(trajectory = traj, stage = stage_labels[i], stage_num = i,
                 omics = "mRNA delta", mean_val = mean(m_v, na.rm=TRUE),
                 se_val = sd(m_v, na.rm=TRUE)/sqrt(n_g), n_genes = n_g),
      data.frame(trajectory = traj, stage = stage_labels[i], stage_num = i,
                 omics = "TE delta", mean_val = mean(t_v, na.rm=TRUE),
                 se_val = sd(t_v, na.rm=TRUE)/sqrt(n_g), n_genes = n_g)
    ))
  }
}
f3b_df <- do.call(rbind, rows)
f3b_df$stage      <- factor(f3b_df$stage, levels = stage_labels)
f3b_df$trajectory <- factor(f3b_df$trajectory, levels = traj_order)
traj_label_map    <- setNames(
  c("C1: TE-Only", "C2: Compensatory", "C5: Mild Clearance",
    "C3: Clearance", "C4: Deep Clearance"), traj_order)
f3b_df$traj_label <- factor(traj_label_map[as.character(f3b_df$trajectory)],
                              levels = traj_label_map[traj_order])
f3b_df$omics      <- factor(f3b_df$omics, levels = c("mRNA delta", "TE delta"))

p3b <- ggplot(f3b_df, aes(x = stage, y = mean_val, group = omics)) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "gray70", linewidth = 0.35) +
  geom_ribbon(aes(ymin = mean_val - se_val, ymax = mean_val + se_val, fill = omics),
              color = NA, alpha = 0.15) +
  geom_line(aes(color = omics, linetype = omics), linewidth = 0.85) +
  geom_point(aes(color = omics, shape = omics), size = 2.0, fill = "white", stroke = 0.9) +
  scale_color_manual(values = omics_colors, name = NULL) +
  scale_fill_manual(values  = omics_colors, name = NULL, guide = "none") +
  scale_linetype_manual(values = omics_lty, name = NULL) +
  scale_shape_manual(values = omics_shape, name = NULL) +
  facet_wrap(~ traj_label, nrow = 1, scales = "fixed") +
  labs(x = "Stage", y = expression(Delta~"from GV (mean ± SE)"),
       title = "Five trajectory mRNA–TE coordination programs") +
  theme_oocyte(base_size = 7) +
  theme(legend.position  = "bottom", legend.key.width = unit(5, "mm"),
        strip.text = element_text(size = 6, face = "bold"),
        axis.text.x = element_text(angle = 35, hjust = 1, size = 5.5),
        panel.spacing = unit(1.5, "mm"))

# ============================================================
# 6. 全图拼接 F3_main_assembly（A+D / B / E+F，不含C热图）
# ============================================================

p3_top    <- (p3a | p3d) + plot_layout(widths = c(1.2, 1))
p3_mid    <- p3b
p3_bottom <- (p3e | p3f) + plot_layout(widths = c(1, 1))

p3_full <- (p3_top / p3_mid / p3_bottom) +
  plot_annotation(
    tag_levels = list(c("A", "B", "C", "D", "E"))
  ) +
  plot_layout(heights = c(1.1, 0.9, 0.9)) &
  theme(plot.subtitle = element_text(size = 6.3, colour = "grey40"),
        plot.margin = margin(2, 4, 2, 3),
        plot.tag = element_text(size = 10, face = "bold"))

save_pub(p3_full, "F3_main_assembly", w = 183, h = 220)
cat("F3_main_assembly saved\n")

cat("\n=== Batch 5 / F3 完成 ===\n")
cat("输出目录:", dir_fig, "\n")
