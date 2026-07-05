# S1_continuous_compensation.R — IU-S1: 连续补偿指标 TBR/CS
# 目的：用连续指标替代 ARI=0.204 的不稳定 TTA 聚类
#
# 定义：
#   Clearance Severity (CS) = -delta_mrna  (>0 表示 mRNA 降解)
#   Translational Buffering Ratio (TBR) = delta_TE / max(|delta_mrna|, 0.1)
#     TBR > 0  : 翻译补偿 (TE 上升抵消 mRNA 下降)
#     TBR > 1  : 超补偿
#     TBR ≈ 0  : 协调清除
#     TBR < -0.5 : 翻译加速下降
#
# 输入：results/normalized/normalized_data.RData, results/tta/cluster_assignments.csv
# 输出：results/supplementary/S1_tbr_cs_matrix.csv
#        results/supplementary/S1_compensation_landscape.pdf
#        results/supplementary/S1_tbr_distribution.pdf
#        results/supplementary/S1_quadrant_summary.csv

source("src/R/00_config.R")
library(ggplot2)
library(dplyr)

set.seed(params$seed)
dir_sup <- file.path(DIR_RESULTS, "supplementary")
dir.create(dir_sup, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# 1. Load data & compute group-mean TE/mRNA per timepoint
# ============================================================

load(file.path(DIR_RESULTS, "normalized", "normalized_data.RData"))
tta <- read.csv(file.path(DIR_RESULTS, "tta", "cluster_assignments.csv"),
                stringsAsFactors = FALSE)

cat("Computing group-mean TE/mRNA...\n")
tp_levels <- c("0h", "3h", "6h", "9h", "12h")

mean_te   <- matrix(NA_real_, nrow = length(gene_ids_filt), ncol = 5,
                    dimnames = list(gene_ids_filt, tp_levels))
mean_mrna <- matrix(NA_real_, nrow = length(gene_ids_filt), ncol = 5,
                    dimnames = list(gene_ids_filt, tp_levels))

for (i in seq_along(tp_levels)) {
  tp      <- tp_levels[i]
  pm      <- pair_map[pair_map$time_point == tp, ]
  m_cols  <- pm$transcriptome_id
  te_cols <- grep(paste0("^", tp, "_rep"), colnames(te_A), value = TRUE)
  mean_mrna[, i] <- rowMeans(norm_trans_A[gene_ids_filt, m_cols,  drop = FALSE], na.rm = TRUE)
  mean_te[, i]   <- rowMeans(te_A[gene_ids_filt,         te_cols, drop = FALSE], na.rm = TRUE)
}

# log1p transform for mRNA (linear scale in norm_trans_A)
log1p_mrna <- log1p(mean_mrna)

# ============================================================
# 2. Compute CS and TBR
# ============================================================

cat("Computing CS and TBR...\n")

delta_mrna <- log1p_mrna[, "12h"] - log1p_mrna[, "0h"]   # GV→MII mRNA change (log1p scale)
delta_te   <- mean_te[, "12h"]    - mean_te[, "0h"]       # GV→MII TE change (log2 ratio scale)

cs  <- -delta_mrna                                          # Clearance Severity: >0 = degraded
tbr <- delta_te / pmax(abs(delta_mrna), 0.1)               # Translational Buffering Ratio

df <- data.frame(
  gene_id         = gene_ids_filt,
  delta_mrna      = delta_mrna,
  delta_te        = delta_te,
  cs              = cs,
  tbr             = tbr,
  stringsAsFactors = FALSE
)

# Merge TTA labels
df <- merge(df, tta[, c("gene_id", "trajectory_type")], by = "gene_id", all.x = TRUE)
df$trajectory_type[is.na(df$trajectory_type)] <- "Unclassified"

cat("  Total genes:", nrow(df), "\n")
cat("  CS range: [", round(min(df$cs, na.rm=TRUE), 2), ",",
    round(max(df$cs, na.rm=TRUE), 2), "]\n")
cat("  TBR range: [", round(quantile(df$tbr, 0.01, na.rm=TRUE), 2), ",",
    round(quantile(df$tbr, 0.99, na.rm=TRUE), 2), "] (1-99 pctile)\n")

# ============================================================
# 3. Quadrant classification (continuous, no clustering needed)
# ============================================================

# Threshold definitions
CS_THRESH  <- 1.0   # mRNA must drop by log1p(1) ≈ 0.69 (modest clearance) — OR use 0.5 (log1p scale)
TBR_HIGH   <- 0.3   # TBR > 0.3: meaningful buffering
TBR_LOW    <- 0.05  # TBR < 0.05: coordinated clearance

df$quadrant <- with(df, case_when(
  cs > CS_THRESH & tbr > TBR_HIGH ~ "Compensatory",
  cs > CS_THRESH & tbr < TBR_LOW  ~ "Coordinated Clearance",
  cs > CS_THRESH & tbr >= TBR_LOW & tbr <= TBR_HIGH ~ "Partial Buffering",
  cs <= CS_THRESH & delta_te > 0.5 ~ "TE-Only Activation",
  TRUE ~ "Stable / Minor Change"
))

quad_summary <- as.data.frame(table(df$quadrant))
colnames(quad_summary) <- c("Quadrant", "n_genes")
quad_summary$pct <- round(quad_summary$n_genes / nrow(df) * 100, 1)
quad_summary <- quad_summary[order(-quad_summary$n_genes), ]

cat("\nQuadrant summary:\n")
print(quad_summary)
write.csv(quad_summary, file.path(dir_sup, "S1_quadrant_summary.csv"), row.names = FALSE)

# Save full matrix
write.csv(df[, c("gene_id", "delta_mrna", "delta_te", "cs", "tbr",
                 "trajectory_type", "quadrant")],
          file.path(dir_sup, "S1_tbr_cs_matrix.csv"), row.names = FALSE)
cat("Matrix saved:", nrow(df), "genes\n")

# ============================================================
# 4. TBR distribution plot
# ============================================================

cat("Generating TBR distribution plot...\n")

tbr_clip <- pmax(pmin(df$tbr, 3), -2)  # clip for display
med_tbr  <- median(df$tbr, na.rm = TRUE)
pct_pos  <- round(mean(df$tbr > 0, na.rm = TRUE) * 100, 1)

p_tbr <- ggplot(df, aes(x = tbr_clip)) +
  geom_histogram(bins = 80, fill = col_omics[["Translatome"]], alpha = 0.8,
                 color = NA) +
  geom_vline(xintercept = 0,        linetype = "solid",  color = "gray30", linewidth = 0.6) +
  geom_vline(xintercept = med_tbr,  linetype = "dashed", color = "gray50", linewidth = 0.5) +
  geom_vline(xintercept = TBR_HIGH, linetype = "dotted", color = col_omics[["Transcriptome"]],
             linewidth = 0.5) +
  annotate("text", x = med_tbr + 0.05, y = Inf, label = paste0("median=", round(med_tbr, 2)),
           hjust = 0, vjust = 1.5, size = 2.5, color = "gray40") +
  annotate("text", x = 2.5, y = Inf,
           label = paste0(pct_pos, "% genes\nTBR > 0"),
           hjust = 1, vjust = 1.5, size = 2.8, color = col_omics[["Translatome"]]) +
  scale_x_continuous(limits = c(-2, 3),
                     breaks = c(-2, -1, 0, TBR_HIGH, 1, 2, 3),
                     labels = c("-2", "-1", "0", "0.3", "1", "2", "3")) +
  labs(x = "Translational Buffering Ratio (TBR)",
       y = "Number of genes",
       title = "Distribution of Translational Buffering Ratio (GV→MII)",
       subtitle = paste0("TBR = ΔTE / |ΔmRNA|  |  n = ", nrow(df), " genes")) +
  theme_oocyte()

ggsave(file.path(dir_sup, "S1_tbr_distribution.pdf"), p_tbr,
       width = 140, height = 100, units = "mm")

# ============================================================
# 5. CS × TBR 2D landscape
# ============================================================

cat("Generating CS-TBR landscape plot...\n")

# Color by TTA trajectory type
tta_colors <- c(
  "TE-Only Activation"          = col_timepoint[["GV"]],
  "Late Compensatory Buffering" = col_timepoint[["GVBD"]],
  "Coordinated Clearance"       = col_timepoint[["MI-6"]],
  "Deep Coordinated Clearance"  = col_timepoint[["MI-9"]],
  "Mild Coordinated Clearance"  = col_timepoint[["MII"]],
  "Unclassified"                = "gray70"
)

# Clip for display
df$cs_clip  <- pmax(pmin(df$cs,  3), -1.5)
df$tbr_clip <- pmax(pmin(df$tbr, 2.5), -2)

# Key gene labels (subset)
key_genes <- c("Ccnb1", "Cnot7", "Cnot8", "Zp3", "Cdk1", "Gdf9")
df_label  <- df[df$gene_id %in% gene_ids_filt, ]  # all genes in df already filtered
# Try to get gene symbols from TTA table
if ("gene_symbol" %in% colnames(tta)) {
  df <- merge(df, tta[, c("gene_id", "gene_symbol")], by = "gene_id", all.x = TRUE)
  df_key <- df[!is.na(df$gene_symbol) & df$gene_symbol %in% key_genes, ]
}

p_land <- ggplot(df, aes(x = cs_clip, y = tbr_clip)) +
  geom_point(aes(color = trajectory_type), size = 0.3, alpha = 0.25) +
  # Reference lines
  geom_hline(yintercept = 0,        linetype = "solid",  color = "gray40", linewidth = 0.4) +
  geom_hline(yintercept = TBR_HIGH, linetype = "dashed", color = "gray60", linewidth = 0.3) +
  geom_vline(xintercept = CS_THRESH, linetype = "dashed", color = "gray60", linewidth = 0.3) +
  # Quadrant annotations
  annotate("text", x = 2.5,   y = 2.3,  label = "Compensatory",
           size = 2.8, color = "gray30", hjust = 1) +
  annotate("text", x = 2.5,   y = -1.8, label = "Coordinated\nClearance",
           size = 2.8, color = "gray30", hjust = 1) +
  annotate("text", x = -1.3,  y = 2.3,  label = "TE-Only\nActivation",
           size = 2.8, color = "gray30", hjust = 0) +
  scale_color_manual(values = tta_colors, name = "TTA Trajectory") +
  scale_x_continuous(breaks = c(-1, 0, 1, 2, 3)) +
  scale_y_continuous(breaks = c(-2, -1, 0, 0.3, 1, 2)) +
  labs(x = "Clearance Severity (CS = -delta_mRNA)",
       y = "Translational Buffering Ratio (TBR)",
       title = "Translational Compensation Landscape (GV→MII)",
       subtitle = "Each point = one gene; TTA trajectory type from IU-5 clustering") +
  theme_oocyte() +
  theme(legend.position = "right",
        legend.key.size = unit(3, "mm"),
        legend.text = element_text(size = 7))

# Add key gene labels if found
if (exists("df_key") && nrow(df_key) > 0) {
  p_land <- p_land +
    ggrepel::geom_text_repel(
      data = df_key,
      aes(label = gene_symbol),
      size = 2.5, color = "black",
      max.overlaps = 20,
      segment.size = 0.3
    )
}

ggsave(file.path(dir_sup, "S1_compensation_landscape.pdf"), p_land,
       width = 180, height = 140, units = "mm")

# ============================================================
# 6. TBR by TTA trajectory type (violin)
# ============================================================

cat("Generating TBR by trajectory type plot...\n")

df_tta <- df[df$trajectory_type != "Unclassified", ]
df_tta$trajectory_type <- factor(df_tta$trajectory_type, levels = c(
  "TE-Only Activation", "Late Compensatory Buffering",
  "Mild Coordinated Clearance", "Coordinated Clearance", "Deep Coordinated Clearance"
))

# Wilcoxon test: each class vs. global median
global_med <- median(df$tbr, na.rm = TRUE)
wtest <- sapply(levels(df_tta$trajectory_type), function(cl) {
  g <- df_tta$tbr[df_tta$trajectory_type == cl]
  wilcox.test(g, mu = global_med)$p.value
})
cat("  Wilcoxon vs global median (BH-corrected):\n")
wtest_adj <- p.adjust(wtest, method = "BH")
print(round(wtest_adj, 4))

p_violin <- ggplot(df_tta, aes(x = trajectory_type, y = pmin(tbr, 2.5))) +
  geom_violin(aes(fill = trajectory_type), alpha = 0.7, color = NA) +
  geom_boxplot(width = 0.12, outlier.shape = NA, color = "gray30", linewidth = 0.4) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.4) +
  scale_fill_manual(values = tta_colors, guide = "none") +
  scale_x_discrete(labels = function(x) gsub(" ", "\n", x)) +
  labs(x = NULL, y = "TBR (clipped at 2.5)",
       title = "TBR Distribution by TTA Trajectory Type") +
  theme_oocyte() +
  theme(axis.text.x = element_text(size = 7))

ggsave(file.path(dir_sup, "S1_tbr_by_trajectory.pdf"), p_violin,
       width = 180, height = 110, units = "mm")

# ============================================================
# 7. Summary
# ============================================================

cat("\n=== IU-S1 Summary ===\n")
cat("Total genes analyzed:", nrow(df), "\n")
cat("\nQuadrant counts (CS threshold=", CS_THRESH, ", TBR threshold=", TBR_HIGH, "):\n")
print(quad_summary)
cat("\nGlobal TBR median:", round(median(df$tbr, na.rm=TRUE), 3), "\n")
cat("Genes with TBR > 0:", sum(df$tbr > 0, na.rm=TRUE),
    sprintf("(%.1f%%)", mean(df$tbr > 0, na.rm=TRUE)*100), "\n")
cat("Genes with CS > 1 AND TBR > 0.3 (Compensatory):",
    sum(df$cs > CS_THRESH & df$tbr > TBR_HIGH, na.rm=TRUE), "\n")
cat("Results saved to:", dir_sup, "\n")

sink(file.path(dir_sup, "S1_sessionInfo.txt")); sessionInfo(); sink()
