# 06_trd.R — Translational Rate Dynamics (IU-7)
#
# Input:  results/normalized/normalized_data.RData, results/tta/cluster_assignments.csv,
#         results/baseline/deseq2_lrt_translatome.csv
# Output: results/trd/
#
# Analyses:
#   1. Per-gene per-interval rate calculation (v_TE, v_M)
#   2. Global rate heatmap (ordered by TTA trajectory type)
#   3. Population-level rate distributions per interval
#   4. MI hidden window analysis (6→9h vs adjacent intervals)
#   5. Acceleration analysis (descriptive)
#   6. Compensation wave detection (TTA C1/C2 genes)
#
# Design note: time intervals are equal (3h each), so rate ∝ difference.
# "Rate" = group-mean ΔTE or ΔM per 3h interval (consistent with TPI's fold-change approach).

source("src/R/00_config.R")

library(matrixStats)
library(ggplot2)
library(reshape2)
library(clusterProfiler)
library(org.Mm.eg.db)

dir_trd <- file.path(DIR_RESULTS, "trd")
dir.create(dir_trd, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# 1. Load data
# ============================================================

load(file.path(DIR_RESULTS, "normalized", "normalized_data.RData"))
tta_assignments <- read.csv(file.path(DIR_RESULTS, "tta", "cluster_assignments.csv"),
                            stringsAsFactors = FALSE)
lrt_translat <- read.csv(file.path(DIR_RESULTS, "baseline", "deseq2_lrt_translatome.csv"),
                         stringsAsFactors = FALSE)
tpi_scores <- read.csv(file.path(DIR_RESULTS, "tpi", "tpi_scores.csv"),
                       stringsAsFactors = FALSE)

cat("Loaded:", length(gene_ids_filt), "genes,", nrow(pair_map), "pairs\n")

tp_levels   <- c("0h", "3h", "6h", "9h", "12h")
stage_labels <- c("GV", "GVBD", "MI-6", "MI-9", "MII")
intervals    <- list(c("0h","3h"), c("3h","6h"), c("6h","9h"), c("9h","12h"))
interval_labels <- c("GV→GVBD", "GVBD→MI6", "MI6→MI9", "MI9→MII")
delta_hours  <- 3  # constant interval width

# ============================================================
# 2. Compute group-mean TE and mRNA at each timepoint
# ============================================================

cat("Computing group means per timepoint...\n")

compute_tp_means <- function(norm_trans, te_mat, pair_map, gene_ids) {
  mean_mRNA <- matrix(NA_real_, nrow = length(gene_ids), ncol = length(tp_levels),
                      dimnames = list(gene_ids, tp_levels))
  mean_TE   <- matrix(NA_real_, nrow = length(gene_ids), ncol = length(tp_levels),
                      dimnames = list(gene_ids, tp_levels))

  for (i in seq_along(tp_levels)) {
    tp <- tp_levels[i]
    pairs_tp <- pair_map[pair_map$time_point == tp, ]

    m_cols <- pairs_tp$transcriptome_id
    mean_mRNA[, i] <- rowMeans(norm_trans[gene_ids, m_cols, drop = FALSE], na.rm = TRUE)

    te_cols <- grep(paste0("^", tp, "_rep"), colnames(te_mat), value = TRUE)
    stopifnot(length(te_cols) == nrow(pairs_tp))
    mean_TE[, i] <- rowMeans(te_mat[gene_ids, te_cols, drop = FALSE], na.rm = TRUE)
  }
  list(mean_mRNA = mean_mRNA, mean_TE = mean_TE)
}

tp_means <- compute_tp_means(norm_trans_A, te_A, pair_map, gene_ids_filt)

# ============================================================
# 3. Rate matrix (v_TE, v_M) — ΔTE/Δt per interval
# ============================================================

cat("Computing rate matrices...\n")

n_genes <- length(gene_ids_filt)
rate_TE <- matrix(NA_real_, nrow = n_genes, ncol = length(intervals),
                  dimnames = list(gene_ids_filt, interval_labels))
rate_mRNA <- matrix(NA_real_, nrow = n_genes, ncol = length(intervals),
                    dimnames = list(gene_ids_filt, interval_labels))

for (iv in seq_along(intervals)) {
  tp_from_idx <- which(tp_levels == intervals[[iv]][1])
  tp_to_idx   <- which(tp_levels == intervals[[iv]][2])
  rate_TE[, iv]   <- (tp_means$mean_TE[, tp_to_idx]   - tp_means$mean_TE[, tp_from_idx])   / delta_hours
  rate_mRNA[, iv] <- (tp_means$mean_mRNA[, tp_to_idx] - tp_means$mean_mRNA[, tp_from_idx]) / delta_hours
}

rate_df <- data.frame(
  gene_id = gene_ids_filt,
  gene_symbol = symbols_filt,
  rate_TE,
  rate_mRNA,
  stringsAsFactors = FALSE
)
colnames(rate_df)[3:6]  <- paste0("vTE_", interval_labels)
colnames(rate_df)[7:10] <- paste0("vM_", interval_labels)
write.csv(rate_df, file.path(dir_trd, "rate_matrix.csv"), row.names = FALSE)
cat("Rate matrix saved.\n")

# ============================================================
# 4. Rate heatmap (ordered by TTA trajectory type)
# ============================================================

cat("Generating rate heatmap...\n")

rate_with_traj <- merge(rate_df,
                        tta_assignments[, c("gene_id", "trajectory_type")],
                        by = "gene_id", all.x = TRUE)
rate_with_traj$trajectory_type[is.na(rate_with_traj$trajectory_type)] <- "Unknown"

traj_order <- c("TE-Only Activation", "Late Compensatory Buffering",
                 "Coordinated Clearance", "Deep Coordinated Clearance",
                 "Mild Coordinated Clearance", "Unknown")
rate_with_traj$trajectory_type <- factor(rate_with_traj$trajectory_type,
                                          levels = traj_order)
rate_with_traj <- rate_with_traj[order(rate_with_traj$trajectory_type), ]

te_cols_hm <- paste0("vTE_", interval_labels)
vTE_mat <- as.matrix(rate_with_traj[, te_cols_hm])

vTE_melt <- melt(cbind(
  data.frame(gene_rank = seq_len(nrow(rate_with_traj)),
             trajectory_type = rate_with_traj$trajectory_type),
  as.data.frame(vTE_mat)
), id.vars = c("gene_rank", "trajectory_type"), variable.name = "interval", value.name = "rate")
vTE_melt$interval <- factor(gsub("vTE_", "", as.character(vTE_melt$interval)),
                              levels = interval_labels)

# Use decile binning on rate for visualization (avoid outlier distortion)
rate_lim <- quantile(abs(vTE_mat), 0.95, na.rm = TRUE)

p_hm <- ggplot(vTE_melt, aes(x = interval, y = gene_rank, fill = rate)) +
  geom_raster() +
  scale_fill_gradient2(
    low = col_heatmap_diverging["low"], mid = col_heatmap_diverging["mid"],
    high = col_heatmap_diverging["high"],
    midpoint = 0, limits = c(-rate_lim, rate_lim), oob = scales::squish,
    name = "v_TE\n(Δlog2 TE/h)"
  ) +
  facet_grid(rows = vars(trajectory_type), scales = "free_y", space = "free_y") +
  labs(x = "Interval", y = "Gene (ordered by TTA type)",
       title = "TE Rate Dynamics by Trajectory Type") +
  theme_oocyte() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1))
ggsave(file.path(dir_trd, "rate_heatmap.pdf"), p_hm,
       width = 120, height = 200, units = "mm")

# ============================================================
# 5. Population-level rate distributions per interval
# ============================================================

cat("Generating interval rate distributions...\n")

vTE_pop <- melt(data.frame(gene_id = gene_ids_filt, as.data.frame(rate_TE),
                           check.names = FALSE, stringsAsFactors = FALSE),
                id.vars = "gene_id", variable.name = "interval", value.name = "rate")
vTE_pop$interval <- factor(as.character(vTE_pop$interval), levels = interval_labels)
vTE_pop$omics <- "TE"

vM_pop <- melt(data.frame(gene_id = gene_ids_filt, as.data.frame(rate_mRNA),
                           check.names = FALSE, stringsAsFactors = FALSE),
               id.vars = "gene_id", variable.name = "interval", value.name = "rate")
vM_pop$interval <- factor(as.character(vM_pop$interval), levels = interval_labels)
vM_pop$omics <- "mRNA"

rate_pop <- rbind(vTE_pop, vM_pop)
rate_pop$omics <- factor(rate_pop$omics, levels = c("mRNA", "TE"))

omics_fill <- c("mRNA" = col_omics[["Transcriptome"]],
                "TE"   = col_omics[["Translatome"]])

# Per-interval stats
rate_stats <- data.frame(
  interval = interval_labels,
  mean_vTE   = colMeans(rate_TE, na.rm = TRUE),
  sd_vTE     = apply(rate_TE, 2, sd, na.rm = TRUE),
  median_vTE = apply(rate_TE, 2, median, na.rm = TRUE),
  mean_vM    = colMeans(rate_mRNA, na.rm = TRUE),
  sd_vM      = apply(rate_mRNA, 2, sd, na.rm = TRUE),
  median_vM  = apply(rate_mRNA, 2, median, na.rm = TRUE),
  stringsAsFactors = FALSE
)
cat("\nInterval rate statistics (TE):\n")
print(rate_stats[, c("interval", "mean_vTE", "sd_vTE", "median_vTE")])
write.csv(rate_stats, file.path(dir_trd, "interval_rate_stats.csv"), row.names = FALSE)

# Violin plot
p_violin <- ggplot(rate_pop, aes(x = interval, y = rate, fill = omics)) +
  geom_violin(alpha = 0.7, trim = TRUE, scale = "width") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  scale_fill_manual(values = omics_fill) +
  coord_cartesian(ylim = c(-rate_lim * 1.2, rate_lim * 1.2)) +
  labs(x = "Interval", y = "Rate (Δlog2 / h)",
       title = "Population-Level Rate Distributions",
       fill = "Omics") +
  facet_wrap(~omics, nrow = 2) +
  theme_oocyte() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        legend.position = "none")
ggsave(file.path(dir_trd, "interval_rate_distributions.pdf"), p_violin,
       width = 170, height = 160, units = "mm")

# ============================================================
# 6. MI hidden window analysis (6→9h)
# ============================================================

cat("\nMI hidden window analysis (MI6→MI9)...\n")

mi_iv_idx <- 3  # index of MI6→MI9 in intervals

abs_rate_TE <- abs(rate_TE)
mi_abs_rate <- abs_rate_TE[, mi_iv_idx]

# LRT-significant translatome genes for threshold and gating
sig_translat_ids <- lrt_translat$gene_id[!is.na(lrt_translat$padj) &
                                          lrt_translat$padj < params$trd$mi_window_fdr]
sig_idx_te <- which(gene_ids_filt %in% sig_translat_ids)

# Wilcoxon: MI6→MI9 absolute rate vs adjacent intervals (complete cases only)
keep_prev <- complete.cases(abs_rate_TE[, mi_iv_idx], abs_rate_TE[, mi_iv_idx - 1])
keep_next <- complete.cases(abs_rate_TE[, mi_iv_idx], abs_rate_TE[, mi_iv_idx + 1])
wt_vs_prev <- wilcox.test(abs_rate_TE[keep_prev, mi_iv_idx],
                          abs_rate_TE[keep_prev, mi_iv_idx - 1],
                          paired = TRUE, alternative = "two.sided", conf.int = TRUE)
wt_vs_next <- wilcox.test(abs_rate_TE[keep_next, mi_iv_idx],
                          abs_rate_TE[keep_next, mi_iv_idx + 1],
                          paired = TRUE, alternative = "two.sided", conf.int = TRUE)

cat("  n paired (vs prev):", sum(keep_prev), "\n")
cat("  n paired (vs next):", sum(keep_next), "\n")
cat("  |vTE(MI6→MI9)| vs |vTE(GVBD→MI6)|: Wilcoxon p =",
    format(wt_vs_prev$p.value, digits = 3),
    ", HL shift =", round(wt_vs_prev$estimate, 4), "\n")
cat("  |vTE(MI6→MI9)| vs |vTE(MI9→MII)|: Wilcoxon p =",
    format(wt_vs_next$p.value, digits = 3),
    ", HL shift =", round(wt_vs_next$estimate, 4), "\n")

# Threshold on LRT-significant genes only (top 10% within sig genes)
mi_threshold <- quantile(mi_abs_rate[sig_idx_te], 0.90, na.rm = TRUE)
mi_window_genes_idx <- sig_idx_te[mi_abs_rate[sig_idx_te] > mi_threshold]
mi_window_gene_ids  <- gene_ids_filt[mi_window_genes_idx]
cat("  MI window genes (LRT sig, |vTE| top 10%):", length(mi_window_gene_ids), "\n")

mi_window_sig <- mi_window_gene_ids  # already LRT-filtered
cat("  MI window genes:", length(mi_window_sig), "\n")

mi_window_df <- data.frame(
  gene_id = mi_window_sig,
  gene_symbol = symbols_filt[match(mi_window_sig, gene_ids_filt)],
  vTE_MI6_MI9 = rate_TE[mi_window_sig, mi_iv_idx],
  direction = ifelse(rate_TE[mi_window_sig, mi_iv_idx] > 0, "increasing", "decreasing"),
  stringsAsFactors = FALSE
)
mi_window_df <- mi_window_df[order(abs(mi_window_df$vTE_MI6_MI9), decreasing = TRUE), ]
write.csv(mi_window_df, file.path(dir_trd, "mi_window_genes.csv"), row.names = FALSE)

cat("  MI window direction split:\n")
print(table(mi_window_df$direction))

# GO enrichment of MI window genes
cat("  GO enrichment of MI window genes...\n")
mi_genes_clean <- strip_ensembl_version(mi_window_sig)
bg_clean       <- strip_ensembl_version(gene_ids_filt)

mi_entrez <- tryCatch(
  clusterProfiler::bitr(mi_genes_clean, fromType = "ENSEMBL",
                        toType = "ENTREZID", OrgDb = org.Mm.eg.db),
  error = function(e) { message("WARN [mi_entrez bitr]: ", conditionMessage(e)); NULL }
)
bg_entrez <- tryCatch(
  clusterProfiler::bitr(bg_clean, fromType = "ENSEMBL",
                        toType = "ENTREZID", OrgDb = org.Mm.eg.db),
  error = function(e) { message("WARN [bg_entrez bitr]: ", conditionMessage(e)); NULL }
)

if (!is.null(mi_entrez) && nrow(mi_entrez) > 0 &&
    !is.null(bg_entrez) && nrow(bg_entrez) > 0) {
  mi_go <- tryCatch(
    clusterProfiler::enrichGO(
      gene          = mi_entrez$ENTREZID,
      universe      = bg_entrez$ENTREZID,
      OrgDb         = org.Mm.eg.db,
      ont           = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.05,
      minGSSize     = 10,
      readable      = TRUE
    ),
    error = function(e) { message("WARN [mi_go enrichGO]: ", conditionMessage(e)); NULL }
  )
  if (!is.null(mi_go) && nrow(as.data.frame(mi_go)) > 0) {
    mi_go_df <- as.data.frame(mi_go)
    write.csv(mi_go_df, file.path(dir_trd, "mi_window_go.csv"), row.names = FALSE)
    cat("  MI window GO terms:", nrow(mi_go_df), "\n")
    cat("  Top 5:\n")
    for (i in seq_len(min(5, nrow(mi_go_df)))) {
      cat("    ", mi_go_df$Description[i],
          " (padj=", format(mi_go_df$p.adjust[i], digits = 3), ")\n")
    }
  } else {
    cat("  No significant GO terms for MI window genes.\n")
    write.csv(data.frame(), file.path(dir_trd, "mi_window_go.csv"), row.names = FALSE)
  }
} else {
  cat("  WARN: Gene ID mapping failed for MI window GO analysis.\n")
  write.csv(data.frame(), file.path(dir_trd, "mi_window_go.csv"), row.names = FALSE)
}

# MI window visualization
mi_vis_df <- data.frame(
  interval = factor(rep(interval_labels, 2), levels = interval_labels),
  omics = rep(c("TE", "mRNA"), each = 4),
  mean_abs_rate = c(
    colMeans(abs_rate_TE, na.rm = TRUE),
    colMeans(abs(rate_mRNA), na.rm = TRUE)
  ),
  sd_rate = c(
    apply(abs_rate_TE, 2, sd, na.rm = TRUE),
    apply(abs(rate_mRNA), 2, sd, na.rm = TRUE)
  ),
  stringsAsFactors = FALSE
)
p_mi <- ggplot(mi_vis_df, aes(x = interval, y = mean_abs_rate, fill = omics,
                               ymin = mean_abs_rate - sd_rate,
                               ymax = mean_abs_rate + sd_rate)) +
  geom_col(position = "dodge", alpha = 0.8) +
  geom_errorbar(position = position_dodge(0.9), width = 0.2, color = "gray30") +
  scale_fill_manual(values = omics_fill) +
  labs(x = "Interval", y = "|Rate| (mean ± SD)",
       title = "MI Hidden Window: Absolute Rate by Interval",
       fill = "Omics") +
  theme_oocyte() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
ggsave(file.path(dir_trd, "mi_window_analysis.pdf"), p_mi,
       width = 170, height = 120, units = "mm")

# ============================================================
# 7. Acceleration analysis (descriptive)
# ============================================================

cat("\nAcceleration analysis (descriptive)...\n")

# a(t_i) = v(t_{i+1}) - v(t_i), 3 values per gene (intervals 2,3,4 relative to 1,2,3)
accel_TE <- matrix(NA_real_, nrow = n_genes, ncol = 3,
                   dimnames = list(gene_ids_filt,
                                   c("GVBD→MI6 accel", "MI6→MI9 accel", "MI9→MII accel")))
for (i in 1:3) {
  accel_TE[, i] <- rate_TE[, i + 1] - rate_TE[, i]
}

accel_stats <- data.frame(
  transition = c("GVBD→MI6 accel", "MI6→MI9 accel", "MI9→MII accel"),
  mean_accel = colMeans(accel_TE, na.rm = TRUE),
  median_accel = apply(accel_TE, 2, median, na.rm = TRUE),
  frac_pos = apply(accel_TE > 0, 2, mean, na.rm = TRUE)
)
cat("  Acceleration statistics (positive = rate increasing):\n")
print(accel_stats)
write.csv(accel_stats, file.path(dir_trd, "acceleration_stats.csv"), row.names = FALSE)

accel_melt <- melt(as.data.frame(accel_TE, check.names = FALSE),
                   variable.name = "transition", value.name = "acceleration")
accel_levels <- c("GVBD→MI6 accel", "MI6→MI9 accel", "MI9→MII accel")
accel_melt$transition <- factor(as.character(accel_melt$transition), levels = accel_levels)
p_accel <- ggplot(accel_melt, aes(x = transition, y = acceleration)) +
  geom_violin(fill = col_omics["Translatome"], alpha = 0.7, trim = TRUE) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  coord_cartesian(ylim = quantile(accel_melt$acceleration, c(0.01, 0.99), na.rm = TRUE)) +
  labs(x = "Transition", y = "TE rate acceleration",
       title = "TE Rate Acceleration (Descriptive)") +
  theme_oocyte() +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))
ggsave(file.path(dir_trd, "acceleration_visualization.pdf"), p_accel,
       width = 120, height = 100, units = "mm")

# ============================================================
# 8. Compensation wave detection
# ============================================================

cat("\nCompensation wave detection...\n")

comp_types <- c("TE-Only Activation", "Late Compensatory Buffering")
comp_genes_df <- tta_assignments[tta_assignments$trajectory_type %in% comp_types, ]
comp_gene_ids <- intersect(comp_genes_df$gene_id, gene_ids_filt)
cat("  Compensation genes (C1+C2):", length(comp_gene_ids), "\n")

# Define wave onset: interval with largest positive TE rate
# Genes with max(vTE) <= 0 are not "compensatory" — mark as Unclassified
rate_TE_comp <- rate_TE[comp_gene_ids, , drop = FALSE]
wave_onset <- apply(rate_TE_comp, 1, function(r) {
  if (all(is.na(r)) || max(r, na.rm = TRUE) <= 0) return(NA_integer_)
  which.max(r)
})
wave_labels <- c("Wave 1 (GV→GVBD)", "Wave 2 (GVBD→MI6)",
                  "Wave 3 (MI6→MI9)", "Wave 4 (MI9→MII)")
wave_assignment <- ifelse(is.na(wave_onset), "Unclassified",
                          wave_labels[wave_onset])

wave_dist <- table(wave_assignment)
cat("\n  Wave distribution:\n"); print(wave_dist)

wave_df <- data.frame(
  gene_id = comp_gene_ids,
  gene_symbol = symbols_filt[match(comp_gene_ids, gene_ids_filt)],
  trajectory_type = comp_genes_df$trajectory_type[
    match(comp_gene_ids, comp_genes_df$gene_id)],
  wave = wave_assignment,
  max_vTE_interval = ifelse(is.na(wave_onset), NA_character_, interval_labels[wave_onset]),
  max_vTE = apply(rate_TE_comp, 1,
                  function(r) if (all(is.na(r))) NA_real_ else max(r, na.rm = TRUE)),
  stringsAsFactors = FALSE
)
write.csv(wave_df, file.path(dir_trd, "compensation_waves.csv"), row.names = FALSE)

# GO enrichment per wave (only waves with >= min_genes)
min_wave_genes <- params$trd$min_genes_per_wave
cat("\n  GO enrichment per compensation wave:\n")

wave_go_list <- list()
bg_entrez_wave <- tryCatch(
  clusterProfiler::bitr(strip_ensembl_version(gene_ids_filt),
                        fromType = "ENSEMBL", toType = "ENTREZID",
                        OrgDb = org.Mm.eg.db),
  error = function(e) { message("WARN [bg_entrez_wave bitr]: ", conditionMessage(e)); NULL }
)

for (wv in wave_labels) {
  wv_genes <- wave_df$gene_id[wave_df$wave == wv]
  if (length(wv_genes) < min_wave_genes) {
    cat("  ", wv, ": n=", length(wv_genes), " (below min", min_wave_genes, ")\n")
    next
  }
  wv_entrez <- tryCatch(
    clusterProfiler::bitr(strip_ensembl_version(wv_genes),
                          fromType = "ENSEMBL", toType = "ENTREZID",
                          OrgDb = org.Mm.eg.db),
    error = function(e) { message("WARN [wv_entrez bitr ", wv, "]: ", conditionMessage(e)); NULL }
  )
  if (is.null(wv_entrez) || nrow(wv_entrez) == 0 ||
      is.null(bg_entrez_wave)) {
    cat("  ", wv, ": gene mapping failed\n"); next
  }
  wv_go <- tryCatch(
    clusterProfiler::enrichGO(
      gene          = wv_entrez$ENTREZID,
      universe      = bg_entrez_wave$ENTREZID,
      OrgDb         = org.Mm.eg.db,
      ont           = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.05,
      minGSSize     = 10,
      readable      = TRUE
    ),
    error = function(e) { message("WARN [wv_go enrichGO ", wv, "]: ", conditionMessage(e)); NULL }
  )
  if (!is.null(wv_go) && nrow(as.data.frame(wv_go)) > 0) {
    wv_go_df <- as.data.frame(wv_go)
    wv_go_df$wave <- wv
    wave_go_list[[wv]] <- wv_go_df
    cat("  ", wv, ": n=", length(wv_genes),
        ", GO terms=", nrow(wv_go_df),
        ", top=", wv_go_df$Description[1], "\n")
  } else {
    cat("  ", wv, ": n=", length(wv_genes), ", no enriched GO terms\n")
  }
}

if (length(wave_go_list) > 0) {
  wave_go_all <- do.call(rbind, wave_go_list)
  write.csv(wave_go_all, file.path(dir_trd, "compensation_wave_go.csv"),
            row.names = FALSE)
} else {
  write.csv(data.frame(), file.path(dir_trd, "compensation_wave_go.csv"),
            row.names = FALSE)
}

# Wave visualization
wave_count_df <- data.frame(
  wave = names(wave_dist),
  count = as.integer(wave_dist),
  stringsAsFactors = FALSE
)
wave_count_df <- wave_count_df[wave_count_df$wave != "Unclassified", ]
wave_count_df$wave <- factor(wave_count_df$wave, levels = wave_labels)

if (nrow(wave_count_df) > 0) {
  wave_fill <- c(
    "Wave 1 (GV→GVBD)"  = col_timepoint[["GVBD"]],
    "Wave 2 (GVBD→MI6)" = col_timepoint[["MI-6"]],
    "Wave 3 (MI6→MI9)"  = col_timepoint[["MI-9"]],
    "Wave 4 (MI9→MII)"  = col_timepoint[["MII"]]
  )
  p_wave <- ggplot(wave_count_df, aes(x = wave, y = count, fill = wave)) +
    geom_col(alpha = 0.85) +
    scale_fill_manual(values = wave_fill) +
    labs(x = "Compensation Wave (onset interval)",
         y = "Number of genes",
         title = "Translational Compensation Wave Distribution") +
    theme_oocyte() +
    theme(axis.text.x = element_text(angle = 25, hjust = 1),
          legend.position = "none")
  ggsave(file.path(dir_trd, "compensation_waves.pdf"), p_wave,
         width = 140, height = 100, units = "mm")
}

# ============================================================
# 9. Control gene validation
# ============================================================

cat("\nControl gene rate validation:\n")
ctrl_genes <- params$control_genes
ctrl_ids <- gene_ids_filt[symbols_filt %in% ctrl_genes]

for (gid in ctrl_ids) {
  gsym <- symbols_filt[gene_ids_filt == gid]
  vTE_vals  <- round(rate_TE[gid, ], 3)
  vM_vals   <- round(rate_mRNA[gid, ], 3)
  cat("  ", gsym, "\n")
  cat("    vTE:", paste(vTE_vals, collapse = " | "), "\n")
  cat("    vM: ", paste(vM_vals, collapse = " | "), "\n")
}

# ============================================================
# 10. Summary
# ============================================================

cat("\n=== TRD Summary ===\n")
cat("Interval with highest mean |vTE|:",
    interval_labels[which.max(colMeans(abs_rate_TE, na.rm = TRUE))], "\n")
cat("Interval with highest vTE variance:",
    interval_labels[which.max(apply(rate_TE, 2, var, na.rm = TRUE))], "\n")
cat("MI6→MI9 vs GVBD→MI6 (|vTE|): Wilcoxon p =",
    format(wt_vs_prev$p.value, digits = 3), "\n")
cat("MI6→MI9 vs MI9→MII (|vTE|): Wilcoxon p =",
    format(wt_vs_next$p.value, digits = 3), "\n")
cat("MI window genes (top 10% + LRT sig):", nrow(mi_window_df), "\n")
cat("Compensation waves with >=", min_wave_genes, "genes:",
    sum(as.integer(wave_dist[names(wave_dist) != "Unclassified"]) >= min_wave_genes), "\n")
cat("Results saved to:", dir_trd, "\n")

sink(file.path(dir_trd, "sessionInfo.txt"))
sessionInfo()
sink()
