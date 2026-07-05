# 03_baseline.R — Baseline Time-Series Analysis (IU-4)
#
# Input:  results/normalized/normalized_data.RData
# Output: results/baseline/ (DESeq2 LRT, soft clustering, endpoint comparison, control genes)
#
# Analyses:
#   1. DESeq2 LRT for time effect (transcriptome + translatome)
#   2. Polynomial time-series regression (maSigPro-style)
#   3. Fuzzy c-means soft clustering (Mfuzz-style via e1071)
#   4. GV vs MII endpoint comparison
#   5. Known gene control panel time courses

source("src/R/00_config.R")

library(DESeq2)
library(e1071)
library(cluster)
library(ggplot2)
library(matrixStats)
library(reshape2)

dir_base <- file.path(DIR_RESULTS, "baseline")
dir.create(dir_base, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# 1. Load normalized data
# ============================================================

load(file.path(DIR_RESULTS, "normalized", "normalized_data.RData"))
load(file.path(DIR_RESULTS, "qc", "qc_filtered_data.RData"))

cat(sprintf("Loaded: %d genes, %d transcriptome, %d translatome, %d pairs\n",
            nrow(counts_trans_filt), ncol(counts_trans_filt),
            ncol(counts_translat_filt), nrow(pair_map)))

# ============================================================
# 2. DESeq2 LRT — time-dependent expression
# ============================================================

run_deseq2_lrt <- function(counts_mat, sample_meta, omics_label) {
  sample_ids <- sample_meta$sample_id[sample_meta$omics_type == omics_label]
  sample_ids <- sample_ids[sample_ids %in% colnames(counts_mat)]
  sub_meta <- sample_meta[match(sample_ids, sample_meta$sample_id), ]
  sub_counts <- counts_mat[, sample_ids]

  dds <- DESeqDataSetFromMatrix(
    countData = sub_counts,
    colData = sub_meta,
    design = ~ stage
  )
  dds <- DESeq(dds, test = "LRT", reduced = ~ 1, quiet = TRUE)

  res <- results(dds, independentFiltering = TRUE)
  res_df <- as.data.frame(res)
  res_df$gene_id <- rownames(res_df)
  res_df$gene_symbol <- symbols_filt[match(res_df$gene_id, gene_ids_filt)]
  res_df
}

cat("Running DESeq2 LRT — Transcriptome...\n")
lrt_trans <- run_deseq2_lrt(counts_trans_filt, meta, "Transcriptome")

cat("Running DESeq2 LRT — Translatome...\n")
lrt_translat <- run_deseq2_lrt(counts_translat_filt, meta, "Translatome")

stopifnot(identical(lrt_trans$gene_id, lrt_translat$gene_id))

fdr_threshold <- params$baseline$deseq2_fdr

n_sig_trans <- sum(lrt_trans$padj < fdr_threshold, na.rm = TRUE)
n_sig_translat <- sum(lrt_translat$padj < fdr_threshold, na.rm = TRUE)

cat(sprintf("DESeq2 LRT significant genes (FDR < %.2f):\n", fdr_threshold))
cat(sprintf("  Transcriptome: %d / %d (%.1f%%)\n",
            n_sig_trans, nrow(lrt_trans), n_sig_trans / nrow(lrt_trans) * 100))
cat(sprintf("  Translatome:   %d / %d (%.1f%%)\n",
            n_sig_translat, nrow(lrt_translat), n_sig_translat / nrow(lrt_translat) * 100))

sig_both <- sum(lrt_trans$padj < fdr_threshold & lrt_translat$padj < fdr_threshold,
                na.rm = TRUE)
sig_trans_only <- sum(lrt_trans$padj < fdr_threshold &
                        (lrt_translat$padj >= fdr_threshold | is.na(lrt_translat$padj)),
                      na.rm = TRUE)
sig_translat_only <- sum((lrt_trans$padj >= fdr_threshold | is.na(lrt_trans$padj)) &
                           lrt_translat$padj < fdr_threshold,
                         na.rm = TRUE)

cat(sprintf("  Both: %d | Transcriptome-only: %d | Translatome-only: %d\n",
            sig_both, sig_trans_only, sig_translat_only))

write.csv(lrt_trans, file.path(dir_base, "deseq2_lrt_transcriptome.csv"),
          row.names = FALSE)
write.csv(lrt_translat, file.path(dir_base, "deseq2_lrt_translatome.csv"),
          row.names = FALSE)

# ============================================================
# 3. Polynomial time-series regression
# ============================================================

run_poly_timeseries <- function(norm_mat, sample_meta, omics_label,
                                degree = 2, fdr_threshold = 0.05) {
  sample_ids <- sample_meta$sample_id[sample_meta$omics_type == omics_label]
  sample_ids <- sample_ids[sample_ids %in% colnames(norm_mat)]
  sub_meta <- sample_meta[match(sample_ids, sample_meta$sample_id), ]
  sub_counts <- norm_mat[, sample_ids]

  time_numeric <- as.numeric(sub(
    "h$", "", as.character(sub_meta$time_point)))

  log_expr <- log2(sub_counts + 1)

  pvals <- numeric(nrow(log_expr))
  r_squared <- numeric(nrow(log_expr))

  for (i in seq_len(nrow(log_expr))) {
    y <- log_expr[i, ]
    fit_full <- lm(y ~ poly(time_numeric, degree))
    fit_null <- lm(y ~ 1)
    anova_res <- anova(fit_null, fit_full)
    pvals[i] <- anova_res$`Pr(>F)`[2]
    r_squared[i] <- summary(fit_full)$r.squared
  }

  padj <- p.adjust(pvals, method = "BH")

  data.frame(
    gene_id = rownames(log_expr),
    gene_symbol = symbols_filt[match(rownames(log_expr), gene_ids_filt)],
    pvalue = pvals,
    padj = padj,
    r_squared = r_squared,
    significant = padj < fdr_threshold,
    stringsAsFactors = FALSE
  )
}

poly_degree <- params$baseline$masigpro_degree
# 5 timepoints with degree=2 leaves only 2 residual df — interpret R² with caution
if (poly_degree >= 2) {
  cat(sprintf("WARNING: degree=%d polynomial on 5 timepoints — risk of overfitting (residual df = %d)\n",
              poly_degree, 5 - poly_degree - 1))
}

cat(sprintf("Running polynomial regression (degree=%d) — Transcriptome...\n", poly_degree))
poly_trans <- run_poly_timeseries(norm_trans_A, meta, "Transcriptome",
                                   degree = poly_degree, fdr_threshold = fdr_threshold)

cat(sprintf("Running polynomial regression (degree=%d) — Translatome...\n", poly_degree))
poly_translat <- run_poly_timeseries(norm_translat_A, meta, "Translatome",
                                      degree = poly_degree, fdr_threshold = fdr_threshold)

n_poly_trans <- sum(poly_trans$significant, na.rm = TRUE)
n_poly_translat <- sum(poly_translat$significant, na.rm = TRUE)

med_r2_trans <- if (n_poly_trans > 0) median(poly_trans$r_squared[poly_trans$significant & !is.na(poly_trans$significant)]) else NA
med_r2_translat <- if (n_poly_translat > 0) median(poly_translat$r_squared[poly_translat$significant & !is.na(poly_translat$significant)]) else NA

cat(sprintf("Polynomial regression significant (FDR < %.2f):\n", fdr_threshold))
cat(sprintf("  Transcriptome: %d (median R²=%.3f for significant)\n",
            n_poly_trans, med_r2_trans))
cat(sprintf("  Translatome:   %d (median R²=%.3f for significant)\n",
            n_poly_translat, med_r2_translat))

write.csv(poly_trans, file.path(dir_base, "poly_regression_transcriptome.csv"),
          row.names = FALSE)
write.csv(poly_translat, file.path(dir_base, "poly_regression_translatome.csv"),
          row.names = FALSE)

# Concordance: LRT vs polynomial
lrt_sig_ids <- lrt_trans$gene_id[lrt_trans$padj < fdr_threshold & !is.na(lrt_trans$padj)]
poly_sig_ids <- poly_trans$gene_id[poly_trans$significant & !is.na(poly_trans$significant)]
union_trans <- union(lrt_sig_ids, poly_sig_ids)
concordance_trans <- if (length(union_trans) == 0) 0 else {
  length(intersect(lrt_sig_ids, poly_sig_ids)) / length(union_trans)
}

lrt_sig_ids_tl <- lrt_translat$gene_id[lrt_translat$padj < fdr_threshold &
                                          !is.na(lrt_translat$padj)]
poly_sig_ids_tl <- poly_translat$gene_id[poly_translat$significant & !is.na(poly_translat$significant)]
union_translat <- union(lrt_sig_ids_tl, poly_sig_ids_tl)
concordance_translat <- if (length(union_translat) == 0) 0 else {
  length(intersect(lrt_sig_ids_tl, poly_sig_ids_tl)) / length(union_translat)
}

cat(sprintf("LRT vs poly (degree=%d) Jaccard concordance: transcriptome=%.3f, translatome=%.3f\n",
            poly_degree, concordance_trans, concordance_translat))

# Sensitivity: degree=1 (linear trend) for comparison
cat("Running polynomial regression (degree=1) as sensitivity check...\n")
poly1_trans <- run_poly_timeseries(norm_trans_A, meta, "Transcriptome",
                                    degree = 1, fdr_threshold = fdr_threshold)
poly1_translat <- run_poly_timeseries(norm_translat_A, meta, "Translatome",
                                       degree = 1, fdr_threshold = fdr_threshold)

n_poly1_trans <- sum(poly1_trans$significant, na.rm = TRUE)
n_poly1_translat <- sum(poly1_translat$significant, na.rm = TRUE)

poly1_sig_ids <- poly1_trans$gene_id[poly1_trans$significant & !is.na(poly1_trans$significant)]
union1_trans <- union(lrt_sig_ids, poly1_sig_ids)
concordance1_trans <- if (length(union1_trans) == 0) 0 else {
  length(intersect(lrt_sig_ids, poly1_sig_ids)) / length(union1_trans)
}

poly1_sig_ids_tl <- poly1_translat$gene_id[poly1_translat$significant & !is.na(poly1_translat$significant)]
union1_translat <- union(lrt_sig_ids_tl, poly1_sig_ids_tl)
concordance1_translat <- if (length(union1_translat) == 0) 0 else {
  length(intersect(lrt_sig_ids_tl, poly1_sig_ids_tl)) / length(union1_translat)
}

cat(sprintf("  degree=1: %d trans, %d translat significant\n", n_poly1_trans, n_poly1_translat))
cat(sprintf("  LRT vs poly (degree=1) Jaccard: transcriptome=%.3f, translatome=%.3f\n",
            concordance1_trans, concordance1_translat))

# degree=1 vs degree=2 concordance
poly12_jaccard_trans <- {
  u <- union(poly_sig_ids, poly1_sig_ids)
  if (length(u) == 0) 0 else length(intersect(poly_sig_ids, poly1_sig_ids)) / length(u)
}
cat(sprintf("  degree=1 vs degree=2 Jaccard (transcriptome): %.3f\n", poly12_jaccard_trans))

# ============================================================
# 4. Fuzzy c-means soft clustering
# ============================================================

prepare_expression_profiles <- function(norm_mat, sample_meta, omics_label,
                                        sig_genes) {
  sample_ids <- sample_meta$sample_id[sample_meta$omics_type == omics_label]
  sample_ids <- sample_ids[sample_ids %in% colnames(norm_mat)]
  sub_meta <- sample_meta[match(sample_ids, sample_meta$sample_id), ]
  sub_counts <- norm_mat[sig_genes, sample_ids]

  log_expr <- log2(sub_counts + 1)

  stages <- levels(sub_meta$stage)
  prof_mat <- matrix(NA, nrow = length(sig_genes), ncol = length(stages),
                     dimnames = list(sig_genes, stages))

  for (s in seq_along(stages)) {
    cols <- sub_meta$sample_id[sub_meta$stage == stages[s]]
    prof_mat[, s] <- rowMeans(log_expr[, cols, drop = FALSE])
  }

  # Standardize per gene (z-score across timepoints)
  row_means <- rowMeans(prof_mat)
  row_sds <- apply(prof_mat, 1, sd)
  keep <- row_sds > 0
  n_dropped <- sum(!keep)
  if (n_dropped > 0) {
    cat(sprintf("  Dropped %d genes with zero variance across timepoints\n", n_dropped))
  }
  prof_mat <- prof_mat[keep, ]
  row_means <- row_means[keep]
  row_sds <- row_sds[keep]
  prof_std <- (prof_mat - row_means) / row_sds

  prof_std
}

run_fuzzy_clustering <- function(prof_std, cluster_range, m = 1.25) {
  results <- list()
  best_score <- -Inf
  best_c <- cluster_range[1]

  for (c_val in cluster_range) {
    set.seed(params$seed)
    cl <- cmeans(prof_std, centers = c_val, m = m,
                 iter.max = 200, dist = "euclidean")

    # Fuzzy partition coefficient (FPC): higher = better separation
    fpc <- sum(cl$membership^2) / nrow(cl$membership)

    # Xie-Beni index: lower = better
    dists <- as.matrix(dist(rbind(prof_std, cl$centers)))
    n <- nrow(prof_std)
    numerator <- sum(cl$membership^2 *
                       dists[1:n, (n + 1):(n + c_val)]^2)
    min_center_dist <- min(dist(cl$centers))^2
    xb <- numerator / (n * max(min_center_dist, .Machine$double.eps))

    results[[as.character(c_val)]] <- list(
      clustering = cl,
      fpc = fpc,
      xie_beni = xb,
      c = c_val
    )

    # Silhouette on hard cluster assignments (subsample for speed if n > 2000)
    n_sil <- nrow(prof_std)
    if (n_sil > 2000) {
      set.seed(params$seed)
      sil_idx <- sample(n_sil, 2000)
      sil <- silhouette(cl$cluster[sil_idx], dist(prof_std[sil_idx, ]))
    } else {
      sil <- silhouette(cl$cluster, dist(prof_std))
    }
    avg_sil <- mean(sil[, "sil_width"])

    results[[as.character(c_val)]]$silhouette <- avg_sil

    cat(sprintf("  c=%d: FPC=%.4f, Xie-Beni=%.2f, Silhouette=%.3f\n",
                c_val, fpc, xb, avg_sil))

    if (fpc > best_score) {
      best_score <- fpc
      best_c <- c_val
    }
  }

  list(results = results, best_c = best_c)
}

# Cluster transcriptome time-dependent genes
sig_trans_genes <- lrt_trans$gene_id[lrt_trans$padj < fdr_threshold &
                                       !is.na(lrt_trans$padj)]
sig_translat_genes <- lrt_translat$gene_id[lrt_translat$padj < fdr_threshold &
                                             !is.na(lrt_translat$padj)]

cat(sprintf("\nFuzzy c-means clustering (m=%.2f):\n", params$baseline$mfuzz_m))

cat("Transcriptome profiles:\n")
prof_trans <- prepare_expression_profiles(norm_trans_A, meta, "Transcriptome",
                                           sig_trans_genes)
cat(sprintf("  %d genes with nonzero variance for clustering\n", nrow(prof_trans)))

fuzz_trans <- run_fuzzy_clustering(prof_trans,
                                    params$baseline$mfuzz_clusters,
                                    params$baseline$mfuzz_m)

cat("Translatome profiles:\n")
prof_translat <- prepare_expression_profiles(norm_translat_A, meta, "Translatome",
                                              sig_translat_genes)
cat(sprintf("  %d genes with nonzero variance for clustering\n", nrow(prof_translat)))

fuzz_translat <- run_fuzzy_clustering(prof_translat,
                                       params$baseline$mfuzz_clusters,
                                       params$baseline$mfuzz_m)

cat(sprintf("Best c: transcriptome=%d, translatome=%d\n",
            fuzz_trans$best_c, fuzz_translat$best_c))

# Extract cluster assignments at best c
best_cl_trans <- fuzz_trans$results[[as.character(fuzz_trans$best_c)]]$clustering
best_cl_translat <- fuzz_translat$results[[as.character(fuzz_translat$best_c)]]$clustering

cluster_trans_df <- data.frame(
  gene_id = names(best_cl_trans$cluster),
  cluster = best_cl_trans$cluster,
  max_membership = apply(best_cl_trans$membership, 1, max),
  stringsAsFactors = FALSE
)
cluster_trans_df$gene_symbol <- symbols_filt[match(cluster_trans_df$gene_id, gene_ids_filt)]

cluster_translat_df <- data.frame(
  gene_id = names(best_cl_translat$cluster),
  cluster = best_cl_translat$cluster,
  max_membership = apply(best_cl_translat$membership, 1, max),
  stringsAsFactors = FALSE
)
cluster_translat_df$gene_symbol <- symbols_filt[match(cluster_translat_df$gene_id, gene_ids_filt)]

write.csv(cluster_trans_df, file.path(dir_base, "mfuzz_clusters_transcriptome.csv"),
          row.names = FALSE)
write.csv(cluster_translat_df, file.path(dir_base, "mfuzz_clusters_translatome.csv"),
          row.names = FALSE)

# Plot cluster centers
plot_cluster_centers <- function(cl_result, title_suffix) {
  centers <- cl_result$centers
  stages <- colnames(centers)
  time_numeric <- c(0, 3, 6, 9, 12)

  plot_df <- data.frame()
  for (k in seq_len(nrow(centers))) {
    n_genes <- sum(cl_result$cluster == k)
    plot_df <- rbind(plot_df, data.frame(
      cluster = paste0("C", k, " (n=", n_genes, ")"),
      time = time_numeric,
      stage = stages,
      expression = centers[k, ],
      stringsAsFactors = FALSE
    ))
  }

  p <- ggplot(plot_df, aes(x = time, y = expression, color = cluster)) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    scale_x_continuous(breaks = time_numeric, labels = stages) +
    labs(title = paste("Fuzzy cluster centers —", title_suffix),
         x = "Stage", y = "Standardized expression (z-score)") +
    theme_oocyte() +
    theme(legend.position = "right")
  p
}

p_cl_trans <- plot_cluster_centers(best_cl_trans, "Transcriptome")
p_cl_translat <- plot_cluster_centers(best_cl_translat, "Translatome")

ggsave(file.path(dir_base, "mfuzz_centers_transcriptome.pdf"),
       p_cl_trans, width = 8, height = 5)
ggsave(file.path(dir_base, "mfuzz_centers_translatome.pdf"),
       p_cl_translat, width = 8, height = 5)

# ============================================================
# 5. GV vs MII endpoint comparison
# ============================================================

cat("\nGV vs MII endpoint comparison (Wald test)...\n")

run_endpoint_test <- function(counts_mat, sample_meta, omics_label) {
  sample_ids <- sample_meta$sample_id[sample_meta$omics_type == omics_label &
                                        sample_meta$stage %in% c("GV", "MII")]
  sample_ids <- sample_ids[sample_ids %in% colnames(counts_mat)]
  sub_meta <- sample_meta[match(sample_ids, sample_meta$sample_id), ]
  sub_meta$stage <- droplevels(sub_meta$stage)
  sub_counts <- counts_mat[, sample_ids]

  dds <- DESeqDataSetFromMatrix(
    countData = sub_counts,
    colData = sub_meta,
    design = ~ stage
  )
  dds <- DESeq(dds, quiet = TRUE)
  res <- results(dds, contrast = c("stage", "MII", "GV"))
  res_df <- as.data.frame(res)
  res_df$gene_id <- rownames(res_df)
  res_df$gene_symbol <- symbols_filt[match(res_df$gene_id, gene_ids_filt)]
  res_df
}

endpoint_trans <- run_endpoint_test(counts_trans_filt, meta, "Transcriptome")
endpoint_translat <- run_endpoint_test(counts_translat_filt, meta, "Translatome")

n_up_trans <- sum(endpoint_trans$padj < fdr_threshold &
                    endpoint_trans$log2FoldChange > 0, na.rm = TRUE)
n_down_trans <- sum(endpoint_trans$padj < fdr_threshold &
                      endpoint_trans$log2FoldChange < 0, na.rm = TRUE)
n_up_translat <- sum(endpoint_translat$padj < fdr_threshold &
                       endpoint_translat$log2FoldChange > 0, na.rm = TRUE)
n_down_translat <- sum(endpoint_translat$padj < fdr_threshold &
                         endpoint_translat$log2FoldChange < 0, na.rm = TRUE)

cat(sprintf("GV→MII significant (FDR < %.2f):\n", fdr_threshold))
cat(sprintf("  Transcriptome: %d up, %d down\n", n_up_trans, n_down_trans))
cat(sprintf("  Translatome:   %d up, %d down\n", n_up_translat, n_down_translat))

# Differential TE: genes where translatome changes differently from transcriptome
common_genes <- intersect(endpoint_trans$gene_id, endpoint_translat$gene_id)
delta_te <- endpoint_translat$log2FoldChange[match(common_genes, endpoint_translat$gene_id)] -
  endpoint_trans$log2FoldChange[match(common_genes, endpoint_trans$gene_id)]
names(delta_te) <- common_genes

# Require both omics to be individually significant for delta_TE to be meaningful
sig_both_endpoint <- endpoint_trans$padj[match(common_genes, endpoint_trans$gene_id)] < fdr_threshold &
  endpoint_translat$padj[match(common_genes, endpoint_translat$gene_id)] < fdr_threshold
sig_both_endpoint[is.na(sig_both_endpoint)] <- FALSE

n_te_up <- sum(abs(delta_te) > 1 & delta_te > 0 & sig_both_endpoint, na.rm = TRUE)
n_te_down <- sum(abs(delta_te) > 1 & delta_te < 0 & sig_both_endpoint, na.rm = TRUE)

cat(sprintf("  Differential TE (|ΔTE| > 1): %d TE-up, %d TE-down\n",
            n_te_up, n_te_down))

endpoint_combined <- data.frame(
  gene_id = common_genes,
  gene_symbol = symbols_filt[match(common_genes, gene_ids_filt)],
  lfc_transcriptome = endpoint_trans$log2FoldChange[match(common_genes, endpoint_trans$gene_id)],
  padj_transcriptome = endpoint_trans$padj[match(common_genes, endpoint_trans$gene_id)],
  lfc_translatome = endpoint_translat$log2FoldChange[match(common_genes, endpoint_translat$gene_id)],
  padj_translatome = endpoint_translat$padj[match(common_genes, endpoint_translat$gene_id)],
  delta_te = delta_te,
  sig_both = sig_both_endpoint,
  stringsAsFactors = FALSE
)

write.csv(endpoint_combined, file.path(dir_base, "gv_vs_mii_endpoint.csv"),
          row.names = FALSE)

# Per-adjacent-transition ΔTE using paired TE matrix
cat("\nPer-transition ΔTE analysis (paired TE matrix):\n")
stages <- c("GV", "GVBD", "MI-6", "MI-9", "MII")
tp_labels <- c("0h", "3h", "6h", "9h", "12h")
transition_stats <- data.frame()

for (t_idx in seq_len(length(stages) - 1)) {
  from_tp <- tp_labels[t_idx]
  to_tp <- tp_labels[t_idx + 1]
  from_cols <- grep(paste0("^", from_tp, "_"), colnames(te_A), value = TRUE)
  to_cols <- grep(paste0("^", to_tp, "_"), colnames(te_A), value = TRUE)

  if (length(from_cols) == 0 || length(to_cols) == 0) next

  te_from <- rowMeans(te_A[, from_cols, drop = FALSE], na.rm = TRUE)
  te_to <- rowMeans(te_A[, to_cols, drop = FALSE], na.rm = TRUE)
  te_change <- te_to - te_from

  n_up <- sum(te_change > 1, na.rm = TRUE)
  n_down <- sum(te_change < -1, na.rm = TRUE)
  med_change <- median(te_change, na.rm = TRUE)

  transition_stats <- rbind(transition_stats, data.frame(
    transition = paste0(stages[t_idx], "→", stages[t_idx + 1]),
    n_te_up = n_up,
    n_te_down = n_down,
    median_delta_te = med_change,
    stringsAsFactors = FALSE
  ))

  cat(sprintf("  %s→%s: %d TE-up, %d TE-down (median ΔTE=%.3f)\n",
              stages[t_idx], stages[t_idx + 1], n_up, n_down, med_change))
}

write.csv(transition_stats, file.path(dir_base, "per_transition_delta_te.csv"),
          row.names = FALSE)

# Scatter plot: transcriptome LFC vs translatome LFC
plot_data_scatter <- endpoint_combined[!is.na(endpoint_combined$lfc_transcriptome) &
                                         !is.na(endpoint_combined$lfc_translatome), ]
cat(sprintf("  Scatter plot: %d genes with valid LFC in both omics (%d excluded with NA)\n",
            nrow(plot_data_scatter), nrow(endpoint_combined) - nrow(plot_data_scatter)))

p_endpoint <- ggplot(plot_data_scatter,
                      aes(x = lfc_transcriptome, y = lfc_translatome)) +
  geom_point(size = 0.5, alpha = 0.3, color = "gray50") +
  geom_point(data = plot_data_scatter[abs(plot_data_scatter$delta_te) > 1 &
                                        plot_data_scatter$sig_both &
                                        !is.na(plot_data_scatter$delta_te), ],
             aes(color = ifelse(delta_te > 0, "TE-up", "TE-down")),
             size = 1, alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray30") +
  scale_color_manual(values = c("TE-up" = unname(col_tpi["TE-leading"]),
                                 "TE-down" = unname(col_tpi["mRNA-leading"])),
                      name = "Differential TE") +
  labs(title = "GV → MII: Transcriptome vs Translatome fold change",
       x = "Transcriptome log2FC (MII/GV)",
       y = "Translatome log2FC (MII/GV)") +
  theme_oocyte()

ggsave(file.path(dir_base, "gv_vs_mii_scatter.pdf"), p_endpoint, width = 7, height = 6)

# ============================================================
# 6. Known gene control panel — dual time courses
# ============================================================

plot_known_gene <- function(gene_name, norm_trans, norm_translat,
                             te_mat, sample_meta) {
  gene_row <- which(symbols_filt == gene_name)
  if (length(gene_row) == 0) return(NULL)

  gene_id <- gene_ids_filt[gene_row[1]]
  stages <- c("GV", "GVBD", "MI-6", "MI-9", "MII")
  time_numeric <- c(0, 3, 6, 9, 12)

  # Transcriptome expression per stage
  plot_data <- data.frame()
  for (s_idx in seq_along(stages)) {
    s <- stages[s_idx]
    trans_cols <- sample_meta$sample_id[sample_meta$stage == s &
                                          sample_meta$omics_type == "Transcriptome"]
    trans_cols <- trans_cols[trans_cols %in% colnames(norm_trans)]
    translat_cols <- sample_meta$sample_id[sample_meta$stage == s &
                                             sample_meta$omics_type == "Translatome"]
    translat_cols <- translat_cols[translat_cols %in% colnames(norm_translat)]

    if (length(trans_cols) > 0) {
      vals <- log2(norm_trans[gene_id, trans_cols] + 1)
      plot_data <- rbind(plot_data, data.frame(
        gene = gene_name, time = time_numeric[s_idx], stage = s,
        value = vals, layer = "Transcriptome", stringsAsFactors = FALSE))
    }
    if (length(translat_cols) > 0) {
      vals <- log2(norm_translat[gene_id, translat_cols] + 1)
      plot_data <- rbind(plot_data, data.frame(
        gene = gene_name, time = time_numeric[s_idx], stage = s,
        value = vals, layer = "Translatome", stringsAsFactors = FALSE))
    }

    # TE from paired matrix
    tp_label <- c("0h", "3h", "6h", "9h", "12h")[s_idx]
    te_cols <- grep(paste0("^", tp_label, "_"), colnames(te_mat))
    if (length(te_cols) > 0 && gene_id %in% rownames(te_mat)) {
      te_vals <- te_mat[gene_id, te_cols]
      plot_data <- rbind(plot_data, data.frame(
        gene = gene_name, time = time_numeric[s_idx], stage = s,
        value = te_vals, layer = "TE", stringsAsFactors = FALSE))
    }
  }

  if (nrow(plot_data) == 0) return(NULL)

  plot_data$panel <- ifelse(plot_data$layer == "TE", "TE", "Expression")

  p <- ggplot(plot_data, aes(x = time, y = value, color = layer)) +
    stat_summary(fun = mean, geom = "line", linewidth = 0.8) +
    stat_summary(fun = mean, geom = "point", size = 2) +
    stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.4, alpha = 0.5) +
    facet_wrap(~ panel, scales = "free_y", ncol = 1) +
    scale_x_continuous(breaks = time_numeric, labels = stages) +
    scale_color_manual(values = c("Transcriptome" = col_omics["Transcriptome"],
                                   "Translatome" = col_omics["Translatome"],
                                   "TE" = "#F0B323")) +
    labs(title = gene_name, x = "Stage", y = NULL) +
    theme_oocyte()
  p
}

control_genes <- params$control_genes
ctrl_plots <- list()
for (g in control_genes) {
  p <- plot_known_gene(g, norm_trans_A, norm_translat_A, te_A, meta)
  if (!is.null(p)) ctrl_plots[[g]] <- p
}

cat(sprintf("Control gene plots generated: %d/%d\n",
            length(ctrl_plots), length(control_genes)))

# Save individual plots
for (g in names(ctrl_plots)) {
  ggsave(file.path(dir_base, paste0("ctrl_", g, ".pdf")),
         ctrl_plots[[g]], width = 5, height = 4)
}

# ============================================================
# 7. Control gene statistics table
# ============================================================

ctrl_stats <- data.frame()
for (g in control_genes) {
  gene_row <- which(symbols_filt == g)
  if (length(gene_row) == 0) next
  gene_id <- gene_ids_filt[gene_row[1]]

  lrt_p_trans <- lrt_trans$padj[lrt_trans$gene_id == gene_id]
  lrt_p_translat <- lrt_translat$padj[lrt_translat$gene_id == gene_id]
  ep_lfc_trans <- endpoint_trans$log2FoldChange[endpoint_trans$gene_id == gene_id]
  ep_lfc_translat <- endpoint_translat$log2FoldChange[endpoint_translat$gene_id == gene_id]
  ep_dte <- delta_te[gene_id]

  # Cluster assignment
  cl_trans <- if (gene_id %in% cluster_trans_df$gene_id) {
    cluster_trans_df$cluster[cluster_trans_df$gene_id == gene_id]
  } else { NA }
  cl_translat <- if (gene_id %in% cluster_translat_df$gene_id) {
    cluster_translat_df$cluster[cluster_translat_df$gene_id == gene_id]
  } else { NA }

  ctrl_stats <- rbind(ctrl_stats, data.frame(
    gene = g,
    gene_id = gene_id,
    lrt_padj_trans = if (length(lrt_p_trans) > 0) lrt_p_trans else NA,
    lrt_padj_translat = if (length(lrt_p_translat) > 0) lrt_p_translat else NA,
    gv_mii_lfc_trans = if (length(ep_lfc_trans) > 0) ep_lfc_trans else NA,
    gv_mii_lfc_translat = if (length(ep_lfc_translat) > 0) ep_lfc_translat else NA,
    delta_te = if (length(ep_dte) > 0) ep_dte else NA,
    cluster_trans = cl_trans,
    cluster_translat = cl_translat,
    stringsAsFactors = FALSE
  ))
}

write.csv(ctrl_stats, file.path(dir_base, "control_gene_stats.csv"),
          row.names = FALSE)

cat("\nControl gene summary:\n")
for (i in seq_len(nrow(ctrl_stats))) {
  cat(sprintf("  %s: LRT trans p=%.1e, translat p=%.1e | GV→MII LFC: trans=%.2f, translat=%.2f | ΔTE=%.2f\n",
              ctrl_stats$gene[i],
              ctrl_stats$lrt_padj_trans[i],
              ctrl_stats$lrt_padj_translat[i],
              ctrl_stats$gv_mii_lfc_trans[i],
              ctrl_stats$gv_mii_lfc_translat[i],
              ctrl_stats$delta_te[i]))
}

# ============================================================
# 8. Summary
# ============================================================

cat(sprintf("\n=== Baseline Analysis Summary ===\n"))
cat(sprintf("DESeq2 LRT (FDR < %.2f): %d transcriptome, %d translatome\n",
            fdr_threshold, n_sig_trans, n_sig_translat))
cat(sprintf("Polynomial regression (degree=%d): %d transcriptome, %d translatome\n",
            poly_degree, n_poly_trans, n_poly_translat))
cat(sprintf("Polynomial regression (degree=1): %d transcriptome, %d translatome\n",
            n_poly1_trans, n_poly1_translat))
cat(sprintf("LRT vs poly Jaccard: degree=%d: trans=%.3f translat=%.3f | degree=1: trans=%.3f translat=%.3f\n",
            poly_degree, concordance_trans, concordance_translat,
            concordance1_trans, concordance1_translat))
cat(sprintf("Fuzzy clusters (best c): transcriptome=%d, translatome=%d\n",
            fuzz_trans$best_c, fuzz_translat$best_c))
for (cv in as.character(params$baseline$mfuzz_clusters)) {
  sil_t <- fuzz_trans$results[[cv]]$silhouette
  sil_tl <- fuzz_translat$results[[cv]]$silhouette
  cat(sprintf("  c=%s: Silhouette trans=%.3f translat=%.3f\n", cv, sil_t, sil_tl))
}
cat(sprintf("GV→MII: %d up + %d down (trans), %d up + %d down (translat)\n",
            n_up_trans, n_down_trans, n_up_translat, n_down_translat))
cat(sprintf("Differential TE (|ΔTE| > 1, both sig): %d up, %d down\n", n_te_up, n_te_down))
cat(sprintf("Results saved to: %s\n", dir_base))

# ============================================================
# 9. Session info
# ============================================================

writeLines(capture.output(sessionInfo()),
           file.path(dir_base, "sessionInfo.txt"))
