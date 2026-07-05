# 02_normalize.R — Normalization & Sensitivity Analysis (IU-3)
#
# Input:  results/qc/qc_filtered_data.RData
# Output: results/normalized/ (normalized counts, TE matrices, sensitivity plots)
#
# Four normalization schemes:
#   A: 147 constGenes reference (primary)
#   B: DESeq2 median-of-ratios
#   C: TMM (edgeR)
#   D: RUVg (unwanted variation removal using constGenes)

source("src/R/00_config.R")

library(DESeq2)
library(edgeR)
library(ggplot2)
library(matrixStats)
library(reshape2)

dir_norm <- file.path(DIR_RESULTS, "normalized")
dir.create(dir_norm, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# 1. Load QC-filtered data
# ============================================================

load(file.path(DIR_RESULTS, "qc", "qc_filtered_data.RData"))
constgene_ids <- load_constgenes()

cat(sprintf("Loaded filtered data: %d genes, %d transcriptome samples, %d translatome samples\n",
            nrow(counts_trans_filt), ncol(counts_trans_filt), ncol(counts_translat_filt)))

pair_map <- build_pair_map(meta)
cat(sprintf("Complete pairs: %d\n", nrow(pair_map)))

# ============================================================
# 2. Map constGenes to filtered gene set
# ============================================================

gene_ids_base_filt <- strip_ensembl_version(gene_ids_filt)
constgene_idx <- which(gene_ids_base_filt %in% constgene_ids)

cat(sprintf("constGenes mapped to filtered set: %d / %d\n",
            length(constgene_idx), length(constgene_ids)))

if (length(constgene_idx) < 50) {
  stop("Too few constGenes in filtered set — check gene ID mapping")
}

# ============================================================
# 3. Validate constGene stability (CV filter)
# ============================================================

# CV on raw counts (not CPM): CPM inflates CV for stable genes when 80% mRNA degrades,
# because the denominator (total library) changes dramatically. Raw count CV, while
# affected by sequencing depth, better reflects absolute expression stability.
compute_constgene_cv <- function(counts_mat, sample_meta, constgene_rows, omics_label) {
  stopifnot(is.factor(sample_meta$stage))

  constgene_counts <- counts_mat[constgene_rows, , drop = FALSE]

  stages <- levels(sample_meta$stage)
  stage_means <- matrix(NA, nrow = nrow(constgene_counts), ncol = length(stages),
                        dimnames = list(rownames(constgene_counts), stages))

  for (s in stages) {
    cols <- sample_meta$sample_id[sample_meta$stage == s &
                                    sample_meta$omics_type == omics_label]
    cols <- cols[cols %in% colnames(constgene_counts)]
    if (length(cols) > 0) {
      stage_means[, s] <- rowMeans(constgene_counts[, cols, drop = FALSE])
    }
  }

  cv_vals <- apply(stage_means, 1, function(x) {
    x <- x[!is.na(x)]
    if (mean(x) == 0) return(Inf)
    sd(x) / mean(x)
  })

  cv_vals
}

meta_trans <- meta[meta$omics_type == "Transcriptome", ]
meta_translat <- meta[meta$omics_type == "Translatome", ]

cv_trans <- compute_constgene_cv(counts_trans_filt, meta, constgene_idx, "Transcriptome")
cv_translat <- compute_constgene_cv(counts_translat_filt, meta, constgene_idx, "Translatome")

cv_threshold <- params$normalization$constgene_cv_threshold

pass_trans <- cv_trans <= cv_threshold
pass_translat <- cv_translat <= cv_threshold
pass_both <- pass_trans & pass_translat

cat(sprintf("constGene CV < %.1f: transcriptome %d/%d, translatome %d/%d, both %d/%d\n",
            cv_threshold,
            sum(pass_trans), length(pass_trans),
            sum(pass_translat), length(pass_translat),
            sum(pass_both), length(pass_both)))

constgene_cv_df <- data.frame(
  gene_id = gene_ids_filt[constgene_idx],
  gene_id_base = gene_ids_base_filt[constgene_idx],
  symbol = symbols_filt[constgene_idx],
  cv_transcriptome = cv_trans,
  cv_translatome = cv_translat,
  pass_transcriptome = pass_trans,
  pass_translatome = pass_translat,
  pass_both = pass_both,
  stringsAsFactors = FALSE
)

write.csv(constgene_cv_df,
          file.path(dir_norm, "constgene_cv_validation.csv"),
          row.names = FALSE)

# Plot CV distribution
cv_plot_df <- data.frame(
  CV = c(cv_trans, cv_translat),
  Omics = rep(c("Transcriptome", "Translatome"), each = length(cv_trans))
)

p_cv <- ggplot(cv_plot_df, aes(x = CV, fill = Omics)) +
  geom_histogram(bins = 30, alpha = 0.6, position = "identity") +
  geom_vline(xintercept = cv_threshold, linetype = "dashed", color = "red") +
  scale_fill_manual(values = col_omics) +
  labs(title = "constGene CV across timepoints",
       x = "Coefficient of Variation", y = "Count") +
  theme_oocyte()

ggsave(file.path(dir_norm, "constgene_cv_distribution.pdf"),
       p_cv, width = 7, height = 4)

# Use genes passing in BOTH omics for primary normalization
constgene_pass_idx <- constgene_idx[pass_both]

if (length(constgene_pass_idx) < 30) {
  # With 80% mRNA degradation, even raw count CV can be high due to variable
  # sequencing depth. Fall back to all mapped constGenes with nonzero expression.
  warning(sprintf(
    "Only %d constGenes pass CV<%.1f in both omics. Using all %d mapped constGenes (externally validated by Di Wu 2022).",
    length(constgene_pass_idx), cv_threshold, length(constgene_idx)))

  # Require nonzero in all samples of at least one omics
  nonzero_trans <- rowSums(counts_trans_filt[constgene_idx, ] == 0) == 0
  nonzero_translat <- rowSums(counts_translat_filt[constgene_idx, ] == 0) == 0
  usable <- nonzero_trans | nonzero_translat
  constgene_pass_idx <- constgene_idx[usable]

  cat(sprintf("Fallback: %d constGenes with nonzero expression in all samples (at least one omics)\n",
              length(constgene_pass_idx)))
}

if (length(constgene_pass_idx) < 20) {
  stop(sprintf("Only %d usable constGenes — insufficient for normalization", length(constgene_pass_idx)))
}

cat(sprintf("Using %d constGenes for normalization\n", length(constgene_pass_idx)))

# ============================================================
# 4. Scheme A — constGenes reference normalization
# ============================================================

compute_constgene_size_factors <- function(counts_mat, ref_rows) {
  ref_counts <- counts_mat[ref_rows, , drop = FALSE]
  nonzero <- rowSums(ref_counts == 0) == 0
  if (sum(nonzero) < 10) {
    stop(sprintf("Only %d constGenes have nonzero counts in all samples", sum(nonzero)))
  }
  ref_counts <- ref_counts[nonzero, , drop = FALSE]
  geo_means <- exp(rowMeans(log(ref_counts)))
  ratios <- sweep(ref_counts, 1, geo_means, "/")
  size_factors <- apply(ratios, 2, median)
  size_factors
}

sf_trans_A <- compute_constgene_size_factors(counts_trans_filt, constgene_pass_idx)
sf_translat_A <- compute_constgene_size_factors(counts_translat_filt, constgene_pass_idx)

norm_trans_A <- t(t(counts_trans_filt) / sf_trans_A)
norm_translat_A <- t(t(counts_translat_filt) / sf_translat_A)

cat("Scheme A (constGenes) size factors:\n")
cat(sprintf("  Transcriptome: range [%.3f, %.3f], CV=%.3f\n",
            min(sf_trans_A), max(sf_trans_A), sd(sf_trans_A) / mean(sf_trans_A)))
cat(sprintf("  Translatome:   range [%.3f, %.3f], CV=%.3f\n",
            min(sf_translat_A), max(sf_translat_A), sd(sf_translat_A) / mean(sf_translat_A)))

# ============================================================
# 5. Scheme B — DESeq2 median-of-ratios
# ============================================================

compute_deseq2_size_factors <- function(counts_mat, sample_meta, omics_label) {
  sample_ids <- sample_meta$sample_id[sample_meta$omics_type == omics_label]
  sample_ids <- sample_ids[sample_ids %in% colnames(counts_mat)]
  sub_meta <- sample_meta[sample_meta$sample_id %in% sample_ids, ]
  sub_meta <- sub_meta[match(sample_ids, sub_meta$sample_id), ]

  dds <- DESeqDataSetFromMatrix(
    countData = counts_mat[, sample_ids],
    colData = sub_meta,
    design = ~ stage
  )
  dds <- estimateSizeFactors(dds)
  sizeFactors(dds)
}

sf_trans_B <- compute_deseq2_size_factors(counts_trans_filt, meta, "Transcriptome")
sf_translat_B <- compute_deseq2_size_factors(counts_translat_filt, meta, "Translatome")

norm_trans_B <- t(t(counts_trans_filt[, names(sf_trans_B)]) / sf_trans_B)
norm_translat_B <- t(t(counts_translat_filt[, names(sf_translat_B)]) / sf_translat_B)

stopifnot(setequal(colnames(norm_trans_B), colnames(norm_trans_A)))
stopifnot(setequal(colnames(norm_translat_B), colnames(norm_translat_A)))

cat("Scheme B (DESeq2) size factors:\n")
cat(sprintf("  Transcriptome: range [%.3f, %.3f], CV=%.3f\n",
            min(sf_trans_B), max(sf_trans_B), sd(sf_trans_B) / mean(sf_trans_B)))
cat(sprintf("  Translatome:   range [%.3f, %.3f], CV=%.3f\n",
            min(sf_translat_B), max(sf_translat_B), sd(sf_translat_B) / mean(sf_translat_B)))

# ============================================================
# 6. Scheme C — TMM (edgeR)
# ============================================================

compute_tmm_size_factors <- function(counts_mat) {
  dge <- DGEList(counts = counts_mat)
  dge <- calcNormFactors(dge, method = "TMM")
  lib_size <- colSums(counts_mat)
  sf <- dge$samples$norm.factors * lib_size / exp(mean(log(dge$samples$norm.factors * lib_size)))
  sf
}

sf_trans_C <- compute_tmm_size_factors(counts_trans_filt)
sf_translat_C <- compute_tmm_size_factors(counts_translat_filt)

norm_trans_C <- t(t(counts_trans_filt) / sf_trans_C)
norm_translat_C <- t(t(counts_translat_filt) / sf_translat_C)

cat("Scheme C (TMM) size factors:\n")
cat(sprintf("  Transcriptome: range [%.3f, %.3f], CV=%.3f\n",
            min(sf_trans_C), max(sf_trans_C), sd(sf_trans_C) / mean(sf_trans_C)))
cat(sprintf("  Translatome:   range [%.3f, %.3f], CV=%.3f\n",
            min(sf_translat_C), max(sf_translat_C), sd(sf_translat_C) / mean(sf_translat_C)))

# ============================================================
# 7. Scheme D — RUVg
# ============================================================

# Manual RUVg implementation (Risso et al. 2014, Nature Biotechnology)
# SVD on log-CPM of negative control genes → top-k singular vectors = W
# Adjusted = original - W %*% alpha (gene-wise regression)
run_ruvg_scheme <- function(counts_mat, sample_meta, omics_label,
                            constgene_rows, k_range) {
  sample_ids <- sample_meta$sample_id[sample_meta$omics_type == omics_label]
  sample_ids <- sample_ids[sample_ids %in% colnames(counts_mat)]
  sub_counts <- counts_mat[, sample_ids]

  lib_sizes <- colSums(sub_counts)
  log_cpm <- log2(t(t(sub_counts) / lib_sizes) * 1e6 + 1)

  ctrl_log_cpm <- log_cpm[constgene_rows, , drop = FALSE]
  ctrl_centered <- t(scale(t(ctrl_log_cpm), center = TRUE, scale = FALSE))
  svd_ctrl <- svd(t(ctrl_centered))

  max_k <- max(k_range)
  stopifnot(max_k <= ncol(svd_ctrl$u))

  results <- list()
  for (k in k_range) {
    W <- svd_ctrl$u[, seq_len(k), drop = FALSE]
    colnames(W) <- paste0("W_", seq_len(k))
    rownames(W) <- sample_ids

    adjusted <- log_cpm
    for (g in seq_len(nrow(log_cpm))) {
      fit <- lm.fit(W, log_cpm[g, ])
      adjusted[g, ] <- log_cpm[g, ] - W %*% fit$coefficients
    }

    norm_counts <- pmax(round(2^adjusted * rep(lib_sizes / 1e6, each = nrow(adjusted)) - 1), 0)

    results[[paste0("k", k)]] <- list(
      normalized = norm_counts,
      W = as.data.frame(W)
    )
  }
  results
}

k_range <- params$normalization$ruvg_k_range

ruvg_trans <- run_ruvg_scheme(counts_trans_filt, meta, "Transcriptome",
                               constgene_pass_idx, k_range)
ruvg_translat <- run_ruvg_scheme(counts_translat_filt, meta, "Translatome",
                                  constgene_pass_idx, k_range)

cat(sprintf("Scheme D (RUVg) computed for k = %s\n",
            paste(k_range, collapse = ", ")))

# Select optimal k by examining RLE spread reduction
compute_rle_iqr <- function(log_counts) {
  med_gene <- rowMedians(log_counts)
  rle <- sweep(log_counts, 1, med_gene)
  median(apply(rle, 2, IQR))
}

rle_results <- data.frame(
  scheme = character(), omics = character(),
  k = integer(), rle_iqr = numeric(),
  stringsAsFactors = FALSE
)

for (k in k_range) {
  key <- paste0("k", k)
  rle_trans <- compute_rle_iqr(log2(pmax(ruvg_trans[[key]]$normalized, 0) + 1))
  rle_translat <- compute_rle_iqr(log2(pmax(ruvg_translat[[key]]$normalized, 0) + 1))
  rle_results <- rbind(rle_results,
    data.frame(scheme = "RUVg", omics = "Transcriptome", k = k, rle_iqr = rle_trans),
    data.frame(scheme = "RUVg", omics = "Translatome", k = k, rle_iqr = rle_translat)
  )
}

# Add baseline RLE for other schemes
rle_raw_trans <- compute_rle_iqr(log2(counts_trans_filt + 1))
rle_raw_translat <- compute_rle_iqr(log2(counts_translat_filt + 1))
rle_A_trans <- compute_rle_iqr(log2(norm_trans_A + 1))
rle_A_translat <- compute_rle_iqr(log2(norm_translat_A + 1))
rle_B_trans <- compute_rle_iqr(log2(norm_trans_B + 1))
rle_B_translat <- compute_rle_iqr(log2(norm_translat_B + 1))
rle_C_trans <- compute_rle_iqr(log2(norm_trans_C + 1))
rle_C_translat <- compute_rle_iqr(log2(norm_translat_C + 1))

rle_comparison <- rbind(
  data.frame(scheme = "Raw", omics = "Transcriptome", k = NA, rle_iqr = rle_raw_trans),
  data.frame(scheme = "Raw", omics = "Translatome", k = NA, rle_iqr = rle_raw_translat),
  data.frame(scheme = "constGenes", omics = "Transcriptome", k = NA, rle_iqr = rle_A_trans),
  data.frame(scheme = "constGenes", omics = "Translatome", k = NA, rle_iqr = rle_A_translat),
  data.frame(scheme = "DESeq2", omics = "Transcriptome", k = NA, rle_iqr = rle_B_trans),
  data.frame(scheme = "DESeq2", omics = "Translatome", k = NA, rle_iqr = rle_B_translat),
  data.frame(scheme = "TMM", omics = "Transcriptome", k = NA, rle_iqr = rle_C_trans),
  data.frame(scheme = "TMM", omics = "Translatome", k = NA, rle_iqr = rle_C_translat),
  rle_results
)

rle_comparison$label <- ifelse(is.na(rle_comparison$k),
                                rle_comparison$scheme,
                                paste0("RUVg_k", rle_comparison$k))

write.csv(rle_comparison, file.path(dir_norm, "rle_comparison.csv"), row.names = FALSE)

cat("RLE IQR comparison:\n")
for (i in seq_len(nrow(rle_comparison))) {
  cat(sprintf("  %s [%s]: %.4f\n",
              rle_comparison$label[i], rle_comparison$omics[i], rle_comparison$rle_iqr[i]))
}

# Select best k: smallest RLE IQR averaged across both omics
avg_rle_by_k <- tapply(rle_results$rle_iqr, rle_results$k, mean)
best_k <- as.integer(names(which.min(avg_rle_by_k)))
cat(sprintf("Best RUVg k = %d (avg RLE IQR = %.4f)\n", best_k, min(avg_rle_by_k)))

norm_trans_D <- ruvg_trans[[paste0("k", best_k)]]$normalized
norm_translat_D <- ruvg_translat[[paste0("k", best_k)]]$normalized
W_trans <- ruvg_trans[[paste0("k", best_k)]]$W
W_translat <- ruvg_translat[[paste0("k", best_k)]]$W

n_neg_trans <- sum(norm_trans_D < 0)
n_neg_translat <- sum(norm_translat_D < 0)
if (n_neg_trans > 0 || n_neg_translat > 0) {
  warning(sprintf("RUVg produced %d negative values in transcriptome, %d in translatome — flooring at 0",
                  n_neg_trans, n_neg_translat))
  norm_trans_D <- pmax(norm_trans_D, 0)
  norm_translat_D <- pmax(norm_translat_D, 0)
}

# ============================================================
# 8. Leave-10%-out cross-validation of constGenes
# ============================================================

n_constgenes <- length(constgene_pass_idx)
leave_out_n <- max(1, round(n_constgenes * params$normalization$leave_out_fraction))
n_iter <- params$normalization$leave_out_iterations

sf_cv_trans <- matrix(NA, nrow = n_iter, ncol = ncol(counts_trans_filt))
sf_cv_translat <- matrix(NA, nrow = n_iter, ncol = ncol(counts_translat_filt))
colnames(sf_cv_trans) <- colnames(counts_trans_filt)
colnames(sf_cv_translat) <- colnames(counts_translat_filt)

set.seed(params$seed)

for (i in seq_len(n_iter)) {
  drop_idx <- sample(seq_len(n_constgenes), leave_out_n)
  subset_idx <- constgene_pass_idx[-drop_idx]
  sf_cv_trans[i, ] <- compute_constgene_size_factors(counts_trans_filt, subset_idx)
  sf_cv_translat[i, ] <- compute_constgene_size_factors(counts_translat_filt, subset_idx)
}

sf_cv_trans_cv <- apply(sf_cv_trans, 2, function(x) sd(x) / mean(x))
sf_cv_translat_cv <- apply(sf_cv_translat, 2, function(x) sd(x) / mean(x))

cat(sprintf("Leave-%.0f%%-out CV of size factors (%d iterations):\n",
            params$normalization$leave_out_fraction * 100, n_iter))
cat(sprintf("  Transcriptome: median CV=%.4f, max=%.4f\n",
            median(sf_cv_trans_cv), max(sf_cv_trans_cv)))
cat(sprintf("  Translatome:   median CV=%.4f, max=%.4f\n",
            median(sf_cv_translat_cv), max(sf_cv_translat_cv)))

cv_stability_df <- data.frame(
  sample_id = c(colnames(counts_trans_filt), colnames(counts_translat_filt)),
  omics = c(rep("Transcriptome", ncol(counts_trans_filt)),
            rep("Translatome", ncol(counts_translat_filt))),
  sf_cv = c(sf_cv_trans_cv, sf_cv_translat_cv)
)

p_stability <- ggplot(cv_stability_df, aes(x = omics, y = sf_cv, fill = omics)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 21) +
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.5) +
  scale_fill_manual(values = col_omics) +
  geom_hline(yintercept = 0.1, linetype = "dashed", color = "red") +
  labs(title = "Size factor stability (leave-10%-out CV)",
       x = NULL, y = "CV of size factors across 100 iterations") +
  theme_oocyte() +
  theme(legend.position = "none")

ggsave(file.path(dir_norm, "constgene_leave_out_stability.pdf"),
       p_stability, width = 5, height = 4)

# ============================================================
# 9. Size factor comparison across schemes
# ============================================================

sf_comparison <- data.frame(
  sample_id = colnames(counts_trans_filt),
  sf_constGenes = sf_trans_A,
  sf_DESeq2 = sf_trans_B[colnames(counts_trans_filt)],
  sf_TMM = sf_trans_C
)

sf_cor_AB <- cor(sf_comparison$sf_constGenes, sf_comparison$sf_DESeq2, method = "spearman")
sf_cor_AC <- cor(sf_comparison$sf_constGenes, sf_comparison$sf_TMM, method = "spearman")
sf_cor_BC <- cor(sf_comparison$sf_DESeq2, sf_comparison$sf_TMM, method = "spearman")

cat(sprintf("Size factor rank correlations (transcriptome):\n"))
cat(sprintf("  constGenes vs DESeq2: rho=%.3f\n", sf_cor_AB))
cat(sprintf("  constGenes vs TMM:    rho=%.3f\n", sf_cor_AC))
cat(sprintf("  DESeq2 vs TMM:        rho=%.3f\n", sf_cor_BC))

p_sf <- ggplot(sf_comparison, aes(x = sf_constGenes, y = sf_DESeq2)) +
  geom_point(size = 2, color = col_omics["Transcriptome"]) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray40") +
  labs(title = sprintf("Size factors: constGenes vs DESeq2 (rho=%.3f)", sf_cor_AB),
       x = "constGenes size factor", y = "DESeq2 size factor") +
  theme_oocyte()

ggsave(file.path(dir_norm, "size_factor_comparison.pdf"), p_sf, width = 5, height = 5)

# ============================================================
# 10. TE calculation (paired design)
# ============================================================

compute_te_matrix <- function(norm_trans, norm_translat, pair_map, pseudo = 1) {
  genes <- rownames(norm_trans)
  n_pairs <- nrow(pair_map)

  te_mat <- matrix(NA, nrow = length(genes), ncol = n_pairs,
                   dimnames = list(genes, paste0(pair_map$time_point, "_rep", pair_map$replicate)))

  for (k in seq_len(n_pairs)) {
    tid <- pair_map$transcriptome_id[k]
    lid <- pair_map$translatome_id[k]

    if (!(tid %in% colnames(norm_trans)) || !(lid %in% colnames(norm_translat))) {
      warning(sprintf("Pair %d: missing sample %s or %s", k, tid, lid))
      next
    }

    te_mat[, k] <- log2((norm_translat[, lid] + pseudo) /
                         (norm_trans[, tid] + pseudo))
  }

  te_mat
}

te_A <- compute_te_matrix(norm_trans_A, norm_translat_A, pair_map)
te_B <- compute_te_matrix(norm_trans_B, norm_translat_B, pair_map)
te_C <- compute_te_matrix(norm_trans_C, norm_translat_C, pair_map)
# RUVg TE is exploratory only: W estimated independently per omics,
# so the ratio does not strictly represent translational efficiency.
te_D <- compute_te_matrix(norm_trans_D, norm_translat_D, pair_map)

cat(sprintf("TE matrices computed: %d genes x %d pairs\n", nrow(te_A), ncol(te_A)))

# ============================================================
# 11. Outlier reassessment post-normalization
# ============================================================

assess_outliers_pca <- function(norm_mat, sample_meta, omics_label,
                                top_n = 2000, threshold_mad = 3, min_group_size = 5) {
  sample_ids <- colnames(norm_mat)
  sub_meta <- sample_meta[match(sample_ids, sample_meta$sample_id), ]

  log_mat <- log2(norm_mat + 1)
  vars <- rowVars(log_mat)
  top_idx <- order(vars, decreasing = TRUE)[seq_len(min(top_n, nrow(log_mat)))]
  pca <- prcomp(t(log_mat[top_idx, ]), scale. = TRUE)

  pc1 <- pca$x[, 1]
  pc2 <- pca$x[, 2]

  outliers <- character(0)
  for (s in levels(sub_meta$stage)) {
    idx <- which(sub_meta$stage == s)
    if (length(idx) < min_group_size) {
      message(sprintf("  %s %s: n=%d, skipping outlier check", omics_label, s, length(idx)))
      next
    }

    for (pc_vals in list(pc1[idx], pc2[idx])) {
      med <- median(pc_vals)
      mad_val <- mad(pc_vals, constant = 1.4826)
      if (mad_val > 0) {
        z <- abs(pc_vals - med) / mad_val
        hits <- names(pc_vals)[z > threshold_mad]
        outliers <- c(outliers, hits)
      }
    }
  }

  list(
    pca = pca,
    pc_scores = data.frame(sample_id = sample_ids,
                           PC1 = pc1, PC2 = pc2,
                           stage = sub_meta$stage,
                           stringsAsFactors = FALSE),
    outliers = unique(outliers),
    var_explained = summary(pca)$importance[2, 1:2]
  )
}

cat("\n=== Outlier reassessment after constGenes normalization ===\n")

outlier_trans_A <- assess_outliers_pca(norm_trans_A, meta, "Transcriptome")
outlier_translat_A <- assess_outliers_pca(norm_translat_A, meta, "Translatome")

cat(sprintf("Post-normalization PCA variance:\n"))
cat(sprintf("  Transcriptome: PC1=%.1f%%, PC2=%.1f%%\n",
            outlier_trans_A$var_explained[1] * 100,
            outlier_trans_A$var_explained[2] * 100))
cat(sprintf("  Translatome:   PC1=%.1f%%, PC2=%.1f%%\n",
            outlier_translat_A$var_explained[1] * 100,
            outlier_translat_A$var_explained[2] * 100))

all_post_outliers <- unique(c(outlier_trans_A$outliers, outlier_translat_A$outliers))

if (length(all_post_outliers) > 0) {
  cat(sprintf("Post-normalization outliers: %s\n",
              paste(all_post_outliers, collapse = ", ")))

  # Cross-reference with QC outliers (loaded from qc_filtered_data.RData)
  qc_outliers <- if (exists("all_outliers")) all_outliers else character(0)
  persistent <- intersect(all_post_outliers, qc_outliers)
  new_outliers <- setdiff(all_post_outliers, qc_outliers)
  resolved <- setdiff(qc_outliers, all_post_outliers)

  cat(sprintf("  Persistent (pre+post): %d — %s\n",
              length(persistent), paste(persistent, collapse = ", ")))
  cat(sprintf("  Resolved by normalization: %d — %s\n",
              length(resolved), paste(resolved, collapse = ", ")))
  if (length(new_outliers) > 0) {
    cat(sprintf("  New post-normalization: %d — %s\n",
                length(new_outliers), paste(new_outliers, collapse = ", ")))
  }
} else {
  cat("No outliers detected post-normalization.\n")
}

# Check flagged pairs specifically
flagged_pairs <- c("F36.bam", "F36R.bam", "F26.bam", "F26R.bam")
cat("\nFlagged pair status post-normalization:\n")
for (sid in flagged_pairs) {
  in_trans <- sid %in% colnames(norm_trans_A)
  in_translat <- sid %in% colnames(norm_translat_A)

  if (in_trans) {
    pc_row <- outlier_trans_A$pc_scores[outlier_trans_A$pc_scores$sample_id == sid, ]
    is_outlier <- sid %in% outlier_trans_A$outliers
    cat(sprintf("  %s (transcriptome): PC1=%.1f, PC2=%.1f, outlier=%s\n",
                sid, pc_row$PC1, pc_row$PC2, is_outlier))
  }
  if (in_translat) {
    pc_row <- outlier_translat_A$pc_scores[outlier_translat_A$pc_scores$sample_id == sid, ]
    is_outlier <- sid %in% outlier_translat_A$outliers
    cat(sprintf("  %s (translatome):   PC1=%.1f, PC2=%.1f, outlier=%s\n",
                sid, pc_row$PC1, pc_row$PC2, is_outlier))
  }
}

# PCA plots post-normalization
plot_pca_norm <- function(pc_scores, var_explained, title_suffix) {
  p <- ggplot(pc_scores, aes(x = PC1, y = PC2, color = stage)) +
    geom_point(size = 3) +
    scale_color_manual(values = col_timepoint) +
    labs(title = paste("PCA —", title_suffix),
         x = sprintf("PC1 (%.1f%%)", var_explained[1] * 100),
         y = sprintf("PC2 (%.1f%%)", var_explained[2] * 100)) +
    theme_oocyte()
  p
}

p_pca_trans <- plot_pca_norm(outlier_trans_A$pc_scores,
                              outlier_trans_A$var_explained,
                              "Transcriptome (constGenes norm)")
p_pca_translat <- plot_pca_norm(outlier_translat_A$pc_scores,
                                 outlier_translat_A$var_explained,
                                 "Translatome (constGenes norm)")

ggsave(file.path(dir_norm, "pca_transcriptome_constgenes.pdf"),
       p_pca_trans, width = 7, height = 5)
ggsave(file.path(dir_norm, "pca_translatome_constgenes.pdf"),
       p_pca_translat, width = 7, height = 5)

# ============================================================
# 12. Gene ranking sensitivity across schemes
# ============================================================

compute_te_fold_change <- function(te_mat) {
  gv_cols <- grep("^0h_", colnames(te_mat))
  mii_cols <- grep("^12h_", colnames(te_mat))

  if (length(gv_cols) == 0 || length(mii_cols) == 0) {
    stop("Cannot find GV (0h) or MII (12h) columns in TE matrix")
  }

  te_gv <- rowMeans(te_mat[, gv_cols, drop = FALSE], na.rm = TRUE)
  te_mii <- rowMeans(te_mat[, mii_cols, drop = FALSE], na.rm = TRUE)

  te_mii - te_gv
}

fc_A <- compute_te_fold_change(te_A)
fc_B <- compute_te_fold_change(te_B)
fc_C <- compute_te_fold_change(te_C)
fc_D <- compute_te_fold_change(te_D)

genes_common <- Reduce(intersect, list(names(fc_A), names(fc_B), names(fc_C), names(fc_D)))
fc_A_common <- fc_A[genes_common]
fc_B_common <- fc_B[genes_common]
fc_C_common <- fc_C[genes_common]
fc_D_common <- fc_D[genes_common]

rank_cor_AB <- cor(fc_A_common, fc_B_common, method = "spearman", use = "complete.obs")
rank_cor_AC <- cor(fc_A_common, fc_C_common, method = "spearman", use = "complete.obs")
rank_cor_AD <- cor(fc_A_common, fc_D_common, method = "spearman", use = "complete.obs")
rank_cor_BC <- cor(fc_B_common, fc_C_common, method = "spearman", use = "complete.obs")

rank_divergence_AB <- 1 - rank_cor_AB

cat(sprintf("\nTE fold-change (GV→MII) rank correlations:\n"))
cat(sprintf("  constGenes vs DESeq2: rho=%.3f (divergence=%.1f%%)\n",
            rank_cor_AB, rank_divergence_AB * 100))
cat(sprintf("  constGenes vs TMM:    rho=%.3f\n", rank_cor_AC))
cat(sprintf("  constGenes vs RUVg:   rho=%.3f\n", rank_cor_AD))
cat(sprintf("  DESeq2 vs TMM:        rho=%.3f\n", rank_cor_BC))

# ============================================================
# 13. Known gene TE trends
# ============================================================

plot_control_gene_te <- function(te_mat, pair_map, gene_symbols_vec, gene_ids_vec,
                                  control_genes, scheme_label) {
  stages <- c("GV", "GVBD", "MI-6", "MI-9", "MII")
  time_numeric <- c(0, 3, 6, 9, 12)

  plot_data <- data.frame()

  for (g in control_genes) {
    gene_rows <- which(gene_symbols_vec == g)
    if (length(gene_rows) == 0) next

    gene_row_name <- gene_ids_vec[gene_rows[1]]
    if (!(gene_row_name %in% rownames(te_mat))) next

    for (tp_idx in seq_along(stages)) {
      tp_label <- c("0h", "3h", "6h", "9h", "12h")[tp_idx]
      cols <- grep(paste0("^", tp_label, "_"), colnames(te_mat))
      if (length(cols) == 0) next

      te_vals <- te_mat[gene_row_name, cols]
      plot_data <- rbind(plot_data, data.frame(
        gene = g,
        time = time_numeric[tp_idx],
        stage = stages[tp_idx],
        te = te_vals,
        stringsAsFactors = FALSE
      ))
    }
  }

  n_found <- length(unique(plot_data$gene))
  cat(sprintf("Control gene TE plot (%s): %d/%d genes found\n",
              scheme_label, n_found, length(control_genes)))

  if (nrow(plot_data) == 0) return(NULL)

  p <- ggplot(plot_data, aes(x = time, y = te, color = gene, group = gene)) +
    stat_summary(fun = mean, geom = "line", linewidth = 0.8) +
    stat_summary(fun = mean, geom = "point", size = 2) +
    stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.5, alpha = 0.5) +
    scale_x_continuous(breaks = time_numeric, labels = stages) +
    labs(title = paste("Control gene TE trends —", scheme_label),
         x = "Stage", y = "log2(TE)") +
    theme_oocyte()
  p
}

p_ctrl_A <- plot_control_gene_te(te_A, pair_map, symbols_filt, gene_ids_filt,
                                  params$control_genes, "constGenes")
p_ctrl_B <- plot_control_gene_te(te_B, pair_map, symbols_filt, gene_ids_filt,
                                  params$control_genes, "DESeq2")

if (!is.null(p_ctrl_A)) {
  ggsave(file.path(dir_norm, "control_gene_te_constgenes.pdf"),
         p_ctrl_A, width = 8, height = 5)
}
if (!is.null(p_ctrl_B)) {
  ggsave(file.path(dir_norm, "control_gene_te_deseq2.pdf"),
         p_ctrl_B, width = 8, height = 5)
}

# ============================================================
# 14. Save all normalized data
# ============================================================

save(
  norm_trans_A, norm_translat_A, sf_trans_A, sf_translat_A,
  norm_trans_B, norm_translat_B, sf_trans_B, sf_translat_B,
  norm_trans_C, norm_translat_C, sf_trans_C, sf_translat_C,
  norm_trans_D, norm_translat_D, W_trans, W_translat, best_k,
  te_A, te_B, te_C, te_D,
  constgene_pass_idx, constgene_cv_df,
  pair_map, meta, gene_ids_filt, symbols_filt,
  file = file.path(dir_norm, "normalized_data.RData")
)

# ============================================================
# 15. Summary
# ============================================================

cat(sprintf("\n=== Normalization Summary ===\n"))
cat(sprintf("constGenes passing CV<%.1f: %d/%d\n",
            cv_threshold, sum(pass_both), length(pass_both)))
cat(sprintf("Leave-10%%-out size factor CV: median=%.4f (threshold <0.1)\n",
            median(c(sf_cv_trans_cv, sf_cv_translat_cv))))
cat(sprintf("Gene ranking divergence (constGenes vs DESeq2): %.1f%% (threshold >20%%)\n",
            rank_divergence_AB * 100))
cat(sprintf("Best RUVg k: %d\n", best_k))
cat(sprintf("Post-normalization outliers: %d\n", length(all_post_outliers)))
cat(sprintf("TE matrices: %d genes x %d pairs x 4 schemes\n",
            nrow(te_A), ncol(te_A)))
cat(sprintf("Results saved to: %s\n", dir_norm))

# ============================================================
# 16. Session info
# ============================================================

writeLines(capture.output(sessionInfo()),
           file.path(dir_norm, "sessionInfo.txt"))
