# 01_qc.R — Data QC & Preprocessing (IU-2)
#
# Input:  rawcount/*.csv, metadata_hcg.csv
# Output: results/qc/ (QC plots + filtered count matrices)

source("src/R/00_config.R")

library(ggplot2)
library(pheatmap)
library(matrixStats)
library(reshape2)

dir_qc <- file.path(DIR_RESULTS, "qc")
dir.create(dir_qc, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# 1. Load data
# ============================================================

cat("Loading metadata...\n")
meta <- load_metadata()

cat("Loading transcriptome counts...\n")
trans <- load_counts(params$paths$rawcount_transcriptome)

cat("Loading translatome counts...\n")
translat <- load_counts(params$paths$rawcount_translatome)

cat(sprintf("Transcriptome: %d genes x %d samples\n",
            nrow(trans$counts), ncol(trans$counts)))
cat(sprintf("Translatome:   %d genes x %d samples\n",
            nrow(translat$counts), ncol(translat$counts)))

stopifnot(identical(trans$gene_ids_full, translat$gene_ids_full))

# Validate sample IDs match between counts and metadata (IMPORTANT-4)
trans_samples <- colnames(trans$counts)
translat_samples <- colnames(translat$counts)
meta_trans_ids <- meta$sample_id[meta$omics_type == "Transcriptome"]
meta_translat_ids <- meta$sample_id[meta$omics_type == "Translatome"]

if (!setequal(trans_samples, meta_trans_ids)) {
  stop(sprintf("Transcriptome sample ID mismatch.\n  In counts not meta: %s\n  In meta not counts: %s",
               paste(setdiff(trans_samples, meta_trans_ids), collapse = ", "),
               paste(setdiff(meta_trans_ids, trans_samples), collapse = ", ")))
}
if (!setequal(translat_samples, meta_translat_ids)) {
  stop(sprintf("Translatome sample ID mismatch.\n  In counts not meta: %s\n  In meta not counts: %s",
               paste(setdiff(translat_samples, meta_translat_ids), collapse = ", "),
               paste(setdiff(meta_translat_ids, translat_samples), collapse = ", ")))
}

# ============================================================
# 2. Gene annotation (mt, ribo)
# ============================================================

gene_info <- trans$gene_info
gene_ids_base <- trans$gene_ids_base

# Map Ensembl IDs to gene symbols for annotation
# Use org.Mm.eg.db if available, otherwise flag mt/ribo by known Ensembl IDs
annotate_genes <- function(gene_ids_base) {
  if (!requireNamespace("org.Mm.eg.db", quietly = TRUE)) {
    stop("org.Mm.eg.db is required for gene annotation. Install with:\n",
         "  BiocManager::install('org.Mm.eg.db')")
  }
  library(org.Mm.eg.db)
  symbols <- AnnotationDbi::mapIds(
    org.Mm.eg.db,
    keys = gene_ids_base,
    column = "SYMBOL",
    keytype = "ENSEMBL",
    multiVals = "first"
  )
  return(symbols)
}

gene_symbols <- annotate_genes(gene_ids_base)

is_mt <- grepl("^chrM", gene_info$Chr)
if (!is.null(gene_symbols)) {
  is_mt <- is_mt | grepl(params$qc$mt_pattern, gene_symbols, ignore.case = TRUE)
}

is_ribo <- rep(FALSE, length(gene_ids_base))
if (!is.null(gene_symbols)) {
  is_ribo <- grepl(params$qc$ribo_pattern, gene_symbols)
}

cat(sprintf("Mitochondrial genes: %d\n", sum(is_mt)))
cat(sprintf("Ribosomal protein genes: %d\n", sum(is_ribo)))

# ============================================================
# 3. Sample-level QC
# ============================================================

compute_sample_qc <- function(counts, is_mt, is_ribo, label) {
  lib_size <- colSums(counts)
  n_detected <- colSums(counts > 0)
  pct_mt <- colSums(counts[is_mt, , drop = FALSE]) / lib_size * 100
  pct_ribo <- colSums(counts[is_ribo, , drop = FALSE]) / lib_size * 100
  pct_top50 <- apply(counts, 2, function(x) {
    total <- sum(x)
    if (total == 0) return(NA_real_)
    sum(sort(x, decreasing = TRUE)[1:min(50, length(x))]) / total * 100
  })

  data.frame(
    sample_id = colnames(counts),
    omics = label,
    lib_size = lib_size,
    n_detected = n_detected,
    pct_mt = pct_mt,
    pct_ribo = pct_ribo,
    pct_top50 = pct_top50,
    stringsAsFactors = FALSE
  )
}

qc_trans <- compute_sample_qc(trans$counts, is_mt, is_ribo, "Transcriptome")
qc_translat <- compute_sample_qc(translat$counts, is_mt, is_ribo, "Translatome")
qc_all <- rbind(qc_trans, qc_translat)

# Add metadata
qc_all <- merge(qc_all, meta[, c("sample_id", "time_point", "stage", "replicate")],
                by = "sample_id", all.x = TRUE)

n_unmatched <- sum(is.na(qc_all$stage))
if (n_unmatched > 0) {
  stop(sprintf("FAIL: %d samples not found in metadata: %s",
               n_unmatched, paste(qc_all$sample_id[is.na(qc_all$stage)], collapse = ", ")))
}

write.csv(qc_all, file.path(dir_qc, "sample_qc_metrics.csv"), row.names = FALSE)
cat("Sample QC metrics saved.\n")

# ============================================================
# 4. QC Plots
# ============================================================

# 4a. Library size barplot
p_libsize <- ggplot(qc_all, aes(x = reorder(sample_id, -lib_size),
                                 y = lib_size / 1e6,
                                 fill = stage)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = col_timepoint) +
  facet_wrap(~omics, scales = "free_x") +
  theme_oocyte() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6)) +
  labs(x = NULL, y = "Library size (millions)", title = "Library size per sample")

ggsave(file.path(dir_qc, "library_size.pdf"), p_libsize,
       width = 12, height = 5)

# 4b. Genes detected
p_detected <- ggplot(qc_all, aes(x = stage, y = n_detected, color = stage)) +
  geom_jitter(width = 0.15, size = 2.5) +
  stat_summary(fun = median, geom = "crossbar", width = 0.5, color = "black") +
  scale_color_manual(values = col_timepoint) +
  facet_wrap(~omics) +
  theme_oocyte() +
  labs(x = NULL, y = "Genes detected (count > 0)",
       title = "Genes detected per sample")

ggsave(file.path(dir_qc, "genes_detected.pdf"), p_detected,
       width = 8, height = 4)

# 4c. Mitochondrial & ribosomal percentage
qc_melt <- melt(qc_all,
                 id.vars = c("sample_id", "omics", "stage"),
                 measure.vars = c("pct_mt", "pct_ribo", "pct_top50"))

p_pct <- ggplot(qc_melt, aes(x = stage, y = value, color = stage)) +
  geom_jitter(width = 0.15, size = 2) +
  stat_summary(fun = median, geom = "crossbar", width = 0.5, color = "black") +
  scale_color_manual(values = col_timepoint) +
  facet_grid(variable ~ omics, scales = "free_y") +
  theme_oocyte() +
  labs(x = NULL, y = "Percentage", title = "QC metrics by stage")

ggsave(file.path(dir_qc, "qc_percentages.pdf"), p_pct,
       width = 8, height = 8)

# ============================================================
# 5. Gene filtering
# ============================================================

min_counts <- params$qc$min_counts
min_samples <- params$qc$min_samples

pass_trans <- rowSums(trans$counts >= min_counts) >= min_samples
pass_translat <- rowSums(translat$counts >= min_counts) >= min_samples
pass_both <- pass_trans & pass_translat
pass_either <- pass_trans | pass_translat
pass_filter <- pass_either & !is_mt

n_before <- length(pass_filter)
n_after <- sum(pass_filter)
n_both <- sum(pass_both & !is_mt)
cat(sprintf("Gene filtering: %d → %d (removed %d; %d mt excluded)\n",
            n_before, n_after, n_before - n_after, sum(is_mt)))
cat(sprintf("  Genes passing in BOTH omics: %d\n", n_both))
cat(sprintf("  Genes passing in EITHER omics: %d (union filter used)\n", n_after))
cat("  Note: downstream TE scripts must handle genes absent in one omics.\n")

counts_trans_filt <- trans$counts[pass_filter, ]
counts_translat_filt <- translat$counts[pass_filter, ]
gene_ids_filt <- trans$gene_ids_full[pass_filter]
gene_ids_base_filt <- trans$gene_ids_base[pass_filter]
gene_info_filt <- gene_info[pass_filter, ]

symbols_filt <- gene_symbols[pass_filter]

write.csv(data.frame(gene_id = gene_ids_filt,
                     gene_id_base = gene_ids_base_filt),
          file.path(dir_qc, "filtered_gene_list.csv"),
          row.names = FALSE)

# ============================================================
# 6. Sample correlation heatmap
# ============================================================

plot_correlation <- function(counts, meta_sub, title, filename) {
  log_counts <- log2(counts + 1)
  cor_mat <- cor(log_counts, method = "spearman")

  ann_col <- data.frame(
    Stage = meta_sub$stage,
    row.names = colnames(counts)
  )

  ann_colors <- list(Stage = col_timepoint)

  pdf(file.path(dir_qc, filename), width = 8, height = 7)
  on.exit(dev.off(), add = TRUE)
  pheatmap(cor_mat,
           annotation_col = ann_col,
           annotation_colors = ann_colors,
           main = title,
           fontsize = 7,
           display_numbers = FALSE)
}

meta_trans <- meta[meta$omics_type == "Transcriptome", ]
meta_trans <- meta_trans[match(colnames(counts_trans_filt), meta_trans$sample_id), ]

meta_translat <- meta[meta$omics_type == "Translatome", ]
meta_translat <- meta_translat[match(colnames(counts_translat_filt), meta_translat$sample_id), ]

plot_correlation(counts_trans_filt, meta_trans,
                 "Spearman correlation — Transcriptome",
                 "correlation_transcriptome.pdf")

plot_correlation(counts_translat_filt, meta_translat,
                 "Spearman correlation — Translatome",
                 "correlation_translatome.pdf")

# ============================================================
# 7. PCA
# ============================================================

run_pca <- function(counts, meta_sub, omics_label) {
  log_counts <- log2(counts + 1)
  vars <- rowVars(log_counts)
  n_top <- params$qc$pca_top_genes
  top_idx <- order(vars, decreasing = TRUE)[1:min(n_top, nrow(log_counts))]
  pca_res <- prcomp(t(log_counts[top_idx, ]), center = TRUE, scale. = TRUE)

  pca_df <- data.frame(
    PC1 = pca_res$x[, 1],
    PC2 = pca_res$x[, 2],
    sample_id = colnames(counts)
  )
  pca_df <- merge(pca_df, meta_sub[, c("sample_id", "stage", "time_point")],
                  by = "sample_id")

  var_explained <- summary(pca_res)$importance[2, 1:2] * 100

  p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = stage)) +
    geom_point(size = 3) +
    scale_color_manual(values = col_timepoint) +
    theme_oocyte() +
    labs(x = sprintf("PC1 (%.1f%%)", var_explained[1]),
         y = sprintf("PC2 (%.1f%%)", var_explained[2]),
         title = paste("PCA —", omics_label))

  list(plot = p, pca = pca_res, var_explained = var_explained)
}

pca_trans <- run_pca(counts_trans_filt, meta_trans, "Transcriptome")
pca_translat <- run_pca(counts_translat_filt, meta_translat, "Translatome")

ggsave(file.path(dir_qc, "pca_transcriptome.pdf"),
       pca_trans$plot, width = 6, height = 5)
ggsave(file.path(dir_qc, "pca_translatome.pdf"),
       pca_translat$plot, width = 6, height = 5)

cat(sprintf("PCA variance explained — Transcriptome: PC1=%.1f%%, PC2=%.1f%%\n",
            pca_trans$var_explained[1], pca_trans$var_explained[2]))
cat(sprintf("PCA variance explained — Translatome:   PC1=%.1f%%, PC2=%.1f%%\n",
            pca_translat$var_explained[1], pca_translat$var_explained[2]))

# ============================================================
# 8. Pair validation
# ============================================================

pair_map <- build_pair_map(meta)
cat(sprintf("Complete pairs: %d\n", nrow(pair_map)))

# Document MI-9 reduced sample size per AGENTS.md requirement
mi9_translat_n <- sum(meta$omics_type == "Translatome" & meta$stage == "MI-9")
if (mi9_translat_n < 5) {
  cat(sprintf("NOTE: MI-9 translatome has n=%d (expected 5). Reduced power documented.\n",
              mi9_translat_n))
}
cat("F67.bam (MI-9 transcriptome rep5): included in transcriptome, excluded from TE pairs.\n")

pair_cors <- numeric(nrow(pair_map))
for (i in seq_len(nrow(pair_map))) {
  tid <- pair_map$transcriptome_id[i]
  rid <- pair_map$translatome_id[i]
  if (tid %in% colnames(counts_trans_filt) &&
      rid %in% colnames(counts_translat_filt)) {
    pair_cors[i] <- cor(log2(counts_trans_filt[, tid] + 1),
                        log2(counts_translat_filt[, rid] + 1),
                        method = "spearman")
  } else {
    pair_cors[i] <- NA
  }
}

pair_map$spearman_cor <- pair_cors
write.csv(pair_map, file.path(dir_qc, "pair_correlation.csv"), row.names = FALSE)

p_pair <- ggplot(pair_map, aes(x = stage, y = spearman_cor, color = stage)) +
  geom_jitter(width = 0.15, size = 3) +
  stat_summary(fun = median, geom = "crossbar", width = 0.5, color = "black") +
  scale_color_manual(values = col_timepoint) +
  theme_oocyte() +
  labs(x = NULL, y = "Spearman correlation (paired samples)",
       title = "Within-pair correlation (transcriptome vs translatome)")

ggsave(file.path(dir_qc, "pair_correlation.pdf"), p_pair,
       width = 6, height = 4)

cat(sprintf("Pair correlation: median = %.3f, range = [%.3f, %.3f]\n",
            median(pair_cors, na.rm = TRUE),
            min(pair_cors, na.rm = TRUE),
            max(pair_cors, na.rm = TRUE)))

# ============================================================
# 9. Known gene control panel
# ============================================================

control_genes <- params$control_genes

in_annotation <- control_genes %in% gene_symbols
in_filtered <- control_genes %in% symbols_filt

ctrl_not_annotated <- control_genes[!in_annotation]
ctrl_filtered_out <- control_genes[in_annotation & !in_filtered]
ctrl_found <- control_genes[in_filtered]

if (length(ctrl_not_annotated) > 0) {
  cat(sprintf("Control genes NOT in annotation: %s\n",
              paste(ctrl_not_annotated, collapse = ", ")))
}
if (length(ctrl_filtered_out) > 0) {
  cat(sprintf("Control genes FILTERED OUT (low expression): %s\n",
              paste(ctrl_filtered_out, collapse = ", ")))
}
cat(sprintf("Control genes retained: %d/%d\n",
            length(ctrl_found), length(control_genes)))

# ============================================================
# 10. Outlier detection
# ============================================================

detect_outliers <- function(qc_df, metric, threshold_mad = 3, min_group_size = 5) {
  by_group <- split(qc_df, list(qc_df$omics, qc_df$stage))
  outliers <- character(0)
  for (g in by_group) {
    if (nrow(g) < min_group_size) {
      message(sprintf("Group %s.%s has n=%d, skipping outlier detection",
                      g$omics[1], g$stage[1], nrow(g)))
      next
    }
    vals <- g[[metric]]
    med <- median(vals)
    mad_val <- mad(vals, constant = 1.4826)
    if (mad_val > 0) {
      z <- abs(vals - med) / mad_val
      outliers <- c(outliers, g$sample_id[z > threshold_mad])
    }
  }
  outliers
}

outlier_libsize <- detect_outliers(qc_all, "lib_size")
outlier_detected <- detect_outliers(qc_all, "n_detected")
all_outliers <- unique(c(outlier_libsize, outlier_detected))

if (length(all_outliers) > 0) {
  cat(sprintf("WARNING: potential outlier samples (>3 SD): %s\n",
              paste(all_outliers, collapse = ", ")))
} else {
  cat("No outlier samples detected (3 SD threshold).\n")
}

# ============================================================
# 11. Remove confirmed outlier samples
# ============================================================

# F36 pair: extreme PCA outlier in both omics (confirmed post-normalization)
# F26 pair: low pair correlation (rho=0.654), transcriptome outlier
# F67.bam: MI-9 unpaired transcriptome, outlier post-normalization
exclude_samples <- c("F36.bam", "F36R.bam", "F26.bam", "F26R.bam", "F67.bam")

cat(sprintf("\nExcluding %d samples: %s\n",
            length(exclude_samples), paste(exclude_samples, collapse = ", ")))

# Remove from transcriptome
trans_exclude <- exclude_samples[exclude_samples %in% colnames(counts_trans_filt)]
counts_trans_filt <- counts_trans_filt[, !colnames(counts_trans_filt) %in% trans_exclude]
cat(sprintf("  Transcriptome: removed %d → %d samples\n",
            length(trans_exclude), ncol(counts_trans_filt)))

# Remove from translatome
translat_exclude <- exclude_samples[exclude_samples %in% colnames(counts_translat_filt)]
counts_translat_filt <- counts_translat_filt[, !colnames(counts_translat_filt) %in% translat_exclude]
cat(sprintf("  Translatome: removed %d → %d samples\n",
            length(translat_exclude), ncol(counts_translat_filt)))

# Update metadata
meta <- meta[!meta$sample_id %in% exclude_samples, ]

# Rebuild pair map after exclusion
pair_map <- build_pair_map(meta)
cat(sprintf("  Complete pairs after exclusion: %d\n", nrow(pair_map)))

# Document per-timepoint sample counts
for (s in levels(meta$stage)) {
  n_trans <- sum(meta$stage == s & meta$omics_type == "Transcriptome")
  n_translat <- sum(meta$stage == s & meta$omics_type == "Translatome")
  cat(sprintf("  %s: n=%d transcriptome, n=%d translatome\n", s, n_trans, n_translat))
}

# ============================================================
# 12. Save filtered counts
# ============================================================

save(counts_trans_filt, counts_translat_filt,
     gene_ids_filt, gene_ids_base_filt, gene_info_filt,
     symbols_filt, meta, pair_map, gene_symbols, pass_filter,
     all_outliers, exclude_samples,
     file = file.path(dir_qc, "qc_filtered_data.RData"))

cat(sprintf("\n=== QC Summary ===\n"))
cat(sprintf("Samples: %d transcriptome + %d translatome\n",
            ncol(counts_trans_filt), ncol(counts_translat_filt)))
cat(sprintf("Genes after filtering: %d\n", nrow(counts_trans_filt)))
cat(sprintf("Complete pairs: %d\n", nrow(pair_map)))
cat(sprintf("Outliers flagged (pre-exclusion): %d\n", length(all_outliers)))
cat(sprintf("Samples excluded: %d\n", length(exclude_samples)))
cat(sprintf("Results saved to: %s\n", dir_qc))

# ============================================================
# 13. Session info
# ============================================================

writeLines(capture.output(sessionInfo()),
           file.path(dir_qc, "sessionInfo.txt"))
