# 05_tpi.R — Temporal Precedence Index (IU-6)
#
# Input:  results/normalized/normalized_data.RData, results/tta/cluster_assignments.csv,
#         results/baseline/ (DESeq2 LRT results)
# Output: results/tpi/
#
# Analyses:
#   1. Group-mean fold-change per interval (mRNA and TE separately)
#   2. First substantial change detection (|log2FC| > threshold)
#   3. TPI calculation (TPI = interval_first_ΔTE - interval_first_ΔM)
#   4. Gene-set level inference (by TTA trajectory type)
#   5. GO pathway-level TPI enrichment (fgsea on continuous TPI)
#   6. Known gene validation
#   7. Bootstrap stability & sensitivity across normalization schemes
#
# Design note: per-interval unpaired t-tests (n=3-5) have insufficient power
# after genome-wide FDR correction. This module uses fold-change magnitudes
# on group means, pre-filtered by DESeq2 LRT (time-dependent genes only).

source("src/R/00_config.R")

library(matrixStats)
library(ggplot2)
library(clusterProfiler)
library(org.Mm.eg.db)

dir_tpi <- file.path(DIR_RESULTS, "tpi")
dir.create(dir_tpi, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# 1. Load data
# ============================================================

load(file.path(DIR_RESULTS, "normalized", "normalized_data.RData"))
tta_assignments <- read.csv(file.path(DIR_RESULTS, "tta", "cluster_assignments.csv"),
                            stringsAsFactors = FALSE)

lrt_trans <- read.csv(file.path(DIR_RESULTS, "baseline", "deseq2_lrt_transcriptome.csv"),
                      stringsAsFactors = FALSE)
lrt_translat <- read.csv(file.path(DIR_RESULTS, "baseline", "deseq2_lrt_translatome.csv"),
                         stringsAsFactors = FALSE)

cat("Loaded:", length(gene_ids_filt), "genes,", nrow(pair_map), "pairs\n")

tp_levels <- c("0h", "3h", "6h", "9h", "12h")
stage_levels <- c("GV", "GVBD", "MI-6", "MI-9", "MII")
intervals <- list(
  c("0h", "3h"), c("3h", "6h"), c("6h", "9h"), c("9h", "12h")
)
interval_labels <- c("GV→GVBD", "GVBD→MI6", "MI6→MI9", "MI9→MII")

# ============================================================
# 2. Pre-filter: DESeq2 LRT time-dependent genes
# ============================================================

sig_trans_ids <- lrt_trans$gene_id[!is.na(lrt_trans$padj) &
                                    lrt_trans$padj < params$tpi$fdr_threshold]
sig_translat_ids <- lrt_translat$gene_id[!is.na(lrt_translat$padj) &
                                          lrt_translat$padj < params$tpi$fdr_threshold]
sig_either <- union(sig_trans_ids, sig_translat_ids)
sig_either <- intersect(sig_either, gene_ids_filt)

cat("DESeq2 LRT time-dependent genes:\n")
cat("  Transcriptome:", length(sig_trans_ids), "\n")
cat("  Translatome:", length(sig_translat_ids), "\n")
cat("  Union (pre-filter):", length(sig_either), "\n")

# ============================================================
# 3. Group-mean fold-change per interval
# ============================================================

compute_interval_fc <- function(norm_trans, te_mat, pair_map, gene_ids) {
  n_genes <- length(gene_ids)
  n_iv <- length(intervals)

  fc_mRNA <- matrix(NA_real_, nrow = n_genes, ncol = n_iv,
                    dimnames = list(gene_ids, interval_labels))
  fc_TE   <- matrix(NA_real_, nrow = n_genes, ncol = n_iv,
                    dimnames = list(gene_ids, interval_labels))

  for (iv in seq_along(intervals)) {
    tp_from <- intervals[[iv]][1]
    tp_to   <- intervals[[iv]][2]

    pairs_from <- pair_map[pair_map$time_point == tp_from, ]
    pairs_to   <- pair_map[pair_map$time_point == tp_to, ]

    m_from <- norm_trans[gene_ids, pairs_from$transcriptome_id, drop = FALSE]
    m_to   <- norm_trans[gene_ids, pairs_to$transcriptome_id, drop = FALSE]

    te_from_cols <- grep(paste0("^", tp_from, "_rep"), colnames(te_mat), value = TRUE)
    te_to_cols   <- grep(paste0("^", tp_to, "_rep"), colnames(te_mat), value = TRUE)
    stopifnot(length(te_from_cols) == nrow(pairs_from),
              length(te_to_cols) == nrow(pairs_to))
    te_from <- te_mat[gene_ids, te_from_cols, drop = FALSE]
    te_to   <- te_mat[gene_ids, te_to_cols, drop = FALSE]

    mean_m_from <- rowMeans(m_from, na.rm = TRUE)
    mean_m_to   <- rowMeans(m_to, na.rm = TRUE)
    mean_te_from <- rowMeans(te_from, na.rm = TRUE)
    mean_te_to   <- rowMeans(te_to, na.rm = TRUE)

    fc_mRNA[, iv] <- log2((mean_m_to + 1) / (mean_m_from + 1))
    fc_TE[, iv]   <- log2((mean_te_to + 1) / (mean_te_from + 1))
  }

  list(fc_mRNA = fc_mRNA, fc_TE = fc_TE)
}

cat("\nComputing group-mean fold changes per interval...\n")
fc_res <- compute_interval_fc(norm_trans_A, te_A, pair_map, gene_ids_filt)

# ============================================================
# 4. First substantial change detection & TPI
# ============================================================

log2fc_threshold <- 0.5

find_first_substantial <- function(fc_row, threshold) {
  substantial <- which(abs(fc_row) >= threshold)
  if (length(substantial) == 0) return(NA_integer_)
  substantial[1]
}

t_first_mRNA_all <- apply(fc_res$fc_mRNA, 1, find_first_substantial,
                          threshold = log2fc_threshold)
t_first_TE_all   <- apply(fc_res$fc_TE, 1, find_first_substantial,
                          threshold = log2fc_threshold)

# Apply pre-filter: only genes with LRT significance in at least one omics
t_first_mRNA <- t_first_mRNA_all
t_first_TE   <- t_first_TE_all
t_first_mRNA[!(gene_ids_filt %in% sig_either)] <- NA
t_first_TE[!(gene_ids_filt %in% sig_either)] <- NA
names(t_first_mRNA) <- gene_ids_filt
names(t_first_TE)   <- gene_ids_filt

has_mRNA_change <- !is.na(t_first_mRNA)
has_TE_change   <- !is.na(t_first_TE)
has_either      <- has_mRNA_change | has_TE_change
has_both        <- has_mRNA_change & has_TE_change

tpi_raw <- rep(NA_real_, length(gene_ids_filt))
names(tpi_raw) <- gene_ids_filt
tpi_raw[has_both] <- t_first_TE[has_both] - t_first_mRNA[has_both]

tpi_category <- rep("No change", length(gene_ids_filt))
names(tpi_category) <- gene_ids_filt
tpi_category[!(gene_ids_filt %in% sig_either)] <- "Not time-dependent"
tpi_category[has_mRNA_change & !has_TE_change] <- "mRNA-only"
tpi_category[!has_mRNA_change & has_TE_change] <- "TE-only"
tpi_category[has_both & tpi_raw < 0]  <- "TE-leading"
tpi_category[has_both & tpi_raw > 0]  <- "mRNA-leading"
tpi_category[has_both & tpi_raw == 0] <- "Simultaneous"

cat("\nTPI category distribution:\n")
print(table(tpi_category))

cat("\n  |log2FC| threshold:", log2fc_threshold, "\n")
cat("  Genes with mRNA substantial change:", sum(has_mRNA_change), "\n")
cat("  Genes with TE substantial change:", sum(has_TE_change), "\n")
cat("  Genes with both (valid TPI):", sum(has_both), "\n")

# Maximum fold change across intervals (for continuous ranking)
max_abs_fc_mRNA <- apply(abs(fc_res$fc_mRNA[gene_ids_filt, , drop = FALSE]), 1,
                         max, na.rm = TRUE)
max_abs_fc_TE   <- apply(abs(fc_res$fc_TE[gene_ids_filt, , drop = FALSE]), 1,
                         max, na.rm = TRUE)

tpi_df <- data.frame(
  gene_id = gene_ids_filt,
  gene_symbol = symbols_filt,
  t_first_mRNA = t_first_mRNA,
  t_first_TE = t_first_TE,
  tpi = tpi_raw,
  tpi_category = tpi_category,
  t_first_mRNA_label = ifelse(is.na(t_first_mRNA), NA,
                               interval_labels[t_first_mRNA]),
  t_first_TE_label = ifelse(is.na(t_first_TE), NA,
                             interval_labels[t_first_TE]),
  max_fc_mRNA = max_abs_fc_mRNA,
  max_fc_TE = max_abs_fc_TE,
  stringsAsFactors = FALSE
)

tta_merged <- merge(tpi_df, tta_assignments[, c("gene_id", "cluster", "trajectory_type")],
                    by = "gene_id", all.x = TRUE)

write.csv(tta_merged, file.path(dir_tpi, "tpi_scores.csv"), row.names = FALSE)

# Save per-interval fold changes
fc_export <- data.frame(
  gene_id = gene_ids_filt,
  gene_symbol = symbols_filt,
  fc_res$fc_mRNA,
  fc_res$fc_TE,
  stringsAsFactors = FALSE
)
colnames(fc_export)[3:6] <- paste0("mRNA_", interval_labels)
colnames(fc_export)[7:10] <- paste0("TE_", interval_labels)
write.csv(fc_export, file.path(dir_tpi, "interval_fold_changes.csv"),
          row.names = FALSE)

cat("TPI scores saved.\n")

# ============================================================
# 5. TPI distribution visualization
# ============================================================

cat("\nGenerating TPI distribution plots...\n")

tpi_valid <- tta_merged[!is.na(tta_merged$tpi), ]
cat("  Genes with valid TPI:", nrow(tpi_valid), "\n")

if (nrow(tpi_valid) > 0) {
  p_hist <- ggplot(tpi_valid, aes(x = tpi)) +
    geom_histogram(binwidth = 1, fill = col_omics["Transcriptome"],
                   color = "white", alpha = 0.8) +
    geom_vline(xintercept = 0, linetype = "dashed", color = col_tpi["Simultaneous"]) +
    scale_x_continuous(breaks = -3:3,
                       labels = c("-3\n(TE leads)", "-2", "-1", "0\n(Simul.)",
                                  "1", "2", "3\n(mRNA leads)")) +
    labs(x = "TPI (interval_first_ΔTE - interval_first_ΔmRNA)",
         y = "Number of genes",
         title = "Genome-wide Temporal Precedence Index",
         subtitle = paste0("n=", nrow(tpi_valid),
                           " genes (LRT pre-filtered, |log2FC|≥", log2fc_threshold, ")")) +
    annotate("text", x = -2, y = Inf, vjust = 2, hjust = 0.5,
             label = paste0("TE-leading: ", sum(tpi_valid$tpi < 0),
                             "\nSimultaneous: ", sum(tpi_valid$tpi == 0),
                             "\nmRNA-leading: ", sum(tpi_valid$tpi > 0)),
             size = 3, family = "Helvetica") +
    theme_oocyte()
  ggsave(file.path(dir_tpi, "tpi_distribution.pdf"), p_hist,
         width = 170, height = 100, units = "mm")
}

cat_counts <- table(tta_merged$tpi_category)
cat_df <- data.frame(
  category = names(cat_counts),
  count = as.integer(cat_counts),
  stringsAsFactors = FALSE
)
cat_df$pct <- round(100 * cat_df$count / sum(cat_df$count), 1)

cat("\nFull TPI category summary:\n")
print(cat_df)

# ============================================================
# 6. TPI by TTA trajectory type
# ============================================================

cat("\nGene-set level TPI inference by TTA trajectory type...\n")

traj_types <- unique(tta_merged$trajectory_type[!is.na(tta_merged$trajectory_type)])
tpi_by_traj <- list()

for (tt in traj_types) {
  sub <- tpi_valid[tpi_valid$trajectory_type == tt, ]
  n_valid <- nrow(sub)
  if (n_valid < 5) {
    cat("  ", tt, ": n=", n_valid, " (too few for testing)\n")
    next
  }

  wt <- tryCatch(
    wilcox.test(sub$tpi, mu = 0, conf.int = TRUE),
    error = function(e) NULL
  )

  set.seed(params$seed + which(traj_types == tt))
  boot_medians <- replicate(params$tpi$bootstrap_n, {
    idx <- sample(nrow(sub), replace = TRUE)
    median(sub$tpi[idx])
  })

  tpi_by_traj[[tt]] <- data.frame(
    trajectory_type = tt,
    n_genes = n_valid,
    n_te_leading = sum(sub$tpi < 0),
    n_simultaneous = sum(sub$tpi == 0),
    n_mrna_leading = sum(sub$tpi > 0),
    median_tpi = median(sub$tpi),
    mean_tpi = mean(sub$tpi),
    boot_ci_lo = quantile(boot_medians, 0.025),
    boot_ci_hi = quantile(boot_medians, 0.975),
    wilcox_p = if (!is.null(wt)) wt$p.value else NA_real_,
    stringsAsFactors = FALSE
  )

  wt_p_str <- if (!is.null(wt)) format(wt$p.value, digits = 3) else "NA"
  cat("  ", tt, ": n=", n_valid,
      ", median TPI=", round(median(sub$tpi), 2),
      ", 95% CI=[", round(quantile(boot_medians, 0.025), 2),
      ",", round(quantile(boot_medians, 0.975), 2),
      "], Wilcoxon p=", wt_p_str, "\n")
}

if (length(tpi_by_traj) > 0) {
  tpi_by_traj_df <- do.call(rbind, tpi_by_traj)
  rownames(tpi_by_traj_df) <- NULL
  tpi_by_traj_df$wilcox_padj <- p.adjust(tpi_by_traj_df$wilcox_p, method = "BH")
  write.csv(tpi_by_traj_df, file.path(dir_tpi, "tpi_by_trajectory_type.csv"),
            row.names = FALSE)
} else {
  cat("  No trajectory types had enough genes for testing.\n")
  write.csv(data.frame(), file.path(dir_tpi, "tpi_by_trajectory_type.csv"),
            row.names = FALSE)
}

# Plot TPI by trajectory type
if (nrow(tpi_valid) > 0 && any(!is.na(tpi_valid$trajectory_type))) {
  tpi_valid$trajectory_type <- factor(tpi_valid$trajectory_type, levels = traj_types)
  p_by_traj <- ggplot(tpi_valid[!is.na(tpi_valid$trajectory_type), ],
                      aes(x = trajectory_type, y = tpi, fill = trajectory_type)) +
    geom_boxplot(outlier.size = 0.5, alpha = 0.7) +
    geom_hline(yintercept = 0, linetype = "dashed", color = col_tpi["Simultaneous"]) +
    scale_y_continuous(breaks = -3:3) +
    labs(x = "TTA Trajectory Type", y = "TPI",
         title = "TPI Distribution by Trajectory Type") +
    theme_oocyte() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1),
          legend.position = "none")
  ggsave(file.path(dir_tpi, "tpi_by_trajectory_type.pdf"), p_by_traj,
         width = 170, height = 120, units = "mm")
}

# ============================================================
# 7. First-change interval distribution
# ============================================================

cat("\nFirst-change interval distribution:\n")
first_change_mRNA_dist <- table(factor(interval_labels[t_first_mRNA[has_mRNA_change]],
                                        levels = interval_labels))
first_change_TE_dist   <- table(factor(interval_labels[t_first_TE[has_TE_change]],
                                        levels = interval_labels))
cat("  mRNA first change:\n"); print(first_change_mRNA_dist)
cat("  TE first change:\n");   print(first_change_TE_dist)

fc_df <- data.frame(
  interval = rep(interval_labels, 2),
  omics = rep(c("mRNA", "TE"), each = 4),
  count = c(as.integer(first_change_mRNA_dist), as.integer(first_change_TE_dist)),
  stringsAsFactors = FALSE
)
fc_df$interval <- factor(fc_df$interval, levels = interval_labels)

p_fc <- ggplot(fc_df, aes(x = interval, y = count, fill = omics)) +
  geom_bar(stat = "identity", position = "dodge", alpha = 0.8) +
  scale_fill_manual(values = c("mRNA" = col_omics["Transcriptome"],
                                "TE"  = col_omics["Translatome"])) +
  labs(x = "Interval", y = "Number of genes (first substantial change)",
       title = "Timing of First Substantial Change (|log2FC| ≥ 0.5)",
       fill = "Omics") +
  theme_oocyte()
ggsave(file.path(dir_tpi, "first_change_intervals.pdf"), p_fc,
       width = 170, height = 100, units = "mm")

# ============================================================
# 8. GO pathway-level TPI enrichment
# ============================================================

cat("\nPathway-level TPI enrichment (GSEA on continuous TPI)...\n")

if (nrow(tpi_valid) >= 50) {
  tpi_continuous <- tpi_valid$tpi + 0.1 * sign(tpi_valid$tpi) * tpi_valid$max_fc_TE
  names(tpi_continuous) <- tpi_valid$gene_id
  tpi_ranked <- sort(tpi_continuous, decreasing = TRUE)

  gene_ids_clean <- strip_ensembl_version(names(tpi_ranked))
  entrez_map <- tryCatch(
    clusterProfiler::bitr(gene_ids_clean, fromType = "ENSEMBL",
                           toType = "ENTREZID", OrgDb = org.Mm.eg.db),
    error = function(e) { cat("  WARN: bitr failed:", e$message, "\n"); NULL }
  )

  if (!is.null(entrez_map) && nrow(entrez_map) > 100) {
    ranked_entrez <- tpi_ranked[match(entrez_map$ENSEMBL,
                                       strip_ensembl_version(names(tpi_ranked)))]
    names(ranked_entrez) <- entrez_map$ENTREZID
    ranked_entrez <- ranked_entrez[!is.na(ranked_entrez)]
    ranked_entrez <- sort(ranked_entrez, decreasing = TRUE)

    gsea_res <- tryCatch(
      clusterProfiler::gseGO(
        geneList = ranked_entrez,
        OrgDb = org.Mm.eg.db,
        ont = "BP",
        minGSSize = params$tpi$gsea_min_set_size,
        maxGSSize = 500,
        pvalueCutoff = 0.05,
        seed = params$seed
      ),
      error = function(e) { cat("  WARN: gseGO failed:", e$message, "\n"); NULL }
    )

    if (!is.null(gsea_res) && nrow(as.data.frame(gsea_res)) > 0) {
      gsea_df <- as.data.frame(gsea_res)
      write.csv(gsea_df, file.path(dir_tpi, "gsea_tpi.csv"), row.names = FALSE)
      cat("  GSEA enriched terms:", nrow(gsea_df), "\n")
      top_terms <- head(gsea_df[order(gsea_df$pvalue), ], 5)
      for (i in seq_len(nrow(top_terms))) {
        cat("    ", top_terms$Description[i],
            " (NES=", round(top_terms$NES[i], 2),
            ", padj=", format(top_terms$p.adjust[i], digits = 3), ")\n")
      }
    } else {
      cat("  No significant GSEA terms found.\n")
      write.csv(data.frame(), file.path(dir_tpi, "gsea_tpi.csv"), row.names = FALSE)
    }
  } else {
    cat("  Insufficient gene ID mappings for GSEA.\n")
    write.csv(data.frame(), file.path(dir_tpi, "gsea_tpi.csv"), row.names = FALSE)
  }
} else {
  cat("  Too few valid TPI genes (", nrow(tpi_valid), ") for GSEA.\n")
  write.csv(data.frame(), file.path(dir_tpi, "gsea_tpi.csv"), row.names = FALSE)
}

# ============================================================
# 9. Known gene validation
# ============================================================

cat("\nKnown gene TPI validation:\n")
ctrl_genes <- params$control_genes
ctrl_tpi <- tta_merged[tta_merged$gene_symbol %in% ctrl_genes, ]
ctrl_tpi <- ctrl_tpi[match(ctrl_genes, ctrl_tpi$gene_symbol), ]
ctrl_tpi <- ctrl_tpi[!is.na(ctrl_tpi$gene_id), ]

for (i in seq_len(nrow(ctrl_tpi))) {
  fc_m <- fc_res$fc_mRNA[ctrl_tpi$gene_id[i], ]
  fc_t <- fc_res$fc_TE[ctrl_tpi$gene_id[i], ]
  cat("  ", ctrl_tpi$gene_symbol[i],
      " → TPI=", ifelse(is.na(ctrl_tpi$tpi[i]), "NA",
                         as.character(ctrl_tpi$tpi[i])),
      " (", ctrl_tpi$tpi_category[i], ")",
      " [mRNA:", ctrl_tpi$t_first_mRNA_label[i],
      ", TE:", ctrl_tpi$t_first_TE_label[i], "]",
      " FC_mRNA:", paste(round(fc_m, 2), collapse = "/"),
      " FC_TE:", paste(round(fc_t, 2), collapse = "/"),
      "\n")
}

write.csv(ctrl_tpi, file.path(dir_tpi, "known_gene_validation.csv"),
          row.names = FALSE)

ctrl_plot <- ctrl_tpi[!is.na(ctrl_tpi$tpi), ]
if (nrow(ctrl_plot) > 0) {
  ctrl_plot$gene_symbol <- factor(ctrl_plot$gene_symbol,
                                   levels = ctrl_plot$gene_symbol[order(ctrl_plot$tpi)])
  ctrl_plot$direction <- ifelse(ctrl_plot$tpi < 0, "TE-leading",
                                 ifelse(ctrl_plot$tpi > 0, "mRNA-leading",
                                        "Simultaneous"))
  p_ctrl <- ggplot(ctrl_plot, aes(x = gene_symbol, y = tpi, fill = direction)) +
    geom_col(alpha = 0.8) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    scale_fill_manual(values = col_tpi) +
    labs(x = "Control Gene", y = "TPI",
         title = "Known Gene TPI Validation") +
    coord_flip() +
    theme_oocyte()
  ggsave(file.path(dir_tpi, "known_gene_validation.pdf"), p_ctrl,
         width = 120, height = 100, units = "mm")
}

# ============================================================
# 10. Sensitivity across normalization schemes
# ============================================================

cat("\nSensitivity analysis across normalization schemes...\n")

run_tpi_for_norm <- function(norm_trans, te_mat, label) {
  fc_alt <- compute_interval_fc(norm_trans, te_mat, pair_map, gene_ids_filt)

  t1_m <- apply(fc_alt$fc_mRNA, 1, find_first_substantial,
                threshold = log2fc_threshold)
  t1_t <- apply(fc_alt$fc_TE, 1, find_first_substantial,
                threshold = log2fc_threshold)

  t1_m[!(gene_ids_filt %in% sig_either)] <- NA
  t1_t[!(gene_ids_filt %in% sig_either)] <- NA

  both <- !is.na(t1_m) & !is.na(t1_t)
  tpi_alt <- rep(NA_real_, length(gene_ids_filt))
  names(tpi_alt) <- gene_ids_filt
  tpi_alt[both] <- t1_t[both] - t1_m[both]

  common <- intersect(names(tpi_alt)[!is.na(tpi_alt)],
                      names(tpi_raw)[!is.na(tpi_raw)])
  if (length(common) < 10) {
    cat("  ", label, ": too few common genes (", length(common), ")\n")
    return(data.frame(normalization = label, n_common = length(common),
                      spearman_rho = NA, concordance_pct = NA,
                      median_tpi = NA, stringsAsFactors = FALSE))
  }

  rho <- cor(tpi_raw[common], tpi_alt[common], method = "spearman")
  concordance <- mean(sign(tpi_raw[common]) == sign(tpi_alt[common]))

  cat("  ", label, ": n_common=", length(common),
      ", rho=", round(rho, 3),
      ", concordance=", round(100 * concordance, 1), "%",
      ", median TPI=", round(median(tpi_alt[common]), 2), "\n")

  data.frame(normalization = label, n_common = length(common),
             spearman_rho = round(rho, 3),
             concordance_pct = round(100 * concordance, 1),
             median_tpi = round(median(tpi_alt[common]), 2),
             stringsAsFactors = FALSE)
}

sens_results <- list()
sens_results[[1]] <- run_tpi_for_norm(norm_trans_B, te_B, "DESeq2")
sens_results[[2]] <- run_tpi_for_norm(norm_trans_C, te_C, "TMM")
sens_results[[3]] <- run_tpi_for_norm(norm_trans_D, te_D, "RUVg")

sens_df <- do.call(rbind, sens_results)
rownames(sens_df) <- NULL
write.csv(sens_df, file.path(dir_tpi, "sensitivity_across_normalizations.csv"),
          row.names = FALSE)

# ============================================================
# 11. Summary
# ============================================================

cat("\n=== TPI Summary ===\n")
cat("Total genes:", length(gene_ids_filt), "\n")
cat("Pre-filter (LRT time-dependent):", length(sig_either), "\n")
cat("|log2FC| threshold:", log2fc_threshold, "\n")
cat("Genes with mRNA substantial change:", sum(has_mRNA_change), "\n")
cat("Genes with TE substantial change:", sum(has_TE_change), "\n")
cat("Genes with both (valid TPI):", sum(has_both), "\n")

if (sum(has_both) > 0) {
  cat("TPI distribution:\n")
  cat("  TE-leading (TPI<0):", sum(tpi_raw[has_both] < 0), "\n")
  cat("  Simultaneous (TPI=0):", sum(tpi_raw[has_both] == 0), "\n")
  cat("  mRNA-leading (TPI>0):", sum(tpi_raw[has_both] > 0), "\n")
  cat("  Median TPI:", median(tpi_raw[has_both]), "\n")
}

cat("Results saved to:", dir_tpi, "\n")

sink(file.path(dir_tpi, "sessionInfo.txt"))
sessionInfo()
sink()
