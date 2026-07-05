# S4_functional_cpe.R — IU-S4: 功能性 CPE 特征工程
# 按 Piqué et al. 2008 (Genes Dev) 定义功能性 CPE：
#   CPE (UUUUAU/TTTTAT) 位于 PAS 上游 <100nt 时才具有激活活性
#
# 新增特征：
#   n_functional_cpe  : CPE 在任意 PAS 上游 <100nt 的数量
#   n_proximal_cpe    : 3'UTR 末端 200nt 内的功能性 CPE
#   min_cpe_pas_dist  : 最近 CPE-PAS 对距离（连续值，越小越可能激活）
#   has_cpe_pas_pair  : 是否有任意功能性 CPE（二元变量）
#   n_are_5prime      : ARE (ATTTA) 在 3'UTR 前半段的数量
#   au_content_3end   : 3'UTR 末端 50nt 的 AU 含量
#
# 输入：results/cis_grammar/utr3_sequences.txt.gz (本地 BioMart FASTA)
#        results/supplementary/S1_tbr_cs_matrix.csv
#        results/tta/cluster_assignments.csv
# 输出：results/supplementary/S4_functional_cpe_features.csv
#        results/supplementary/S4_cpe_enrichment.csv
#        results/supplementary/S4_cpe_enrichment.pdf
#        results/supplementary/S4_compensation_by_cpe.pdf

source("src/R/00_config.R")
library(Biostrings)
library(ggplot2)
library(dplyr)

set.seed(params$seed)
dir_sup   <- file.path(DIR_RESULTS, "supplementary")
dir_cis   <- file.path(DIR_RESULTS, "cis_grammar")
dir.create(dir_sup, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# 1. Load 3'UTR sequences (same parsing as 07_cis_grammar.R)
# ============================================================

local_fasta <- file.path(dir_cis, "utr3_sequences.txt.gz")
stopifnot(file.exists(local_fasta))

cat("Parsing local 3'UTR FASTA...\n")
con   <- gzfile(local_fasta, "rt")
lines <- readLines(con); close(con)

h_idx    <- grep("^>", lines)
rec_list <- vector("list", length(h_idx))
for (i in seq_along(h_idx)) {
  h_line  <- lines[h_idx[i]]
  end_idx <- if (i < length(h_idx)) h_idx[i + 1L] - 1L else length(lines)
  seq_str <- paste(lines[(h_idx[i] + 1L):end_idx], collapse = "")
  gene_id <- strsplit(sub("^>", "", h_line), "\\|")[[1]][1]
  rec_list[[i]] <- list(g = gene_id, s = seq_str, n = nchar(seq_str))
}

rec_df <- data.frame(
  gene_id = sapply(rec_list, `[[`, "g"),
  seq     = sapply(rec_list, `[[`, "s"),
  len     = sapply(rec_list, `[[`, "n"),
  stringsAsFactors = FALSE
)
valid  <- grepl("^[ATCGNatcgn]+$", rec_df$seq) & rec_df$len >= 10
rec_df <- rec_df[valid, ]

# Keep longest UTR per gene
rec_best <- do.call(rbind, lapply(split(rec_df, rec_df$gene_id), function(d)
  d[which.max(d$len), , drop = FALSE]))
rownames(rec_best) <- NULL
cat("  Unique genes with valid 3'UTR:", nrow(rec_best), "\n")

# ============================================================
# 2. Compute functional CPE features
# ============================================================

cat("Computing functional CPE features...\n")

# Helper functions
count_exact <- function(s, p) tryCatch(
  countPattern(DNAString(p), DNAString(s), fixed = TRUE),
  error = function(e) 0L)

get_positions <- function(s, p) {
  tryCatch({
    m <- matchPattern(DNAString(p), DNAString(s), fixed = TRUE)
    if (length(m) == 0) integer(0) else start(m)
  }, error = function(e) integer(0))
}

compute_functional_cpe <- function(seq_str) {
  n <- nchar(seq_str)

  # CPE positions: TTTTAT (DNA for UUUUAU)
  cpe_pos <- get_positions(seq_str, "TTTTAT")
  # PAS positions: AATAAA or ATTAAA
  pas_pos  <- c(get_positions(seq_str, "AATAAA"), get_positions(seq_str, "ATTAAA"))

  # --- Functional CPE: CPE upstream of PAS (CPE < PAS), distance < 100nt ---
  n_functional <- 0L
  min_dist     <- NA_real_

  if (length(cpe_pos) > 0 && length(pas_pos) > 0) {
    dists <- outer(cpe_pos, pas_pos, FUN = function(c, p) p - c)  # PAS - CPE position
    # Valid pairs: CPE is upstream of PAS (dist > 0) and within 100nt
    valid_pairs <- which(dists > 0 & dists < 100, arr.ind = TRUE)
    n_functional <- nrow(valid_pairs)
    if (n_functional > 0) {
      min_dist <- min(dists[valid_pairs], na.rm = TRUE)
    }
  }

  # --- Proximal functional CPE: functional CPE in last 200nt ---
  n_proximal <- 0L
  if (n_functional > 0 && length(pas_pos) > 0) {
    # CPE in the 3' end 200nt AND has a nearby PAS
    cpe_near_end <- cpe_pos[cpe_pos > max(0, n - 200)]
    pas_near_end <- pas_pos[pas_pos > max(0, n - 200)]
    if (length(cpe_near_end) > 0 && length(pas_near_end) > 0) {
      dists_prox <- outer(cpe_near_end, pas_near_end, FUN = function(c, p) p - c)
      n_proximal <- sum(dists_prox > 0 & dists_prox < 100, na.rm = TRUE)
    }
  }

  # --- ARE in 5' half of UTR ---
  n_are <- count_exact(seq_str, "ATTTA")
  half  <- ceiling(n / 2)
  n_are_5p <- if (half >= 5) count_exact(substr(seq_str, 1, half), "ATTTA") else 0L

  # --- AU content in last 50nt ---
  last50 <- if (n >= 50) substr(seq_str, n - 49, n) else seq_str
  au_cnt <- nchar(gsub("[^ATat]", "", last50))
  au_content_3end <- au_cnt / nchar(last50)

  list(
    n_functional_cpe  = n_functional,
    n_proximal_cpe    = n_proximal,
    has_cpe_pas_pair  = as.integer(n_functional > 0),
    min_cpe_pas_dist  = if (is.na(min_dist)) 9999L else as.integer(min_dist),
    n_are             = n_are,
    n_are_5prime      = n_are_5p,
    au_content_3end   = au_content_3end,
    utr3_len          = n
  )
}

feat_list <- vector("list", nrow(rec_best))
for (i in seq_len(nrow(rec_best))) {
  if (i %% 2000 == 0) cat("  Processing", i, "/", nrow(rec_best), "\r")
  feat_list[[i]] <- compute_functional_cpe(rec_best$seq[i])
}
cat("\n")

feat_df <- as.data.frame(do.call(rbind, lapply(feat_list, as.data.frame)))
feat_df$gene_id_clean <- rec_best$gene_id
feat_df[] <- lapply(feat_df, unlist)

cat("  Genes with >= 1 functional CPE:", sum(feat_df$n_functional_cpe > 0),
    sprintf("(%.1f%%)", mean(feat_df$n_functional_cpe > 0) * 100), "\n")
cat("  Genes with proximal CPE:", sum(feat_df$n_proximal_cpe > 0), "\n")
cat("  Min CPE-PAS distance (where functional):",
    round(mean(feat_df$min_cpe_pas_dist[feat_df$n_functional_cpe > 0]), 1), "nt mean\n")

write.csv(feat_df, file.path(dir_sup, "S4_functional_cpe_features.csv"), row.names = FALSE)

# ============================================================
# 3. Enrichment test: functional CPE in compensatory vs. clearance genes
# ============================================================

cat("\nEnrichment analysis...\n")

s1_df <- read.csv(file.path(dir_sup, "S1_tbr_cs_matrix.csv"), stringsAsFactors = FALSE)

# trajectory_type already in s1_df; clean gene IDs for merge with feat_df
s1_df$gene_id_clean <- strip_ensembl_version(s1_df$gene_id)
merged <- merge(s1_df, feat_df, by = "gene_id_clean")
cat("  Genes with UTR features + TBR/CS:", nrow(merged), "\n")

# Define groups
merged$comp_group <- with(merged, case_when(
  cs > 1 & tbr > 0.3  ~ "Compensatory",
  cs > 1 & tbr < 0.05 ~ "Coordinated Clearance",
  TRUE                 ~ "Other"
))

# Fisher exact test: has_cpe_pas_pair ~ comp_group
grp_comp <- merged[merged$comp_group == "Compensatory", ]
grp_clea <- merged[merged$comp_group == "Coordinated Clearance", ]

tbl <- matrix(c(
  sum(grp_comp$has_cpe_pas_pair == 1), sum(grp_comp$has_cpe_pas_pair == 0),
  sum(grp_clea$has_cpe_pas_pair == 1), sum(grp_clea$has_cpe_pas_pair == 0)
), nrow = 2, byrow = TRUE,
  dimnames = list(c("Compensatory", "Clearance"), c("Has_fCPE", "No_fCPE")))

fisher_res <- fisher.test(tbl)
cat("  Fisher exact test (Compensatory vs Clearance, functional CPE):\n")
cat("    OR =", round(fisher_res$estimate, 3),
    "  95% CI [", round(fisher_res$conf.int[1], 3), ",",
    round(fisher_res$conf.int[2], 3), "]\n")
cat("    p =", format(fisher_res$p.value, scientific = TRUE), "\n")
print(tbl)

# Wilcoxon: n_functional_cpe by TTA trajectory type
traj_levels <- c("TE-Only Activation", "Late Compensatory Buffering",
                 "Mild Coordinated Clearance", "Coordinated Clearance",
                 "Deep Coordinated Clearance")
merged_tta <- merged[merged$trajectory_type %in% traj_levels, ]
merged_tta$trajectory_type <- factor(merged_tta$trajectory_type, levels = traj_levels)

wtest_traj <- sapply(traj_levels, function(traj) {
  g1 <- merged$n_functional_cpe[merged$comp_group == "Compensatory"]
  g2 <- merged_tta$n_functional_cpe[merged_tta$trajectory_type == traj]
  if (length(g2) < 5) return(NA_real_)
  wilcox.test(g2, g1, alternative = "two.sided")$p.value
})
wtest_adj <- p.adjust(wtest_traj, method = "BH")

enrich_summary <- data.frame(
  group     = c("Compensatory", "Coordinated Clearance"),
  n         = c(nrow(grp_comp), nrow(grp_clea)),
  pct_has_fcpe = round(c(mean(grp_comp$has_cpe_pas_pair)*100,
                         mean(grp_clea$has_cpe_pas_pair)*100), 1),
  mean_n_fcpe  = round(c(mean(grp_comp$n_functional_cpe),
                         mean(grp_clea$n_functional_cpe)), 3),
  fisher_OR = c(round(fisher_res$estimate, 3), NA),
  fisher_p  = c(format(fisher_res$p.value, scientific=TRUE), NA),
  stringsAsFactors = FALSE
)
write.csv(enrich_summary, file.path(dir_sup, "S4_cpe_enrichment.csv"), row.names = FALSE)

# ============================================================
# 4. Visualizations
# ============================================================

cat("Generating plots...\n")

# Plot 1: % with functional CPE by trajectory type
traj_fcpe <- merged_tta %>%
  group_by(trajectory_type) %>%
  summarise(
    n             = n(),
    pct_has_fcpe  = mean(has_cpe_pas_pair == 1) * 100,
    mean_n_fcpe   = mean(n_functional_cpe),
    .groups = "drop"
  )

traj_colors <- c(
  "TE-Only Activation"          = col_timepoint[["GV"]],
  "Late Compensatory Buffering" = col_timepoint[["GVBD"]],
  "Mild Coordinated Clearance"  = col_timepoint[["MII"]],
  "Coordinated Clearance"       = col_timepoint[["MI-6"]],
  "Deep Coordinated Clearance"  = col_timepoint[["MI-9"]]
)

p_enrich <- ggplot(traj_fcpe, aes(x = trajectory_type, y = pct_has_fcpe,
                                   fill = trajectory_type)) +
  geom_col(alpha = 0.85, width = 0.6) +
  geom_text(aes(label = sprintf("%.1f%%\n(n=%d)", pct_has_fcpe, n)),
            vjust = -0.3, size = 2.5, color = "gray20") +
  scale_fill_manual(values = traj_colors, guide = "none") +
  scale_x_discrete(labels = function(x) gsub(" ", "\n", x)) +
  scale_y_continuous(limits = c(0, max(traj_fcpe$pct_has_fcpe) * 1.2)) +
  labs(x = NULL, y = "% genes with functional CPE\n(CPE within 100nt upstream of PAS)",
       title = "Functional CPE Enrichment by TTA Trajectory Type",
       subtitle = sprintf("Fisher OR (Compensatory vs Clearance) = %.2f, p = %s",
                          fisher_res$estimate,
                          format(fisher_res$p.value, digits=2, scientific=TRUE))) +
  theme_oocyte() +
  theme(axis.text.x = element_text(size = 7))

ggsave(file.path(dir_sup, "S4_cpe_enrichment.pdf"), p_enrich,
       width = 180, height = 120, units = "mm")

# Plot 2: n_functional_cpe distribution (violin) by compensation quadrant
merged_quad <- merged[merged$comp_group %in% c("Compensatory", "Coordinated Clearance"), ]
merged_quad$comp_group <- factor(merged_quad$comp_group,
                                 levels = c("Compensatory", "Coordinated Clearance"))

p_violin <- ggplot(merged_quad,
                   aes(x = comp_group, y = pmin(n_functional_cpe, 5))) +
  geom_violin(aes(fill = comp_group), alpha = 0.7, color = NA) +
  geom_boxplot(width = 0.1, outlier.shape = NA, color = "gray30", linewidth = 0.4) +
  scale_fill_manual(values = c("Compensatory"          = col_omics[["Translatome"]],
                               "Coordinated Clearance" = col_omics[["Transcriptome"]]),
                    guide = "none") +
  scale_y_continuous(breaks = 0:5) +
  labs(x = NULL, y = "n functional CPE (capped at 5)",
       title = "Functional CPE Count: Compensatory vs Clearance Genes",
       subtitle = sprintf("n = %d compensatory, %d clearance | Fisher p = %s",
                          nrow(grp_comp), nrow(grp_clea),
                          format(fisher_res$p.value, digits=2, scientific=TRUE))) +
  theme_oocyte()

ggsave(file.path(dir_sup, "S4_compensation_by_cpe.pdf"), p_violin,
       width = 130, height = 110, units = "mm")

# ============================================================
# 5. Summary
# ============================================================

cat("\n=== IU-S4 Summary ===\n")
cat("Genes analyzed:", nrow(feat_df), "\n")
cat("Genes with >= 1 functional CPE:",
    sum(feat_df$n_functional_cpe > 0),
    sprintf("(%.1f%%)\n", mean(feat_df$n_functional_cpe > 0) * 100))
cat("\nEnrichment (Compensatory vs Clearance):\n")
cat("  % with functional CPE — Compensatory:", tbl["Compensatory","Has_fCPE"],
    "/", nrow(grp_comp),
    sprintf("= %.1f%%\n", tbl["Compensatory","Has_fCPE"]/nrow(grp_comp)*100))
cat("  % with functional CPE — Clearance:", tbl["Clearance","Has_fCPE"],
    "/", nrow(grp_clea),
    sprintf("= %.1f%%\n", tbl["Clearance","Has_fCPE"]/nrow(grp_clea)*100))
cat("  Odds Ratio:", round(fisher_res$estimate, 3), "\n")
cat("  p-value:", format(fisher_res$p.value, scientific = TRUE), "\n")
cat("  Interpretation:",
    ifelse(fisher_res$p.value < 0.05,
           "Functional CPE IS enriched in compensatory genes -> cis regulation matters",
           "No significant enrichment -> motif counting (even functional) insufficient"),
    "\n")
cat("Results saved to:", dir_sup, "\n")

sink(file.path(dir_sup, "S4_sessionInfo.txt")); sessionInfo(); sink()
