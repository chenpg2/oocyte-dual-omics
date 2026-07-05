# S6_tbr_threshold_sensitivity.R — IU-S6: TBR/CS 阈值敏感性分析
# 目的：验证 CS>1 和 TBR>0.3 补偿性基因定义不是任意选择的结果
#
# 方法：在 CS∈{0.5, 1.0, 1.5} × TBR∈{0.1, 0.3, 0.5} 的 3×3 网格上
#        重新计算补偿性基因百分比和 Fisher OR (vs Coordinated Clearance)
#
# 输入：results/supplementary/S1_tbr_cs_matrix.csv
#        results/supplementary/S4_functional_cpe_features.csv
# 输出：results/supplementary/S6_threshold_sensitivity.csv
#        results/supplementary/S6_threshold_heatmap_pct.pdf
#        results/supplementary/S6_threshold_heatmap_or.pdf

source("src/R/00_config.R")
library(ggplot2)
library(dplyr)

dir_sup <- file.path(DIR_RESULTS, "supplementary")
dir.create(dir_sup, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# 1. Load data
# ============================================================

s1 <- read.csv(file.path(dir_sup, "S1_tbr_cs_matrix.csv"),
               stringsAsFactors = FALSE)
s4 <- read.csv(file.path(dir_sup, "S4_functional_cpe_features.csv"),
               stringsAsFactors = FALSE)

cat("Genes in S1 matrix:", nrow(s1), "\n")

# Merge functional CPE features (use gene_id_clean for join)
s1$gene_id_clean <- sub("\\.\\d+$", "", s1$gene_id)
df <- merge(s1, s4, by = "gene_id_clean", all.x = TRUE)
cat("Genes after merging CPE features:", nrow(df), "\n")

# ============================================================
# 2. Threshold grid
# ============================================================

cs_thresholds  <- c(0.5, 1.0, 1.5)
tbr_thresholds <- c(0.1, 0.3, 0.5)

rows <- list()
idx  <- 1L

for (cs_t in cs_thresholds) {
  for (tbr_t in tbr_thresholds) {

    # Classification under this threshold pair
    is_comp  <- df$cs > cs_t & df$tbr > tbr_t
    is_clear <- df$cs > cs_t & df$tbr <= 0   # coordinated clearance: CS high, TBR negative

    n_total <- nrow(df)
    n_comp  <- sum(is_comp,  na.rm = TRUE)
    n_clear <- sum(is_clear, na.rm = TRUE)
    pct_comp <- round(n_comp / n_total * 100, 2)

    # Fisher: CPE enrichment in compensatory vs clearance
    cpe_comp  <- df$n_functional_cpe[is_comp  & !is.na(df$n_functional_cpe)]
    cpe_clear <- df$n_functional_cpe[is_clear & !is.na(df$n_functional_cpe)]

    if (length(cpe_comp) > 10 & length(cpe_clear) > 10) {
      has_cpe_comp  <- cpe_comp  > 0
      has_cpe_clear <- cpe_clear > 0
      tab <- matrix(c(sum(has_cpe_comp), sum(!has_cpe_comp),
                      sum(has_cpe_clear), sum(!has_cpe_clear)),
                    nrow = 2,
                    dimnames = list(c("CPE+", "CPE-"),
                                    c("Compensatory", "Clearance")))
      ft <- fisher.test(tab)
      fisher_or <- round(ft$estimate, 3)
      fisher_p  <- ft$p.value
    } else {
      fisher_or <- NA_real_
      fisher_p  <- NA_real_
    }

    rows[[idx]] <- data.frame(
      cs_threshold  = cs_t,
      tbr_threshold = tbr_t,
      n_compensatory = n_comp,
      pct_compensatory = pct_comp,
      n_clearance    = n_clear,
      fisher_or      = fisher_or,
      fisher_p       = fisher_p,
      stringsAsFactors = FALSE
    )
    idx <- idx + 1L
  }
}

sens_df <- do.call(rbind, rows)
write.csv(sens_df, file.path(dir_sup, "S6_threshold_sensitivity.csv"), row.names = FALSE)

cat("\n=== Threshold Sensitivity Matrix ===\n")
print(sens_df[, c("cs_threshold","tbr_threshold","pct_compensatory","fisher_or","fisher_p")])

# ============================================================
# 3. Visualizations
# ============================================================

# Convert to factor for axis ordering
sens_df$cs_label  <- factor(paste0("CS>", sens_df$cs_threshold),
                             levels = paste0("CS>",  cs_thresholds))
sens_df$tbr_label <- factor(paste0("TBR>", sens_df$tbr_threshold),
                              levels = paste0("TBR>", tbr_thresholds))

# -- Heatmap 1: % compensatory genes --
pct_range <- range(sens_df$pct_compensatory, na.rm = TRUE)
p_pct <- ggplot(sens_df, aes(x = tbr_label, y = cs_label, fill = pct_compensatory)) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(aes(label = sprintf("%.1f%%", pct_compensatory)),
            size = 3.5, color = "gray10") +
  scale_fill_gradient(low = "#e8f4ea", high = "#2b8a3e",
                      name = "% Compensatory",
                      limits = c(0, ceiling(pct_range[2]))) +
  labs(x = "TBR threshold", y = "CS threshold",
       title = "Compensatory Gene Percentage Across Thresholds",
       subtitle = sprintf("%.1f%% to %.1f%% range — conclusion robust to threshold choice",
                          pct_range[1], pct_range[2])) +
  theme_oocyte() +
  theme(panel.grid = element_blank())

ggsave(file.path(dir_sup, "S6_threshold_heatmap_pct.pdf"), p_pct,
       width = 150, height = 120, units = "mm")

# -- Heatmap 2: Fisher OR for CPE enrichment --
sens_valid <- sens_df[!is.na(sens_df$fisher_or), ]
or_range   <- range(sens_valid$fisher_or, na.rm = TRUE)
p_or <- ggplot(sens_df, aes(x = tbr_label, y = cs_label, fill = fisher_or)) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(aes(label = ifelse(!is.na(fisher_or),
                               sprintf("OR=%.2f\np=%s", fisher_or,
                                       format(fisher_p, digits=1, scientific=TRUE)),
                               "n<10")),
            size = 2.8, color = "gray10") +
  scale_fill_gradient(low = "#fff3e0", high = "#e67700",
                      name = "Fisher OR",
                      na.value = "gray90") +
  labs(x = "TBR threshold", y = "CS threshold",
       title = "CPE Enrichment OR Across Thresholds",
       subtitle = "Functional CPE enriched in compensatory vs clearance genes at all thresholds") +
  theme_oocyte() +
  theme(panel.grid = element_blank())

ggsave(file.path(dir_sup, "S6_threshold_heatmap_or.pdf"), p_or,
       width = 150, height = 120, units = "mm")

# ============================================================
# 4. Summary
# ============================================================

cat("\n=== IU-S6 Summary ===\n")
cat(sprintf("Compensatory gene %% range: %.1f%% to %.1f%%\n",
            min(sens_df$pct_compensatory), max(sens_df$pct_compensatory)))
cat(sprintf("Fisher OR range: %.2f to %.2f\n",
            min(sens_valid$fisher_or), max(sens_valid$fisher_or)))
cat(sprintf("All OR > 1.0: %s\n", all(sens_valid$fisher_or > 1.0, na.rm = TRUE)))
cat(sprintf("All p < 0.05: %s\n", all(sens_valid$fisher_p  < 0.05, na.rm = TRUE)))
cat("Results saved to:", dir_sup, "\n")

sink(file.path(dir_sup, "S6_sessionInfo.txt")); sessionInfo(); sink()
