# S7_functional_cpe_go.R — IU-S7: 功能性CPE基因 GO 富集分析
# 目的：阐明功能性CPE+ 补偿性基因的生物学功能，揭示翻译补偿的分子基础
#
# 三个富集分析：
#   分析1: 功能性CPE+ vs CPE- 全基因组 GO 富集 (BP+MF)
#   分析2: 高TBR vs 低TBR 分组 GO 富集 (补偿性发生的机制)
#   分析3: 清除型 vs 补偿型 差异 GO 富集 (两种调控命运对比)
#
# 输入：results/supplementary/S4_functional_cpe_features.csv
#        results/supplementary/S1_tbr_cs_matrix.csv
# 输出：results/supplementary/S7_*.{csv,pdf}

source("src/R/00_config.R")
library(ggplot2)
library(dplyr)
library(clusterProfiler)
library(org.Mm.eg.db)
library(AnnotationDbi)

set.seed(params$seed)
dir_sup <- file.path(DIR_RESULTS, "supplementary")
dir.create(dir_sup, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# 1. Load and merge
# ============================================================

s1 <- read.csv(file.path(dir_sup, "S1_tbr_cs_matrix.csv"),
               stringsAsFactors = FALSE)
s4 <- read.csv(file.path(dir_sup, "S4_functional_cpe_features.csv"),
               stringsAsFactors = FALSE)

s1$gene_id_clean <- sub("\\.\\d+$", "", s1$gene_id)
df <- merge(s1, s4, by = "gene_id_clean", all.x = TRUE)
cat("Merged genes:", nrow(df), "\n")

# ============================================================
# 2. Gene ID → Entrez ID conversion
# ============================================================

entrez_map <- AnnotationDbi::select(
  org.Mm.eg.db,
  keys    = df$gene_id_clean,
  columns = c("ENTREZID", "SYMBOL"),
  keytype = "ENSEMBL"
)
entrez_map <- entrez_map[!is.na(entrez_map$ENTREZID), ]
entrez_map <- entrez_map[!duplicated(entrez_map$ENSEMBL), ]

df <- merge(df, entrez_map, by.x = "gene_id_clean", by.y = "ENSEMBL", all.x = TRUE)
cat("Genes with Entrez ID:", sum(!is.na(df$ENTREZID)), "/", nrow(df), "\n")

# Background = all genes with Entrez IDs
universe <- df$ENTREZID[!is.na(df$ENTREZID)]

# ============================================================
# 3. Helper: run GO enrichment + dotplot
# ============================================================

run_go <- function(gene_entrez, universe, label, outprefix) {
  gene_entrez <- unique(gene_entrez[!is.na(gene_entrez)])
  cat(sprintf("  %s: %d genes\n", label, length(gene_entrez)))
  if (length(gene_entrez) < 20) {
    cat("  Too few genes — skipping\n")
    return(NULL)
  }

  ego <- enrichGO(
    gene          = gene_entrez,
    universe      = universe,
    OrgDb         = org.Mm.eg.db,
    ont           = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.20,
    readable      = TRUE
  )
  if (is.null(ego) || nrow(ego@result) == 0) {
    cat("  No significant GO terms\n")
    return(NULL)
  }

  res <- ego@result[ego@result$p.adjust < 0.05, ]
  res <- res[order(res$p.adjust), ]
  write.csv(res, paste0(outprefix, ".csv"), row.names = FALSE)
  cat(sprintf("  Significant GO BP terms: %d\n", nrow(res)))

  if (nrow(res) == 0) {
    cat("  No significant terms — skipping plot\n")
    return(ego)
  }

  top_n <- min(20L, nrow(res))
  p <- dotplot(ego, showCategory = top_n) +
    labs(title = paste("GO BP Enrichment:", label)) +
    theme_oocyte() +
    theme(axis.text.y = element_text(size = 7))
  ggsave(paste0(outprefix, ".pdf"), p, width = 180, height = 150, units = "mm")
  invisible(ego)
}

# ============================================================
# 4. Analysis 1: CPE+ vs universe
# ============================================================

cat("\n--- Analysis 1: Functional CPE+ genes (n_functional_cpe >= 1) ---\n")
cpe_pos <- df$ENTREZID[!is.na(df$n_functional_cpe) & df$n_functional_cpe >= 1]
ego_cpe <- run_go(cpe_pos, universe, "CPE+ genes",
                  file.path(dir_sup, "S7_go_cpe_positive"))

# ============================================================
# 5. Analysis 2: High TBR vs Low TBR
# ============================================================

cat("\n--- Analysis 2: High TBR (top quartile) vs Low TBR (bottom quartile) ---\n")
tbr_q <- quantile(df$tbr, probs = c(0.25, 0.75), na.rm = TRUE)
high_tbr <- df$ENTREZID[!is.na(df$tbr) & df$tbr >= tbr_q[2]]
low_tbr  <- df$ENTREZID[!is.na(df$tbr) & df$tbr <= tbr_q[1]]
cat(sprintf("  High TBR (>Q3=%.2f): %d genes\n", tbr_q[2], length(high_tbr[!is.na(high_tbr)])))
cat(sprintf("  Low TBR  (<Q1=%.2f): %d genes\n", tbr_q[1], length(low_tbr[!is.na(low_tbr)])))

ego_high <- run_go(high_tbr, universe, "High TBR (compensatory)",
                   file.path(dir_sup, "S7_go_high_tbr"))

# ============================================================
# 6. Analysis 3: Compensatory vs Clearance (CS>1)
# ============================================================

cat("\n--- Analysis 3: Compensatory (CS>1, TBR>0.3) vs Clearance (CS>1, TBR<0) ---\n")
is_comp  <- !is.na(df$cs) & !is.na(df$tbr) & df$cs > 1 & df$tbr > 0.3
is_clear <- !is.na(df$cs) & !is.na(df$tbr) & df$cs > 1 & df$tbr < 0
comp_ent  <- df$ENTREZID[is_comp]
clear_ent <- df$ENTREZID[is_clear]
cat(sprintf("  Compensatory: %d genes | Clearance: %d genes\n",
            sum(is_comp), sum(is_clear)))

# Compare compensatory against clearance background
ego_comp_vs_clear <- run_go(comp_ent, clear_ent,
                             "Compensatory vs Clearance",
                             file.path(dir_sup, "S7_go_comp_vs_clear"))

# ============================================================
# 7. Summary table
# ============================================================

summary_rows <- list(
  data.frame(analysis = "CPE+ genes (universe: all)",
             n_genes = sum(!is.na(cpe_pos)),
             n_sig_terms = if (!is.null(ego_cpe)) nrow(ego_cpe@result[ego_cpe@result$p.adjust<0.05,]) else 0,
             stringsAsFactors = FALSE),
  data.frame(analysis = "High TBR Q4 (universe: all)",
             n_genes = sum(!is.na(high_tbr)),
             n_sig_terms = if (!is.null(ego_high)) nrow(ego_high@result[ego_high@result$p.adjust<0.05,]) else 0,
             stringsAsFactors = FALSE),
  data.frame(analysis = "Compensatory vs Clearance",
             n_genes = sum(is_comp),
             n_sig_terms = if (!is.null(ego_comp_vs_clear)) nrow(ego_comp_vs_clear@result[ego_comp_vs_clear@result$p.adjust<0.05,]) else 0,
             stringsAsFactors = FALSE)
)
summary_df <- do.call(rbind, summary_rows)
write.csv(summary_df, file.path(dir_sup, "S7_go_summary.csv"), row.names = FALSE)

cat("\n=== IU-S7 Summary ===\n")
print(summary_df)
cat("Results saved to:", dir_sup, "\n")

sink(file.path(dir_sup, "S7_sessionInfo.txt")); sessionInfo(); sink()
