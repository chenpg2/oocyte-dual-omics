# S8_cpe_multivariate.R — IU-S8: 多变量逻辑回归（CPE 效应控制混杂因子）
# 目的：验证 S4 的 OR=1.619 (功能性CPE在补偿型基因中富集) 不是
#       UTR长度/AU含量/基础表达量混杂造成的
#
# 结果变量：is_compensatory (quadrant == "Compensatory" vs "Coordinated Clearance")
#            ← 与 S4 Fisher test 相同的分组定义
# 主要预测变量：has_cpe_pas_pair (0/1) ← 与 S4 Fisher test 相同
# 备选预测变量：n_functional_cpe (连续)
# 混杂变量：utr3_len, au_content_3end, mrna_gv
#
# 四个模型：
#   M0: 单变量 (has_cpe_pas_pair)  ← 重现 S4 Fisher OR
#   M1: + utr3_len + au_content_3end
#   M2: + mrna_gv (加入基础表达量)
#   M3: n_functional_cpe (连续版) + 全混杂
#
# 输入：results/supplementary/S4_functional_cpe_features.csv
#        results/supplementary/S1_tbr_cs_matrix.csv
#        results/normalized/normalized_data.RData
# 输出：results/supplementary/S8_*.{csv,pdf}

source("src/R/00_config.R")
library(ggplot2)
library(dplyr)

set.seed(params$seed)
dir_sup <- file.path(DIR_RESULTS, "supplementary")
dir.create(dir_sup, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# 1. Load data
# ============================================================

s1 <- read.csv(file.path(dir_sup, "S1_tbr_cs_matrix.csv"),
               stringsAsFactors = FALSE)
s4 <- read.csv(file.path(dir_sup, "S4_functional_cpe_features.csv"),
               stringsAsFactors = FALSE)

# GV-state expression (mrna_gv)
load(file.path(DIR_RESULTS, "normalized", "normalized_data.RData"))
pm_gv   <- pair_map[pair_map$time_point == "0h", ]
m_gv    <- pm_gv$transcriptome_id
mrna_gv_vec <- rowMeans(norm_trans_A[gene_ids_filt, m_gv, drop = FALSE], na.rm = TRUE)
gv_df <- data.frame(
  gene_id_clean = sub("\\.\\d+$", "", gene_ids_filt),
  mrna_gv       = log1p(mrna_gv_vec),
  stringsAsFactors = FALSE
)

# ============================================================
# 2. Merge — using same grouping as S4 (quadrant from S1)
# ============================================================

s1$gene_id_clean <- sub("\\.\\d+$", "", s1$gene_id)
df <- merge(s1, s4, by = "gene_id_clean", all.x = TRUE)
df <- merge(df, gv_df, by = "gene_id_clean", all.x = TRUE)

# Restrict to Compensatory vs Coordinated Clearance (matching S4 Fisher test)
df_rc <- df[df$quadrant %in% c("Compensatory", "Coordinated Clearance") &
              !is.na(df$has_cpe_pas_pair), ]
df_rc$is_compensatory <- as.integer(df_rc$quadrant == "Compensatory")

# Filter to genes with complete confounders
df_m <- df_rc[!is.na(df_rc$utr3_len) &
                !is.na(df_rc$au_content_3end) &
                !is.na(df_rc$mrna_gv), ]

cat(sprintf("Analysis dataset: %d genes (Compensatory=%d, Clearance=%d)\n",
            nrow(df_m), sum(df_m$is_compensatory), sum(1 - df_m$is_compensatory)))

# Scale continuous confounders
df_m$utr3_len_z   <- scale(df_m$utr3_len)[, 1]
df_m$au_content_z <- scale(df_m$au_content_3end)[, 1]
df_m$mrna_gv_z    <- scale(df_m$mrna_gv)[, 1]

# ============================================================
# 3. Logistic regression models
# ============================================================

fit_logistic <- function(formula_str, data, model_name, predictor_name) {
  fit <- glm(as.formula(formula_str), data = data, family = binomial(link = "logit"))
  s   <- summary(fit)$coefficients
  cis <- confint.default(fit)

  r <- predictor_name
  data.frame(
    model        = model_name,
    predictor    = r,
    beta         = s[r, "Estimate"],
    se           = s[r, "Std. Error"],
    z_val        = s[r, "z value"],
    p_value      = s[r, "Pr(>|z|)"],
    OR           = exp(s[r, "Estimate"]),
    OR_ci_low    = exp(cis[r, 1]),
    OR_ci_high   = exp(cis[r, 2]),
    n_genes      = nrow(data),
    n_comp       = sum(data$is_compensatory),
    AIC          = AIC(fit),
    stringsAsFactors = FALSE
  )
}

cat("\n--- Binary predictor: has_cpe_pas_pair ---\n")
m0 <- fit_logistic("is_compensatory ~ has_cpe_pas_pair",
                   df_m, "M0: Unadjusted", "has_cpe_pas_pair")
m1 <- fit_logistic("is_compensatory ~ has_cpe_pas_pair + utr3_len_z + au_content_z",
                   df_m, "M1: + UTR + AU content", "has_cpe_pas_pair")
m2 <- fit_logistic("is_compensatory ~ has_cpe_pas_pair + utr3_len_z + au_content_z + mrna_gv_z",
                   df_m, "M2: + GV expression", "has_cpe_pas_pair")

cat("\n--- Continuous predictor: n_functional_cpe ---\n")
m3 <- fit_logistic("is_compensatory ~ n_functional_cpe + utr3_len_z + au_content_z + mrna_gv_z",
                   df_m, "M3: n_func_cpe + all confounders", "n_functional_cpe")

all_res <- rbind(m0, m1, m2, m3)
all_res$OR         <- round(all_res$OR, 3)
all_res$OR_ci_low  <- round(all_res$OR_ci_low, 3)
all_res$OR_ci_high <- round(all_res$OR_ci_high, 3)
all_res$p_fmt      <- format(all_res$p_value, digits = 2, scientific = TRUE)

cat("\nOR summary:\n")
print(all_res[, c("model","predictor","OR","OR_ci_low","OR_ci_high","p_fmt","AIC")])
write.csv(all_res, file.path(dir_sup, "S8_logistic_or.csv"), row.names = FALSE)

# ============================================================
# 4. Continuous TBR linear regression
# ============================================================

cat("\n--- Linear regression: TBR ~ n_functional_cpe + confounders (all genes) ---\n")
df_all_m <- df[!is.na(df$n_functional_cpe) & !is.na(df$utr3_len) &
                 !is.na(df$au_content_3end) & !is.na(df$mrna_gv), ]
df_all_m$utr3_len_z   <- scale(df_all_m$utr3_len)[, 1]
df_all_m$au_content_z <- scale(df_all_m$au_content_3end)[, 1]
df_all_m$mrna_gv_z    <- scale(df_all_m$mrna_gv)[, 1]

fit_tbr <- lm(tbr ~ n_functional_cpe + utr3_len_z + au_content_z + mrna_gv_z, data = df_all_m)
tbr_s   <- summary(fit_tbr)$coefficients
r       <- "n_functional_cpe"
tbr_res <- data.frame(
  predictor = r,
  beta      = tbr_s[r, "Estimate"],
  se        = tbr_s[r, "Std. Error"],
  t_val     = tbr_s[r, "t value"],
  p_value   = tbr_s[r, "Pr(>|t|)"],
  R2_model  = summary(fit_tbr)$r.squared,
  n_genes   = nrow(df_all_m),
  stringsAsFactors = FALSE
)
cat(sprintf("TBR ~ n_functional_cpe (adjusted): beta=%.4f, p=%s\n",
            tbr_res$beta, format(tbr_res$p_value, scientific = TRUE, digits = 2)))
write.csv(tbr_res, file.path(dir_sup, "S8_tbr_linear.csv"), row.names = FALSE)

# ============================================================
# 5. Forest plot
# ============================================================

forest_df <- all_res[, c("model","OR","OR_ci_low","OR_ci_high","p_fmt")]
forest_df$model <- factor(forest_df$model, levels = rev(forest_df$model))
x_max <- min(max(forest_df$OR_ci_high) * 1.4, 3.5)

p_forest <- ggplot(forest_df, aes(x = OR, y = model)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray50", linewidth = 0.5) +
  geom_errorbarh(aes(xmin = OR_ci_low, xmax = OR_ci_high),
                 height = 0.2, color = col_omics[["Translatome"]], linewidth = 0.9) +
  geom_point(size = 3.5, color = col_omics[["Translatome"]]) +
  geom_text(aes(label = sprintf("OR=%.3f [%.3f, %.3f]  p=%s",
                                OR, OR_ci_low, OR_ci_high, p_fmt)),
            hjust = -0.05, size = 2.5, color = "gray20") +
  scale_x_continuous(limits = c(0.7, x_max)) +
  labs(x = "Odds Ratio (compensatory vs clearance)",
       y = NULL,
       title = "Functional CPE Enrichment After Confounder Adjustment",
       subtitle = "Outcome: Compensatory quadrant vs Coordinated Clearance (same as S4 Fisher test)") +
  theme_oocyte()

ggsave(file.path(dir_sup, "S8_forest_plot.pdf"), p_forest,
       width = 200, height = 110, units = "mm")

# ============================================================
# 6. Summary
# ============================================================

cat("\n=== IU-S8 Summary ===\n")
or_attenuation <- if (m0$OR > 1)
  (m0$OR - m2$OR) / (m0$OR - 1) * 100 else NA

cat(sprintf("M0 (unadjusted)       OR = %.3f [%.3f–%.3f]  p=%s\n",
            m0$OR, m0$OR_ci_low, m0$OR_ci_high, m0$p_fmt))
cat(sprintf("M1 (+ UTR/AU)         OR = %.3f [%.3f–%.3f]  p=%s\n",
            m1$OR, m1$OR_ci_low, m1$OR_ci_high, m1$p_fmt))
cat(sprintf("M2 (+ GV expression)  OR = %.3f [%.3f–%.3f]  p=%s\n",
            m2$OR, m2$OR_ci_low, m2$OR_ci_high, m2$p_fmt))

if (!is.na(or_attenuation)) {
  cat(sprintf("OR attenuation M0→M2: %.1f%%\n", or_attenuation))
}

if (m2$OR > 1 && m2$p_value < 0.05) {
  cat("CONCLUSION: Functional CPE enrichment is ROBUST to confounder adjustment\n")
} else {
  cat("CONCLUSION: Functional CPE effect needs careful interpretation after adjustment\n")
}

cat(sprintf("\nLinear TBR model: n_functional_cpe beta=%.4f, p=%s (all genes adjusted)\n",
            tbr_res$beta, format(tbr_res$p_value, scientific = TRUE, digits = 2)))

cat("\nResults saved to:", dir_sup, "\n")
sink(file.path(dir_sup, "S8_sessionInfo.txt")); sessionInfo(); sink()
