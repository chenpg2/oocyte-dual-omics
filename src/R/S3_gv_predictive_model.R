# S3_gv_predictive_model.R — IU-S3: GV 基态预测模型（消除循环性）
# 问题：IU-9 的 XGBoost AUC=0.927 存在循环性（特征和标签来自同一数据）
# 解决：仅用 GV 时间点特征预测后续轨迹——完全消除循环性
#
# 三种模型对比：
#   M1: GV-only omics (te_gv, mrna_gv)
#   M2: GV-only omics + sequence (te_gv, mrna_gv + UTR特征含功能性CPE)
#   M3: Sequence-only (UTR特征 only, 无任何组学信息)
#   M_full: IU-9 full model (基准对照, AUC=0.927)
#
# 输入：results/normalized/normalized_data.RData
#        results/tta/cluster_assignments.csv
#        results/supplementary/S4_functional_cpe_features.csv
# 输出：results/supplementary/S3_model_comparison.csv + S3_auc_comparison.pdf
#        results/supplementary/S3_gv_feature_importance.pdf

source("src/R/00_config.R")
library(xgboost)
library(pROC)
library(ggplot2)
library(dplyr)

set.seed(params$seed)
dir_sup <- file.path(DIR_RESULTS, "supplementary")
dir.create(dir_sup, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# 1. Load data
# ============================================================

load(file.path(DIR_RESULTS, "normalized", "normalized_data.RData"))
tta  <- read.csv(file.path(DIR_RESULTS, "tta", "cluster_assignments.csv"),
                 stringsAsFactors = FALSE)
cpe  <- read.csv(file.path(dir_sup, "S4_functional_cpe_features.csv"),
                 stringsAsFactors = FALSE)

tp_levels <- c("0h", "3h", "6h", "9h", "12h")

# GV group means (timepoint 0h only)
pm_gv   <- pair_map[pair_map$time_point == "0h", ]
m_gv    <- pm_gv$transcriptome_id
te_gv   <- grep("^0h_rep", colnames(te_A), value = TRUE)
mrna_gv <- rowMeans(norm_trans_A[gene_ids_filt, m_gv,  drop = FALSE], na.rm = TRUE)
te_gv_m <- rowMeans(te_A[gene_ids_filt,         te_gv, drop = FALSE], na.rm = TRUE)

gv_df <- data.frame(
  gene_id    = gene_ids_filt,
  mrna_gv    = log1p(mrna_gv),    # log1p because norm_trans_A is linear
  te_gv      = te_gv_m,
  stringsAsFactors = FALSE
)

# Merge sequence features
gv_df$gene_id_clean <- strip_ensembl_version(gv_df$gene_id)
gv_df <- merge(gv_df, cpe, by = "gene_id_clean", all.x = TRUE)

# Merge TTA labels
gv_df <- merge(gv_df, tta[, c("gene_id", "trajectory_type")], by = "gene_id")
gv_df <- gv_df[!is.na(gv_df$trajectory_type), ]

# Collapse small classes (same threshold as IU-9)
min_cls     <- params$cis_grammar$min_class_size
class_sizes <- table(gv_df$trajectory_type)
small_cls   <- names(class_sizes)[class_sizes < min_cls]
if (length(small_cls) > 0) {
  gv_df$trajectory_type[gv_df$trajectory_type %in% small_cls] <- "Other"
}
class_levels <- sort(unique(gv_df$trajectory_type))
n_classes    <- length(class_levels)
gv_df$label  <- as.integer(factor(gv_df$trajectory_type, levels = class_levels)) - 1L

cat("Genes with GV data + TTA labels:", nrow(gv_df), "\n")
cat("Classes:", paste(class_levels, collapse = " | "), "\n")
cat("Genes with sequence features (has_cpe_pas_pair not NA):",
    sum(!is.na(gv_df$n_functional_cpe)), "\n")

# ============================================================
# 2. Feature sets
# ============================================================

seq_cols <- c("n_functional_cpe", "n_proximal_cpe", "has_cpe_pas_pair",
              "n_are_5prime", "au_content_3end", "utr3_len",
              "min_cpe_pas_dist")
# Clip min_cpe_pas_dist for non-functional genes (9999 = no CPE)
gv_df$min_cpe_pas_dist[gv_df$min_cpe_pas_dist == 9999 | is.na(gv_df$min_cpe_pas_dist)] <- 200L

feat_sets <- list(
  M1_gv_omics  = c("mrna_gv", "te_gv"),
  M2_gv_plus_seq = c("mrna_gv", "te_gv", seq_cols),
  M3_seq_only  = seq_cols
)

# ============================================================
# 3. XGBoost params (same as IU-9)
# ============================================================

xgb_params <- list(
  objective        = "multi:softprob",
  num_class        = n_classes,
  eta              = 0.05,
  max_depth        = 5,
  subsample        = 0.8,
  colsample_bytree = 0.8,
  min_child_weight = 10,
  eval_metric      = "mlogloss"
)
k <- params$cis_grammar$cv_folds

multiclass_auc <- function(y_true, probs, classes) {
  sapply(seq_along(classes), function(ci)
    suppressMessages(as.numeric(auc(y_true == (ci - 1L), probs[, ci]))))
}

train_and_eval <- function(feat_cols, df, tag) {
  cat("\n--- Model:", tag, "---\n")
  df_m <- df[, c(feat_cols, "trajectory_type", "label"), drop = FALSE]
  # Impute with train-set medians (within CV folds — simplified: use full set medians)
  for (fc in feat_cols) {
    vals <- df_m[[fc]]
    if (!is.numeric(vals)) next
    med <- median(vals[is.finite(vals)], na.rm = TRUE)
    if (is.na(med) || !is.finite(med)) med <- 0
    df_m[[fc]][!is.finite(df_m[[fc]])] <- med
  }

  # Stratified 80/20 split
  set.seed(params$seed)
  train_idx <- unlist(lapply(class_levels, function(cl) {
    idx_cl  <- which(df_m$trajectory_type == cl)
    sample(idx_cl, ceiling(length(idx_cl) * params$cis_grammar$train_fraction))
  }))
  test_idx <- setdiff(seq_len(nrow(df_m)), train_idx)
  train_df <- df_m[train_idx, ]
  test_df  <- df_m[test_idx, ]

  X_train <- as.matrix(train_df[, feat_cols])
  y_train <- train_df$label
  X_test  <- as.matrix(test_df[, feat_cols])
  y_test  <- test_df$label

  # 5-fold CV
  fold_ids <- integer(nrow(train_df))
  for (cl in class_levels) {
    idx_cl <- sample(which(train_df$trajectory_type == cl))
    fold_ids[idx_cl] <- ((seq_along(idx_cl) - 1L) %% k) + 1L
  }

  cv_aucs    <- matrix(NA_real_, k, n_classes)
  best_iters <- integer(k)
  for (fold in seq_len(k)) {
    val_idx  <- which(fold_ids == fold)
    tr_idx   <- which(fold_ids != fold)
    dtrain_f <- xgb.DMatrix(X_train[tr_idx, ],  label = y_train[tr_idx])
    dval_f   <- xgb.DMatrix(X_train[val_idx, ], label = y_train[val_idx])
    m_fold   <- xgb.train(params = xgb_params, data = dtrain_f, nrounds = 500,
                          watchlist = list(val = dval_f),
                          early_stopping_rounds = 30, verbose = 0)
    best_iters[fold] <- m_fold$best_iteration
    probs_val <- matrix(predict(m_fold, dval_f), ncol = n_classes, byrow = TRUE)
    cv_aucs[fold, ] <- multiclass_auc(y_train[val_idx], probs_val, class_levels)
  }
  cv_macro <- mean(rowMeans(cv_aucs))
  cv_sd    <- sd(rowMeans(cv_aucs))

  # Final model
  final_nrounds <- max(50L, round(median(best_iters) * 1.1))
  dtrain_full   <- xgb.DMatrix(X_train, label = y_train)
  m_final       <- xgb.train(params = xgb_params, data = dtrain_full,
                              nrounds = final_nrounds, verbose = 0)
  dtest         <- xgb.DMatrix(X_test)
  probs_test    <- matrix(predict(m_final, dtest), ncol = n_classes, byrow = TRUE)
  test_aucs     <- multiclass_auc(y_test, probs_test, class_levels)
  macro_test    <- mean(test_aucs)

  cat(sprintf("  Train: %d | Test: %d | Features: %d\n",
              nrow(train_df), nrow(test_df), length(feat_cols)))
  cat(sprintf("  CV macro AUC: %.3f +/- %.3f\n", cv_macro, cv_sd))
  cat(sprintf("  Test macro AUC: %.3f\n", macro_test))

  imp <- xgb.importance(model = m_final)

  list(model = m_final, cv_macro = cv_macro, cv_sd = cv_sd,
       test_macro = macro_test, test_aucs = test_aucs,
       importance = imp, n_train = nrow(train_df), n_test = nrow(test_df))
}

# ============================================================
# 4. Run all models
# ============================================================

results <- list()
for (nm in names(feat_sets)) {
  cols <- feat_sets[[nm]]
  # For seq models, only use genes with sequence data
  if ("n_functional_cpe" %in% cols) {
    df_use <- gv_df[!is.na(gv_df$n_functional_cpe), ]
  } else {
    df_use <- gv_df
  }
  results[[nm]] <- train_and_eval(cols, df_use, nm)
}

# Reference: IU-9 full model AUC
full_model_auc <- 0.927

# ============================================================
# 5. Comparison table
# ============================================================

comp_df <- data.frame(
  model = c(names(results), "M_full (IU-9 reference)"),
  description = c(
    "GV omics only (mrna_gv + te_gv)",
    "GV omics + sequence features (incl. functional CPE)",
    "Sequence features only (no omics)",
    "Full IU-9 model (all 17 features, all timepoints)"
  ),
  n_features = c(sapply(feat_sets, length), 17L),
  cv_macro_auc = c(sapply(results, `[[`, "cv_macro"), NA),
  cv_sd        = c(sapply(results, `[[`, "cv_sd"), NA),
  test_macro_auc = c(sapply(results, `[[`, "test_macro"), full_model_auc),
  stringsAsFactors = FALSE
)
comp_df$cv_macro_auc <- round(comp_df$cv_macro_auc, 3)
comp_df$test_macro_auc <- round(comp_df$test_macro_auc, 3)

cat("\n=== Model Comparison ===\n")
print(comp_df[, c("model","n_features","cv_macro_auc","test_macro_auc")])
write.csv(comp_df, file.path(dir_sup, "S3_model_comparison.csv"), row.names = FALSE)

# ============================================================
# 6. Feature importance (M2 = GV omics + seq)
# ============================================================

imp_m2 <- results[["M2_gv_plus_seq"]]$importance
write.csv(imp_m2, file.path(dir_sup, "S3_gv_feature_importance.csv"), row.names = FALSE)

# ============================================================
# 7. Visualizations
# ============================================================

cat("Generating plots...\n")

# AUC comparison bar chart
plot_df <- comp_df
plot_df$model_label <- c("M1: GV omics\n(2 features)",
                         "M2: GV omics\n+ sequence",
                         "M3: Sequence\nonly",
                         "M_full: IU-9\n(17 features)")
plot_df$model_label <- factor(plot_df$model_label, levels = plot_df$model_label)
plot_df$fill_col <- c(col_omics[["Translatome"]], "steelblue3",
                      col_omics[["Transcriptome"]], "gray60")
plot_df$circularity <- c("No circularity", "No circularity",
                          "No circularity", "Circular (same data)")

p_comp <- ggplot(plot_df, aes(x = model_label, y = test_macro_auc)) +
  geom_col(fill = plot_df$fill_col, alpha = 0.85, width = 0.6) +
  geom_hline(yintercept = params$cis_grammar$auc_threshold,
             linetype = "dashed", color = "gray40", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.3f", test_macro_auc)),
            vjust = -0.4, size = 3, color = "gray20") +
  annotate("text", x = 4, y = full_model_auc * 0.97,
           label = "Circular\n(reference)", hjust = 0.5, size = 2.3, color = "gray50") +
  scale_y_continuous(limits = c(0, 1.0), breaks = seq(0, 1, 0.2)) +
  labs(x = NULL, y = "Test Macro AUC",
       title = "GV-State Predictive Models vs Full IU-9 Reference",
       subtitle = "M1-M3 use only GV timepoint features (no circularity)") +
  theme_oocyte() +
  theme(axis.text.x = element_text(size = 8))

ggsave(file.path(dir_sup, "S3_auc_comparison.pdf"), p_comp,
       width = 160, height = 110, units = "mm")

# Feature importance M2
top_n <- min(12L, nrow(imp_m2))
imp_plot <- imp_m2[seq_len(top_n), ]
imp_plot$Feature <- factor(imp_plot$Feature, levels = rev(imp_plot$Feature))
p_imp <- ggplot(imp_plot, aes(x = Gain, y = Feature)) +
  geom_col(fill = "steelblue3", alpha = 0.85) +
  labs(x = "XGBoost Gain", y = NULL,
       title = "Feature Importance: GV + Sequence Model (M2)") +
  theme_oocyte()

ggsave(file.path(dir_sup, "S3_gv_feature_importance.pdf"), p_imp,
       width = 140, height = 100, units = "mm")

# ============================================================
# 8. Summary
# ============================================================

cat("\n=== IU-S3 Summary ===\n")
m1_auc <- results[["M1_gv_omics"]]$test_macro
m2_auc <- results[["M2_gv_plus_seq"]]$test_macro
m3_auc <- results[["M3_seq_only"]]$test_macro

cat(sprintf("M1 (GV omics only):     Test AUC = %.3f\n", m1_auc))
cat(sprintf("M2 (GV + sequence):     Test AUC = %.3f\n", m2_auc))
cat(sprintf("M3 (sequence only):     Test AUC = %.3f\n", m3_auc))
cat(sprintf("M_full (IU-9 ref):      Test AUC = %.3f (circular)\n", full_model_auc))

cat("\nInterpretation:\n")
if (m1_auc > 0.70) {
  cat("  GV-only omics AUC > 0.70: Trajectory fate is largely determined at GV stage (STRONG)\n")
} else if (m1_auc > 0.60) {
  cat("  GV-only omics AUC 0.60-0.70: Partial early determination, more emerges during maturation\n")
} else {
  cat("  GV-only omics AUC < 0.60: Fate is primarily determined dynamically, not at GV baseline\n")
}

delta_seq <- m2_auc - m1_auc
cat(sprintf("  Adding sequence to GV omics: dAUC = +%.3f", delta_seq))
if (abs(delta_seq) > 0.02) {
  cat(" -> Sequence features add meaningful information\n")
} else {
  cat(" -> Sequence features add negligible information\n")
}
cat("Results saved to:", dir_sup, "\n")

sink(file.path(dir_sup, "S3_sessionInfo.txt")); sessionInfo(); sink()
