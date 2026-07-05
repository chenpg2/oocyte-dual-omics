# 07_cis_grammar.R — Cis-Regulatory Grammar (IU-9)
# H4: 3'UTR architecture + omics features predict mRNA fate trajectory (TTA type)
#
# Tier 1 features (omics-derived, always available):
#   te_gv, delta_te, te_var, mrna_gv, delta_mrna, mrna_var,
#   tpi, max_fc_te, max_fc_mrna, max_abs_vte
# Tier 2 features (sequence-derived, cached via biomaRt):
#   utr3_len_log, cds_len_log, utr3_gc, n_cpe, n_pas, n_are, n_m6a, n_rg4
#
# Model: XGBoost multiclass; 5-fold stratified CV + 20% held-out test
# Success: macro AUC >= params$cis_grammar$auc_threshold (0.65)
#
# Input:  results/normalized/normalized_data.RData,
#         results/tta/cluster_assignments.csv,
#         results/tpi/tpi_scores.csv,
#         results/trd/rate_matrix.csv
# Output: results/cis_grammar/

source("src/R/00_config.R")

library(xgboost)
library(pROC)
library(biomaRt)
library(Biostrings)
library(ggplot2)
library(reshape2)

set.seed(params$seed)
dir_cis <- file.path(DIR_RESULTS, "cis_grammar")
dir.create(dir_cis, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# 1. Load data
# ============================================================

load(file.path(DIR_RESULTS, "normalized", "normalized_data.RData"))
tta   <- read.csv(file.path(DIR_RESULTS, "tta",  "cluster_assignments.csv"),
                  stringsAsFactors = FALSE)
tpi   <- read.csv(file.path(DIR_RESULTS, "tpi",  "tpi_scores.csv"),
                  stringsAsFactors = FALSE)
rates <- read.csv(file.path(DIR_RESULTS, "trd",  "rate_matrix.csv"),
                  stringsAsFactors = FALSE, check.names = FALSE)

cat("Loaded:", length(gene_ids_filt), "genes,", nrow(tta), "TTA labels\n")

tp_levels <- c("0h","3h","6h","9h","12h")

# ============================================================
# 2. Compute group-mean TE and mRNA per timepoint
# ============================================================

cat("Computing group-mean TE/mRNA per timepoint...\n")

mean_te   <- matrix(NA_real_, nrow = length(gene_ids_filt), ncol = 5,
                    dimnames = list(gene_ids_filt, tp_levels))
mean_mrna <- matrix(NA_real_, nrow = length(gene_ids_filt), ncol = 5,
                    dimnames = list(gene_ids_filt, tp_levels))

for (i in seq_along(tp_levels)) {
  tp       <- tp_levels[i]
  pairs_tp <- pair_map[pair_map$time_point == tp, ]
  m_cols   <- pairs_tp$transcriptome_id
  te_cols  <- grep(paste0("^", tp, "_rep"), colnames(te_A), value = TRUE)
  stopifnot(length(te_cols) == nrow(pairs_tp))
  mean_mrna[, i] <- rowMeans(norm_trans_A[gene_ids_filt, m_cols,   drop = FALSE], na.rm = TRUE)
  mean_te[, i]   <- rowMeans(te_A[gene_ids_filt,         te_cols,  drop = FALSE], na.rm = TRUE)
}

# ============================================================
# 3. Tier 1 features
# ============================================================

cat("Building Tier 1 features...\n")

# mRNA in norm_trans_A is in linear (not log2) scale → log1p transform
log1p_mrna <- log1p(mean_mrna)

tier1 <- data.frame(
  gene_id    = gene_ids_filt,
  te_gv      = mean_te[, "0h"],
  delta_te   = mean_te[, "12h"] - mean_te[, "0h"],
  te_var     = apply(mean_te,   1, var, na.rm = TRUE),
  mrna_gv    = log1p_mrna[, "0h"],
  delta_mrna = log1p_mrna[, "12h"] - log1p_mrna[, "0h"],
  mrna_var   = apply(log1p_mrna, 1, var, na.rm = TRUE),
  stringsAsFactors = FALSE
)

# TPI and max-FC features from tpi_scores.csv
tpi_sub <- tpi[, c("gene_id", "tpi", "max_fc_mRNA", "max_fc_TE")]
names(tpi_sub)[3:4] <- c("max_fc_mrna", "max_fc_te")
tier1 <- merge(tier1, tpi_sub, by = "gene_id", all.x = TRUE)
tier1$tpi[is.na(tier1$tpi)] <- 0  # genes with no substantial change → TPI = 0

# Max absolute vTE from rate matrix (IU-7)
vte_cols <- paste0("vTE_", c("GV→GVBD","GVBD→MI6","MI6→MI9","MI9→MII"))
rate_sub  <- rates[, c("gene_id", vte_cols)]
tier1     <- merge(tier1, rate_sub, by = "gene_id", all.x = TRUE)
tier1$max_abs_vte <- apply(abs(tier1[, vte_cols]), 1,
                            function(r) if (all(is.na(r))) NA_real_ else max(r, na.rm = TRUE))
tier1 <- tier1[, setdiff(colnames(tier1), vte_cols)]  # drop raw rate columns

cat("  Tier 1 features:", ncol(tier1) - 1, "\n")

# ============================================================
# 4. Tier 2 features (sequence-derived, cached)
# ============================================================

seq_cache   <- file.path(dir_cis, "seq_feat_cache.RData")
local_fasta <- file.path(dir_cis, "utr3_sequences.txt.gz")

# Helper motif functions (available in all branches below)
count_exact <- function(s, p) tryCatch(
  countPattern(DNAString(p), DNAString(s), fixed = TRUE), error = function(e) 0L)
count_regex <- function(s, p) { m <- gregexpr(p, s, perl = TRUE)[[1]]; sum(m > 0) }

compute_motifs <- function(gene_ids_vec, seqs_vec) {
  # gene_ids_vec and seqs_vec are parallel character vectors (one per gene, longest UTR)
  cat("  Computing motifs for", length(seqs_vec), "genes...\n")
  data.frame(
    gene_id_clean = gene_ids_vec,
    utr3_len_log  = log1p(nchar(seqs_vec)),
    utr3_gc       = sapply(seqs_vec, function(s)
                      as.numeric(letterFrequency(DNAString(s), letters = "GC")) / nchar(s)),
    n_cpe         = sapply(seqs_vec, count_exact, p = "TTTTAT"),
    n_pas         = sapply(seqs_vec, function(s)
                      count_exact(s, "AATAAA") + count_exact(s, "ATTAAA")),
    n_are         = sapply(seqs_vec, count_exact, p = "ATTTA"),
    n_m6a         = sapply(seqs_vec, count_regex, p = "[AGT][AG]AC[ACT]"),
    n_rg4         = sapply(seqs_vec, count_regex,
                           p = "GGG[ACGT]{1,7}GGG[ACGT]{1,7}GGG[ACGT]{1,7}GGG"),
    stringsAsFactors = FALSE
  )
}

if (file.exists(local_fasta)) {
  # --- Parse local BioMart FASTA export (priority over cache) ---
  cat("Parsing local 3'UTR FASTA:", local_fasta, "\n")

  con   <- gzfile(local_fasta, "rt")
  lines <- readLines(con)
  close(con)

  h_idx <- grep("^>", lines)
  cat("  FASTA records:", length(h_idx), "\n")

  rec_list <- vector("list", length(h_idx))
  for (i in seq_along(h_idx)) {
    h_line  <- lines[h_idx[i]]
    end_idx <- if (i < length(h_idx)) h_idx[i + 1L] - 1L else length(lines)
    seq_str <- paste(lines[(h_idx[i] + 1L):end_idx], collapse = "")
    # Header: >ENSMUSG...(no ver)|ENSMUSG...(ver)|ENSMUST...|ENSMUST...
    gene_id <- strsplit(sub("^>", "", h_line), "\\|")[[1]][1]
    rec_list[[i]] <- list(g = gene_id, s = seq_str, n = nchar(seq_str))
  }

  rec_df <- data.frame(
    gene_id = sapply(rec_list, `[[`, "g"),
    seq     = sapply(rec_list, `[[`, "s"),
    len     = sapply(rec_list, `[[`, "n"),
    stringsAsFactors = FALSE
  )

  # Filter invalid sequences (unavailable / non-DNA / too short)
  valid <- grepl("^[ATCGNatcgn]+$", rec_df$seq) & rec_df$len >= 10
  rec_df <- rec_df[valid, ]
  cat("  Valid UTR sequences:", nrow(rec_df), "\n")

  # Keep longest UTR per gene
  rec_best <- do.call(rbind, lapply(split(rec_df, rec_df$gene_id), function(d)
    d[which.max(d$len), , drop = FALSE]))
  rownames(rec_best) <- NULL

  seq_feat <- compute_motifs(rec_best$gene_id, rec_best$seq)
  save(seq_feat, file = seq_cache)
  cat("  Sequence features cached for", nrow(seq_feat), "genes.\n")

} else if (file.exists(seq_cache)) {
  cat("Loading cached sequence features...\n")
  load(seq_cache)  # loads: seq_feat (data.frame or NULL)

} else {
  # --- Fallback: fetch from biomaRt ---
  cat("Fetching sequence features from biomaRt...\n")
  seq_feat <- NULL

  gene_ids_clean <- strip_ensembl_version(gene_ids_filt)
  mart <- tryCatch(
    useEnsembl("ensembl", dataset = "mmusculus_gene_ensembl"),
    error = function(e) { message("WARN [biomaRt connect]: ", conditionMessage(e)); NULL }
  )

  if (!is.null(mart)) {
    batch_sz <- 300
    batches  <- split(gene_ids_clean, ceiling(seq_along(gene_ids_clean) / batch_sz))
    t0       <- proc.time()[3]
    seqs_raw <- list()

    for (i in seq_along(batches)) {
      if (proc.time()[3] - t0 > 300) {
        cat("  Sequence fetch time budget exceeded at batch", i, "\n"); break
      }
      cat("  Batch", i, "/", length(batches), "\r")
      s <- tryCatch(
        getSequence(id = batches[[i]], type = "ensembl_gene_id",
                    seqType = "3utr", mart = mart),
        error = function(e) { message("WARN [biomaRt seq ", i, "]: ", conditionMessage(e)); NULL }
      )
      if (!is.null(s) && nrow(s) > 0) seqs_raw[[i]] <- s
    }
    cat("\n")
    utr_df <- if (length(seqs_raw) > 0) do.call(rbind, seqs_raw) else NULL

    if (!is.null(utr_df) && nrow(utr_df) > 0) {
      utr_df   <- utr_df[nchar(utr_df$`3utr`) >= 10, ]
      utr_df$n <- nchar(utr_df$`3utr`)
      utr_best <- do.call(rbind, lapply(split(utr_df, utr_df$ensembl_gene_id), function(d)
        d[which.max(d$n), , drop = FALSE]))
      seq_feat <- compute_motifs(utr_best$ensembl_gene_id, utr_best$`3utr`)
    }
  }

  save(seq_feat, file = seq_cache)
  cat("Sequence features cached.\n")
}

# ============================================================
# 5. Merge all features + TTA labels
# ============================================================

cat("\nMerging features...\n")

feat_df <- tier1

if (!is.null(seq_feat) && nrow(seq_feat) > 0) {
  feat_df$gene_id_clean <- strip_ensembl_version(feat_df$gene_id)
  feat_df <- merge(feat_df, seq_feat, by = "gene_id_clean", all.x = TRUE)
  feat_df$gene_id_clean <- NULL
  cat("  Tier 2 features added:", sum(colnames(feat_df) %in% c(
    "utr3_len_log","cds_len_log","utr3_gc","n_cpe","n_pas","n_are","n_m6a","n_rg4")), "\n")
}

feat_df <- merge(feat_df, tta[, c("gene_id","trajectory_type")], by = "gene_id")
feat_df <- feat_df[!is.na(feat_df$trajectory_type), ]

# Collapse small classes
min_cls     <- params$cis_grammar$min_class_size
class_sizes <- table(feat_df$trajectory_type)
cat("Class distribution:\n"); print(class_sizes)
small_cls <- names(class_sizes)[class_sizes < min_cls]
if (length(small_cls) > 0) {
  feat_df$trajectory_type[feat_df$trajectory_type %in% small_cls] <- "Other"
  cat("Collapsed to 'Other':", paste(small_cls, collapse = ", "), "\n")
}

class_levels <- sort(unique(feat_df$trajectory_type))
n_classes    <- length(class_levels)
feat_df$label <- as.integer(factor(feat_df$trajectory_type, levels = class_levels)) - 1L
cat("Classes (", n_classes, "):", paste(class_levels, collapse = " | "), "\n")

feat_cols <- setdiff(colnames(feat_df), c("gene_id","trajectory_type","label"))
cat("Total features:", length(feat_cols), "\n")

# ============================================================
# 6. Stratified 80/20 train/test split
# ============================================================

set.seed(params$seed)
train_idx <- unlist(lapply(class_levels, function(cl) {
  idx_cl  <- which(feat_df$trajectory_type == cl)
  n_train <- ceiling(length(idx_cl) * params$cis_grammar$train_fraction)
  sample(idx_cl, n_train)
}))
test_idx <- setdiff(seq_len(nrow(feat_df)), train_idx)

train_df <- feat_df[train_idx, ]
test_df  <- feat_df[test_idx, ]
cat("\nTrain:", nrow(train_df), "| Test:", nrow(test_df), "\n")

# Impute NA / Inf / -Inf / NaN with train-set medians (no leakage)
train_medians <- sapply(feat_cols, function(fc) {
  vals <- train_df[[fc]]
  med  <- median(vals[is.finite(vals)], na.rm = TRUE)
  if (is.na(med) || !is.finite(med)) 0 else med
})
for (fc in feat_cols) {
  train_df[[fc]][!is.finite(train_df[[fc]])] <- train_medians[[fc]]
  test_df[[fc]][!is.finite(test_df[[fc]])]   <- train_medians[[fc]]
}

X_train <- as.matrix(train_df[, feat_cols])
y_train <- train_df$label
X_test  <- as.matrix(test_df[, feat_cols])
y_test  <- test_df$label

# ============================================================
# 7. 5-fold stratified cross-validation
# ============================================================

cat("\n5-fold CV...\n")

set.seed(params$seed)
k        <- params$cis_grammar$cv_folds
fold_ids <- integer(nrow(train_df))
for (cl in class_levels) {
  idx_cl <- sample(which(train_df$trajectory_type == cl))
  fold_ids[idx_cl] <- ((seq_along(idx_cl) - 1L) %% k) + 1L
}

xgb_params <- list(
  objective        = "multi:softprob",
  num_class        = n_classes,
  eta              = 0.05,
  max_depth        = 5,
  subsample        = 0.8,
  colsample_bytree = 0.8,
  min_child_weight = 10,
  eval_metric      = "mlogloss"
  # seed controlled via set.seed() above; XGBoost R package ignores params$seed
)

cv_aucs    <- matrix(NA_real_, k, n_classes, dimnames = list(paste0("fold", 1:k), class_levels))
best_iters <- integer(k)

multiclass_auc <- function(y_true, probs, classes) {
  sapply(seq_along(classes), function(ci)
    suppressMessages(as.numeric(auc(y_true == (ci - 1L), probs[, ci]))))
}

for (fold in seq_len(k)) {
  val_idx <- which(fold_ids == fold)
  tr_idx  <- which(fold_ids != fold)
  dtrain_f <- xgb.DMatrix(X_train[tr_idx, ],  label = y_train[tr_idx])
  dval_f   <- xgb.DMatrix(X_train[val_idx, ], label = y_train[val_idx])

  m_fold <- xgb.train(
    params                = xgb_params,
    data                  = dtrain_f,
    nrounds               = 500,
    watchlist             = list(val = dval_f),
    early_stopping_rounds = 30,
    verbose               = 0
  )
  best_iters[fold] <- m_fold$best_iteration

  probs_val        <- matrix(predict(m_fold, dval_f), ncol = n_classes, byrow = TRUE)
  cv_aucs[fold, ]  <- multiclass_auc(y_train[val_idx], probs_val, class_levels)
  cat("  Fold", fold, "macro AUC:", round(mean(cv_aucs[fold, ]), 3),
      "(best iter:", best_iters[fold], ")\n")
}

cv_macro <- rowMeans(cv_aucs)
cat("  Mean CV macro AUC:", round(mean(cv_macro), 3), "±", round(sd(cv_macro), 3), "\n")
write.csv(data.frame(fold = rownames(cv_aucs), macro_auc = cv_macro, cv_aucs),
          file.path(dir_cis, "cv_results.csv"), row.names = FALSE)

# ============================================================
# 8. Final model (train on all training data)
# ============================================================

cat("\nTraining final model...\n")
final_nrounds <- max(50L, round(median(best_iters) * 1.1))
cat("  nrounds:", final_nrounds, "\n")

dtrain_full <- xgb.DMatrix(X_train, label = y_train)
m_final <- xgb.train(params = xgb_params, data = dtrain_full,
                     nrounds = final_nrounds, verbose = 0)
xgb.save(m_final, file.path(dir_cis, "xgb_model.bin"))

# ============================================================
# 9. Test set evaluation
# ============================================================

cat("\nTest evaluation...\n")
dtest       <- xgb.DMatrix(X_test)
probs_test  <- matrix(predict(m_final, dtest), ncol = n_classes, byrow = TRUE)
test_aucs   <- multiclass_auc(y_test, probs_test, class_levels)
macro_test  <- mean(test_aucs)
pass_h4     <- macro_test >= params$cis_grammar$auc_threshold

cat("  Test macro AUC:", round(macro_test, 3),
    "—", if (pass_h4) "PASS (>= 0.65)" else "FAIL (< 0.65)", "\n")
for (ci in seq_len(n_classes))
  cat("  ", class_levels[ci], "AUC:", round(test_aucs[ci], 3), "\n")

auc_df <- data.frame(class = c(class_levels, "MACRO"),
                     auc   = c(test_aucs, macro_test),
                     stringsAsFactors = FALSE)
write.csv(auc_df, file.path(dir_cis, "test_auc.csv"), row.names = FALSE)

pred_class <- class_levels[apply(probs_test, 1, which.max)]
conf_mat   <- table(Predicted = pred_class, Actual = test_df$trajectory_type)
cat("\n  Confusion matrix:\n"); print(conf_mat)
write.csv(as.data.frame.matrix(conf_mat), file.path(dir_cis, "confusion_matrix.csv"))

# ============================================================
# 10. Feature importance
# ============================================================

cat("\nFeature importance...\n")
imp_mat <- xgb.importance(model = m_final)
cat("  Top 10 by gain:\n"); print(head(imp_mat, 10))
write.csv(imp_mat, file.path(dir_cis, "feature_importance.csv"), row.names = FALSE)

# SHAP values (class 0 for illustration; suppress if package unavailable)
shap_plot <- tryCatch({
  library(shapviz)
  n_shap <- min(2000L, nrow(X_train))
  # mshapviz for multi-class: returns one shapviz object per class
  sv_multi <- mshapviz(m_final, X_pred = X_train[seq_len(n_shap), , drop = FALSE])
  sv_importance(sv_multi, kind = "bar", max_display = 15L)
}, error = function(e) { message("WARN [shapviz]: ", conditionMessage(e)); NULL })

# ============================================================
# 11. Visualizations
# ============================================================

cat("\nGenerating plots...\n")

# Feature importance bar chart (top 15)
top_n  <- min(15L, nrow(imp_mat))
imp_df <- imp_mat[seq_len(top_n), ]
imp_df$Feature <- factor(imp_df$Feature, levels = rev(imp_df$Feature))

p_imp <- ggplot(imp_df, aes(x = Gain, y = Feature)) +
  geom_col(fill = col_omics[["Translatome"]], alpha = 0.85) +
  labs(x = "XGBoost Gain", y = NULL,
       title = "Feature Importance — Cis-Regulatory Grammar") +
  theme_oocyte()
ggsave(file.path(dir_cis, "feature_importance.pdf"), p_imp,
       width = 120, height = 100, units = "mm")

# Per-class AUC bar chart
auc_plot_df <- auc_df[auc_df$class != "MACRO", ]
p_auc <- ggplot(auc_plot_df, aes(x = reorder(class, auc), y = auc)) +
  geom_col(fill = col_omics[["Transcriptome"]], alpha = 0.85) +
  geom_hline(yintercept = params$cis_grammar$auc_threshold,
             linetype = "dashed", color = "gray40") +
  annotate("text", x = 0.7, y = params$cis_grammar$auc_threshold + 0.01,
           label = paste0("threshold=", params$cis_grammar$auc_threshold),
           hjust = 0, size = 2.5, color = "gray40") +
  coord_flip() +
  scale_y_continuous(limits = c(0.5, 1)) +
  labs(x = NULL, y = "AUC (one-vs-rest)",
       title = paste0("Test AUC by Class (macro=", round(macro_test, 3), ")")) +
  theme_oocyte()
ggsave(file.path(dir_cis, "test_auc_by_class.pdf"), p_auc,
       width = 140, height = 100, units = "mm")

# CV fold macro AUC
p_cv <- ggplot(data.frame(fold = seq_len(k), macro = cv_macro),
               aes(x = fold, y = macro)) +
  geom_line(color = col_omics[["Transcriptome"]]) +
  geom_point(color = col_omics[["Transcriptome"]], size = 2) +
  geom_hline(yintercept = params$cis_grammar$auc_threshold,
             linetype = "dashed", color = "gray40") +
  scale_x_continuous(breaks = seq_len(k)) +
  labs(x = "CV Fold", y = "Macro AUC", title = "5-Fold CV Macro AUC") +
  theme_oocyte()
ggsave(file.path(dir_cis, "cv_auc.pdf"), p_cv,
       width = 100, height = 80, units = "mm")

# SHAP beeswarm (class 0 = first class alphabetically)
if (!is.null(shap_plot))
  ggsave(file.path(dir_cis, "shap_importance.pdf"), shap_plot,
         width = 150, height = 120, units = "mm")

# ============================================================
# 12. Summary
# ============================================================

cat("\n=== Cis-Grammar Summary ===\n")
cat("Features used:", length(feat_cols), "\n")
cat("Train:", nrow(train_df), "| Test:", nrow(test_df), "\n")
cat("CV macro AUC:", round(mean(cv_macro), 3), "±", round(sd(cv_macro), 3), "\n")
cat("Test macro AUC:", round(macro_test, 3), "\n")
cat("H4 (AUC >= 0.65):", if (pass_h4) "PASS" else "FAIL", "\n")
cat("Top 3 features by gain:\n")
for (i in seq_len(min(3L, nrow(imp_mat))))
  cat("  ", imp_mat$Feature[i], "(gain=", round(imp_mat$Gain[i], 4), ")\n")
cat("Results saved to:", dir_cis, "\n")

sink(file.path(dir_cis, "sessionInfo.txt"))
sessionInfo()
sink()
