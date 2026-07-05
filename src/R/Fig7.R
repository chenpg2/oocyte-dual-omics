# Fig7.R — assemble Figure 7: GV-state predictability of trajectory
# Run from repo root under UTF-8 locale. Reads results/ produced by 01-07 + S1-S8.
source("src/R/00_config.R")
suppressPackageStartupMessages({
  library(ggplot2); library(patchwork); library(dplyr); library(scales)
  library(ragg); library(ggridges)
})

save_full <- function(p, path_stem, w, h) {
  grDevices::cairo_pdf(paste0(path_stem, ".pdf"),
                       width = w / 25.4, height = h / 25.4, family = "Helvetica")
  print(p); dev.off()
  svglite::svglite(paste0(path_stem, ".svg"), width = w / 25.4, height = h / 25.4)
  print(p); dev.off()
  ragg::agg_tiff(paste0(path_stem, ".tiff"),
                 width = w / 25.4, height = h / 25.4, units = "in", res = 600)
  print(p); dev.off()
  invisible(NULL)
}

{
  dir_f7  <- file.path(DIR_RESULTS, "figures", "F7_gv_prediction")
  dir_sup <- file.path(DIR_RESULTS, "supplementary")
  dir_cis <- file.path(DIR_RESULTS, "cis_grammar")

  mc  <- read.csv(file.path(dir_sup, "S3_model_comparison.csv"), stringsAsFactors = FALSE)
  fi  <- read.csv(file.path(dir_sup, "S3_gv_feature_importance.csv"), stringsAsFactors = FALSE)
  rob <- read.csv(file.path(dir_sup, "S5_robustness_matrix.csv"), stringsAsFactors = FALSE)

  # F7A — model AUC bar
  mc$label <- c("M1\n(GV omics)", "M2\n(+ seq)", "M3\n(seq only)", "M_full\n(MI-9 ref)")
  mc$is_target <- mc$model %in% c("M2_gv_plus_seq")
  mc$auc_use <- ifelse(is.na(mc$test_macro_auc), mc$cv_macro_auc, mc$test_macro_auc)
  mc$label <- factor(mc$label, levels = mc$label)

  p7a2 <- ggplot(mc, aes(x = label, y = auc_use, fill = is_target)) +
    geom_col(width = 0.65, alpha = 0.88) +
    geom_text(aes(label = sprintf("%.3f", auc_use), y = auc_use + 0.02),
              size = 2.4, family = "Helvetica") +
    scale_fill_manual(values = c("TRUE" = col_omics[["Translatome"]], "FALSE" = "#A7A9AC"),
                       guide = "none") +
    scale_y_continuous(limits = c(0, 1.08), expand = c(0, 0)) +
    labs(x = NULL, y = "Macro AUC", title = "GV-state trajectory\nprediction AUC") +
    theme_oocyte(base_size = 8) +
    theme(axis.text.x = element_text(size = 6, lineheight = 0.85))

  # F7B — feature importance
  fi <- fi[order(fi$Gain, decreasing = TRUE), ]
  fi$Feature_label <- gsub("_", " ", fi$Feature)
  fi$Feature_label <- factor(fi$Feature_label, levels = fi$Feature_label)
  col_feat <- ifelse(grepl("cpe|cpe_pas|proximal", fi$Feature, ignore.case = TRUE),
                     col_omics[["Translatome"]],
                     ifelse(grepl("mrna|te_gv", fi$Feature, ignore.case = TRUE),
                            col_omics[["Transcriptome"]], "#6E6E6E"))

  p7b2 <- ggplot(fi, aes(x = reorder(Feature_label, Gain), y = Gain,
                           fill = Feature_label)) +
    geom_col(fill = col_feat[order(fi$Gain)], alpha = 0.88, width = 0.65) +
    coord_flip() +
    scale_y_continuous(expand = expansion(mult = c(0, 0.06))) +
    labs(x = NULL, y = "Feature importance (Gain)",
         title = "Feature importance\n(M2: GV + sequence)") +
    theme_oocyte(base_size = 8)

  # F7C — per-class AUC
  auc_df2 <- read.csv(file.path(dir_cis, "test_auc.csv"), stringsAsFactors = FALSE)
  auc_df2 <- auc_df2[auc_df2$class != "MACRO", ]
  col_traj5 <- c(
    "Coordinated Clearance"      = unname(col_omics["Transcriptome"]),
    "Deep Coordinated Clearance" = "#4A4A4A",
    "Late Compensatory Buffering"= "#5B8DB8",
    "Mild Coordinated Clearance" = "#BCBCBC",
    "TE-Only Activation"         = unname(col_omics["Translatome"])
  )
  p7c2 <- ggplot(auc_df2, aes(x = reorder(class, auc), y = auc, fill = class)) +
    geom_col(width = 0.65, alpha = 0.88) +
    geom_hline(yintercept = 0.927, linetype = "dashed", colour = "grey40", linewidth = 0.5) +
    geom_text(aes(label = sprintf("%.3f", auc), y = auc - 0.015),
              hjust = 1.1, size = 2.5, colour = "white", fontface = "bold") +
    coord_flip() +
    scale_x_discrete(labels = function(x) sub("Late Compensatory Buffering", "Late Compensatory", x)) +
    scale_y_continuous(limits = c(0, 1.02), expand = c(0, 0)) +
    scale_fill_manual(values = col_traj5, guide = "none") +
    labs(x = NULL, y = "AUC (one-vs-rest)",
         title = "Per-class AUC\n(full model, circular ref.)") +
    theme_oocyte(base_size = 8) +
    theme(axis.text.y = element_text(size = 6.5))

  # F7D — confusion matrix (check.names=FALSE keeps class names with spaces intact)
  cm_raw2 <- read.csv(file.path(dir_cis, "confusion_matrix.csv"),
                       row.names = 1, stringsAsFactors = FALSE, check.names = FALSE)
  short_lbl <- c(
    "Coordinated Clearance"      = "CC",
    "Deep Coordinated Clearance" = "DCC",
    "Late Compensatory Buffering"= "LCB",
    "Mild Coordinated Clearance" = "MCC",
    "TE-Only Activation"         = "TEA"
  )
  cm_mat2   <- as.matrix(cm_raw2)
  cm_norm2  <- sweep(cm_mat2, 1, rowSums(cm_mat2), "/")
  cm_long2  <- expand.grid(True = rownames(cm_norm2), Predicted = colnames(cm_norm2),
                             stringsAsFactors = FALSE)
  cm_long2$value <- as.vector(cm_norm2)
  cm_long2$count <- as.vector(cm_mat2)
  cm_long2$True_s <- factor(short_lbl[cm_long2$True],      levels = rev(short_lbl))
  cm_long2$Pred_s <- factor(short_lbl[cm_long2$Predicted], levels = short_lbl)

  p7d2 <- ggplot(cm_long2, aes(x = Pred_s, y = True_s, fill = value)) +
    geom_tile(colour = "white", linewidth = 0.4) +
    geom_text(aes(label = ifelse(count > 0, count, "")),
              size = 2.2, family = "Helvetica",
              colour = ifelse(cm_long2$value > 0.55, "white", "grey25")) +
    scale_fill_gradient2(low = "white", mid = "#C8D8F0",
                          high = col_omics[["Transcriptome"]],
                          midpoint = 0.4, limits = c(0, 1),
                          name = "Recall",
                          labels = scales::percent_format(accuracy = 1)) +
    labs(x = "Predicted", y = "True",
         title = "Confusion matrix\n(full model, circular ref.)") +
    theme_oocyte(base_size = 8) +
    theme(axis.text.x = element_text(size = 7),
          axis.text.y = element_text(size = 7),
          legend.key.height = unit(0.8, "cm"))

  # F7E — robustness heatmap
  if (nrow(rob) > 0 && "conclusion" %in% names(rob)) {
    rob_long <- rob
    norms <- c("norm_A", "norm_B", "norm_C", "norm_D")
    norms_present <- norms[norms %in% names(rob)]
    if (length(norms_present) > 0) {
      rob_l <- tidyr::pivot_longer(rob,
                                    cols = all_of(norms_present),
                                    names_to = "norm", values_to = "result")
      rob_l$consistent_num <- ifelse(rob_l$consistent == "High", 3,
                                      ifelse(rob_l$consistent == "Medium", 2, 1))
      rob_l$conclusion <- sub("translation buffering", "translational compensation", rob_l$conclusion)
      p7e2 <- ggplot(rob_l, aes(x = norm, y = conclusion, fill = consistent_num)) +
        geom_tile(colour = "white") +
        scale_fill_gradient(low = "#F0F0F0", high = col_omics[["Translatome"]],
                             name = "Consistency",
                             breaks = 1:3, labels = c("Low", "Med", "High")) +
        scale_x_discrete(labels = function(x) sub("norm_", "", x)) +
        labs(x = "Normalization", y = NULL, title = "Robustness across\nnormalizations") +
        theme_oocyte(base_size = 8) +
        theme(axis.text.y = element_text(size = 5.5),
              legend.key.size = unit(3.2, "mm"),
              legend.title = element_text(size = 6))
    } else {
      p7e2 <- ggplot() + annotate("text", x=0.5, y=0.5, label="F7E") + theme_void()
    }
  } else {
    p7e2 <- ggplot() + annotate("text", x=0.5, y=0.5, label="F7E") + theme_void()
  }

  p7_full <- (p7a2 | p7b2 | p7c2) /
    (p7d2 | p7e2 | plot_spacer()) +
    plot_layout(heights = c(1, 1.1)) +
    plot_annotation(
      tag_levels = list(c("A", "B", "C", "D", "E")),
      theme = theme(plot.tag = element_text(size = 10, face = "bold",
                                             family = "Helvetica"))
    ) &
    theme(plot.title = element_text(size = 8, face = "bold"),
          plot.subtitle = element_text(size = 6.2, colour = "grey40"),
          plot.margin = margin(3, 5, 3, 4))

  save_full(p7_full, file.path(dir_f7, "F7_main_assembly"), w = 183, h = 155)
  cat("F7 main assembly (A–E) saved\n")
}
