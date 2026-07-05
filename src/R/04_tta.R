# 04_tta.R — Translational Trajectory Analysis (IU-5)
#
# Input:  results/normalized/normalized_data.RData, results/baseline/
# Output: results/tta/ (phase-space trajectories, Euclidean clustering, GO enrichment)
#
# Analyses:
#   1. Phase-space construction (mRNA × TE per gene across 5 timepoints)
#   2. Trajectory feature extraction (angles, step sizes, rotation, curvature)
#   3. Euclidean-based hierarchical clustering (all genes, no subsampling)
#   4. Trajectory type annotation
#   5. Bootstrap stability validation
#   6. Functional enrichment (GO/KEGG)
#   7. Sensitivity across normalization schemes

source("src/R/00_config.R")

library(cluster)
library(matrixStats)
library(clusterProfiler)
library(org.Mm.eg.db)
library(ggplot2)
library(reshape2)

dir_tta <- file.path(DIR_RESULTS, "tta")
dir.create(dir_tta, recursive = TRUE, showWarnings = FALSE)

compute_ari <- function(labels1, labels2) {
  n <- length(labels1)
  if (n < 2) return(NA_real_)
  ct <- table(labels1, labels2)
  a <- sum(choose(ct, 2))
  b1 <- sum(choose(rowSums(ct), 2))
  b2 <- sum(choose(colSums(ct), 2))
  expected <- b1 * b2 / choose(n, 2)
  max_index <- (b1 + b2) / 2
  if (max_index == expected) 1 else (a - expected) / (max_index - expected)
}

# ============================================================
# 1. Load data
# ============================================================

load(file.path(DIR_RESULTS, "normalized", "normalized_data.RData"))
load(file.path(DIR_RESULTS, "qc", "qc_filtered_data.RData"))

stopifnot(exists("te_A"), exists("norm_trans_A"), exists("norm_translat_A"))
stopifnot(exists("meta"), exists("pair_map"))
stopifnot(exists("gene_ids_filt"), exists("symbols_filt"))
stopifnot(identical(rownames(te_A), rownames(norm_trans_A)))

cat(sprintf("Loaded: %d genes, %d pairs\n", nrow(te_A), nrow(pair_map)))

stages <- c("GV", "GVBD", "MI-6", "MI-9", "MII")
tp_labels <- c("0h", "3h", "6h", "9h", "12h")
time_numeric <- c(0, 3, 6, 9, 12)

# ============================================================
# 2. Phase-space construction
# ============================================================

compute_phase_space <- function(norm_trans, norm_translat, te_mat,
                                 sample_meta, tp_labels) {
  n_genes <- nrow(norm_trans)
  n_tp <- length(tp_labels)

  M_mat <- matrix(NA, nrow = n_genes, ncol = n_tp,
                  dimnames = list(rownames(norm_trans), tp_labels))
  TE_mat <- matrix(NA, nrow = n_genes, ncol = n_tp,
                   dimnames = list(rownames(te_mat), tp_labels))

  for (t_idx in seq_len(n_tp)) {
    tp <- tp_labels[t_idx]
    te_cols <- grep(paste0("^", tp, "_"), colnames(te_mat), value = TRUE)

    trans_cols <- sample_meta$sample_id[
      sample_meta$omics_type == "Transcriptome" &
        sample_meta$time_point == tp]
    trans_cols <- trans_cols[trans_cols %in% colnames(norm_trans)]

    M_mat[, t_idx] <- rowMeans(log2(norm_trans[, trans_cols, drop = FALSE] + 1))

    if (length(te_cols) > 0) {
      TE_mat[, t_idx] <- rowMeans(te_mat[, te_cols, drop = FALSE], na.rm = TRUE)
    }
  }

  n_all_na <- sum(apply(TE_mat, 1, function(x) any(is.nan(x))))
  if (n_all_na > 0) {
    warning(sprintf("%d genes have NaN TE at one or more timepoints", n_all_na))
  }

  list(M = M_mat, TE = TE_mat)
}

cat("Computing phase-space coordinates (group means)...\n")
ps <- compute_phase_space(norm_trans_A, norm_translat_A, te_A, meta, tp_labels)

# Center per gene (subtract row mean) but preserve magnitude across genes
# Z-scoring would destroy magnitude differences — a gene dropping 4 log2
# would look identical to one dropping 0.01 after per-gene scaling
M_centered <- ps$M - rowMeans(ps$M, na.rm = TRUE)
TE_centered <- ps$TE - rowMeans(ps$TE, na.rm = TRUE)

# Scale both dimensions by a shared global SD so M and TE are comparable
global_sd_M <- sd(as.vector(M_centered), na.rm = TRUE)
global_sd_TE <- sd(as.vector(TE_centered), na.rm = TRUE)
stopifnot(global_sd_M > 0, global_sd_TE > 0)
M_z <- M_centered / global_sd_M
TE_z <- TE_centered / global_sd_TE

# Remove genes with zero variance (constant across timepoints)
keep_m <- apply(ps$M, 1, sd, na.rm = TRUE) > 0
keep_te <- apply(ps$TE, 1, sd, na.rm = TRUE) > 0
keep_finite <- apply(M_z, 1, function(x) all(is.finite(x))) &
  apply(TE_z, 1, function(x) all(is.finite(x)))
keep <- keep_m & keep_te & keep_finite

n_dropped <- sum(!keep)
cat(sprintf("Phase space: %d genes retained, %d dropped (zero variance or non-finite)\n",
            sum(keep), n_dropped))

M_z <- M_z[keep, ]
TE_z <- TE_z[keep, ]
gene_ids_tta <- rownames(M_z)

# ============================================================
# 3. Trajectory feature extraction
# ============================================================

extract_trajectory_features <- function(M_z, TE_z, gene_ids) {
  n_genes <- nrow(M_z)
  n_tp <- ncol(M_z)
  n_steps <- n_tp - 1

  angles <- matrix(NA, n_genes, n_steps)
  step_sizes <- matrix(NA, n_genes, n_steps)
  cross_products <- matrix(NA, n_genes, n_steps - 1)

  for (i in seq_len(n_steps)) {
    dm <- M_z[, i + 1] - M_z[, i]
    dte <- TE_z[, i + 1] - TE_z[, i]
    angles[, i] <- atan2(dte, dm)
    step_sizes[, i] <- sqrt(dm^2 + dte^2)
  }

  # Cross products for consecutive vectors (rotation direction)
  for (i in seq_len(n_steps - 1)) {
    dm1 <- M_z[, i + 1] - M_z[, i]
    dte1 <- TE_z[, i + 1] - TE_z[, i]
    dm2 <- M_z[, i + 2] - M_z[, i + 1]
    dte2 <- TE_z[, i + 2] - TE_z[, i + 1]
    cross_products[, i] <- dm1 * dte2 - dte1 * dm2
  }

  rotation_sign <- sign(rowSums(cross_products))

  # Net displacement angle (start to end)
  net_dm <- M_z[, n_tp] - M_z[, 1]
  net_dte <- TE_z[, n_tp] - TE_z[, 1]
  net_angle <- atan2(net_dte, net_dm)
  net_disp <- sqrt(net_dm^2 + net_dte^2)

  # Curvature: mean unsigned turning angle
  turning_angles <- matrix(NA, n_genes, n_steps - 1)
  for (i in seq_len(n_steps - 1)) {
    turning_angles[, i] <- abs(angles[, i + 1] - angles[, i])
    turning_angles[, i] <- pmin(turning_angles[, i], 2 * pi - turning_angles[, i])
  }
  curvature <- rowMeans(turning_angles)

  # Signed area (shoelace formula, closed polygon)
  signed_area <- numeric(n_genes)
  for (i in seq_len(n_steps)) {
    signed_area <- signed_area +
      (M_z[, i] * TE_z[, i + 1] - M_z[, i + 1] * TE_z[, i])
  }
  signed_area <- signed_area +
    (M_z[, n_tp] * TE_z[, 1] - M_z[, 1] * TE_z[, n_tp])
  signed_area <- signed_area / 2

  feat_df <- data.frame(
    gene_id = gene_ids,
    angle_1 = angles[, 1], angle_2 = angles[, 2],
    angle_3 = angles[, 3], angle_4 = angles[, 4],
    step_1 = step_sizes[, 1], step_2 = step_sizes[, 2],
    step_3 = step_sizes[, 3], step_4 = step_sizes[, 4],
    rotation = rotation_sign,
    net_angle = net_angle,
    net_displacement = net_disp,
    curvature = curvature,
    signed_area = signed_area,
    stringsAsFactors = FALSE
  )
  feat_df
}

cat("Extracting trajectory features...\n")
traj_features <- extract_trajectory_features(M_z, TE_z, gene_ids_tta)
cat(sprintf("Feature matrix: %d genes × %d features\n",
            nrow(traj_features), ncol(traj_features) - 1))

write.csv(traj_features, file.path(dir_tta, "trajectory_features.csv"),
          row.names = FALSE)

# ============================================================
# 4. Euclidean hierarchical clustering
# ============================================================

# Euclidean distance on concatenated [M_z, TE_z] vectors (10D per gene)
# DTW warping is meaningless on 5 fixed, aligned timepoints
cat("\nEuclidean clustering (all genes):\n")
phase_mat <- cbind(M_z, TE_z)
colnames(phase_mat) <- c(paste0("M_", tp_labels), paste0("TE_", tp_labels))
cat(sprintf("  Phase matrix: %d genes × %d dimensions\n", nrow(phase_mat), ncol(phase_mat)))

euc_dist <- dist(phase_mat, method = "euclidean")
hc <- hclust(euc_dist, method = "ward.D2")
dtw_genes <- rownames(phase_mat)  # keep variable name for downstream compatibility

# Test cluster range
cluster_range <- params$tta$cluster_range
cat("Testing cluster counts:\n")

best_k <- cluster_range[1]
best_sil <- -Inf
sil_results <- data.frame()

for (k in cluster_range) {
  cl <- cutree(hc, k = k)
  sil <- silhouette(cl, euc_dist)
  avg_sil <- mean(sil[, "sil_width"])

  sil_results <- rbind(sil_results, data.frame(
    k = k, avg_silhouette = avg_sil, stringsAsFactors = FALSE
  ))

  cat(sprintf("  k=%d: avg silhouette=%.3f\n", k, avg_sil))

  if (avg_sil > best_sil) {
    best_sil <- avg_sil
    best_k <- k
  }
}

cat(sprintf("Best k=%d (silhouette=%.3f)\n", best_k, best_sil))
write.csv(sil_results, file.path(dir_tta, "silhouette_scores.csv"), row.names = FALSE)

# Final cluster assignments
cl_final <- cutree(hc, k = best_k)

rm(euc_dist, hc)
gc(verbose = FALSE)

# ============================================================
# 5. All genes clustered directly (no subsampling needed with Euclidean)
# ============================================================

all_assignments <- cl_final
names(all_assignments) <- dtw_genes

# ============================================================
# 6. Trajectory type annotation
# ============================================================

# Compute median trajectory per cluster for annotation
annotate_clusters <- function(M_z, TE_z, assignments) {
  k <- max(assignments)
  medians <- list()

  for (cl in seq_len(k)) {
    members <- names(assignments)[assignments == cl]
    members <- members[members %in% rownames(M_z)]
    if (length(members) == 0) next

    med_m <- colMedians(M_z[members, , drop = FALSE])
    med_te <- colMedians(TE_z[members, , drop = FALSE])

    # Net mRNA change
    delta_m <- med_m[5] - med_m[1]
    # Net TE change
    delta_te <- med_te[5] - med_te[1]
    # Signed area (closed polygon)
    sa <- 0
    for (i in seq_len(4)) {
      sa <- sa + (med_m[i] * med_te[i + 1] - med_m[i + 1] * med_te[i])
    }
    sa <- sa + (med_m[5] * med_te[1] - med_m[1] * med_te[5])
    sa <- sa / 2

    # Temporal sub-classification: early = GV→MI-6, late = MI-9→MII
    delta_m_early <- med_m[3] - med_m[1]
    delta_m_late  <- med_m[5] - med_m[3]

    traj_type <- if (delta_m < -0.3 && delta_te > 0.3) {
      if (abs(delta_m_early) > abs(delta_m_late) * 1.5) {
        "Early Compensatory Buffering"
      } else if (abs(delta_m_late) > abs(delta_m_early) * 1.5) {
        "Late Compensatory Buffering"
      } else {
        "Compensatory Buffering"
      }
    } else if (delta_m < -0.3 && delta_te < -0.3) {
      if (delta_m < -3.0) {
        "Deep Coordinated Clearance"
      } else if (delta_m < -1.0) {
        "Coordinated Clearance"
      } else {
        "Mild Coordinated Clearance"
      }
    } else if (delta_m > 0.3 && delta_te > 0.3) {
      "Selective Activation"
    } else if (abs(delta_m) <= 0.3 && delta_te > 0.3) {
      "TE-Only Activation"
    } else if (abs(delta_m) <= 0.3 && delta_te < -0.3) {
      "Dissociated TE Decline"
    } else if (abs(delta_m) <= 0.3 && abs(delta_te) <= 0.3) {
      if (abs(sa) > 0.3) "Transient Utilization" else "Dormant Storage"
    } else if (delta_m < -0.3 && abs(delta_te) <= 0.3) {
      "Translate-to-Degrade"
    } else if (delta_m > 0.3 && delta_te < -0.3) {
      "Dissociated Accumulation"
    } else if (delta_m > 0.3 && abs(delta_te) <= 0.3) {
      "Transcription-Driven"
    } else {
      "Unclassified"
    }

    medians[[cl]] <- list(
      cluster = cl,
      n_genes = length(members),
      med_m = med_m,
      med_te = med_te,
      delta_m = delta_m,
      delta_te = delta_te,
      signed_area = sa,
      type = traj_type
    )
  }
  medians
}

cluster_info <- annotate_clusters(M_z, TE_z, all_assignments)

cat("\nTrajectory type annotation:\n")
type_map <- character(max(all_assignments))
for (cl in seq_along(cluster_info)) {
  info <- cluster_info[[cl]]
  type_map[info$cluster] <- info$type
  cat(sprintf("  C%d (%s): n=%d, Δm=%.2f, Δte=%.2f, area=%.2f\n",
              info$cluster, info$type, info$n_genes,
              info$delta_m, info$delta_te, info$signed_area))
}

empty_types <- which(type_map == "")
if (length(empty_types) > 0) {
  stop(sprintf("Clusters with no annotated members: %s",
               paste(empty_types, collapse = ", ")))
}

# Build assignment table
cluster_df <- data.frame(
  gene_id = names(all_assignments),
  cluster = all_assignments,
  trajectory_type = type_map[all_assignments],
  stringsAsFactors = FALSE
)
cluster_df$gene_symbol <- symbols_filt[match(cluster_df$gene_id, gene_ids_filt)]

write.csv(cluster_df, file.path(dir_tta, "cluster_assignments.csv"),
          row.names = FALSE)

# ============================================================
# 7. Median trajectory plots
# ============================================================

plot_median_trajectories <- function(cluster_info, stages) {
  plot_df <- data.frame()
  for (info in cluster_info) {
    plot_df <- rbind(plot_df, data.frame(
      cluster = paste0("C", info$cluster, ": ", info$type,
                       " (n=", info$n_genes, ")"),
      M = info$med_m,
      TE = info$med_te,
      stage = stages,
      time_idx = seq_along(stages),
      stringsAsFactors = FALSE
    ))
  }

  # Phase space plot
  p_phase <- ggplot(plot_df, aes(x = M, y = TE, color = cluster)) +
    geom_path(linewidth = 1, arrow = arrow(length = unit(0.15, "cm"),
                                            type = "closed")) +
    geom_point(data = plot_df[plot_df$time_idx == 1, ],
               shape = 16, size = 3) +
    geom_text(data = plot_df,
              aes(label = stage), size = 2.5, vjust = -0.8, show.legend = FALSE) +
    labs(title = "TTA median trajectories — mRNA × TE phase space",
         x = "mRNA (centered, global-scaled)", y = "TE (centered, global-scaled)") +
    theme_oocyte() +
    theme(legend.position = "right")

  # Dual time course
  plot_long <- rbind(
    data.frame(plot_df[, c("cluster", "stage", "time_idx")],
               value = plot_df$M, metric = "mRNA", stringsAsFactors = FALSE),
    data.frame(plot_df[, c("cluster", "stage", "time_idx")],
               value = plot_df$TE, metric = "TE", stringsAsFactors = FALSE)
  )

  p_time <- ggplot(plot_long, aes(x = time_idx, y = value, color = cluster)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2) +
    facet_wrap(~ metric, ncol = 1, scales = "free_y") +
    scale_x_continuous(breaks = seq_along(stages), labels = stages) +
    labs(title = "TTA cluster time courses",
         x = "Stage", y = "Centered value (global-scaled)") +
    theme_oocyte() +
    theme(legend.position = "right")

  list(phase = p_phase, time = p_time)
}

cat("Generating trajectory plots...\n")
plots <- plot_median_trajectories(cluster_info, stages)

ggsave(file.path(dir_tta, "cluster_medians_phase.pdf"),
       plots$phase, width = 9, height = 7)
ggsave(file.path(dir_tta, "cluster_medians_time.pdf"),
       plots$time, width = 9, height = 7)

# ============================================================
# 8. Bootstrap stability validation
# ============================================================

# dist() on ~19k genes is O(n²) per iteration; cap at 200 for feasibility
boot_n <- min(params$tta$bootstrap_n, 200)
if (boot_n < params$tta$bootstrap_n) {
  cat(sprintf("Bootstrap capped: %d → %d (Euclidean dist on %d genes per iteration)\n",
              params$tta$bootstrap_n, boot_n, length(gene_ids_tta)))
}
cat(sprintf("\nBootstrap stability (n=%d, all %d genes, Euclidean)...\n",
            boot_n, length(gene_ids_tta)))

bootstrap_ari_euc <- function(te_mat, norm_trans, norm_translat,
                               sample_meta, tp_labels, original_cl,
                               all_genes, k, n_boot) {
  ari_vals <- numeric(n_boot)

  for (b in seq_len(n_boot)) {
    set.seed(params$seed + b)

    n_tp <- length(tp_labels)
    M_boot <- matrix(NA, length(all_genes), n_tp)
    TE_boot <- matrix(NA, length(all_genes), n_tp)
    rownames(M_boot) <- all_genes
    rownames(TE_boot) <- all_genes

    for (t_idx in seq_len(n_tp)) {
      tp <- tp_labels[t_idx]
      te_cols <- grep(paste0("^", tp, "_"), colnames(te_mat), value = TRUE)
      trans_cols <- sample_meta$sample_id[
        sample_meta$omics_type == "Transcriptome" &
          sample_meta$time_point == tp]
      trans_cols <- trans_cols[trans_cols %in% colnames(norm_trans)]

      boot_te_cols <- sample(te_cols, length(te_cols), replace = TRUE)
      boot_trans_cols <- sample(trans_cols, length(trans_cols), replace = TRUE)

      M_boot[, t_idx] <- rowMeans(log2(
        norm_trans[all_genes, boot_trans_cols, drop = FALSE] + 1))
      TE_boot[, t_idx] <- rowMeans(
        te_mat[all_genes, boot_te_cols, drop = FALSE], na.rm = TRUE)
    }

    # Same centering as primary
    M_bc <- M_boot - rowMeans(M_boot, na.rm = TRUE)
    TE_bc <- TE_boot - rowMeans(TE_boot, na.rm = TRUE)
    sd_M_b <- sd(as.vector(M_bc), na.rm = TRUE)
    sd_TE_b <- sd(as.vector(TE_bc), na.rm = TRUE)
    if (sd_M_b == 0 || sd_TE_b == 0) { ari_vals[b] <- NA; next }
    M_bz <- M_bc / sd_M_b
    TE_bz <- TE_bc / sd_TE_b

    ok <- apply(M_bz, 1, function(x) all(is.finite(x))) &
      apply(TE_bz, 1, function(x) all(is.finite(x)))
    if (sum(ok) < k * 5) { ari_vals[b] <- NA; next }

    phase_boot <- cbind(M_bz[ok, ], TE_bz[ok, ])
    d_boot <- dist(phase_boot, method = "euclidean")
    hc_boot <- hclust(d_boot, method = "ward.D2")
    cl_boot <- cutree(hc_boot, k = k)
    names(cl_boot) <- rownames(phase_boot)

    common <- intersect(names(cl_boot), names(original_cl))
    if (length(common) < k * 5) { ari_vals[b] <- NA; next }
    ari_vals[b] <- compute_ari(original_cl[common], cl_boot[common])

    if (b %% 100 == 0) cat(sprintf("  bootstrap %d/%d, ARI=%.3f\n",
                                     b, n_boot, ari_vals[b]))
  }
  ari_vals
}

stopifnot(all(gene_ids_tta %in% rownames(te_A)))
stopifnot(all(gene_ids_tta %in% rownames(norm_trans_A)))

ari_vals <- bootstrap_ari_euc(
  te_A, norm_trans_A, norm_translat_A,
  meta, tp_labels, all_assignments,
  gene_ids_tta, best_k, boot_n
)
n_boot_genes <- length(gene_ids_tta)

ari_valid <- ari_vals[!is.na(ari_vals)]
ari_mean <- mean(ari_valid)
ari_ci <- quantile(ari_valid, c(0.025, 0.975))
ari_pass <- ari_mean > params$tta$ari_threshold

cat(sprintf("\nBootstrap ARI: mean=%.3f, 95%% CI=[%.3f, %.3f]\n",
            ari_mean, ari_ci[1], ari_ci[2]))
cat(sprintf("ARI threshold: %.2f — %s\n",
            params$tta$ari_threshold,
            if (ari_pass) "PASS" else "FAIL (labels are exploratory)"))

stability_df <- data.frame(
  metric = c("mean_ari", "ci_lower", "ci_upper", "n_bootstrap",
             "n_genes_tested", "threshold", "pass"),
  value = c(ari_mean, ari_ci[1], ari_ci[2], boot_n,
            n_boot_genes, params$tta$ari_threshold, as.numeric(ari_pass)),
  stringsAsFactors = FALSE
)
write.csv(stability_df, file.path(dir_tta, "stability_report.csv"),
          row.names = FALSE)

# ============================================================
# 9. Functional enrichment
# ============================================================

cat("\nFunctional enrichment (GO BP) per trajectory type...\n")

run_go_enrichment <- function(cluster_df, background_genes) {
  types <- unique(cluster_df$trajectory_type)
  all_results <- list()

  bg_entrez <- tryCatch(
    bitr(background_genes, fromType = "ENSEMBL", toType = "ENTREZID",
         OrgDb = org.Mm.eg.db)$ENTREZID,
    error = function(e) NULL
  )

  for (tt in types) {
    genes <- cluster_df$gene_id[cluster_df$trajectory_type == tt]
    genes_clean <- sub("\\.\\d+$", "", genes)

    entrez <- tryCatch(
      bitr(genes_clean, fromType = "ENSEMBL", toType = "ENTREZID",
           OrgDb = org.Mm.eg.db),
      error = function(e) NULL
    )

    if (is.null(entrez) || nrow(entrez) < 5) {
      cat(sprintf("  %s: <5 mapped genes, skipping\n", tt))
      next
    }

    ego <- tryCatch(
      enrichGO(gene = entrez$ENTREZID,
               universe = bg_entrez,
               OrgDb = org.Mm.eg.db,
               ont = "BP",
               pAdjustMethod = "BH",
               pvalueCutoff = 0.05,
               qvalueCutoff = 0.1,
               readable = TRUE),
      error = function(e) { cat(sprintf("  %s: GO error: %s\n", tt, e$message)); NULL }
    )

    if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
      res <- as.data.frame(ego)
      res$trajectory_type <- tt
      all_results[[tt]] <- head(res, 10)
      cat(sprintf("  %s: %d enriched terms (top: %s)\n",
                  tt, nrow(res), res$Description[1]))
    } else {
      cat(sprintf("  %s: no enriched terms\n", tt))
    }
  }

  if (length(all_results) > 0) do.call(rbind, all_results) else data.frame()
}

bg_genes <- sub("\\.\\d+$", "", gene_ids_tta)
go_results <- run_go_enrichment(cluster_df, bg_genes)

if (nrow(go_results) > 0) {
  write.csv(go_results, file.path(dir_tta, "go_enrichment.csv"),
            row.names = FALSE)
}

# ============================================================
# 10. Sensitivity across normalization schemes
# ============================================================

cat("\nSensitivity: re-clustering under alternative normalizations...\n")

run_sensitivity <- function(te_scheme_label, te_mat_alt,
                             norm_trans_alt, norm_translat_alt,
                             sample_meta, tp_labels, all_genes,
                             orig_assignments, k) {
  n_tp <- length(tp_labels)
  n_g <- length(all_genes)
  M_alt <- matrix(NA, n_g, n_tp, dimnames = list(all_genes, tp_labels))
  TE_alt <- matrix(NA, n_g, n_tp, dimnames = list(all_genes, tp_labels))

  for (t_idx in seq_len(n_tp)) {
    tp <- tp_labels[t_idx]
    te_cols <- grep(paste0("^", tp, "_"), colnames(te_mat_alt), value = TRUE)
    trans_cols <- sample_meta$sample_id[
      sample_meta$omics_type == "Transcriptome" &
        sample_meta$time_point == tp]
    trans_cols <- trans_cols[trans_cols %in% colnames(norm_trans_alt)]

    valid_genes <- all_genes[all_genes %in% rownames(norm_trans_alt)]
    M_alt[valid_genes, t_idx] <- rowMeans(log2(
      norm_trans_alt[valid_genes, trans_cols, drop = FALSE] + 1))
    te_cols_valid <- te_cols[te_cols %in% colnames(te_mat_alt)]
    if (length(te_cols_valid) > 0) {
      valid_te <- valid_genes[valid_genes %in% rownames(te_mat_alt)]
      TE_alt[valid_te, t_idx] <- rowMeans(
        te_mat_alt[valid_te, te_cols_valid, drop = FALSE], na.rm = TRUE)
    }
  }

  # Same centering as primary (not z-scoring)
  M_alt_c <- M_alt - rowMeans(M_alt, na.rm = TRUE)
  TE_alt_c <- TE_alt - rowMeans(TE_alt, na.rm = TRUE)
  sd_M_alt <- sd(as.vector(M_alt_c), na.rm = TRUE)
  sd_TE_alt <- sd(as.vector(TE_alt_c), na.rm = TRUE)
  if (sd_M_alt == 0 || sd_TE_alt == 0) return(NA)
  M_alt_z <- M_alt_c / sd_M_alt
  TE_alt_z <- TE_alt_c / sd_TE_alt

  ok <- apply(M_alt_z, 1, function(x) all(is.finite(x))) &
    apply(TE_alt_z, 1, function(x) all(is.finite(x)))

  if (sum(ok) < k * 5) {
    cat(sprintf("  %s: too few valid genes (%d)\n", te_scheme_label, sum(ok)))
    return(NA)
  }

  phase_alt <- cbind(M_alt_z[ok, ], TE_alt_z[ok, ])
  d_alt <- dist(phase_alt, method = "euclidean")
  hc_alt <- hclust(d_alt, method = "ward.D2")
  cl_alt <- cutree(hc_alt, k = k)
  names(cl_alt) <- rownames(phase_alt)

  # ARI with original (named vector indexing — no positional match)
  common <- intersect(names(cl_alt), names(orig_assignments))
  if (length(common) < k * 5) return(NA)
  compute_ari(orig_assignments[common], cl_alt[common])
}

sensitivity_results <- data.frame()

te_schemes <- list()
if (exists("te_B")) te_schemes[["DESeq2"]] <- list(te = te_B, trans = norm_trans_B, translat = norm_translat_B)
if (exists("te_C")) te_schemes[["TMM"]] <- list(te = te_C, trans = norm_trans_C, translat = norm_translat_C)

for (scheme_name in names(te_schemes)) {
  scheme <- te_schemes[[scheme_name]]
  cat(sprintf("  Testing %s normalization...\n", scheme_name))

  ari_sens <- run_sensitivity(
    scheme_name, scheme$te, scheme$trans, scheme$translat,
    meta, tp_labels, gene_ids_tta, all_assignments, best_k
  )

  sensitivity_results <- rbind(sensitivity_results, data.frame(
    scheme = scheme_name, ari_vs_constgenes = ari_sens, stringsAsFactors = FALSE
  ))
  cat(sprintf("  %s vs constGenes ARI: %.3f\n", scheme_name, ari_sens))
}

if (nrow(sensitivity_results) > 0) {
  write.csv(sensitivity_results,
            file.path(dir_tta, "sensitivity_across_normalizations.csv"),
            row.names = FALSE)
}

# ============================================================
# 11. Control gene trajectory check
# ============================================================

cat("\nControl gene trajectory assignments:\n")
for (g in params$control_genes) {
  gene_row <- which(symbols_filt == g)
  if (length(gene_row) == 0) next
  gid <- gene_ids_filt[gene_row[1]]
  if (gid %in% cluster_df$gene_id) {
    info <- cluster_df[cluster_df$gene_id == gid, ]
    cat(sprintf("  %s → C%d (%s)\n", g, info$cluster, info$trajectory_type))
  } else {
    cat(sprintf("  %s → not in TTA gene set\n", g))
  }
}

# ============================================================
# 12. Comparison with Mfuzz clusters (IU-4)
# ============================================================

mfuzz_file <- file.path(DIR_RESULTS, "baseline", "mfuzz_clusters_transcriptome.csv")
if (file.exists(mfuzz_file)) {
  mfuzz_cl <- read.csv(mfuzz_file, stringsAsFactors = FALSE)

  common_mfuzz <- intersect(cluster_df$gene_id, mfuzz_cl$gene_id)
  if (length(common_mfuzz) > 0) {
    tta_cl_common <- cluster_df$cluster[match(common_mfuzz, cluster_df$gene_id)]
    mfuzz_cl_common <- mfuzz_cl$cluster[match(common_mfuzz, mfuzz_cl$gene_id)]

    tta_mfuzz_ari <- compute_ari(tta_cl_common, mfuzz_cl_common)

    cat(sprintf("\nTTA vs Mfuzz (transcriptome) ARI: %.3f (n=%d common genes)\n",
                tta_mfuzz_ari, length(common_mfuzz)))
  }
}

# ============================================================
# 13. Summary
# ============================================================

cat(sprintf("\n=== TTA Summary ===\n"))
cat(sprintf("Genes in phase space: %d\n", length(gene_ids_tta)))
cat(sprintf("Euclidean clustering: k=%d (silhouette=%.3f)\n", best_k, best_sil))
cat("Trajectory types:\n")
for (info in cluster_info) {
  cat(sprintf("  C%d %s: %d genes (%.1f%%)\n",
              info$cluster, info$type, info$n_genes,
              info$n_genes / length(gene_ids_tta) * 100))
}
cat(sprintf("Bootstrap ARI: %.3f [%.3f, %.3f] — %s\n",
            ari_mean, ari_ci[1], ari_ci[2],
            if (ari_pass) "PASS" else "FAIL"))
if (nrow(go_results) > 0) {
  cat(sprintf("GO enrichment: %d terms across %d trajectory types\n",
              nrow(go_results), length(unique(go_results$trajectory_type))))
}
cat(sprintf("Results saved to: %s\n", dir_tta))

# ============================================================
# 14. Session info
# ============================================================

writeLines(capture.output(sessionInfo()),
           file.path(dir_tta, "sessionInfo.txt"))
