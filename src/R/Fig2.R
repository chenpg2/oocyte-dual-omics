# Fig2_tpi.R — Figure 2: Translational Temporal Precedence Index (TPI)
#
# Core claim: Translational changes systematically precede mRNA changes
#             genome-wide during oocyte maturation.
#
# Panels:
#   A — TPI concept schematic (4 intervals, first-change arrows for TE vs mRNA)
#   B — Genome-wide TPI distribution histogram (hero panel, n=8,818)
#   C — TPI category proportions (ring chart)
#   D — TPI by TTA trajectory type (violin + median + bootstrap CI)
#   E — Known gene TPI lollipop chart
#   F — Normalization sensitivity (concordance %)
#
# Input:  results/tpi/tpi_scores.csv
#         results/tpi/tpi_by_trajectory_type.csv
#         results/tpi/known_gene_validation.csv
#         results/tpi/sensitivity_across_normalizations.csv
# Output: results/figures/F2_tpi/

source("src/R/00_config.R")

library(ggplot2)
library(patchwork)
library(dplyr)
library(scales)

dir_fig <- file.path(DIR_RESULTS, "figures", "F2_tpi")
dir.create(dir_fig, recursive = TRUE, showWarnings = FALSE)

# ---- export helper -------------------------------------------------------

save_pub <- function(p, name, w = 88, h = 70) {
  grDevices::cairo_pdf(file.path(dir_fig, paste0(name, ".pdf")),
                       width = w / 25.4, height = h / 25.4, family = "Helvetica")
  print(p); dev.off()
  svglite::svglite(file.path(dir_fig, paste0(name, ".svg")),
                   width = w / 25.4, height = h / 25.4)
  print(p); dev.off()
}

# ============================================================
# 0. Load data
# ============================================================

tpi_df      <- read.csv(file.path(DIR_RESULTS, "tpi", "tpi_scores.csv"),
                         stringsAsFactors = FALSE)
traj_df     <- read.csv(file.path(DIR_RESULTS, "tpi", "tpi_by_trajectory_type.csv"),
                         stringsAsFactors = FALSE)
ctrl_df     <- read.csv(file.path(DIR_RESULTS, "tpi", "known_gene_validation.csv"),
                         stringsAsFactors = FALSE)
sens_df     <- read.csv(file.path(DIR_RESULTS, "tpi", "sensitivity_across_normalizations.csv"),
                         stringsAsFactors = FALSE)

# Validated TPI genes
tpi_valid <- tpi_df[!is.na(tpi_df$tpi), ]
n_total   <- nrow(tpi_valid)

# Category counts
n_te   <- sum(tpi_valid$tpi < 0)
n_sim  <- sum(tpi_valid$tpi == 0)
n_mrna <- sum(tpi_valid$tpi > 0)

cat("TPI valid genes:", n_total, "\n")
cat("  TE-leading:", n_te, "(", round(100*n_te/n_total, 1), "%)\n")
cat("  Simultaneous:", n_sim, "(", round(100*n_sim/n_total, 1), "%)\n")
cat("  mRNA-leading:", n_mrna, "(", round(100*n_mrna/n_total, 1), "%)\n")

# TTA order for plots
traj_order <- c("TE-Only Activation", "Late Compensatory Buffering",
                 "Deep Coordinated Clearance", "Coordinated Clearance",
                 "Mild Coordinated Clearance")
traj_short <- c("TE-Only\nActivation", "Late\nCompensatory",
                 "Deep Coord.\nClearance", "Coord.\nClearance",
                 "Mild Coord.\nClearance")

# ============================================================
# Panel A — TPI Concept Schematic
# ============================================================
# Shows 4 meiotic intervals; TE changes in interval 1 (GV→GVBD),
# mRNA changes in interval 2 (GVBD→MI6) → TPI = 1-2 = -1

interval_labs <- c("GV→GVBD", "GVBD→MI6", "MI6→MI9", "MI9→MII")

schema_df <- data.frame(
  molecule = c("mRNA", "mRNA", "TE", "TE"),
  xstart   = c(2, 1, 1, 1),
  xend     = c(4, 2, 2, 4),
  y        = c(2, 2, 1, 1),
  type     = c("flat", "change", "change", "flat")
)

# Segment data for trajectories
seg_mrna <- data.frame(
  x = c(1,  2,  3,  4),
  y = c(0.1, 0.1, -0.4, -0.8),
  omics = "mRNA"
)
seg_te <- data.frame(
  x = c(1,  2,  3,  4),
  y = c(0.1, 0.5, 0.7, 0.6),
  omics = "TE"
)

schematic_pts <- rbind(seg_mrna, seg_te)
schematic_pts$omics <- factor(schematic_pts$omics, levels = c("mRNA", "TE"))

arrow_df <- data.frame(
  x1 = c(1.5, 2.5), y1 = c(0.48, -0.2),
  x2 = c(1.5, 2.5), y2 = c(0.12, -0.38),
  label = c("TE first\n(interval 1)", "mRNA first\n(interval 2)")
)

p_schema_bg <- ggplot() +
  # shaded intervals
  annotate("rect", xmin = 1, xmax = 2, ymin = -1.1, ymax = 1.2,
           fill = "#F5F5F5", alpha = 0.7) +
  annotate("rect", xmin = 2, xmax = 3, ymin = -1.1, ymax = 1.2,
           fill = "#EBEBEB", alpha = 0.6) +
  # trajectory lines
  geom_line(data = schematic_pts,
            aes(x = x, y = y, colour = omics, linetype = omics),
            linewidth = 0.9) +
  geom_point(data = schematic_pts,
             aes(x = x, y = y, colour = omics), size = 2.2) +
  # first-change arrows
  annotate("segment",
           x = 1.5, xend = 1.5, y = 0.47, yend = 0.13,
           arrow = arrow(length = unit(0.12, "cm"), type = "closed"),
           colour = col_omics["Translatome"], linewidth = 0.6) +
  annotate("segment",
           x = 2.5, xend = 2.5, y = -0.19, yend = -0.37,
           arrow = arrow(length = unit(0.12, "cm"), type = "closed"),
           colour = col_omics["Transcriptome"], linewidth = 0.6) +
  annotate("text", x = 1.6, y = 0.30, label = "TE first\n(int. 1)",
           hjust = 0, size = 2.6, colour = col_omics["Translatome"],
           family = "Helvetica", lineheight = 0.9) +
  annotate("text", x = 2.6, y = -0.28, label = "mRNA first\n(int. 2)",
           hjust = 0, size = 2.6, colour = col_omics["Transcriptome"],
           family = "Helvetica", lineheight = 0.9) +
  # TPI label
  annotate("text", x = 2.5, y = -0.82,
           label = "TPI = t₁(TE) − t₁(mRNA) = 1 − 2 = −1",
           hjust = 0.5, size = 2.3, fontface = "italic",
           colour = "grey30", family = "Helvetica") +
  # x axis labels
  scale_x_continuous(breaks = 1:4, labels = interval_labs,
                     expand = c(0.1, 0.1)) +
  scale_colour_manual(
    values = c("mRNA" = unname(col_omics["Transcriptome"]),
               "TE"   = unname(col_omics["Translatome"])),
    labels = c("mRNA (transcriptome)", "TE (translatome)")
  ) +
  scale_linetype_manual(
    values = c("mRNA" = "dashed", "TE" = "solid"),
    labels = c("mRNA (transcriptome)", "TE (translatome)")
  ) +
  labs(x = "Meiotic interval", y = "Relative expression change",
       title = "TPI concept",
       colour = NULL, linetype = NULL) +
  theme_oocyte(base_size = 8) +
  theme(legend.position = "inside",
        legend.position.inside = c(0.70, 0.92),
        legend.key.width = unit(0.7, "cm"),
        legend.text = element_text(size = 6),
        legend.background = element_blank(),
        axis.text.x = element_text(size = 6, angle = 18, hjust = 1),
        plot.title = element_text(size = 8, face = "bold"))

save_pub(p_schema_bg, "F2A_tpi_schematic", w = 62, h = 60)
cat("F2A saved.\n")

# ============================================================
# Panel B — Genome-wide TPI distribution (HERO panel)
# ============================================================

tpi_valid$tpi_clamp <- pmax(-3, pmin(3, tpi_valid$tpi))
tpi_valid$direction <- ifelse(tpi_valid$tpi < 0, "TE-leading",
                               ifelse(tpi_valid$tpi > 0, "mRNA-leading",
                                      "Simultaneous"))
tpi_valid$direction <- factor(tpi_valid$direction,
                               levels = c("TE-leading", "Simultaneous", "mRNA-leading"))

# Per-bin counts for annotation
bin_counts <- tpi_valid %>%
  mutate(tpi_bin = tpi_clamp) %>%
  count(tpi_bin)

pct_te   <- round(100 * n_te   / n_total, 1)
pct_sim  <- round(100 * n_sim  / n_total, 1)
pct_mrna <- round(100 * n_mrna / n_total, 1)

p_hist <- ggplot(tpi_valid, aes(x = tpi_clamp, fill = direction)) +
  geom_histogram(binwidth = 1, colour = "white", linewidth = 0.3, alpha = 0.88) +
  geom_vline(xintercept = 0, linetype = "dashed",
             colour = col_tpi["Simultaneous"], linewidth = 0.5) +
  # direction labels in the free top corners, above the short outer bars
  annotate("text", x = -3.4, y = Inf, vjust = 1.4, hjust = 0,
           label = paste0("TE-leading\n", format(n_te, big.mark = ","),
                          " (", pct_te, "%)"),
           size = 2.7, colour = col_tpi["TE-leading"],
           family = "Helvetica", fontface = "bold", lineheight = 0.95) +
  annotate("text", x = 3.4, y = Inf, vjust = 1.4, hjust = 1,
           label = paste0("mRNA-leading\n", format(n_mrna, big.mark = ","),
                          " (", pct_mrna, "%)"),
           size = 2.7, colour = col_tpi["mRNA-leading"],
           family = "Helvetica", fontface = "bold", lineheight = 0.95) +
  annotate("text", x = 0, y = Inf, vjust = 1.4, hjust = 0.5,
           label = paste0("Simultaneous ", pct_sim, "%"),
           size = 2.4, colour = "grey45",
           family = "Helvetica", lineheight = 0.95) +
  scale_x_continuous(
    breaks = -3:3,
    labels = c("-3\n(TE leads\nearlier)", "-2", "-1",
               "0\n(Sim.)",
               "1", "2", "3\n(mRNA leads\nearlier)"),
    expand = c(0.02, 0.02)
  ) +
  scale_y_continuous(labels = scales::comma, expand = expansion(mult = c(0, 0.20))) +
  scale_fill_manual(values = col_tpi, guide = "none") +
  labs(x = "TPI  (interval of first ΔTE − interval of first ΔmRNA)",
       y = "Number of genes",
       title = "Genome-wide Temporal Precedence Index",
       subtitle = paste0("n = ", format(n_total, big.mark = ","),
                          " genes · |log₂FC| ≥ 0.5 · DESeq2 LRT pre-filtered"),
       fill = NULL) +
  theme_oocyte(base_size = 8.5) +
  theme(plot.subtitle = element_text(size = 7, colour = "grey40"))

save_pub(p_hist, "F2B_tpi_histogram", w = 118, h = 70)
cat("F2B saved.\n")

# ============================================================
# Panel C — TPI category ring chart
# ============================================================

ring_df <- data.frame(
  category = factor(c("TE-leading", "Simultaneous", "mRNA-leading"),
                    levels = c("TE-leading", "Simultaneous", "mRNA-leading")),
  n = c(n_te, n_sim, n_mrna)
)
ring_df$pct <- round(100 * ring_df$n / sum(ring_df$n), 1)
ring_short <- c("TE-leading" = "TE-lead.", "Simultaneous" = "Sim.",
                "mRNA-leading" = "mRNA-lead.")

# Compute cumulative positions for labels
ring_df <- ring_df %>%
  arrange(desc(category)) %>%
  mutate(
    ymax = cumsum(n),
    ymin = ymax - n,
    label_pos = (ymin + ymax) / 2,
    label = paste0(ring_short[as.character(category)], "\n", pct, "%")
  )

p_ring <- ggplot(ring_df, aes(ymax = ymax, ymin = ymin,
                               xmax = 4, xmin = 2.4,
                               fill = category)) +
  geom_rect(colour = "white", linewidth = 0.3) +
  geom_text(aes(x = 5.1, y = label_pos, label = label, colour = category),
            size = 2.2, lineheight = 0.9, family = "Helvetica", fontface = "bold") +
  annotate("text", x = 0, y = 0,
           label = paste0(round(n_total / 1000, 2), "k\ngenes"),
           size = 2.6, hjust = 0.5, vjust = 0.5,
           colour = "grey30", family = "Helvetica") +
  coord_polar(theta = "y") +
  scale_fill_manual(values = col_tpi) +
  scale_colour_manual(values = col_tpi) +
  xlim(0, 6.5) +
  theme_void() +
  theme(legend.position = "none",
        plot.title = element_text(size = 8, face = "bold",
                                  family = "Helvetica", hjust = 0.5)) +
  labs(title = "TPI categories")

save_pub(p_ring, "F2C_tpi_ring", w = 52, h = 52)
cat("F2C saved.\n")

# ============================================================
# Panel D — TPI by TTA trajectory type
# ============================================================

tpi_traj <- tpi_valid[!is.na(tpi_valid$trajectory_type) &
                        tpi_valid$trajectory_type %in% traj_order, ]
tpi_traj$trajectory_type <- factor(tpi_traj$trajectory_type, levels = traj_order)

# Summary for overlay
traj_summary <- traj_df[traj_df$trajectory_type %in% traj_order, ]
traj_summary$trajectory_type <- factor(traj_summary$trajectory_type, levels = traj_order)

# Significance asterisks
traj_summary$sig_label <- ifelse(
  is.na(traj_summary$wilcox_padj), "",
  ifelse(traj_summary$wilcox_padj < 0.001, "***",
         ifelse(traj_summary$wilcox_padj < 0.01, "**",
                ifelse(traj_summary$wilcox_padj < 0.05, "*", "ns")))
)

# Short labels for x-axis
traj_label_map <- setNames(traj_short, traj_order)

p_violin <- ggplot(tpi_traj, aes(x = trajectory_type, y = tpi_clamp,
                                   fill = trajectory_type)) +
  geom_violin(trim = TRUE, scale = "width", alpha = 0.65,
              colour = NA, bw = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed",
             colour = col_tpi["Simultaneous"], linewidth = 0.4) +
  # bootstrap CI bar
  geom_errorbar(data = traj_summary,
                aes(x = trajectory_type, ymin = boot_ci_lo, ymax = boot_ci_hi),
                colour = "grey20", linewidth = 0.55, width = 0.15,
                inherit.aes = FALSE) +
  # median dot
  geom_point(data = traj_summary,
             aes(x = trajectory_type, y = median_tpi),
             colour = "white", fill = "grey20",
             shape = 21, size = 2.2, stroke = 0.5,
             inherit.aes = FALSE) +
  # significance
  geom_text(data = traj_summary,
            aes(x = trajectory_type, y = 3.3, label = sig_label),
            size = 3.5, family = "Helvetica", inherit.aes = FALSE) +
  scale_x_discrete(labels = traj_label_map) +
  scale_y_continuous(breaks = -3:3, limits = c(-3.8, 3.8)) +
  scale_fill_manual(values = c(
    "TE-Only Activation"       = col_omics["Translatome"],
    "Late Compensatory Buffering" = "#5B8DB8",
    "Deep Coordinated Clearance"  = "#6E6E6E",
    "Coordinated Clearance"       = "#9A9A9A",
    "Mild Coordinated Clearance"  = "#BCBCBC"
  )) +
  labs(x = NULL, y = "TPI",
       title = "TPI by trajectory type",
       subtitle = "Dot = median · bar = 95% bootstrap CI") +
  theme_oocyte(base_size = 8) +
  theme(axis.text.x = element_text(size = 6, lineheight = 0.85),
        legend.position = "none",
        plot.subtitle = element_text(size = 6.5, colour = "grey40"))

save_pub(p_violin, "F2D_tpi_by_trajectory", w = 88, h = 65)
cat("F2D saved.\n")

# ============================================================
# Panel E — Known gene TPI lollipop
# ============================================================

ctrl_plot <- ctrl_df[!is.na(ctrl_df$tpi), ]
ctrl_plot <- ctrl_plot[order(ctrl_plot$tpi), ]

# Choose informative genes: keep all with TPI data
ctrl_plot$gene_symbol <- factor(ctrl_plot$gene_symbol,
                                 levels = ctrl_plot$gene_symbol)
ctrl_plot$direction <- ifelse(ctrl_plot$tpi < 0, "TE-leading",
                               ifelse(ctrl_plot$tpi > 0, "mRNA-leading",
                                      "Simultaneous"))
ctrl_plot$direction <- factor(ctrl_plot$direction,
                               levels = c("TE-leading", "Simultaneous", "mRNA-leading"))
# short numeric labels placed just OUTSIDE the dot, away from the gene-name axis
ctrl_plot$tpi_label <- sprintf("%+d", ctrl_plot$tpi)

p_ctrl <- ggplot(ctrl_plot, aes(x = gene_symbol, y = tpi, colour = direction)) +
  geom_hline(yintercept = 0, linetype = "dashed",
             colour = col_tpi["Simultaneous"], linewidth = 0.4) +
  geom_segment(aes(xend = gene_symbol, y = 0, yend = tpi),
               linewidth = 0.7) +
  geom_point(size = 3.0, alpha = 0.9) +
  geom_text(aes(label = tpi_label,
                hjust = ifelse(tpi < 0, 1.5, -0.5)),
            size = 2.4, family = "Helvetica", colour = "grey30") +
  scale_y_continuous(breaks = -3:3, limits = c(-4.2, 4.2)) +
  scale_colour_manual(values = col_tpi) +
  coord_flip() +
  labs(x = NULL, y = "TPI",
       title = "Validated marker genes",
       colour = NULL) +
  theme_oocyte(base_size = 8) +
  theme(legend.position = "none",
        axis.text.y = element_text(face = "italic", size = 7.5))

save_pub(p_ctrl, "F2E_known_genes", w = 62, h = 55)
cat("F2E saved.\n")

# ============================================================
# Panel F — Normalization sensitivity
# ============================================================

# Add constGenes reference row
const_row <- data.frame(
  normalization = "constGenes\n(primary)",
  n_common      = n_total,
  spearman_rho  = 1.0,
  concordance_pct = 100.0,
  median_tpi    = -1
)

sens_df$normalization <- c("DESeq2", "TMM", "RUVg")
sens_plot <- rbind(const_row, sens_df)
sens_plot$normalization <- factor(
  sens_plot$normalization,
  levels = c("constGenes\n(primary)", "DESeq2", "TMM", "RUVg")
)
sens_plot$is_primary <- sens_plot$normalization == "constGenes\n(primary)"

p_sens <- ggplot(sens_plot, aes(x = normalization, y = concordance_pct,
                                 fill = is_primary)) +
  geom_col(width = 0.65, alpha = 0.88) +
  geom_text(aes(label = paste0(concordance_pct, "%"),
                y = concordance_pct + 1.5),
            size = 2.8, family = "Helvetica", vjust = 0) +
  geom_text(aes(label = paste0("ρ = ", ifelse(is_primary, "—", spearman_rho))),
            y = 5, size = 2.5, colour = "white", family = "Helvetica") +
  scale_fill_manual(values = c("TRUE" = col_omics["Translatome"],
                                "FALSE" = "#A7A9AC")) +
  scale_y_continuous(limits = c(0, 110), expand = c(0, 0)) +
  labs(x = "Normalization", y = "TPI direction concordance (%)",
       title = "Normalization robustness") +
  theme_oocyte(base_size = 8) +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 7.5, lineheight = 0.9))

save_pub(p_sens, "F2F_norm_sensitivity", w = 52, h = 52)
cat("F2F saved.\n")

# ============================================================
# Full Figure Assembly — F2 main (183 × 150 mm)
# ============================================================

# Row 1: A (schematic 62mm) + B (histogram 121mm)
# Row 2: C (ring 52mm) + D (violin 80mm) + E (known genes 51mm)
# — total width = 183mm, heights 70 / 65 mm

row_top <- (p_schema_bg | p_hist) + plot_layout(widths = c(1, 1.9))
row_bot <- (p_ring | p_violin | p_ctrl) + plot_layout(widths = c(0.85, 1.5, 0.95))

p_assembled <- (row_top / row_bot) +
  plot_layout(heights = c(1.05, 1)) +
  plot_annotation(
    tag_levels = list(c("A", "B", "C", "D", "E")),
    theme = theme(plot.tag = element_text(size = 10, face = "bold",
                                          family = "Helvetica"))
  ) &
  theme(plot.title = element_text(size = 8, face = "bold"),
        plot.subtitle = element_text(size = 6.3, colour = "grey40"),
        plot.margin = margin(3, 5, 3, 4))

# --- save main assembly ---
grDevices::cairo_pdf(
  file.path(dir_fig, "F2_main_assembly.pdf"),
  width  = 183 / 25.4,
  height = 150 / 25.4,
  family = "Helvetica"
)
print(p_assembled); dev.off()

svglite::svglite(
  file.path(dir_fig, "F2_main_assembly.svg"),
  width  = 183 / 25.4,
  height = 150 / 25.4
)
print(p_assembled); dev.off()

ragg::agg_tiff(
  file.path(dir_fig, "F2_main_assembly.tiff"),
  width  = 183 / 25.4, height = 150 / 25.4,
  units  = "in", res = 600
)
print(p_assembled); dev.off()

cat("\nFigure 2 assembly saved.\n")
cat("Output:", dir_fig, "\n")
