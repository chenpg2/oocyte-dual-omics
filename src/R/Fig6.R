# Fig6.R — assemble Figure 6: CPE cis-regulatory basis
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
  dir_f6  <- file.path(DIR_RESULTS, "figures", "F6_cpe_mechanism")
  dir_sup <- file.path(DIR_RESULTS, "supplementary")
  dir_cis <- file.path(DIR_RESULTS, "cis_grammar")

  s4  <- read.csv(file.path(dir_sup, "S4_cpe_enrichment.csv"), stringsAsFactors = FALSE)
  s4f <- read.csv(file.path(dir_sup, "S4_functional_cpe_features.csv"), stringsAsFactors = FALSE)
  s8  <- read.csv(file.path(dir_sup, "S8_logistic_or.csv"), stringsAsFactors = FALSE)
  # F6E uses CPE+ gene GO enrichment (S7_go_comp_vs_clear.csv has 0 significant terms)
  s7  <- read.csv(file.path(dir_sup, "S7_go_cpe_positive.csv"), stringsAsFactors = FALSE)
  s1  <- read.csv(file.path(DIR_RESULTS, "supplementary", "S1_tbr_cs_matrix.csv"),
                   stringsAsFactors = FALSE)
  # Version-stripped id for merging with functional-CPE features (keyed on gene_id_clean)
  s1$gene_id_clean <- sub("\\.\\d+$", "", s1$gene_id)

  # F6A schematic (re-use the standalone version saved above)
  # Rebuild in-memory for patchwork
  p6a_pw <- ggplot() +
    annotate("rect", xmin=0,xmax=3.2,ymin=0.82,ymax=1.18,
             fill="#D9E8F5", colour=col_omics[["Transcriptome"]], linewidth=0.5) +
    annotate("text", x=1.6,y=1.0, label="CDS", size=2.8,
             colour=col_omics[["Transcriptome"]], fontface="bold", family="Helvetica") +
    annotate("rect", xmin=3.2,xmax=7.8,ymin=0.82,ymax=1.18,
             fill="#FFF0D9", colour="#E07B39", linewidth=0.5) +
    annotate("text", x=3.6,y=1.0, label="3'UTR", size=2.3,
             colour="#E07B39", fontface="bold", family="Helvetica", hjust=0.5) +
    annotate("rect", xmin=7.8,xmax=9.0,ymin=0.82,ymax=1.18,
             fill="#F0F0F0", colour="grey60", linewidth=0.5) +
    annotate("text", x=8.4,y=1.0, label="AAAA", size=2.3,
             colour="grey55", family="Helvetica") +
    annotate("rect", xmin=4.3,xmax=5.7,ymin=0.82,ymax=1.18,
             fill=col_omics[["Translatome"]], colour=col_omics[["Translatome"]],
             linewidth=0.7, alpha=0.88) +
    annotate("text", x=5.0,y=1.0, label="CPE\nUUUUAU", size=2.2,
             colour="white", fontface="bold", lineheight=0.9, family="Helvetica") +
    annotate("rect", xmin=6.9,xmax=7.75,ymin=0.82,ymax=1.18,
             fill="#2D4B8E", colour="#2D4B8E", linewidth=0.7, alpha=0.82) +
    annotate("text", x=7.32,y=1.0, label="PAS\nAATAAA", size=2.2,
             colour="white", fontface="bold", lineheight=0.9, family="Helvetica") +
    annotate("segment", x=5.72,xend=6.88,y=1.26,yend=1.26,
             colour="grey40", linewidth=0.4,
             arrow=arrow(ends="both",length=unit(0.07,"cm"),type="open")) +
    annotate("text", x=6.3,y=1.36, label="≤100 nt", size=2.1, colour="grey40") +
    annotate("rect", xmin=4.5,xmax=5.5,ymin=0.44,ymax=0.76,
             fill="#F0B323", colour="#F0B323", linewidth=0.4, alpha=0.88) +
    annotate("text", x=5.0,y=0.60, label="CPEB1", size=2.3,
             colour="white", fontface="bold", family="Helvetica") +
    annotate("segment", x=5.0,xend=5.0,y=0.78,yend=0.82,
             arrow=arrow(length=unit(0.09,"cm"),type="closed"),
             colour="#F0B323", linewidth=0.6) +
    annotate("text", x=5.0,y=0.24,
             label="GV: short poly(A) → low TE",
             size=2.2, colour=col_omics[["Transcriptome"]], family="Helvetica") +
    annotate("text", x=5.0,y=0.10,
             label="GVBD: long poly(A) → high TE",
             size=2.2, colour=col_omics[["Translatome"]], family="Helvetica") +
    coord_cartesian(xlim=c(-0.5,9.6), ylim=c(-0.05,1.55)) +
    labs(title="CPE mechanism") +
    theme_void(base_family="Helvetica") +
    theme(legend.position="none",
          plot.title=element_text(size=8,face="bold",hjust=0))

  # F6B — CPE+ vs CPE- TBR violin (merge on version-stripped gene id)
  df_b <- merge(s1, s4f, by = "gene_id_clean", all.x = TRUE)
  df_b <- df_b[!is.na(df_b$n_functional_cpe) & is.finite(df_b$tbr), ]
  df_b$cpe_group <- ifelse(df_b$n_functional_cpe > 0, "CPE+", "CPE-")
  df_b$cpe_group <- factor(df_b$cpe_group, levels = c("CPE-", "CPE+"))
  n_pos <- sum(df_b$cpe_group == "CPE+"); n_neg <- sum(df_b$cpe_group == "CPE-")
  med_b <- df_b %>% group_by(cpe_group) %>%
    summarise(med = median(tbr, na.rm = TRUE), .groups = "drop")

  p6b_pw <- ggplot(df_b, aes(x = cpe_group, y = pmax(-5, pmin(5, tbr)),
                               fill = cpe_group)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey55", linewidth = 0.4) +
    geom_violin(trim = TRUE, alpha = 0.55, colour = NA, scale = "width") +
    geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white",
                 colour = "grey30", linewidth = 0.5) +
    geom_point(data = med_b, aes(x = cpe_group, y = med),
               inherit.aes = FALSE, size = 2.2, shape = 18, colour = "grey20") +
    scale_fill_manual(values = c("CPE-" = "#A7A9AC",
                                  "CPE+" = unname(col_omics["Translatome"])),
                       guide = "none") +
    scale_x_discrete(labels = c("CPE-" = sprintf("CPE−\n(n=%s)", format(n_neg, big.mark = ",")),
                                "CPE+" = sprintf("CPE+\n(n=%s)", format(n_pos, big.mark = ",")))) +
    coord_cartesian(ylim = c(-3, 5)) +
    labs(x = NULL, y = "TBR (log₂)",
         title = "Functional CPE → higher TBR") +
    theme_oocyte(base_size = 8)

  # F6C — enrichment bar (S4 columns: group / pct_has_fcpe / fisher_OR)
  s4_sub <- s4[s4$group %in% c("Compensatory", "Coordinated Clearance"), ]
  s4_sub$group <- factor(s4_sub$group,
                         levels = c("Compensatory", "Coordinated Clearance"))
  or_val <- s4$fisher_OR[s4$group == "Compensatory"]
  p6c_pw <- ggplot(s4_sub, aes(x = group, y = pct_has_fcpe, fill = group)) +
    geom_col(width = 0.6, alpha = 0.88) +
    geom_text(aes(label = sprintf("%.1f%%", pct_has_fcpe)),
              vjust = -0.4, size = 2.3, fontface = "bold", colour = "grey20") +
    annotate("text", x = 1.5, y = max(s4_sub$pct_has_fcpe) * 1.28,
             label = sprintf("OR = %.2f", or_val), size = 2.3,
             fontface = "bold", colour = "grey30") +
    scale_fill_manual(
      values = c("Compensatory" = "#00857C",
                 "Coordinated Clearance" = unname(col_omics["Transcriptome"])),
      guide = "none") +
    scale_x_discrete(labels = c("Compensatory" = "Compen.",
                                "Coordinated Clearance" = "Clearance")) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
    labs(x = NULL, y = "% with functional CPE",
         title = "CPE enrichment\nby quadrant") +
    theme_oocyte(base_size = 8)

  # F6D — forest plot (S8 columns: model / predictor / OR / OR_ci_low / OR_ci_high)
  s8_sub <- s8[s8$predictor %in% c("has_cpe_pas_pair", "n_functional_cpe"), ]
  if (nrow(s8_sub) == 0) s8_sub <- head(s8, 4)
  s8_sub$model_short <- sub(":.*", "", s8_sub$model)          # "M0".."M3"
  s8_sub$model_short <- factor(s8_sub$model_short, levels = rev(s8_sub$model_short))
  p6d_pw <- ggplot(s8_sub, aes(x = OR, y = model_short)) +
    geom_vline(xintercept = 1, linetype = "dashed", colour = "grey55", linewidth = 0.4) +
    geom_errorbarh(aes(xmin = OR_ci_low, xmax = OR_ci_high),
                   height = 0.22, colour = "grey30", linewidth = 0.6) +
    geom_point(size = 2.8, colour = col_omics[["Translatome"]], shape = 18) +
    geom_text(aes(label = sprintf("%.2f", OR)), vjust = -1.05,
              size = 2.0, colour = "grey25") +
    scale_x_continuous(expand = expansion(mult = c(0.05, 0.12))) +
    labs(x = "Odds Ratio (95% CI)", y = NULL,
         title = "CPE effect after\nconfounder adjustment") +
    theme_oocyte(base_size = 8)

  # F6E — GO bubble (CPE+ genes, biological process; top 10 by adj. p)
  if (nrow(s7) > 0 && "Description" %in% names(s7)) {
    s7_top <- head(s7[order(s7$p.adjust), ], 10)
    # wrap long GO term names so they fit the narrow panel
    s7_top$Description_wrap <- vapply(
      s7_top$Description,
      function(x) paste(strwrap(x, width = 26), collapse = "\n"),
      character(1))
    s7_top$Description_wrap <- factor(s7_top$Description_wrap,
                                       levels = rev(s7_top$Description_wrap))
    p6e_pw <- ggplot(s7_top,
                      aes(x = -log10(p.adjust),
                          y = Description_wrap,
                          size = Count,
                          colour = -log10(p.adjust))) +
      geom_point(alpha = 0.85) +
      scale_colour_gradient(low = "#9BBBD4", high = col_omics[["Transcriptome"]],
                             guide = "none") +
      scale_size_continuous(range = c(2, 6), name = "Gene\ncount") +
      scale_x_continuous(expand = expansion(mult = c(0.06, 0.10))) +
      labs(x = "−log₁₀(adj. p)", y = NULL,
           title = "GO enrichment (CPE+ genes)") +
      theme_oocyte(base_size = 8) +
      theme(axis.text.y = element_text(size = 5.4, lineheight = 0.82),
            legend.key.size = unit(3, "mm"),
            legend.title = element_text(size = 6),
            legend.text = element_text(size = 5.5))
  } else {
    p6e_pw <- ggplot() +
      annotate("text",x=0.5,y=0.5,label="F6E: see F6E_go_bubble.pdf",size=3) +
      theme_void()
  }

  # F6F — dose response (merge on version-stripped gene id)
  df_f6 <- merge(s1, s4f, by = "gene_id_clean", all.x = TRUE)
  df_f6 <- df_f6[!is.na(df_f6$n_functional_cpe) & is.finite(df_f6$tbr), ]
  df_f6$cpe_g <- cut(df_f6$n_functional_cpe, breaks=c(-Inf,0,1,2,Inf),
                      labels=c("0","1","2","≥3"))
  dose6 <- df_f6 %>% group_by(cpe_g) %>%
    summarise(mean_tbr=mean(tbr,na.rm=TRUE),
              se_tbr=sd(tbr,na.rm=TRUE)/sqrt(n()),
              n=n(), .groups="drop")
  p6f_pw <- ggplot(dose6, aes(x=cpe_g,y=mean_tbr,group=1)) +
    geom_hline(yintercept=0,linetype="dashed",colour="grey55",linewidth=0.4) +
    geom_ribbon(aes(ymin=mean_tbr-se_tbr,ymax=mean_tbr+se_tbr),
                fill=col_omics[["Translatome"]],alpha=0.18,colour=NA) +
    geom_line(colour=col_omics[["Translatome"]],linewidth=1.0) +
    geom_point(size=2.5,colour=col_omics[["Translatome"]],
               shape=21,fill="white",stroke=1.2) +
    labs(x="Functional CPEs",y="Mean TBR (± SE)",
         title="CPE–TBR\ndose response") +
    theme_oocyte(base_size=8)

  row1 <- (p6a_pw | p6b_pw | p6c_pw) + plot_layout(widths = c(1.7, 1, 1))
  row2 <- (p6d_pw | p6e_pw | p6f_pw) + plot_layout(widths = c(1, 1.7, 1))
  p6_full <- (row1 / row2) +
    plot_layout(heights = c(0.9, 1)) +
    plot_annotation(
      tag_levels = list(c("A","B","C","D","E","F")),
      theme = theme(plot.tag = element_text(size=10,face="bold",family="Helvetica"))
    ) &
    theme(plot.title = element_text(size = 8, face = "bold"),
          plot.subtitle = element_text(size = 6.2, colour = "grey40"),
          plot.margin = margin(3, 4, 3, 4))

  save_full(p6_full, file.path(dir_f6, "F6_main_assembly"), w=183, h=155)
  cat("F6 main assembly (A–F) saved\n")
}
