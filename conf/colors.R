# Pantone Color Palette — Oocyte Dual-Omics Project
# All visualization scripts source this file. Never hardcode colors elsewhere.

# --- Timepoint / Meiotic Stage Colors ---
# Sequential palette: cool → warm progression through maturation
col_timepoint <- c(
  "GV"   = "#2D4B8E",  # Pantone 2758 C (deep blue) — quiescent GV
  "GVBD" = "#00857C",  # Pantone 3285 C (teal) — GVBD transition
  "MI-6" = "#F0B323",  # Pantone 7406 C (golden) — early MI
  "MI-9" = "#E07E3C",  # Pantone 7413 C (amber) — late MI / hidden window
  "MII"  = "#B42F37"   # Pantone 7621 C (crimson) — MII completion
)

col_timepoint_label <- c(
  "0h\nGV", "3h\nGVBD", "6h\nMI", "9h\nMI", "12h\nMII"
)

# --- Trajectory Type Colors (TTA) ---
col_trajectory <- c(
  "Translate-to-Degrade"   = "#B42F37",  # Pantone 7621 C (crimson)
  "Compensatory Buffering"  = "#2D4B8E",  # Pantone 2758 C (deep blue)
  "Coordinated Clearance"   = "#6E6E6E",  # Pantone 424 C (neutral gray)
  "Selective Activation"    = "#00857C",  # Pantone 3285 C (teal)
  "Dormant Storage"         = "#A7A9AC",  # Pantone 422 C (light gray)
  "Transient Utilization"   = "#F0B323"   # Pantone 7406 C (golden)
)

# --- Omics Type Colors ---
col_omics <- c(
  "Transcriptome" = "#2D4B8E",  # Pantone 2758 C
  "Translatome"   = "#B42F37"   # Pantone 7621 C
)

# --- TPI Direction Colors ---
col_tpi <- c(
  "TE-leading"    = "#B42F37",  # Pantone 7621 C — TE changes first
  "mRNA-leading"  = "#2D4B8E",  # Pantone 2758 C — mRNA changes first
  "Simultaneous"  = "#6E6E6E"   # Pantone 424 C
)

# --- Heatmap Palette ---
col_heatmap_diverging <- c(
  low  = "#2D4B8E",  # Pantone 2758 C
  mid  = "#FFFFFF",
  high = "#B42F37"   # Pantone 7621 C
)

col_heatmap_sequential <- c(
  low  = "#F7F7F7",
  high = "#2D4B8E"   # Pantone 2758 C
)

# --- Significance / Annotation ---
col_sig <- c(
  "significant"     = "#B42F37",  # Pantone 7621 C
  "not significant" = "#A7A9AC"   # Pantone 422 C
)

# --- ggplot2 Theme ---
theme_oocyte <- function(base_size = 10) {
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      text = ggplot2::element_text(family = "Helvetica"),
      plot.title = ggplot2::element_text(size = base_size + 2, face = "bold"),
      axis.title = ggplot2::element_text(size = base_size),
      axis.text = ggplot2::element_text(size = base_size - 1, color = "black"),
      legend.title = ggplot2::element_text(size = base_size),
      legend.text = ggplot2::element_text(size = base_size - 1),
      strip.text = ggplot2::element_text(size = base_size, face = "bold"),
      strip.background = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank()
    )
}
