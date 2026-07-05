# 00_config.R — Project-wide configuration
# Source this file at the top of every analysis script.

library(yaml)
library(readxl)
library(dplyr)
library(tibble)

# --- Load parameters ---
params <- yaml::read_yaml("conf/params.yaml")
set.seed(params$seed)

# --- Load color palette ---
source("conf/colors.R")

# --- Paths ---
DIR_RESULTS <- params$paths$results
dir.create(DIR_RESULTS, recursive = TRUE, showWarnings = FALSE)

# --- Metadata ---
load_metadata <- function() {
  meta <- read.csv(params$paths$metadata, stringsAsFactors = FALSE)
  meta$time_point <- factor(meta$time_point,
                            levels = c("0h", "3h", "6h", "9h", "12h"))
  meta$stage <- factor(meta$stage,
                       levels = c("GV", "GVBD", "MI-6", "MI-9", "MII"))
  meta$omics_type <- factor(meta$omics_type,
                            levels = c("Transcriptome", "Translatome"))
  meta
}

# --- Load raw counts ---
load_counts <- function(path) {
  raw <- read.csv(path, check.names = FALSE, fileEncoding = "UTF-8-BOM")
  gene_info <- raw[, 1:6]
  counts <- raw[, 7:ncol(raw)]

  gene_ids_full <- raw$Geneid
  gene_ids_base <- sub("\\.\\d+$", "", gene_ids_full)

  rownames(counts) <- gene_ids_full
  list(
    counts = as.matrix(counts),
    gene_info = gene_info,
    gene_ids_full = gene_ids_full,
    gene_ids_base = gene_ids_base
  )
}

# --- Load 147 constGenes ---
load_constgenes <- function() {
  tbl <- readxl::read_xlsx(params$paths$constgenes, skip = 3)
  gene_col <- tbl[[1]]
  gene_ids <- gene_col[grepl("^ENSMUSG", gene_col)]
  unique(gene_ids)
}

# --- Pair mapping ---
build_pair_map <- function(meta) {
  trans <- meta[meta$omics_type == "Transcriptome", ]
  pairs <- trans[!is.na(trans$paired) & trans$paired != "NA",
                 c("sample_id", "paired", "time_point", "stage", "replicate")]
  colnames(pairs) <- c("transcriptome_id", "translatome_id",
                        "time_point", "stage", "replicate")
  pairs
}

# --- Gene symbol mapping ---
strip_ensembl_version <- function(ids) {
  sub("\\.\\d+$", "", ids)
}
