# Oocyte dual-omics: translational engagement precedes maternal mRNA clearance

Analysis and visualization code for a time-resolved transcriptome and translatome
(T&T-seq) study of mouse oocyte meiotic maturation (GV to MII). The project asks whether
changes in translational engagement precede or follow changes in transcript abundance at
each meiotic stage, and whether a transcript's fate is encoded in its 3'UTR sequence.

Private repository. Every script lives in `src/R/`. Configuration is in `conf/`
(`params.yaml`, colour palette, ggplot theme); `metadata_hcg.csv` is the sample sheet.

## Analysis pipeline

Run in order; each step reads outputs written by earlier steps into `results/`.

| Script | Analysis |
|--------|----------|
| `src/R/00_config.R` | Shared paths, parameters, and data loaders (sourced by every script) |
| `src/R/01_qc.R` | Quality control, pairing, and gene filtering |
| `src/R/02_normalize.R` | constGenes normalization and cross-scheme comparison (DESeq2, TMM, RUVg) |
| `src/R/03_baseline.R` | Time-dependent gene detection and endpoint (GV to MII) differential engagement |
| `src/R/04_tta.R` | mRNA-TE trajectory analysis; five-program phase-space clustering |
| `src/R/05_tpi.R` | Temporal precedence index (TPI) |
| `src/R/06_trd.R` | Translational rate dynamics (vTE), biphasic oscillation, compensation waves |
| `src/R/07_cis_grammar.R` | 3'UTR cis-features and XGBoost trajectory prediction |

## Supplementary analyses

| Script | Analysis |
|--------|----------|
| `src/R/S1_continuous_compensation.R` | Clearance severity (CS) and translational balance ratio (TBR) matrix |
| `src/R/S2_tpi_null_model.R` | TPI significance against a permutation null |
| `src/R/S3_gv_predictive_model.R` | Anti-circular GV-state prediction (models M1-M3 vs the full model) |
| `src/R/S4_functional_cpe.R` | Functional CPE annotation and enrichment |
| `src/R/S5_normalization_robustness.R` | Robustness of each conclusion across normalization schemes |
| `src/R/S6_tbr_threshold_sensitivity.R` | CS/TBR threshold sensitivity grid |
| `src/R/S7_functional_cpe_go.R` | Gene Ontology enrichment of CPE-positive genes |
| `src/R/S8_cpe_multivariate.R` | Multivariate logistic adjustment of the CPE effect |

## Figures

One script per figure. Each reads the tables written by the analysis and supplementary
steps and writes an editable SVG/PDF plus a 600 dpi TIFF to `results/figures/`.

| Script | Figure |
|--------|--------|
| `src/R/Fig1.R` | Figure 1 — global transcriptome remodelling and mRNA-TE decoupling |
| `src/R/Fig2.R` | Figure 2 — translational engagement precedes transcript change (TPI) |
| `src/R/Fig3.R` | Figure 3 — five mRNA-TE coordination programs |
| `src/R/Fig4.R` | Figure 4 — translational-engagement compensation (CS and TBR) |
| `src/R/Fig5.R` | Figure 5 — biphasic vTE oscillation and compensation waves |
| `src/R/Fig6.R` | Figure 6 — functional 3'UTR CPE cis-regulation |
| `src/R/Fig7.R` | Figure 7 — GV-state predictability of trajectory |

## Requirements

- R >= 4.3
- CRAN: `ggplot2`, `patchwork`, `dplyr`, `tidyr`, `ggrepel`, `ggridges`, `ggplotify`,
  `svglite`, `ragg`, `yaml`, `readxl`, `scales`, `xgboost`, `pROC`
- Bioconductor: `DESeq2`, `edgeR`, `RUVSeq`, `ComplexHeatmap`, `circlize`,
  `clusterProfiler`, `org.Mm.eg.db`, `Mfuzz`

## Running

Run every script from the repository root, **under a UTF-8 locale** (the config and the
figure labels use characters such as Delta, arrow, and subscripts; under a `C` locale
`conf/params.yaml` fails to parse and labels render as mojibake):

```bash
export LANG=en_US.UTF-8 LC_ALL=en_US.UTF-8

for s in 01_qc 02_normalize 03_baseline 04_tta 05_tpi 06_trd 07_cis_grammar; do
  Rscript src/R/$s.R
done
for s in src/R/S[1-8]_*.R; do Rscript "$s"; done
for n in 1 2 3 4 5 6 7; do Rscript src/R/Fig$n.R; done
```

## Data availability

Raw counts, the constGenes reference panel, and all generated outputs (`results/`,
high-resolution TIFF) are not tracked here. Sequencing data will be deposited in a public
archive on publication; contact the corresponding author for access before then.
