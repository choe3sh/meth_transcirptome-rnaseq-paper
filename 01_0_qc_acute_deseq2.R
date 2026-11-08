# =============================================================================
# 02_0_qc_acute_deseq2.R
# Acute QC
#   1. Library size / Detected genes / lib vs detected
#   2. VST PCA (AgeGroup × timepoint)
# =============================================================================
suppressPackageStartupMessages({
  library(DESeq2)
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(readr)
})

source("/storage/home/csh91/lab2/RNAseq/RES_04_DrugMETH_cyno/final/final_fi/code_acute/00_utils.R")

# -----------------------------------------------------------------------------
# 경로 설정
# -----------------------------------------------------------------------------
PROJECT_DIR       <- "/storage/home/csh91/lab2/RNAseq/RES_04_DrugMETH_cyno"
COUNT_MATRIX_FILE <- file.path(PROJECT_DIR, "counts_matrix.tsv.gz")
META_FILE         <- file.path(PROJECT_DIR, "sample_sheet.txt")
GMT_FILE          <- file.path(PROJECT_DIR, "h.all.v2023.2.Hs.symbols.gmt")
ORTHOLOG_FILE     <- file.path(PROJECT_DIR, "macaque_human_ortholog.tsv")
OUT_DIR           <- file.path(PROJECT_DIR, "results/qc/acute")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

AGE_COLORS <- c(G1 = "#4393C3", G2 = "#74C476", G3 = "#FD8D3C", G4 = "#D6604D")

# -----------------------------------------------------------------------------
# 1. 데이터 로드
# -----------------------------------------------------------------------------
log_msg("Loading counts & metadata...")
counts_all <- read_count_matrix_tsv(COUNT_MATRIX_FILE)
meta_all   <- read_metadata(META_FILE)

meta    <- meta_all %>% filter(administration == "acute") %>% droplevels()
aligned <- align_counts_meta(counts_all, meta, strict = FALSE)
counts  <- aligned$counts
meta    <- aligned$meta

log_msg(sprintf("Acute: %d genes x %d samples", nrow(counts), ncol(counts)))
print(table(meta$AgeGroup, meta$timepoint))

# -----------------------------------------------------------------------------
# 2. QC 기본 지표
# -----------------------------------------------------------------------------
log_msg("Building QC summary...")
qc_tbl <- qc_summary_tbl(counts, meta)
write_tsv_gz(qc_tbl, file.path(OUT_DIR, "qc_summary.tsv.gz"))

save_gg(plot_qc_libsize(qc_tbl, color_by = "AgeGroup"),
        file.path(OUT_DIR, "01_libsize_by_AgeGroup.png"))
save_gg(plot_qc_detected(qc_tbl, color_by = "AgeGroup"),
        file.path(OUT_DIR, "02_detected_genes_by_AgeGroup.png"))
save_gg(plot_qc_lib_vs_detected(qc_tbl, color_by = "AgeGroup",
                                label_points = TRUE),
        file.path(OUT_DIR, "03_libsize_vs_detected.png"))

# -----------------------------------------------------------------------------
# 3. acute QC용 VST 생성, acute 전체 (이후 PCA, ssGSEA, shift distance 공용)
# -----------------------------------------------------------------------------
log_msg("VST normalization...")
dds_acute_qc <- DESeqDataSetFromMatrix(
  countData = counts,
  colData   = meta,
  design    = ~ monkey + timepoint
)
dds_acute_qc <- dds_acute_qc[rowSums(counts(dds_acute_qc) >= 10) >= 3, ]
dds_acute_qc <- estimateSizeFactors(dds_acute_qc)
vsd_acute <- vst(dds_acute_qc, blind = TRUE)

# -----------------------------------------------------------------------------
# 4. VST PCA
# -----------------------------------------------------------------------------
log_msg("VST PCA...")
pca    <- prcomp(t(assay(vsd_acute)), scale. = FALSE)
pct    <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
pca_df <- as.data.frame(pca$x[, 1:2]) %>%
  tibble::rownames_to_column("sample_id") %>%
  left_join(meta, by = "sample_id")

p_pca <- ggplot(pca_df,
                aes(PC1, PC2, color = AgeGroup, shape = timepoint,
                    label = sample_id)) +
  geom_point(size = 3, alpha = 0.9) +
  geom_text_repel(size = 2, show.legend = FALSE,
                  max.overlaps = 20, segment.color = NA) +
  scale_shape_manual(values = c(Base = 16, `12h` = 17)) +
  scale_color_manual(values = AGE_COLORS) +
  labs(title = "Acute QC: VST PCA",
       x = paste0("PC1 (", pct[1], "%)"),
       y = paste0("PC2 (", pct[2], "%)")) +
  theme_bw(base_size = 12)

save_gg(p_pca, file.path(OUT_DIR, "04_vst_PCA.png"), w = 8, h = 6)



log_msg("02_0_qc_acute_deseq2.R 완료 ->", OUT_DIR)


