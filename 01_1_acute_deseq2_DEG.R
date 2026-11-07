# =============================================================================
# 01_1_acute_deseq2_DEG.R
#
# Acute age-group-specific DEG analysis
#   - 12h vs Base within each AgeGroup
#   - Paired design within each AgeGroup: ~ monkey + timepoint
#   - Cynomolgus Ensembl gene ID -> human one-to-one ortholog annotation
#   - DEG count barplot
#   - DEG overlap / UpSet plot
#
# IMPORTANT
#   This script is self-contained and does not require objects created by
#   01_0_qc_acute_deseq2.R.
#
#   This analysis describes within-age-group exposure-associated changes.

# =============================================================================

suppressPackageStartupMessages({
  library(DESeq2)
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(readr)
  library(ComplexUpset)
})

# -----------------------------------------------------------------------------
# 0. Project configuration
# -----------------------------------------------------------------------------

PROJECT_DIR <- Sys.getenv(
  "METH_PROJECT_DIR",
  unset = "/storage/home/csh91/lab2/RNAseq/RES_04_DrugMETH_cyno"
)

CODE_DIR <- file.path(PROJECT_DIR, "final/final_fi/code_acute")

# Common I/O / metadata utility functions only
source(file.path(CODE_DIR, "00_utils.R"))

COUNT_MATRIX_FILE <- file.path(PROJECT_DIR, "counts_matrix.tsv.gz")
META_FILE         <- file.path(PROJECT_DIR, "sample_sheet.txt")
ORTHOLOG_FILE     <- file.path(PROJECT_DIR, "macaque_human_ortholog.tsv")

# DEG-specific output directory (separate from QC output)
OUT_DIR <- file.path(PROJECT_DIR, "results/acute/DEG_by_AgeGroup")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Analysis parameters
MIN_COUNT   <- 10L
MIN_SAMPLES <- 3L
FDR_CUTOFF  <- 0.05

# Local plotting helper: no dependency on 01_qc_plot_helpers.R
save_plot <- function(p, path, w = 7, h = 5) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(
    filename = path,
    plot = p,
    width = w,
    height = h,
    dpi = 300
  )
}

# -----------------------------------------------------------------------------
# 1. Load counts and metadata
# -----------------------------------------------------------------------------

log_msg("Loading counts and metadata...")

counts_all <- read_count_matrix_tsv(COUNT_MATRIX_FILE)
meta_all   <- read_metadata(META_FILE)

meta <- meta_all %>%
  dplyr::filter(administration == "acute") %>%
  droplevels()

aligned <- align_counts_meta(counts_all, meta, strict = FALSE)
counts  <- aligned$counts
meta    <- aligned$meta

stopifnot(ncol(counts) == nrow(meta))
stopifnot(identical(colnames(counts), as.character(meta$sample_id)))

# Explicit factor levels
meta$AgeGroup  <- factor(meta$AgeGroup, levels = c("G1", "G2", "G3", "G4"))
meta$timepoint <- factor(meta$timepoint, levels = c("Base", "12h"))

log_msg(sprintf(
  "Acute cohort loaded: %d genes x %d samples",
  nrow(counts), ncol(counts)
))
print(table(meta$AgeGroup, meta$timepoint))

# -----------------------------------------------------------------------------
# 2. Human ortholog annotation
# -----------------------------------------------------------------------------

log_msg("Loading cynomolgus-human one-to-one ortholog annotation...")

ortholog <- readr::read_tsv(
  ORTHOLOG_FILE,
  show_col_types = FALSE
)

orth_clean <- ortholog %>%
  dplyr::rename(
    gene         = ensembl_gene_id,
    human_symbol = hsapiens_homolog_associated_gene_name,
    orth_type    = hsapiens_homolog_orthology_type
  ) %>%
  dplyr::filter(
    orth_type == "ortholog_one2one",
    !is.na(human_symbol),
    human_symbol != ""
  ) %>%
  dplyr::distinct(
    gene,
    .keep_all = TRUE
  )

log_msg(
  sprintf(
    "One-to-one human orthologs loaded: %d genes",
    nrow(orth_clean)
  )
)


# -----------------------------------------------------------------------------
# Annotation coverage
# -----------------------------------------------------------------------------

annotation_coverage <- mean(
  rownames(counts) %in% orth_clean$gene
)

log_msg(
  sprintf(
    "Human ortholog annotation coverage: %.1f%%",
    annotation_coverage * 100
  )
)

# =============================================================================
# 1. AgeGroup-specific paired DEG analysis
# =============================================================================

# timepoint factor order fixation (Base가 reference)
meta$timepoint <- factor(meta$timepoint, levels = c("Base", "12h"))

age_groups <- c("G1", "G2", "G3", "G4")
deg_list   <- list()

for (ag in age_groups) {
  
  log_msg("Processing", ag)
  
  # ---------------------------------------------------------------------------
  # Metadata subset
  # ---------------------------------------------------------------------------
  
  meta_sub <- meta %>%
    filter(AgeGroup == ag) %>%
    droplevels()
  
  meta_sub$monkey    <- droplevels(factor(meta_sub$monkey))
  meta_sub$timepoint <- factor(meta_sub$timepoint, levels = c("Base", "12h"))
  
  # ---------------------------------------------------------------------------
  # Sample matching
  # ---------------------------------------------------------------------------
  
  samples <- intersect(
    meta_sub$sample_id,
    colnames(counts)
  )
  
  meta_sub <- meta_sub %>%
    filter(sample_id %in% samples) %>%
    arrange(match(sample_id, samples))
  
  rownames(meta_sub) <- meta_sub$sample_id
  
  counts_sub <- counts[, samples, drop = FALSE]
  
  
  # sample order confirmed of count matrix and metadata
  stopifnot(
    identical(
      colnames(counts_sub),
      rownames(meta_sub)
    )
  )
  
  # ---------------------------------------------------------------------------
  # Paired sample check
  # ---------------------------------------------------------------------------
  
  pair_check <- table(
    meta_sub$monkey,
    meta_sub$timepoint
  )

  if (!all(c("Base", "12h") %in% colnames(pair_check))) {
    stop(paste0("Required Base/12h timepoints are missing in ", ag))
  }
  
  if (
    !all(
      pair_check[, "Base"] == 1 &
      pair_check[, "12h"] == 1
    )
  ) {
    
    stop(
      paste0(
        "Incomplete Base-12h pairing detected in ",
        ag
      )
    )
  }
  
  
  log_msg(
    sprintf(
      "%s → monkeys: %d, samples: %d",
      ag,
      length(unique(meta_sub$monkey)),
      nrow(meta_sub)
    )
  )
  
  # ---------------------------------------------------------------------------
  # DESeq2 dataset
  # ---------------------------------------------------------------------------
  
  dds_sub <- DESeqDataSetFromMatrix(
    countData = counts_sub,
    colData   = meta_sub,
    design    = ~ monkey + timepoint
  )
  
  
  # ---------------------------------------------------------------------------
  # Low-count filtering
  #
  # Keep genes with count >= 10 in at least 3 samples
  # ---------------------------------------------------------------------------
  keep <- rowSums(
    counts(dds_sub) >= MIN_COUNT
  ) >= MIN_SAMPLES
  
  dds_sub <- dds_sub[
    keep,
  ]

  if (nrow(dds_sub) == 0L) {
    stop(paste0("No genes passed low-count filtering in ", ag))
  }
  
  
  log_msg(
    sprintf(
      "%s → genes after filtering: %d",
      ag,
      nrow(dds_sub)
    )
  )
  
  # ---------------------------------------------------------------------------
  # Run DESeq2
  # ---------------------------------------------------------------------------
  
  dds_sub <- DESeq(dds_sub)
  
  # ---------------------------------------------------------------------------
  # Contrast: 12h vs Base
  #
  # log2FC > 0 : increased at 12h
  # log2FC < 0 : decreased at 12h
  # ---------------------------------------------------------------------------
  
  res <- results(dds_sub, contrast = c("timepoint", "12h", "Base"))
  
  # ---------------------------------------------------------------------------
  # Add human ortholog annotation
  # ---------------------------------------------------------------------------
  
  res_df <- as.data.frame(res) %>%
    tibble::rownames_to_column(
      "gene"
    ) %>%
    
    dplyr::left_join(
      orth_clean %>%
        dplyr::select(
          gene,
          human_symbol
        ),
      by = "gene"
    ) %>%
    
    dplyr::mutate(
      
      AgeGroup = ag,
      
      significant =
        !is.na(padj) &
        padj < FDR_CUTOFF,
      
      Regulation = dplyr::case_when(
        
        significant &
          log2FoldChange > 0 ~ "up",
        
        significant &
          log2FoldChange < 0 ~ "down",
        
        TRUE ~ "ns"
      )
    ) %>%
    
    dplyr::select(
      gene,
      human_symbol,
      AgeGroup,
      baseMean,
      log2FoldChange,
      lfcSE,
      stat,
      pvalue,
      padj,
      significant,
      Regulation
    )
  
  # ---------------------------------------------------------------------------
  # Save DEG result
  # ---------------------------------------------------------------------------
  
  write.csv(res_df,
            file.path(OUT_DIR, paste0("DEG_", ag, ".csv")),
            row.names = FALSE)
  
  # Save to list
  deg_list[[ag]] <- res_df
  
  # ---------------------------------------------------------------------------
  # DEG summary
  # ---------------------------------------------------------------------------
  
  n_up   <- sum(res_df$Regulation == "up")
  n_down <- sum(res_df$Regulation == "down")
  log_msg(sprintf("%s → UP: %d, DOWN: %d", ag, n_up, n_down))
}

# =============================================================================
# 2. Combine all DEG results
# =============================================================================

log_msg("Combining DEG results...")

deg_all <- dplyr::bind_rows(
  deg_list
)


write.csv(
  deg_all,
  file.path(
    OUT_DIR,
    "DEG_all_groups_combined.csv"
  ),
  row.names = FALSE
)


saveRDS(
  deg_list,
  file.path(
    OUT_DIR,
    "DEG_results_by_age.rds"
  )
)


log_msg("AgeGroup-specific DEG analysis completed.")

# =============================================================================
# DEG visualization and overlap analysis
# =============================================================================
# -----------------------------------------------------------------------------
# 1. UP/DOWN barplot
# -----------------------------------------------------------------------------
log_msg("DEG UP/DOWN barplot...")

deg_ud <- bind_rows(
  lapply(names(deg_list), function(ag) {
    deg_list[[ag]] %>%
      as.data.frame() %>%
      dplyr::filter(!is.na(padj), padj < FDR_CUTOFF) %>%
      mutate(
        direction = dplyr::case_when(
          log2FoldChange > 0 ~ "Up",
          log2FoldChange < 0 ~ "Down",
          TRUE ~ NA_character_
        ),
        AgeGroup  = ag
      ) %>%
      dplyr::filter(!is.na(direction)) %>%
      dplyr::group_by(AgeGroup, direction) %>%
      summarise(n = n(), .groups = "drop")
  })
) %>%
  dplyr::mutate(
    n_plot   = ifelse(direction == "Down", -n, n),
    AgeGroup = factor(AgeGroup, levels = c("G1","G2","G3","G4"))
  )

p_ud <- ggplot(deg_ud,
               aes(x = AgeGroup, y = n_plot, fill = direction)) +
  geom_col(width = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.2) +
  geom_text(aes(label = abs(n_plot),
                vjust = ifelse(n_plot > 0, -0.3, 1.2)),
            size = 3) +
  annotate("text", x = 1, y = max(deg_ud$n_plot) * 0.9,
           label = "Up",   size = 4, fontface = "bold") +
  annotate("text", x = 1, y = min(deg_ud$n_plot) * 0.9,
           label = "Down", size = 4, fontface = "bold") +
  scale_fill_manual(values = c(Up = "grey2", Down = "grey")) +
  labs(x = "AgeGroup", y = "Number of DEGs") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none")

save_plot(p_ud, file.path(OUT_DIR, "deg_updown.png"), w = 4, h = 5)

# -----------------------------------------------------------------------------
# 2. DEG union define
# -----------------------------------------------------------------------------
log_msg("DEG set operations...")

get_deg_genes <- function(df) {
  df %>%
    as.data.frame() %>%
    filter(!is.na(padj), padj < FDR_CUTOFF) %>%
    pull(gene) %>%
    unique()
}

g1_deg <- get_deg_genes(deg_list[["G1"]])
g2_deg <- get_deg_genes(deg_list[["G2"]])
g3_deg <- get_deg_genes(deg_list[["G3"]])
g4_deg <- get_deg_genes(deg_list[["G4"]])

# G2-G4 shared set INCLUDING genes also significant in G1.
# In the current analysis this corresponds to the 321-gene G2-G4-convergent set.
common_G234_inclusive <- Reduce(intersect, list(g2_deg, g3_deg, g4_deg))

# G2-G4 shared set EXCLUDING genes also significant in G1.
# This corresponds to the G2∩G3∩G4-only intersection shown separately in the UpSet plot.
common_G234_exclusive <- setdiff(common_G234_inclusive, g1_deg)

# Genes significant in all four groups.
common_all <- Reduce(intersect, list(g1_deg, g2_deg, g3_deg, g4_deg))

common_23 <- setdiff(
  intersect(g2_deg, g3_deg),
  union(g1_deg, g4_deg)
)

g4_specific <- setdiff(
  g4_deg,
  union(g1_deg, union(g2_deg, g3_deg))
)

log_msg("G2-G4 shared (inclusive of G1 overlap): ", length(common_G234_inclusive))
log_msg("G2-G4-only (excluding G1):             ", length(common_G234_exclusive))
log_msg("All four groups:                       ", length(common_all))
log_msg("G2-G3 specific:                        ", length(common_23))
log_msg("G4 specific:                           ", length(g4_specific))

# Save gene lists
write.table(
  data.frame(gene = common_G234_inclusive),
  file.path(OUT_DIR, "genelist_common_G234_inclusive.txt"),
  quote = FALSE,
  row.names = FALSE
)

write.table(
  data.frame(gene = common_G234_exclusive),
  file.path(OUT_DIR, "genelist_common_G234_exclusive.txt"),
  quote = FALSE,
  row.names = FALSE
)
write.table(data.frame(gene = common_23),
            file.path(OUT_DIR, "genelist_common_G23_specific.txt"),
            quote = FALSE, row.names = FALSE)
write.table(data.frame(gene = g4_specific),
            file.path(OUT_DIR, "genelist_G4_specific.txt"),
            quote = FALSE, row.names = FALSE)

# -----------------------------------------------------------------------------
# 3. UpSet plot
# -----------------------------------------------------------------------------
log_msg("UpSet plot...")

sets     <- list(G1 = g1_deg, G2 = g2_deg, G3 = g3_deg, G4 = g4_deg)
all_genes <- unique(unlist(sets))

upset_df <- data.frame(gene = all_genes)
for (nm in names(sets)) {
  upset_df[[nm]] <- upset_df$gene %in% sets[[nm]]
}

p_upset <- upset(
  upset_df,
  intersect   = c("G1","G2","G3","G4"),
  name        = "DEG overlap",
  sort_sets   = FALSE,
  set_sizes = upset_set_size(position = "right") + scale_x_reverse(),
  stripes = "white",
  themes = upset_default_themes(
    plot.background = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA)
  ),
  queries = list(
    # 314 : G2·G3·G4 common (not included G1) — concordance core without G1-deviate part
    upset_query(intersect = c("G2","G3","G4"),
                color = "#0072B2",
                only_components = c("intersections_matrix", "intersection_size")),
    # 7 : G1·G2·G3·G4 common (all age groups) — distinguished color with 321
    upset_query(intersect = c("G1","G2","G3","G4"),
                color = "#E69F00",
                only_components = c("intersections_matrix", "intersection_size"))
    
  )
)

save_plot(p_upset, file.path(OUT_DIR, "UpSet_DEG_overlap.png"), w = 7, h = 5)

# Save session information for reproducibility
capture.output(
  sessionInfo(),
  file = file.path(OUT_DIR, "sessionInfo_01_1_acute_DEG.txt")
)

log_msg("01_1_acute_deseq2_DEG.R completed.")
log_msg(paste0("Output directory: ", OUT_DIR))
