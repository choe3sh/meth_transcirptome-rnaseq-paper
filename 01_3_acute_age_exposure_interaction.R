# =============================================================================
# 01_3_acute_age_exposure_interaction.R
#
# age-group-by-exposure interaction analysis
#
#  G1-versus-G2-G4 contrast."
#
# Primary planned contrast:
#
#   (12h - Base in G1)
#      -
#   mean[(12h - Base in G2),
#        (12h - Base in G3),
#        (12h - Base in G4)]
#
# Model:
#   individual (monkey) fixed effects
#   + age-group-specific exposure effects
#
#   design matrix = ~ 0 + monkey + AgeGroup:timepoint12h
#
# AgeGroup main effects are NOT included because AgeGroup is constant within
# each monkey and is therefore absorbed by the monkey fixed effects.
#
# This parameterization:
#   - preserves Base/12h within-animal pairing,
#   - estimates one exposure effect for each age group,
#   - is explicitly checked for full column rank,
#   - allows the planned G1-vs-G2-G4 interaction contrast.
#
# IMPORTANT:
# This formal interaction analysis is independent of the 321-gene G2-G4
# shared-DEG selection used for descriptive/convergence analyses in 01_1/01_2.
# Genome-wide multiple-testing correction is performed across all tested genes.
#
# Outputs:
#   interaction_G1_vs_G234_all_genes.csv
#   interaction_G1_vs_G234_significant.csv
#   interaction_principal_genes.csv
#   interaction_ranked_genes_for_pathway.tsv
#   interaction_model_coefficients.csv
#   interaction_model_matrix.csv
#   interaction_summary.txt
#   sessionInfo_01_3_acute_age_exposure_interaction.txt
#
# Optional outputs when 01_1 gene lists are available:
#   interaction_G234_321_genes.csv
#
# NOTE ON PATHWAYS:
# This script creates the genome-wide interaction ranking required for a
# pathway-level interaction analysis. The actual FGSEA should use this ranking
# in the dedicated pathway script so that gene-level modeling and pathway
# enrichment remain modular and reproducible.
# =============================================================================

suppressPackageStartupMessages({
  library(DESeq2)
  library(dplyr)
  library(tibble)
  library(readr)
})

# -----------------------------------------------------------------------------
# 0. Project configuration
# -----------------------------------------------------------------------------

PROJECT_DIR <- Sys.getenv(
  "METH_PROJECT_DIR",
  unset = "/storage/home/csh91/lab2/RNAseq/RES_04_DrugMETH_cyno"
)

CODE_DIR <- file.path(
  PROJECT_DIR,
  "final/final_fi/code_acute"
)

source(file.path(CODE_DIR, "00_utils.R"))

COUNT_MATRIX_FILE <- file.path(
  PROJECT_DIR,
  "counts_matrix.tsv.gz"
)

META_FILE <- file.path(
  PROJECT_DIR,
  "sample_sheet.txt"
)

ORTHOLOG_FILE <- file.path(
  PROJECT_DIR,
  "macaque_human_ortholog.tsv"
)

OUT_DIR <- file.path(
  PROJECT_DIR,
  "results/acute/age_exposure_interaction"
)

dir.create(
  OUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

MIN_COUNT   <- 10L
MIN_SAMPLES <- 3L
FDR_CUTOFF  <- 0.05

AGE_LEVELS  <- c("G1", "G2", "G3", "G4")
TIME_LEVELS <- c("Base", "12h")

PRINCIPAL_GENES <- c(
  "S100A9",
  "CD177",
  "VNN2",
  "CLEC4D",
  "GCA",
  "CCL5"
)

log_msg <- function(...) {
  message(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", ...)
}

# -----------------------------------------------------------------------------
# 1. Load acute raw counts and metadata
# -----------------------------------------------------------------------------

log_msg("Loading count matrix and metadata...")

counts_all <- read_count_matrix_tsv(
  COUNT_MATRIX_FILE
)

meta_all <- read_metadata(
  META_FILE
)

meta <- meta_all %>%
  filter(administration == "acute") %>%
  droplevels()

if (nrow(meta) == 0L) {
  stop("No acute samples found in metadata.")
}

aligned <- align_counts_meta(
  counts_all,
  meta,
  strict = FALSE
)

counts <- aligned$counts
meta   <- aligned$meta

stopifnot(
  ncol(counts) == nrow(meta)
)

stopifnot(
  identical(
    colnames(counts),
    as.character(meta$sample_id)
  )
)

# Explicit factor levels.
meta$AgeGroup <- factor(
  as.character(meta$AgeGroup),
  levels = AGE_LEVELS
)

meta$timepoint <- factor(
  as.character(meta$timepoint),
  levels = TIME_LEVELS
)

meta$monkey <- factor(
  as.character(meta$monkey)
)

if (anyNA(meta$AgeGroup)) {
  stop("Unexpected or missing AgeGroup values.")
}

if (anyNA(meta$timepoint)) {
  stop("Unexpected or missing acute timepoint values.")
}

# -----------------------------------------------------------------------------
# 2. Verify the paired acute design
# -----------------------------------------------------------------------------

log_msg("Checking Base/12h pairing within animals...")

pair_check <- meta %>%
  count(
    monkey,
    AgeGroup,
    timepoint,
    name = "n"
  ) %>%
  tidyr::complete(
    monkey,
    timepoint = factor(
      TIME_LEVELS,
      levels = TIME_LEVELS
    ),
    fill = list(n = 0L)
  )

# Each monkey must have exactly one Base and one 12h sample.
pair_wide <- meta %>%
  count(
    monkey,
    AgeGroup,
    timepoint,
    name = "n"
  ) %>%
  tidyr::pivot_wider(
    names_from = timepoint,
    values_from = n,
    values_fill = 0L
  )

if (!all(c("Base", "12h") %in% colnames(pair_wide))) {
  stop("Base and/or 12h columns are missing from paired-design check.")
}

if (any(pair_wide$Base != 1L | pair_wide$`12h` != 1L)) {
  print(pair_wide)
  stop(
    "Incomplete or duplicated Base/12h pairing detected. ",
    "Interaction analysis stopped."
  )
}

# AgeGroup must be constant within monkey.
age_per_monkey <- meta %>%
  distinct(monkey, AgeGroup) %>%
  count(monkey, name = "n_age")

if (any(age_per_monkey$n_age != 1L)) {
  stop("At least one monkey is assigned to more than one AgeGroup.")
}

group_n <- meta %>%
  distinct(monkey, AgeGroup) %>%
  count(AgeGroup, name = "n_animals")

print(group_n)

# -----------------------------------------------------------------------------
# 3. Human one-to-one ortholog annotation
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

# 4. Filter low-count genes
# -----------------------------------------------------------------------------

log_msg(
  sprintf(
    "Filtering genes: count >= %d in >= %d samples...",
    MIN_COUNT,
    MIN_SAMPLES
  )
)

keep <- rowSums(
  counts >= MIN_COUNT
) >= MIN_SAMPLES

counts_filt <- counts[
  keep,
  ,
  drop = FALSE
]

if (nrow(counts_filt) == 0L) {
  stop("No genes remain after low-count filtering.")
}

log_msg(
  sprintf(
    "Genes retained: %d / %d",
    nrow(counts_filt),
    nrow(counts)
  )
)

# -----------------------------------------------------------------------------
# 5. Construct the full-rank paired interaction design matrix
#
# We create an explicit indicator for 12h:
#   exposure12h = 0 at Base
#   exposure12h = 1 at 12h
#
# Then:
#   ~ 0 + monkey + AgeGroup:exposure12h
#
# The monkey columns absorb all time-invariant between-animal differences,
# including the AgeGroup main effect.
#
# The four interaction columns estimate:
#   beta_G1 = 12h - Base in G1
#   beta_G2 = 12h - Base in G2
#   beta_G3 = 12h - Base in G3
#   beta_G4 = 12h - Base in G4
# -----------------------------------------------------------------------------

meta$exposure12h <- as.numeric(
  meta$timepoint == "12h"
)

design_formula <- ~ 0 + monkey + AgeGroup:exposure12h

X <- model.matrix(
  design_formula,
  data = meta
)

# Make coefficient names syntactically safe and deterministic.
colnames(X) <- make.names(
  colnames(X),
  unique = TRUE
)

rank_X <- qr(X)$rank
ncoef_X <- ncol(X)

log_msg(
  sprintf(
    "Design matrix rank: %d / %d columns",
    rank_X,
    ncoef_X
  )
)

if (rank_X != ncoef_X) {
  # Save the problematic matrix before stopping.
  write.csv(
    cbind(
      sample_id = meta$sample_id,
      as.data.frame(X)
    ),
    file.path(
      OUT_DIR,
      "interaction_model_matrix_RANK_DEFICIENT.csv"
    ),
    row.names = FALSE
  )

  stop(
    "Interaction design matrix is not full rank: rank=",
    rank_X,
    ", columns=",
    ncoef_X,
    ". Analysis stopped."
  )
}

write.csv(
  cbind(
    sample_id = meta$sample_id,
    monkey = as.character(meta$monkey),
    AgeGroup = as.character(meta$AgeGroup),
    timepoint = as.character(meta$timepoint),
    as.data.frame(X)
  ),
  file.path(
    OUT_DIR,
    "interaction_model_matrix.csv"
  ),
  row.names = FALSE
)

# Identify the four age-specific exposure coefficients robustly.
interaction_cols <- grep(
  "AgeGroup.*exposure12h|exposure12h.*AgeGroup",
  colnames(X),
  value = TRUE
)

if (length(interaction_cols) != 4L) {
  stop(
    "Expected four AgeGroup-specific exposure coefficients, found: ",
    paste(interaction_cols, collapse = ", ")
  )
}

coef_lookup <- tibble(
  coefficient = interaction_cols
) %>%
  mutate(
    AgeGroup = case_when(
      grepl("G1", coefficient, fixed = TRUE) ~ "G1",
      grepl("G2", coefficient, fixed = TRUE) ~ "G2",
      grepl("G3", coefficient, fixed = TRUE) ~ "G3",
      grepl("G4", coefficient, fixed = TRUE) ~ "G4",
      TRUE ~ NA_character_
    )
  )

if (
  anyNA(coef_lookup$AgeGroup) ||
  anyDuplicated(coef_lookup$AgeGroup)
) {
  print(coef_lookup)
  stop(
    "Could not uniquely map interaction coefficients to G1-G4."
  )
}

coef_lookup <- coef_lookup %>%
  arrange(
    match(AgeGroup, AGE_LEVELS)
  )

print(coef_lookup)

write.csv(
  coef_lookup,
  file.path(
    OUT_DIR,
    "interaction_model_coefficients.csv"
  ),
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# 6. Fit DESeq2 using the explicit full-rank model matrix
# -----------------------------------------------------------------------------

log_msg("Fitting DESeq2 full-rank paired interaction model...")

dds <- DESeqDataSetFromMatrix(
  countData = round(counts_filt),
  colData = meta,
  design = ~ 1
)

dds <- DESeq(
  dds,
  full = X,
  betaPrior = FALSE,
  quiet = FALSE
)

rn <- resultsNames(dds)

log_msg(
  "DESeq2 result coefficients: ",
  paste(rn, collapse = " | ")
)

# DESeq2 may internally sanitize coefficient names. Map by age group.
rn_interaction <- grep(
  "AgeGroup.*exposure12h|exposure12h.*AgeGroup",
  rn,
  value = TRUE
)

if (length(rn_interaction) != 4L) {
  stop(
    "Could not identify four age-specific exposure coefficients in resultsNames(dds): ",
    paste(rn, collapse = ", ")
  )
}

rn_lookup <- tibble(
  coefficient = rn_interaction
) %>%
  mutate(
    AgeGroup = case_when(
      grepl("G1", coefficient, fixed = TRUE) ~ "G1",
      grepl("G2", coefficient, fixed = TRUE) ~ "G2",
      grepl("G3", coefficient, fixed = TRUE) ~ "G3",
      grepl("G4", coefficient, fixed = TRUE) ~ "G4",
      TRUE ~ NA_character_
    )
  )

if (
  anyNA(rn_lookup$AgeGroup) ||
  anyDuplicated(rn_lookup$AgeGroup)
) {
  print(rn_lookup)
  stop(
    "Could not uniquely map DESeq2 result coefficients to G1-G4."
  )
}

rn_lookup <- rn_lookup %>%
  arrange(
    match(AgeGroup, AGE_LEVELS)
  )

# -----------------------------------------------------------------------------
# 7. Planned G1-vs-G2-G4 interaction contrast
#
# Contrast:
#   beta_G1 - (beta_G2 + beta_G3 + beta_G4)/3
#
# Positive:
#   Base->12h response is more positive / less negative in G1 than in G2-G4.
#
# Negative:
#   Base->12h response is more negative / less positive in G1 than in G2-G4.
#
# This is a difference-in-differences interaction effect on the log2 scale.
# -----------------------------------------------------------------------------

contrast_vector <- rep(
  0,
  length(rn)
)

names(contrast_vector) <- rn

contrast_vector[
  rn_lookup$coefficient[
    rn_lookup$AgeGroup == "G1"
  ]
] <- 1

for (ag in c("G2", "G3", "G4")) {

  contrast_vector[
    rn_lookup$coefficient[
      rn_lookup$AgeGroup == ag
    ]
  ] <- -1 / 3
}

if (
  abs(sum(contrast_vector)) > 1e-12
) {
  stop("Internal contrast-vector check failed.")
}

log_msg(
  "Planned contrast: G1 response - mean(G2,G3,G4 responses)"
)

res_interaction <- results(
  dds,
  contrast = contrast_vector,
  alpha = FDR_CUTOFF,
  independentFiltering = TRUE
)

interaction_tbl <- as.data.frame(
  res_interaction
) %>%
  rownames_to_column("gene") %>%
  left_join(
    orth_clean,
    by = "gene"
  ) %>%
  rename(
    interaction_log2FC = log2FoldChange,
    interaction_SE = lfcSE,
    interaction_stat = stat,
    interaction_pvalue = pvalue,
    interaction_padj = padj
  ) %>%
  mutate(
    interaction_direction = case_when(
      is.na(interaction_log2FC) ~ NA_character_,
      interaction_log2FC > 0 ~ "G1_more_positive",
      interaction_log2FC < 0 ~ "G1_more_negative",
      TRUE ~ "zero"
    )
  )

# -----------------------------------------------------------------------------
# 8. Extract age-specific Base->12h effects from the SAME paired model
#
# These columns make the interaction result directly interpretable and allow
# principal genes to be reported with both:
#   - age-specific exposure effects
#   - formal G1-vs-G2-G4 interaction effect and FDR
# -----------------------------------------------------------------------------

extract_age_effect <- function(ag) {

  coef_name <- rn_lookup$coefficient[
    rn_lookup$AgeGroup == ag
  ]

  if (length(coef_name) != 1L) {
    stop("Could not identify coefficient for ", ag)
  }

  x <- results(
    dds,
    name = coef_name,
    alpha = FDR_CUTOFF,
    independentFiltering = TRUE
  )

  as.data.frame(x) %>%
    rownames_to_column("gene") %>%
    transmute(
      gene,
      !!paste0("log2FC_", ag) := log2FoldChange,
      !!paste0("SE_", ag) := lfcSE,
      !!paste0("pvalue_", ag) := pvalue,
      !!paste0("padj_", ag) := padj
    )
}

age_effects <- extract_age_effect("G1") %>%
  left_join(
    extract_age_effect("G2"),
    by = "gene"
  ) %>%
  left_join(
    extract_age_effect("G3"),
    by = "gene"
  ) %>%
  left_join(
    extract_age_effect("G4"),
    by = "gene"
  )

interaction_tbl <- interaction_tbl %>%
  left_join(
    age_effects,
    by = "gene"
  ) %>%
  mutate(
    mean_log2FC_G234 = rowMeans(
      cbind(
        log2FC_G2,
        log2FC_G3,
        log2FC_G4
      ),
      na.rm = FALSE
    ),
    interaction_effect_check =
      log2FC_G1 - mean_log2FC_G234,
    interaction_effect_difference =
      interaction_log2FC - interaction_effect_check
  )

# Numerical consistency check.
effect_diff <- interaction_tbl$interaction_effect_difference
effect_diff <- effect_diff[is.finite(effect_diff)]

if (
  length(effect_diff) > 0L &&
  max(abs(effect_diff)) > 1e-6
) {
  warning(
    "Age-specific coefficient reconstruction differs from the formal ",
    "interaction contrast by >1e-6 for at least one gene."
  )
}

# Sort by adjusted P, then raw P.
interaction_tbl <- interaction_tbl %>%
  arrange(
    is.na(interaction_padj),
    interaction_padj,
    interaction_pvalue
  )

write.csv(
  interaction_tbl,
  file.path(
    OUT_DIR,
    "interaction_G1_vs_G234_all_genes.csv"
  ),
  row.names = FALSE
)

interaction_sig <- interaction_tbl %>%
  filter(
    !is.na(interaction_padj),
    interaction_padj < FDR_CUTOFF
  )

write.csv(
  interaction_sig,
  file.path(
    OUT_DIR,
    "interaction_G1_vs_G234_significant.csv"
  ),
  row.names = FALSE
)

log_msg(
  sprintf(
    "Significant G1-vs-G2-G4 interaction genes at FDR < %.2f: %d",
    FDR_CUTOFF,
    nrow(interaction_sig)
  )
)

# -----------------------------------------------------------------------------
# 9. Principal genes requested for manuscript/reviewer reporting
# -----------------------------------------------------------------------------

principal_tbl <- interaction_tbl %>%
  filter(
    human_symbol %in% PRINCIPAL_GENES
  ) %>%
  mutate(
    human_symbol = factor(
      human_symbol,
      levels = PRINCIPAL_GENES
    )
  ) %>%
  arrange(human_symbol) %>%
  mutate(
    human_symbol = as.character(human_symbol)
  ) %>%
  select(
    gene,
    human_symbol,
    log2FC_G1,
    log2FC_G2,
    log2FC_G3,
    log2FC_G4,
    mean_log2FC_G234,
    interaction_log2FC,
    interaction_SE,
    interaction_stat,
    interaction_pvalue,
    interaction_padj,
    interaction_direction
  )

write.csv(
  principal_tbl,
  file.path(
    OUT_DIR,
    "interaction_principal_genes.csv"
  ),
  row.names = FALSE
)

print(principal_tbl)

missing_principal <- setdiff(
  PRINCIPAL_GENES,
  principal_tbl$human_symbol
)

if (length(missing_principal) > 0L) {
  warning(
    "Principal genes not found after filtering/annotation: ",
    paste(missing_principal, collapse = ", ")
  )
}

# -----------------------------------------------------------------------------
# 10. Optional: annotate the previously defined G2-G4 shared 321 genes
#
# This is NOT used to calculate interaction FDR.
# It is only a convenient descriptive subset of the genome-wide interaction
# result for linking the formal test back to Figure 2.
# -----------------------------------------------------------------------------

G234_FILE <- file.path(
  PROJECT_DIR,
  "results/acute/DEG_by_AgeGroup",
  "genelist_common_G234_inclusive.txt"
)

if (file.exists(G234_FILE)) {

  g234_genes <- read.table(
    G234_FILE,
    header = TRUE,
    stringsAsFactors = FALSE
  )$gene %>%
    unique()

  interaction_321 <- interaction_tbl %>%
    filter(
      gene %in% g234_genes
    )

  write.csv(
    interaction_321,
    file.path(
      OUT_DIR,
      "interaction_G234_321_genes.csv"
    ),
    row.names = FALSE
  )

  log_msg(
    sprintf(
      "Interaction results linked to G2-G4 shared set: %d / %d genes",
      nrow(interaction_321),
      length(g234_genes)
    )
  )
}

# -----------------------------------------------------------------------------
# 11. Genome-wide ranked list for pathway-level interaction analysis
#
# Ranking statistic = DESeq2 Wald statistic for the formal interaction contrast.
#
# Use human one-to-one ortholog symbols. If more than one macaque gene maps to
# the same human symbol, retain the row with the largest absolute interaction
# statistic so that FGSEA receives unique gene names.
# -----------------------------------------------------------------------------

ranked_interaction <- interaction_tbl %>%
  filter(
    !is.na(human_symbol),
    human_symbol != "",
    is.finite(interaction_stat)
  ) %>%
  arrange(
    desc(abs(interaction_stat))
  ) %>%
  distinct(
    human_symbol,
    .keep_all = TRUE
  ) %>%
  arrange(
    desc(interaction_stat)
  ) %>%
  transmute(
    human_symbol,
    interaction_stat,
    interaction_log2FC,
    interaction_pvalue,
    interaction_padj
  )

if (anyDuplicated(ranked_interaction$human_symbol)) {
  stop("Duplicated human symbols remain in pathway ranking.")
}

write_tsv(
  ranked_interaction,
  file.path(
    OUT_DIR,
    "interaction_ranked_genes_for_pathway.tsv"
  )
)

# -----------------------------------------------------------------------------
# 12. Summary text
# -----------------------------------------------------------------------------

n_tested <- sum(
  !is.na(interaction_tbl$interaction_pvalue)
)

n_fdr <- nrow(interaction_sig)

summary_lines <- c(
  "Acute age-group-by-exposure interaction analysis",
  "================================================",
  "",
  paste0("Samples: ", nrow(meta)),
  paste0("Animals: ", nlevels(meta$monkey)),
  paste0(
    "Animals by age group: ",
    paste(
      paste0(
        group_n$AgeGroup,
        "=",
        group_n$n_animals
      ),
      collapse = ", "
    )
  ),
  "",
  paste0(
    "Low-count filter: count >= ",
    MIN_COUNT,
    " in >= ",
    MIN_SAMPLES,
    " samples"
  ),
  paste0(
    "Genes retained after count filtering: ",
    nrow(counts_filt)
  ),
  paste0(
    "Genes with interaction P value: ",
    n_tested
  ),
  "",
  paste0(
    "Design: ~ 0 + monkey + AgeGroup:exposure12h"
  ),
  paste0(
    "Design matrix rank: ",
    rank_X,
    "/",
    ncoef_X
  ),
  "",
  paste0(
    "Planned interaction contrast: ",
    "(12h-Base in G1) - mean(12h-Base in G2,G3,G4)"
  ),
  paste0(
    "Genome-wide FDR threshold: ",
    FDR_CUTOFF
  ),
  paste0(
    "Significant interaction genes: ",
    n_fdr
  ),
  "",
  "Positive interaction_log2FC = G1 response is more positive / less negative",
  "than the average response across G2-G4.",
  "Negative interaction_log2FC = G1 response is more negative / less positive",
  "than the average response across G2-G4.",
  "",
  "Pathway ranking file:",
  "interaction_ranked_genes_for_pathway.tsv"
)

writeLines(
  summary_lines,
  con = file.path(
    OUT_DIR,
    "interaction_summary.txt"
  )
)

# -----------------------------------------------------------------------------
# 13. Reproducibility information
# -----------------------------------------------------------------------------

capture.output(
  sessionInfo(),
  file = file.path(
    OUT_DIR,
    "sessionInfo_01_3_acute_age_exposure_interaction.txt"
  )
)

log_msg("01_3 acute age-by-exposure interaction analysis completed.")
log_msg("Output directory: ", OUT_DIR)
