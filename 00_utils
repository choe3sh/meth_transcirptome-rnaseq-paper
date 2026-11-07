# ====================================================================
# 00_utils.R
# Utility functions only - packages are loaded in individual analysis
# ====================================================================

log_msg <- function(...) {
  cat(sprintf("[%s] %s\n",
              format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
              paste(..., collapse = " ")))
}

stop_if_not_exists <- function(path, what = "file") {
  if (!file.exists(path)) stop(what, " not found: ", path)
}

assert_no_duplicates <- function(x, what = "values") {
  if (anyDuplicated(x)) {
    dup <- unique(x[duplicated(x)])
    stop("Duplicate ", what, ": ", paste(dup, collapse = ", "))
  }
}

write_tsv_gz <- function(df, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  gz <- gzfile(path, "wt")
  on.exit(close(gz), add = TRUE)
  write.table(df, file = gz, sep = "\t", quote = FALSE, row.names = FALSE)
}

read_count_matrix_tsv <- function(path) {
  stop_if_not_exists(path, "Count matrix file")
  
  # Note: counts_matrix.tsv.gz is a plain (uncompressed) TSV with CRLF line endings, despite the .gz extension

  df <- readr::read_tsv(path, show_col_types = FALSE) %>%
    as.data.frame()

  if (!("gene_id" %in% colnames(df)))
    stop("Count matrix must contain a 'gene_id' column: ", path)

  rownames(df) <- df$gene_id
  mat <- as.matrix(df[, setdiff(colnames(df), "gene_id"), drop = FALSE])

  if (ncol(mat) < 2)  stop("Count matrix has <2 samples: ", path)
  if (anyNA(mat))     stop("Count matrix contains NA values: ", path)
  suppressWarnings(storage.mode(mat) <- "integer")
  if (any(mat < 0))   stop("Count matrix contains negative values: ", path)

  mat
}

align_counts_meta <- function(counts, meta, strict = TRUE, verbose = TRUE) {
  if (!("sample_id" %in% colnames(meta)))
    stop("Metadata must contain 'sample_id' column.")

  meta$sample_id <- as.character(meta$sample_id)
  in_counts_not_meta <- setdiff(colnames(counts), meta$sample_id)
  in_meta_not_counts <- setdiff(meta$sample_id, colnames(counts))

  if (length(in_counts_not_meta) > 0) {
    msg <- paste0("sample_id present in counts but not in metadata: ",
                  paste(in_counts_not_meta, collapse = ", "))
    if (strict) stop(msg) else if (verbose) log_msg("WARN", msg)
  }
  if (length(in_meta_not_counts) > 0) {
    msg <- paste0("sample_id present in metadata but not in counts: ",
                  paste(in_meta_not_counts, collapse = ", "))
    if (strict) stop(msg) else if (verbose) log_msg("WARN", msg)
  }

  meta2   <- meta[meta$sample_id %in% colnames(counts), , drop = FALSE]
  counts2 <- counts[, meta2$sample_id, drop = FALSE]

  if (!identical(colnames(counts2), meta2$sample_id))
    stop("Counts columns are not aligned to metadata sample_id order after subsetting.")

  list(counts = counts2, meta = meta2)
}

read_metadata <- function(path) {
  stop_if_not_exists(path, "Metadata file")
  meta <- readr::read_tsv(path, na = c("NA", "", "NaN"), show_col_types = FALSE)

  required <- c("sample_id", "administration", "timepoint", "monkey")
  missing  <- setdiff(required, colnames(meta))
  if (length(missing) > 0)
    stop("metadata missing columns: ", paste(missing, collapse = ", "),
         "\nColumns found: ", paste(colnames(meta), collapse = ", "))

  # AgeGroup is taken as-is from the sample sheet
  if (!("AgeGroup" %in% colnames(meta))) {
    log_msg("WARN: AgeGroup column not found - setting to NA")
    meta$AgeGroup <- NA_character_
  }

  if (!("time_label" %in% colnames(meta)))
    meta <- meta %>%
      dplyr::mutate(time_label = paste(administration, timepoint, sep = "_"))

  meta <- meta %>%
    dplyr::mutate(
      sample_id      = as.character(sample_id),
      administration = factor(administration, levels = c("acute", "chronic")),
      # Timepoint values: Base, 12h, W2, W5, W9
      timepoint      = factor(timepoint, levels = c("Base", "12h", "W2", "W5", "W9")),
      time_label     = factor(time_label),
      AgeGroup       = factor(AgeGroup, levels = c("G1", "G2", "G3", "G4")),
      monkey         = factor(monkey)
    )

  if ("sex" %in% colnames(meta)) meta$sex <- factor(meta$sex)

  assert_no_duplicates(meta$sample_id, what = "meta$sample_id")
  meta
}

qc_summary_tbl <- function(counts, meta) {
  tibble::tibble(
    sample_id      = colnames(counts),
    lib_size       = as.numeric(colSums(counts)),
    detected_genes = as.integer(colSums(counts > 0)),
    zero_rate      = round(colMeans(counts == 0), 3)
  ) %>% dplyr::left_join(meta, by = "sample_id")
}
