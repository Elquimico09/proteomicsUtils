#' convertFormat
#' @description
#' Converts ANOVA or Kruskal-Wallis results to the format expected by makeVolcano. Extracts
#' the p-value and log2 fold change for a specific comparison and renames them to pvalue,
#' nlog10p, and log2fc. Automatically detects whether input is from performANOVA
#' (tukey_pvalue) or performKW (dunn_pvalue). Prefers `log2fc_*` columns when available
#' and falls back to deriving from `foldchange_*` for backward compatibility.
#'
#' @param results Data frame output from performANOVA/performKW (legacy) or
#'   list output containing `data` and `cohort_columns` (current).
#' @param comparison Character vector of length 2 specifying the comparison,
#'   e.g., c("control", "treated"). Order doesn't matter.
#'
#' @return A data frame with added columns:
#'   - pvalue: the Tukey/Dunn p-value for the specified comparison
#'   - nlog10p: -log10 transformed p-value
#'   - log2fc: log2 fold change for the specified comparison
#' @export
#'
#' @examples
#' # From ANOVA:
#' anova_res <- performANOVA(control, treated, placebo)
#' volcano_data <- convertFormat(anova_res, c("control", "treated"))
#' makeVolcano(volcano_data)
#'
#' # From Kruskal-Wallis:
#' kw_res <- performKW(control, treated, placebo)
#' volcano_data <- convertFormat(kw_res, c("control", "treated"))
#' makeVolcano(volcano_data)
convertFormat <- function(results, comparison) {
  # Backward/forward compatibility:
  # - legacy: results is a data frame
  # - current: results is list(data=..., cohort_columns=...)
  if (is.list(results) && !is.data.frame(results) && "data" %in% names(results)) {
    results <- results$data
  }

  if (length(comparison) != 2) {
    stop("comparison must be a character vector of length 2")
  }

  # Try both orderings of the comparison to find the matching column
  pair_name1 <- paste0(comparison[1], "_", comparison[2])
  pair_name2 <- paste0(comparison[2], "_", comparison[1])

  # Check for tukey_pvalue (from performANOVA) or dunn_pvalue (from performKW)
  tukey_col1 <- paste0("tukey_pvalue_", pair_name1)
  tukey_col2 <- paste0("tukey_pvalue_", pair_name2)
  dunn_col1 <- paste0("dunn_pvalue_", pair_name1)
  dunn_col2 <- paste0("dunn_pvalue_", pair_name2)
  log2fc_col1 <- paste0("log2fc_", pair_name1)
  log2fc_col2 <- paste0("log2fc_", pair_name2)
  fc_col1 <- paste0("foldchange_", pair_name1)
  fc_col2 <- paste0("foldchange_", pair_name2)

  # Find which column exists
  pvalue_col <- NULL
  log2fc_col <- NULL
  fc_col <- NULL
  flip_sign <- FALSE

  if (tukey_col1 %in% names(results)) {
    pvalue_col <- tukey_col1
    log2fc_col <- log2fc_col1
    fc_col <- fc_col1
  } else if (tukey_col2 %in% names(results)) {
    pvalue_col <- tukey_col2
    log2fc_col <- log2fc_col2
    fc_col <- fc_col2
    flip_sign <- TRUE
  } else if (dunn_col1 %in% names(results)) {
    pvalue_col <- dunn_col1
    log2fc_col <- log2fc_col1
    fc_col <- fc_col1
  } else if (dunn_col2 %in% names(results)) {
    pvalue_col <- dunn_col2
    log2fc_col <- log2fc_col2
    fc_col <- fc_col2
    flip_sign <- TRUE
  } else {
    # List available comparisons from either tukey or dunn columns
    available <- grep("^(tukey_pvalue_|dunn_pvalue_)", names(results), value = TRUE)
    available <- gsub("^(tukey_pvalue_|dunn_pvalue_)", "", available)
    stop(paste0("Comparison '", comparison[1], "' vs '", comparison[2], "' not found.\n",
                "Available comparisons: ", paste(available, collapse = ", ")))
  }

  # Create output with renamed columns
  result <- results
  result$pvalue <- result[[pvalue_col]]
  if (log2fc_col %in% names(results)) {
    result$log2fc <- if (flip_sign) -result[[log2fc_col]] else result[[log2fc_col]]
  } else if (fc_col %in% names(results)) {
    # Backward-compatible fallback for older outputs without explicit log2fc columns.
    result$log2fc <- if (flip_sign) -result[[fc_col]] else result[[fc_col]]
  } else {
    stop(paste0("Could not find log2fc/foldchange column for comparison '",
                comparison[1], "' vs '", comparison[2], "'."))
  }
  result$nlog10p <- -log10(result$pvalue)

  return(result)
}
