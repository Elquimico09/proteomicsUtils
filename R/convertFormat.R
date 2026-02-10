#' convertFormat
#' @description
#' Converts ANOVA or Kruskal-Wallis results to the format expected by makeVolcano. Extracts
#' the p-value and fold change for a specific comparison and renames them to pvalue, nlog10p,
#' and log2fc. Automatically detects whether input is from performANOVA (tukey_pvalue) or
#' performKW (dunn_pvalue).
#'
#' @param results Data frame output from performANOVA or performKW
#' @param comparison Character vector of length 2 specifying the comparison,
#'   e.g., c("control", "treated"). Order doesn't matter.
#'
#' @return A data frame with added columns:
#'   - pvalue: the Tukey/Dunn p-value for the specified comparison
#'   - nlog10p: -log10 transformed p-value
#'   - log2fc: the fold change for the specified comparison
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
  fc_col1 <- paste0("foldchange_", pair_name1)
  fc_col2 <- paste0("foldchange_", pair_name2)

  # Find which column exists
  pvalue_col <- NULL
  fc_col <- NULL
  flip_sign <- FALSE

  if (tukey_col1 %in% names(results)) {
    pvalue_col <- tukey_col1
    fc_col <- fc_col1
  } else if (tukey_col2 %in% names(results)) {
    pvalue_col <- tukey_col2
    fc_col <- fc_col2
    flip_sign <- TRUE
  } else if (dunn_col1 %in% names(results)) {
    pvalue_col <- dunn_col1
    fc_col <- fc_col1
  } else if (dunn_col2 %in% names(results)) {
    pvalue_col <- dunn_col2
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
  result$log2fc <- if (flip_sign) -result[[fc_col]] else result[[fc_col]]
  result$nlog10p <- -log10(result$pvalue)

  return(result)
}