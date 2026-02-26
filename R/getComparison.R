#' getComparison
#' @description
#' Extracts cohort-specific columns and comparison metrics from `performANOVA()`
#' or `performKW()` output. Keeps only the original columns for the two requested cohorts, then
#' adds:
#' - `pvalue` (from `tukey_pvalue_*` or `dunn_pvalue_*`)
#' - `fold change` (from `foldchange_*`)
#'
#' @param results List output from `performANOVA()` or `performKW()`
#'   containing `data` and `cohort_columns`.
#' @param comparison Character vector of length 2 specifying the comparison,
#'   e.g., `c("control", "treated")`. Order does not matter for matching.
#'
#' @return A data frame with:
#' - Original columns from only the two requested cohorts
#' - `pvalue`
#' - `fold change`
#' @export
#'
#' @examples
#' anova_res <- performANOVA(control, treated, placebo)
#' cmp1 <- getComparison(anova_res, c("control", "treated"))
#' kw_res <- performKW(control, treated, placebo)
#' cmp2 <- getComparison(kw_res, c("control", "treated"))
getComparison <- function(results, comparison) {
  if (!is.list(results) || !all(c("data", "cohort_columns") %in% names(results))) {
    stop("results must be the list output from performANOVA()/performKW() with `data` and `cohort_columns`.")
  }

  if (length(comparison) != 2) {
    stop("comparison must be a character vector of length 2")
  }

  data <- results$data
  cohort_columns <- results$cohort_columns

  missing_groups <- setdiff(comparison, names(cohort_columns))
  if (length(missing_groups) > 0) {
    stop(paste0(
      "Comparison groups not found in cohort_columns: ",
      paste(missing_groups, collapse = ", "),
      ". Available groups: ",
      paste(names(cohort_columns), collapse = ", ")
    ))
  }

  pair_name1 <- paste0(comparison[1], "_", comparison[2])
  pair_name2 <- paste0(comparison[2], "_", comparison[1])

  tukey_col1 <- paste0("tukey_pvalue_", pair_name1)
  tukey_col2 <- paste0("tukey_pvalue_", pair_name2)
  dunn_col1 <- paste0("dunn_pvalue_", pair_name1)
  dunn_col2 <- paste0("dunn_pvalue_", pair_name2)
  fc_col1 <- paste0("foldchange_", pair_name1)
  fc_col2 <- paste0("foldchange_", pair_name2)

  pvalue_col <- NULL
  fc_col <- NULL

  if (tukey_col1 %in% names(data) && fc_col1 %in% names(data)) {
    pvalue_col <- tukey_col1
    fc_col <- fc_col1
  } else if (tukey_col2 %in% names(data) && fc_col2 %in% names(data)) {
    pvalue_col <- tukey_col2
    fc_col <- fc_col2
  } else if (dunn_col1 %in% names(data) && fc_col1 %in% names(data)) {
    pvalue_col <- dunn_col1
    fc_col <- fc_col1
  } else if (dunn_col2 %in% names(data) && fc_col2 %in% names(data)) {
    pvalue_col <- dunn_col2
    fc_col <- fc_col2
  } else {
    available <- grep("^(tukey_pvalue_|dunn_pvalue_)", names(data), value = TRUE)
    available <- gsub("^(tukey_pvalue_|dunn_pvalue_)", "", available)
    stop(paste0(
      "Comparison '", comparison[1], "' vs '", comparison[2], "' not found.\n",
      "Available comparisons: ", paste(available, collapse = ", ")
    ))
  }

  cohort_cols <- unique(c(
    cohort_columns[[comparison[1]]],
    cohort_columns[[comparison[2]]]
  ))
  cohort_cols <- cohort_cols[cohort_cols %in% names(data)]

  out <- data[, cohort_cols, drop = FALSE]
  out$pvalue <- data[[pvalue_col]]
  out[["fold change"]] <- data[[fc_col]]

  return(out)
}
