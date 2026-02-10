#' Perform Mann-Whitney U Test on Two Cohorts
#'
#' Takes two dataframes and performs a non-parametric Mann-Whitney U test on each row.
#' This is the non-parametric alternative to the t-test for comparing two groups.
#' Computes p-values, fold changes, and log2 fold changes with optional p-value adjustment.
#'
#' @param df1 Dataframe containing the first cohort. Must have the same number of rows as df2.
#' @param df2 Dataframe containing the second cohort. Must have the same number of rows as df1.
#' @param adjust Logical, whether to adjust p-values for multiple testing (default: FALSE).
#' @param p_adjust Method for p-value adjustment if adjust = TRUE. Options: "BH" (default),
#'   "bonferroni", "holm", "hochberg", "hommel", "BY", "none".
#' @param log_scale Logical, whether data is already log-transformed (default: TRUE).
#'   If TRUE: log2fc = mean(df2) - mean(df1) (difference of log values).
#'   If FALSE: foldchange = mean(df2) / mean(df1), log2fc = log2(foldchange).
#'
#' @return A dataframe containing the combined original data plus computed columns:
#'   - pvalue: raw p-value from Mann-Whitney U test
#'   - adjusted_pvalue: adjusted p-value (if adjust = TRUE)
#'   - log2fc: log2 fold change
#'   - foldchange: fold change (if log_scale = FALSE)
#'
#' @export
#'
#' @examples
#' # Basic Mann-Whitney test without adjustment
#' result <- performMWTest(control_data, treatment_data)
#'
#' # Mann-Whitney test with BH adjustment
#' result <- performMWTest(control_data, treatment_data, adjust = TRUE)
#'
#' # Test on non-log data
#' result <- performMWTest(control_data, treatment_data, log_scale = FALSE)
performMWTest <- function(df1, df2, adjust = FALSE, p_adjust = "BH", log_scale = TRUE) {
  # Ensure dataframes have the same number of rows
  if (nrow(df1) != nrow(df2)) {
    stop("df1 and df2 must have the same number of rows")
  }

  # Warn if adjust = TRUE but p_adjust not explicitly specified
  if (adjust && missing(p_adjust)) {
    warning("Adjustment method not specified, defaulting to BH")
  }

  # Print log scale status
  if (log_scale) {
    message("log_scale = TRUE: Data assumed to be log-transformed. Fold change will NOT be calculated.")
  } else {
    message("log_scale = FALSE: Data assumed to be linear scale. Fold change WILL be calculated.")
  }

  # Initialize result vectors
  pvalues <- numeric(nrow(df1))
  foldchange <- numeric(nrow(df1))
  log2fc <- numeric(nrow(df1))

  for (i in 1:nrow(df1)) {
    values_df1 <- as.numeric(df1[i, ])
    values_df1 <- values_df1[!is.na(values_df1)]

    values_df2 <- as.numeric(df2[i, ])
    values_df2 <- values_df2[!is.na(values_df2)]

    # Check if both rows have at least 3 non-NA values
    if (length(values_df1) >= 3 && length(values_df2) >= 3) {
      test_result <- tryCatch(
        wilcox.test(values_df1, values_df2),
        error = function(e) list(p.value = 1)
      )

      pvalues[i] <- test_result$p.value

      mean1 <- mean(values_df1, na.rm = TRUE)
      mean2 <- mean(values_df2, na.rm = TRUE)

      if (log_scale) {
        # Data is log-transformed: difference of means = log2 fold change
        log2fc[i] <- mean2 - mean1
      } else {
        # Data is linear scale: calculate ratio then take log2
        if (mean1 != 0) {
          foldchange[i] <- mean2 / mean1
          log2fc[i] <- log2(foldchange[i])
        } else {
          foldchange[i] <- NA
          log2fc[i] <- NA
        }
      }
    } else {
      pvalues[i] <- 1
      foldchange[i] <- ifelse(log_scale, NA, 0)
      log2fc[i] <- 0
    }
  }

  # Adjust p-values if requested
  if (adjust) {
    pvalues <- p.adjust(pvalues, method = p_adjust)
  }

  nlog10p <- -log10(pvalues)

  # Combine results with original data
  c_df <- cbind(df1, df2)
  c_df$pvalue <- pvalues
  if (!log_scale) {
    c_df$foldchange <- foldchange
  }
  c_df$nlog10p <- nlog10p
  c_df$log2fc <- log2fc

  return(c_df)
}