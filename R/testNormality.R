#' Test Normality of Data Distribution
#'
#' Tests whether data is normally distributed across rows (e.g., proteins, genes) using
#' the Shapiro-Wilk test. Helps decide whether to use parametric (t-test, ANOVA) or
#' non-parametric (Mann-Whitney, Kruskal-Wallis) statistical tests. Provides a
#' recommendation based on the proportion of rows passing the normality test.
#'
#' @param ... Data frames, each containing data for one cohort. All data frames must
#'   have the same number of rows.
#' @param alpha Numeric, significance level for normality test (default: 0.05).
#'   Rows with Shapiro-Wilk p-value > alpha are considered normally distributed.
#' @param threshold Numeric, proportion of rows that must pass normality to recommend
#'   parametric tests (default: 0.8, i.e., 80%).
#'
#' @return Invisibly returns a list containing:
#'   - group_results: detailed statistics for each group (total_rows, valid_tests,
#'     normal_count, non_normal_count, normal_proportion, pvalues)
#'   - overall_normal_proportion: mean proportion of normal rows across all groups
#'   - recommendation: "PARAMETRIC" or "NON-PARAMETRIC"
#'   - alpha: the significance level used
#'   - threshold: the threshold used for recommendation
#'
#' @export
#'
#' @examples
#' # Test normality for two groups
#' testNormality(control_data, treatment_data)
#'
#' # Test with custom threshold
#' testNormality(control_data, treatment_data, alpha = 0.01, threshold = 0.9)
testNormality <- function(..., alpha = 0.05, threshold = 0.8) {
  # Capture group names
  args <- as.list(match.call(expand.dots = TRUE))[-1]
  args$alpha <- NULL
  args$threshold <- NULL
  group_names <- as.character(args)

  df_list <- list(...)
  if (!is.null(names(df_list))) {
    df_list <- df_list[!names(df_list) %in% c("alpha", "threshold")]
  }

  n_groups <- length(df_list)
  n_rows <- nrow(df_list[[1]])

  # Store results for each group
  group_results <- list()

  for (g in seq_len(n_groups)) {
    df <- df_list[[g]]
    shapiro_pvalues <- numeric(n_rows)

    for (i in seq_len(n_rows)) {
      values <- as.numeric(df[i, ])
      values <- values[!is.na(values)]

      # Shapiro-Wilk requires 3-5000 observations
      if (length(values) >= 3 && length(values) <= 5000) {
        test_result <- tryCatch(
          shapiro.test(values),
          error = function(e) list(p.value = NA)
        )
        shapiro_pvalues[i] <- test_result$p.value
      } else {
        shapiro_pvalues[i] <- NA
      }
    }

    # Calculate statistics
    valid_tests <- sum(!is.na(shapiro_pvalues))
    normal_count <- sum(shapiro_pvalues > alpha, na.rm = TRUE)
    normal_proportion <- normal_count / valid_tests

    group_results[[group_names[g]]] <- list(
      total_rows = n_rows,
      valid_tests = valid_tests,
      normal_count = normal_count,
      non_normal_count = valid_tests - normal_count,
      normal_proportion = normal_proportion,
      pvalues = shapiro_pvalues
    )
  }

  # Calculate overall statistics
  all_normal_props <- sapply(group_results, function(x) x$normal_proportion)
  overall_normal_prop <- mean(all_normal_props, na.rm = TRUE)

  # Generate recommendation
  if (overall_normal_prop >= threshold) {
    recommendation <- "PARAMETRIC"
    rec_message <- paste0(
      "Recommendation: Use PARAMETRIC tests (t-test, ANOVA)\n",
      sprintf("%.1f%% of rows pass normality test (threshold: %.0f%%)",
              overall_normal_prop * 100, threshold * 100)
    )
  } else {
    recommendation <- "NON-PARAMETRIC"
    rec_message <- paste0(
      "Recommendation: Use NON-PARAMETRIC tests (Mann-Whitney, Kruskal-Wallis)\n",
      sprintf("Only %.1f%% of rows pass normality test (threshold: %.0f%%)",
              overall_normal_prop * 100, threshold * 100)
    )
  }

  # Print summary
  message("\n===== Normality Test Summary =====\n")
  for (name in names(group_results)) {
    res <- group_results[[name]]
    message(sprintf("%s: %.1f%% normal (%d/%d rows)",
                    name, res$normal_proportion * 100,
                    res$normal_count, res$valid_tests))
  }
  message(sprintf("\nOverall: %.1f%% of rows pass normality", overall_normal_prop * 100))
  message(paste0("\n", rec_message, "\n"))

  # Return results invisibly
  invisible(list(
    group_results = group_results,
    overall_normal_proportion = overall_normal_prop,
    recommendation = recommendation,
    alpha = alpha,
    threshold = threshold
  ))
}
