#' Perform ANOVA on Multiple Cohorts
#'
#' Takes data frames as arguments, each belonging to a cohort, and performs an ANOVA on each row.
#' Works dynamically with 3 or more groups. The function returns a data frame with the p-value
#' from the ANOVA, Tukey's HSD p-values for pairwise comparisons, and fold changes for each
#' pairwise comparison. Column names are automatically derived from the variable names passed in.
#'
#' @param ... Data frames, each containing data for one cohort.
#'   All data frames must have the same number of rows.
#' @param comparisons Which pairwise comparisons to include in output. Options:
#'   - "all" (default): include all pairwise comparisons
#'   - A list of character vectors specifying pairs, e.g.,
#'     list(c("control", "treated"), c("control", "placebo")).
#'     Group names must match the variable names passed in.
#' @param log_scale Logical; whether input data are already on a log scale.
#'   - TRUE (default): fold change is computed as mean(group2) - mean(group1)
#'   - FALSE: fold change is computed as mean(group2) / mean(group1), and log2fc
#'     is computed as log2(fold change)
#'
#' @return A data frame with the combined original data plus computed columns:
#'   - anova_pvalue: overall ANOVA p-value for each row
#'   - tukey_pvalue_X_Y: Tukey's HSD adjusted p-values for each pairwise comparison
#'   - foldchange_X_Y: fold changes for each pairwise comparison (difference if
#'     `log_scale = TRUE`, ratio if `log_scale = FALSE`)
#'   - log2fc_X_Y: log2 fold changes for each pairwise comparison
#'
#' @export
#'
#' @examples
#' # For 3 groups (all pairwise comparisons)
#' result <- performANOVA(control, treated, placebo)
#'
#' # Only specific comparisons
#' result <- performANOVA(control, treated, placebo,
#'                        comparisons = list(c("control", "treated"),
#'                                           c("control", "placebo")))
performANOVA <- function(..., comparisons = "all", log_scale = TRUE) {
  require(dplyr)
  require(stats)

  # Capture the call to extract variable names
  args <- as.list(match.call(expand.dots = TRUE))[-1]
  # Remove 'comparisons' from args if present
  args$comparisons <- NULL
  args$log_scale <- NULL
  arg_names <- as.character(args)

  # Get the actual data frames
  df_list <- list(...)
  # Remove comparisons from df_list if it got included
  if (!is.null(names(df_list))) {
    df_list <- df_list[names(df_list) != "comparisons"]
  }

  n_groups <- length(df_list)
  if (n_groups < 3) stop("Need at least 3 groups for ANOVA")

  # Use variable names as group names
  group_names <- arg_names

  # Number of rows (assuming all data frames have the same number of rows)
  n_rows <- nrow(df_list[[1]])

  # Generate all pairwise combinations (needed for Tukey regardless of output filter)
  all_pairs <- combn(n_groups, 2)
  all_n_pairs <- ncol(all_pairs)
  all_pair_names <- apply(all_pairs, 2, function(p) paste0(group_names[p[1]], "_", group_names[p[2]]))

  # Determine which pairs to include in output
  if (identical(comparisons, "all")) {
    output_pair_indices <- seq_len(all_n_pairs)
  } else {
    # comparisons is a list of pairs like list(c("control", "treated"), ...)
    output_pair_indices <- integer(0)
    for (comp in comparisons) {
      if (length(comp) != 2) stop("Each comparison must be a vector of 2 group names")
      # Find this pair in all_pairs (order doesn't matter)
      found <- FALSE
      for (j in seq_len(all_n_pairs)) {
        g1_name <- group_names[all_pairs[1, j]]
        g2_name <- group_names[all_pairs[2, j]]
        if ((comp[1] == g1_name && comp[2] == g2_name) ||
            (comp[1] == g2_name && comp[2] == g1_name)) {
          output_pair_indices <- c(output_pair_indices, j)
          found <- TRUE
          break
        }
      }
      if (!found) {
        stop(paste0("Comparison not found: ", comp[1], " vs ", comp[2],
                    ". Available groups: ", paste(group_names, collapse = ", ")))
      }
    }
  }

  # Initialize result storage for ALL pairs (Tukey computes all anyway)
  pvalues <- numeric(n_rows)
  tukey_pvalues <- matrix(1, nrow = n_rows, ncol = all_n_pairs)
  foldchanges <- matrix(0, nrow = n_rows, ncol = all_n_pairs)
  log2fcs <- matrix(0, nrow = n_rows, ncol = all_n_pairs)

  # Loop over each row
  for (i in seq_len(n_rows)) {
    tryCatch({
      # Extract values for each group
      values_list <- lapply(df_list, function(df) as.numeric(df[i, ]))

      # Check if all groups have >= 3 non-NA values
      enough_data <- all(sapply(values_list, function(v) sum(!is.na(v)) >= 3))

      if (enough_data) {
        # Combine values and create group factor
        values <- unlist(values_list)
        groups <- factor(rep(group_names, sapply(values_list, length)),
                         levels = group_names)

        # Perform ANOVA
        anova_result <- aov(values ~ groups)
        pvalues[i] <- summary(anova_result)[[1]][["Pr(>F)"]][1]

        if (!is.na(pvalues[i])) {
          tukey_result <- TukeyHSD(anova_result)
          comp <- tukey_result$groups

          # Extract Tukey p-values and fold changes for each pair
          for (j in seq_len(all_n_pairs)) {
            g1 <- all_pairs[1, j]
            g2 <- all_pairs[2, j]
            comp_name <- paste0(group_names[g2], "-", group_names[g1])
            mean_g1 <- mean(values_list[[g1]], na.rm = TRUE)
            mean_g2 <- mean(values_list[[g2]], na.rm = TRUE)

            tukey_pvalues[i, j] <- comp[comp_name, "p adj"]
            if (isTRUE(log_scale)) {
              foldchanges[i, j] <- mean_g2 - mean_g1
              log2fcs[i, j] <- foldchanges[i, j]
            } else {
              if (is.na(mean_g1) || mean_g1 == 0) {
                foldchanges[i, j] <- NA_real_
                log2fcs[i, j] <- NA_real_
              } else {
                foldchanges[i, j] <- mean_g2 / mean_g1
                log2fcs[i, j] <- ifelse(foldchanges[i, j] > 0,
                                        log2(foldchanges[i, j]),
                                        NA_real_)
              }
            }
          }
        }
      } else {
        # If not enough data for ANOVA, set p-value = 1
        pvalues[i] <- 1
        # tukey_pvalues and foldchanges already initialized to 1 and 0
      }
    }, error = function(e) {
      # If an error occurs, use safe defaults (already initialized)
      pvalues[i] <<- 1
    })
  }

  # Build result data frame
  c_df <- dplyr::bind_cols(df_list)
  c_df$anova_pvalue <- pvalues

  # Add Tukey p-values and fold changes ONLY for requested comparisons
  for (j in output_pair_indices) {
    c_df[[paste0("tukey_pvalue_", all_pair_names[j])]] <- tukey_pvalues[, j]
    c_df[[paste0("foldchange_", all_pair_names[j])]] <- foldchanges[, j]
    c_df[[paste0("log2fc_", all_pair_names[j])]] <- log2fcs[, j]
  }

  return(c_df)
}
