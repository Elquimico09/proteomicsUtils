#' Perform 2-Way ANOVA with Automatic Grouping
#'
#' Performs 2-way ANOVA on each row of a dataframe with automatic grouping based on two factors.
#' Takes a single dataframe where columns represent different combinations of two factors.
#' For example, if you have male control, female control, male disease, and female disease samples,
#' you can specify factor1 = c("Male", "Female") and factor2 = c("control", "disease"), and the
#' function will automatically group the columns accordingly.
#'
#' @param df A dataframe where each row is a variable and columns represent different samples.
#' @param factor1 A character vector of levels for the first factor (e.g., c("Male", "Female")).
#' @param factor2 A character vector of levels for the second factor (e.g., c("control", "disease")).
#' @param comparisons Either "all" for all pairwise comparisons, or a list of specific comparisons
#'   for Tukey's post-hoc test. Each comparison should be a vector of 2 group names formed by
#'   combining factor levels with a space (e.g., "Male control").
#'   Example: list(c("Male control", "Male disease"), c("Female control", "Female disease"))
#' @param group_sizes Optional integer vector of replicate counts per group (in the same
#'   order as factor1 x factor2). Use this to support uneven replicates when columns are
#'   ordered by group. If NULL (default), equal replicates per group are assumed.
#'
#' @return A dataframe with ANOVA results, including main effects, interaction, and Tukey post-hoc tests:
#'   - anova2way_pvalue_overall: minimum p-value across all effects
#'   - anova2way_pvalue_factor1: main effect p-value for first factor
#'   - anova2way_pvalue_factor2: main effect p-value for second factor
#'   - anova2way_pvalue_interaction: interaction effect p-value
#'   - tukey_pvalue_X_Y: Tukey's HSD p-values for specified comparisons
#'   - foldchange_X_Y: fold changes (mean differences) for specified comparisons
#'
#' @export
#'
#' @examples
#' # Example with male/female and control/disease
#' result <- perform2WayANOVA(data,
#'                            factor1 = c("Male", "Female"),
#'                            factor2 = c("control", "disease"),
#'                            comparisons = list(c("Male control", "Male disease"),
#'                                             c("Female control", "Female disease")))
perform2WayANOVA <- function(df, factor1, factor2, comparisons = "all", group_sizes = NULL) {
  require(dplyr)
  require(stats)

  # Create all combinations of the two factors
  combinations <- expand.grid(factor1 = factor1, factor2 = factor2, stringsAsFactors = FALSE)
  group_names <- paste(combinations$factor1, combinations$factor2)
  n_groups <- length(group_names)

  if (n_groups < 2) stop("Need at least 2 groups for 2-way ANOVA")

  # Verify that we have the right number of columns in the dataframe
  # Assuming columns are ordered by group, allow uneven replicates via group_sizes
  n_cols <- ncol(df)
  if (is.null(group_sizes)) {
    if (n_cols %% n_groups != 0) {
      stop("Uneven replicates detected. Provide group_sizes to specify counts per group.")
    }
    n_replicates <- n_cols / n_groups
    group_sizes <- rep(n_replicates, n_groups)
  } else {
    if (length(group_sizes) != n_groups) {
      stop("group_sizes length must equal number of groups (length(factor1) * length(factor2))")
    }
    if (sum(group_sizes) != n_cols) {
      stop("Sum of group_sizes must equal ncol(df)")
    }
    if (any(group_sizes < 1) || any(group_sizes %% 1 != 0)) {
      stop("group_sizes must contain positive integers")
    }
  }
  n_rows <- nrow(df)

  # Create group assignment vector (assumes columns are ordered by group)
  group_vector <- rep(group_names, times = group_sizes)
  factor1_vector <- rep(combinations$factor1, times = group_sizes)
  factor2_vector <- rep(combinations$factor2, times = group_sizes)

  # Generate all pairwise combinations for Tukey
  all_pairs <- combn(n_groups, 2)
  all_n_pairs <- ncol(all_pairs)
  all_pair_names <- apply(all_pairs, 2, function(p) paste0(group_names[p[1]], "_", group_names[p[2]]))

  # Determine which pairs to include in output
  if (identical(comparisons, "all")) {
    output_pair_indices <- seq_len(all_n_pairs)
  } else {
    output_pair_indices <- integer(0)
    for (comp in comparisons) {
      if (length(comp) != 2) stop("Each comparison must be a vector of 2 group names")
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

  # Initialize result storage
  pvalues_factor1 <- numeric(n_rows)
  pvalues_factor2 <- numeric(n_rows)
  pvalues_interaction <- numeric(n_rows)
  pvalues_overall <- numeric(n_rows)
  tukey_pvalues <- matrix(1, nrow = n_rows, ncol = all_n_pairs)
  foldchanges <- matrix(0, nrow = n_rows, ncol = all_n_pairs)

  # Loop over each row
  for (i in seq_len(n_rows)) {
    tryCatch({
      # Extract values for this row
      values <- as.numeric(df[i, ])

      # Check if we have enough non-NA values per group
      valid_counts <- sapply(group_names, function(g) {
        sum(!is.na(values[group_vector == g]))
      })

      if (all(valid_counts >= 2)) {  # At least 2 replicates per group minimum
        # Create dataframe for ANOVA
        anova_data <- data.frame(
          value = values,
          factor1 = factor(factor1_vector, levels = factor1),
          factor2 = factor(factor2_vector, levels = factor2),
          group = factor(group_vector, levels = group_names)
        )

        # Remove NA values
        anova_data <- anova_data[!is.na(anova_data$value), ]

        # Perform 2-way ANOVA
        anova_result <- aov(value ~ factor1 * factor2, data = anova_data)
        anova_summary <- summary(anova_result)[[1]]

        # Extract p-values for main effects and interaction
        pvalues_factor1[i] <- anova_summary[["Pr(>F)"]][1]
        pvalues_factor2[i] <- anova_summary[["Pr(>F)"]][2]
        pvalues_interaction[i] <- anova_summary[["Pr(>F)"]][3]

        # Overall model p-value (minimum of the three)
        pvalues_overall[i] <- min(pvalues_factor1[i], pvalues_factor2[i],
                                  pvalues_interaction[i], na.rm = TRUE)

        # Perform Tukey's HSD on the group factor for post-hoc comparisons
        anova_group <- aov(value ~ group, data = anova_data)

        if (!is.na(pvalues_overall[i])) {
          tukey_result <- TukeyHSD(anova_group)
          comp <- tukey_result$group

          # Calculate group means for fold changes
          group_means <- sapply(seq_len(n_groups), function(g) {
            group_data <- values[group_vector == group_names[g]]
            mean(group_data, na.rm = TRUE)
          })

          # Extract Tukey p-values and fold changes for each pair
          for (j in seq_len(all_n_pairs)) {
            g1 <- all_pairs[1, j]
            g2 <- all_pairs[2, j]
            comp_name <- paste0(group_names[g2], "-", group_names[g1])

            tukey_pvalues[i, j] <- comp[comp_name, "p adj"]
            foldchanges[i, j] <- group_means[g2] - group_means[g1]
          }
        }
      } else {
        # Not enough data
        pvalues_factor1[i] <- 1
        pvalues_factor2[i] <- 1
        pvalues_interaction[i] <- 1
        pvalues_overall[i] <- 1
      }
    }, error = function(e) {
      # If an error occurs, use safe defaults
      pvalues_factor1[i] <<- 1
      pvalues_factor2[i] <<- 1
      pvalues_interaction[i] <<- 1
      pvalues_overall[i] <<- 1
    })
  }

  # Build result data frame
  result_df <- df
  result_df$anova2way_pvalue_overall <- pvalues_overall
  result_df[[paste0("anova2way_pvalue_", deparse(substitute(factor1)))]] <- pvalues_factor1
  result_df[[paste0("anova2way_pvalue_", deparse(substitute(factor2)))]] <- pvalues_factor2
  result_df$anova2way_pvalue_interaction <- pvalues_interaction

  # Add Tukey p-values and fold changes ONLY for requested comparisons
  for (j in output_pair_indices) {
    result_df[[paste0("tukey_pvalue_", all_pair_names[j])]] <- tukey_pvalues[, j]
    result_df[[paste0("foldchange_", all_pair_names[j])]] <- foldchanges[, j]
  }

  return(result_df)
}