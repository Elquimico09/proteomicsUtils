#' Perform Kruskal-Wallis Test on Multiple Cohorts
#'
#' Takes data frames as arguments, each belonging to a cohort, and performs a non-parametric
#' Kruskal-Wallis test on each row. This is the non-parametric alternative to ANOVA for comparing
#' three or more groups. Works dynamically with 3 or more groups. The function returns a data frame
#' with the p-value from the Kruskal-Wallis test, Dunn's test p-values for pairwise comparisons,
#' and fold changes for each pairwise comparison. Column names are automatically derived from
#' the variable names passed in.
#'
#' @param ... Data frames, each containing data for one cohort.
#'   All data frames must have the same number of rows.
#' @param comparisons Which pairwise comparisons to include in output. Options:
#'   - "all" (default): include all pairwise comparisons
#'   - A list of character vectors specifying pairs, e.g.,
#'     list(c("control", "treated"), c("control", "placebo")).
#'     Group names must match the variable names passed in.
#' @param p_adjust Method for p-value adjustment in Dunn's test. Default is "BH" (Benjamini-Hochberg).
#'   Other options: "bonferroni", "holm", "hochberg", "hommel", "BY", "none".
#' @param log_scale Logical; whether input data are already on a log scale.
#'   - TRUE (default): fold change is computed as mean(group2) - mean(group1)
#'   - FALSE: fold change is computed as mean(group2) / mean(group1), and log2fc
#'     is computed as log2(fold change)
#'
#' @return A data frame with the combined original data plus computed columns:
#'   - kw_pvalue: overall Kruskal-Wallis p-value for each row
#'   - dunn_pvalue_X_Y: Dunn's test adjusted p-values for each pairwise comparison
#'   - foldchange_X_Y: fold changes for each pairwise comparison (difference if
#'     `log_scale = TRUE`, ratio if `log_scale = FALSE`)
#'   - log2fc_X_Y: log2 fold changes for each pairwise comparison
#'
#' @export
#'
#' @examples
#' # For 3 groups (all pairwise comparisons)
#' result <- performKW(control, treated, placebo)
#'
#' # Only specific comparisons
#' result <- performKW(control, treated, placebo,
#'                     comparisons = list(c("control", "treated"),
#'                                        c("control", "placebo")))
performKW <- function(..., comparisons = "all", p_adjust = "BH", log_scale = TRUE) {
  require(dplyr)
  require(stats)

  # Check if dunn.test is installed
  if (!requireNamespace("dunn.test", quietly = TRUE)) {
    stop("Package 'dunn.test' is required for Dunn's post-hoc test.\n",
         "Please install it with: install.packages('dunn.test')")
  }

  # Capture the call to extract variable names
  args <- as.list(match.call(expand.dots = TRUE))[-1]
  # Remove named arguments from args
  args$comparisons <- NULL
  args$p_adjust <- NULL
  args$log_scale <- NULL
  arg_names <- as.character(args)

  # Get the actual data frames
  df_list <- list(...)
  # Remove named arguments from df_list if they got included
  if (!is.null(names(df_list))) {
    df_list <- df_list[!names(df_list) %in% c("comparisons", "p_adjust")]
  }

  n_groups <- length(df_list)
  if (n_groups < 3) stop("Need at least 3 groups for Kruskal-Wallis")

  # Use variable names as group names
  group_names <- arg_names

  # Number of rows (assuming all data frames have the same number of rows)
  n_rows <- nrow(df_list[[1]])

  # Generate all pairwise combinations
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
  pvalues <- numeric(n_rows)
  dunn_pvalues <- matrix(1, nrow = n_rows, ncol = all_n_pairs)
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

        # Remove NAs
        valid_idx <- !is.na(values)
        values <- values[valid_idx]
        groups <- groups[valid_idx]

        # Perform Kruskal-Wallis test
        kw_result <- kruskal.test(values ~ groups)
        pvalues[i] <- kw_result$p.value

        if (!is.na(pvalues[i])) {
          # Perform Dunn's test (suppress all printed output)
          sink(tempfile())
          dunn_result <- tryCatch(
            dunn.test::dunn.test(values, groups, method = p_adjust, kw = FALSE, table = FALSE, list = FALSE),
            finally = sink()
          )

          # Extract Dunn's p-values for each pair
          # dunn.test returns comparisons in format "group1 - group2"
          for (j in seq_len(all_n_pairs)) {
            g1 <- all_pairs[1, j]
            g2 <- all_pairs[2, j]

            # Try both orderings to find the comparison
            comp_name1 <- paste0(group_names[g1], " - ", group_names[g2])
            comp_name2 <- paste0(group_names[g2], " - ", group_names[g1])

            idx <- which(dunn_result$comparisons == comp_name1 | dunn_result$comparisons == comp_name2)

            if (length(idx) > 0) {
              dunn_pvalues[i, j] <- dunn_result$P.adjusted[idx[1]]
            }

            mean_g1 <- mean(values_list[[g1]], na.rm = TRUE)
            mean_g2 <- mean(values_list[[g2]], na.rm = TRUE)

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
        pvalues[i] <- 1
      }
    }, error = function(e) {
      pvalues[i] <<- 1
    })
  }

  # Build result data frame
  c_df <- dplyr::bind_cols(df_list)
  c_df$kw_pvalue <- pvalues

  # Add Dunn's p-values and fold changes ONLY for requested comparisons
  for (j in output_pair_indices) {
    c_df[[paste0("dunn_pvalue_", all_pair_names[j])]] <- dunn_pvalues[, j]
    c_df[[paste0("foldchange_", all_pair_names[j])]] <- foldchanges[, j]
    c_df[[paste0("log2fc_", all_pair_names[j])]] <- log2fcs[, j]
  }

  return(c_df)
}
