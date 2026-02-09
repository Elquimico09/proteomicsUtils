#' Perform Row-wise ANOVA, Extract Specific Tukey HSD Comparisons, and Calculate Log2 Fold Change
#'
#' This function iterates through rows of a combined data matrix, performs
#' ANOVA across groups defined by a grouping factor, runs Tukey's HSD post-hoc
#' test, calculates the log2 fold change (difference of means, assuming log-transformed input)
#' for user-specified comparisons, and extracts adjusted p-values.
#'
#' @param data_list A named list of numeric data frames or matrices. Each element
#'   represents a cohort (e.g., list(Q = Q_df, G = G_df, QZ = QZ_df, ...)).
#'   All data frames must have the same number of rows, corresponding to the
#'   same features (proteins) in the same order. Columns are replicates.
#'   Input data is assumed to be log-transformed (e.g., log2).
#' @param comparisons_list A character vector listing the desired comparisons
#'   in the format "Group1-Group2" (e.g., c("Q-G", "Q-QZ", "QZ-QAE")).
#'   The fold change will be calculated as mean(Group1) - mean(Group2).
#' @param feature_ids An optional vector of feature identifiers (e.g., protein IDs,
#'   gene names) corresponding to the rows of the data frames in `data_list`.
#'   If provided, these will be included in the output data frame. Must have
#'   the same length as the number of rows in the data frames.
#'
#' @return A data frame containing the results. If `feature_ids` were provided,
#'   the first column will be `FeatureID`. Subsequent columns correspond to the
#'   requested comparisons, with pairs of columns for log2 fold change (`_log2FC`)
#'   and adjusted p-value (`_p.adj`). Rows where analysis failed will have NA.
#' @export
#'
#' @examples
#' # --- Create Sample Data (simulate user's dataframes on log2 scale) ---
#' set.seed(42)
#' n_proteins <- 50
#' Q <- as.data.frame(matrix(rnorm(n_proteins * 4, mean = 10, sd = 1), nrow = n_proteins))
#' G <- as.data.frame(matrix(rnorm(n_proteins * 5, mean = 11, sd = 1), nrow = n_proteins)) # G higher mean (log2FC ~1)
#' QZ <- as.data.frame(matrix(rnorm(n_proteins * 4, mean = 10, sd = 1), nrow = n_proteins))
#' QAE <- as.data.frame(matrix(rnorm(n_proteins * 3, mean = 9, sd = 1), nrow = n_proteins)) # QAE lower mean than QZ (log2FC ~ -1)
#' GX <- as.data.frame(matrix(rnorm(n_proteins * 5, mean = 11, sd = 1), nrow = n_proteins)) # Same as G
#' GYM <- as.data.frame(matrix(rnorm(n_proteins * 4, mean = 12, sd = 1), nrow = n_proteins)) # GYM higher mean than GX (log2FC ~ -1 for GX-GYM)
#'
#' # Assign meaningful column names (optional but good practice)
#' colnames(Q) <- paste0("Q_", 1:ncol(Q))
#' colnames(G) <- paste0("G_", 1:ncol(G))
#' colnames(QZ) <- paste0("QZ_", 1:ncol(QZ))
#' colnames(QAE) <- paste0("QAE_", 1:ncol(QAE))
#' colnames(GX) <- paste0("GX_", 1:ncol(GX))
#' colnames(GYM) <- paste0("GYM_", 1:ncol(GYM))
#'
#' # Create the input list
#' data_for_analysis <- list(
#'   Q = Q, G = G, QZ = QZ, QAE = QAE, GX = GX, GYM = GYM
#' )
#'
#' # Define desired comparisons
#' my_comparisons <- c("Q-G", "Q-QZ", "QZ-QAE", "G-GX", "GX-GYM")
#'
#' # Generate dummy protein IDs
#' protein_ids <- paste0("Prot_", 1:n_proteins)
#'
#' # --- Run the function ---
#' analysis_results <- perform_rowwise_anova_tukey_fc(
#'   data_list = data_for_analysis,
#'   comparisons_list = my_comparisons,
#'   feature_ids = protein_ids
#' )
#'
#' # --- View Results ---
#' print(head(analysis_results))
#' # Note columns like Q-G_log2FC and Q-G_p.adj
#' summary(analysis_results)
#'
#' # Example: Filter for significant results with absolute log2FC > 1
#' # library(dplyr)
#' # significant_fc_Q_G <- analysis_results %>%
#' #   filter(`Q-G_p.adj` < 0.05 & abs(`Q-G_log2FC`) > 1)
#' # print(head(significant_fc_Q_G))

perform_rowwise_anova_tukey_fc <- function(data_list, comparisons_list, feature_ids = NULL) {
  
  # --- Input Validation ---
  if (!is.list(data_list) || is.data.frame(data_list)) {
    stop("Error: 'data_list' must be a named list of data frames or matrices.")
  }
  if (length(names(data_list)) != length(data_list) || any(names(data_list) == "")) {
    stop("Error: 'data_list' must be a *named* list.")
  }
  n_rows <- unique(sapply(data_list, nrow))
  if (length(n_rows) != 1) {
    stop("Error: All data frames in 'data_list' must have the same number of rows.")
  }
  # Check if data appears numeric (basic check on first data frame)
  if(!all(sapply(data_list[[1]], is.numeric))) {
    warning("Warning: Input data in the first element of 'data_list' does not appear fully numeric. Ensure all data frames contain only numeric values.")
  }
  if (!is.null(feature_ids) && length(feature_ids) != n_rows) {
    stop("Error: 'feature_ids' length must match the number of rows in the data frames.")
  }
  if(!is.character(comparisons_list)) {
    stop("Error: 'comparisons_list' must be a character vector.")
  }
  # Check comparison format
  if(!all(grepl("^[^-]+-[^-]+$", comparisons_list))) {
    stop("Error: 'comparisons_list' elements must be in the format 'Group1-Group2'.")
  }
  
  
  # --- Combine Data and Create Grouping Factor ---
  group_names <- names(data_list)
  combined_data <- do.call(cbind, data_list)
  # Ensure unique column names, necessary if original dfs had overlapping names
  colnames(combined_data) <- make.unique(unlist(sapply(data_list, colnames, USE.NAMES = FALSE)))
  
  group_vector <- factor(rep(group_names, times = sapply(data_list, ncol)))
  
  # --- Helper function to standardize comparison names (alphabetical order) ---
  # Used only for matching Tukey output, not for FC calculation direction
  standardize_comparison <- function(comp_str) {
    parts <- sort(strsplit(comp_str, "-", fixed = TRUE)[[1]])
    return(paste(parts, collapse = "-"))
  }
  
  # --- Function to process a single row ---
  process_row <- function(row_values, grouping_factor, target_comparisons_original) {
    
    # Create data frame for aov
    df_row <- data.frame(intensity = as.numeric(row_values), group = grouping_factor)
    df_row <- na.omit(df_row) # Remove NAs for this specific row
    
    # Check if enough data remains for ANOVA (at least 2 groups with >1 data point)
    group_counts <- table(df_row$group)
    valid_groups <- names(group_counts[group_counts >= 1]) # Need at least 1 point per group for mean calc
    valid_groups_anova <- names(group_counts[group_counts > 1]) # Need >1 point per group for ANOVA variance
    
    # Initialize result vector (interleaved p-val, fc)
    n_comp <- length(target_comparisons_original)
    result_vector <- rep(NA_real_, n_comp * 2)
    names(result_vector) <- character(n_comp * 2)
    
    # Calculate group means if possible
    group_means <- tapply(df_row$intensity, df_row$group, mean) # Will have NA for groups not present
    
    # --- Calculate Fold Changes ---
    fc_idx <- 2 # Start index for FC results in the interleaved vector
    for (comp in target_comparisons_original) {
      groups_to_compare <- strsplit(comp, "-", fixed = TRUE)[[1]]
      group1 <- groups_to_compare[1]
      group2 <- groups_to_compare[2]
      
      mean1 <- group_means[group1] # Returns NA if group1 not present in df_row
      mean2 <- group_means[group2] # Returns NA if group2 not present in df_row
      
      # Calculate FC only if both means are available
      if (!is.na(mean1) && !is.na(mean2)) {
        result_vector[fc_idx] <- mean1 - mean2
      } # else remains NA
      
      # Set names for FC column
      names(result_vector)[fc_idx] <- paste0(comp, "_log2FC")
      fc_idx <- fc_idx + 2 # Move to next FC slot
    }
    
    
    # --- Perform ANOVA and Tukey HSD (if possible) ---
    p_val_idx <- 1 # Start index for p-value results
    if (length(valid_groups_anova) < 2) {
      # Cannot perform ANOVA/Tukey, fill p-values with NA
      for (comp in target_comparisons_original) {
        names(result_vector)[p_val_idx] <- paste0(comp, "_p.adj")
        # result_vector[p_val_idx] is already NA
        p_val_idx <- p_val_idx + 2
      }
      return(result_vector) # Return vector with NAs for p-values and calculated FCs
    }
    
    # Try ANOVA/Tukey
    tukey_results <- tryCatch({
      aov_res <- aov(intensity ~ group, data = df_row)
      TukeyHSD(aov_res)
    }, error = function(e) {
      NULL # Return NULL if error occurs
    })
    
    # Extract Tukey p-values if successful
    if (!is.null(tukey_results) && "group" %in% names(tukey_results)) {
      tukey_table <- as.data.frame(tukey_results$group)
      tukey_comparisons_std <- sapply(rownames(tukey_table), standardize_comparison, USE.NAMES = FALSE)
      
      for (comp in target_comparisons_original) {
        req_std_comp <- standardize_comparison(comp) # Standardize requested comparison for matching
        match_idx <- which(tukey_comparisons_std == req_std_comp)
        
        if (length(match_idx) == 1) {
          result_vector[p_val_idx] <- tukey_table[match_idx, "p adj"] # Extract adjusted p-value
        } # else remains NA
        
        names(result_vector)[p_val_idx] <- paste0(comp, "_p.adj")
        p_val_idx <- p_val_idx + 2
      }
    } else {
      # ANOVA/Tukey failed, fill p-values with NA
      for (comp in target_comparisons_original) {
        names(result_vector)[p_val_idx] <- paste0(comp, "_p.adj")
        # result_vector[p_val_idx] is already NA
        p_val_idx <- p_val_idx + 2
      }
    }
    
    return(result_vector) # Return interleaved vector
    
  } # End process_row function
  
  # --- Apply process_row to each row ---
  results_matrix <- t(apply(combined_data, 1, process_row,
                            grouping_factor = group_vector,
                            target_comparisons_original = comparisons_list))
  results_df <- as.data.frame(results_matrix)
  
  # --- Add Feature IDs if provided ---
  if (!is.null(feature_ids)) {
    results_df <- cbind(FeatureID = feature_ids, results_df)
  }
  
  return(results_df)
  
} # End perform_rowwise_anova_tukey_fc function
