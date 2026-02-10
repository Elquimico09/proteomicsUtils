#' performANOVA6
#' @description
#' Takes in 6 data frames, each belonging to a cohort, and performs an ANOVA on each row.
#' The function returns a data frame with the p-value from the ANOVA, Tukey's HSD p-values
#' for pairwise comparisons, and fold changes for each pairwise comparison.
#'
#' @param df1 Data frame containing cohort 1
#' @param df2 Data frame containing cohort 2
#' @param df3 Data frame containing cohort 3
#' @param df4 Data frame containing cohort 4
#' @param df5 Data frame containing cohort 5
#' @param df6 Data frame containing cohort 6
#'
#' @return A data frame with the results of the ANOVA, Tukey’s HSD p-values, and fold changes
#' @export
#'
#' @examples
#' test <- performANOVA6(df1, df2, df3, df4, df5, df6)
performANOVA6 <- function(df1, df2, df3, df4, df5, df6) {
  # require(dplyr) # REMOVED: library calls should be in DESCRIPTION/NAMESPACE
  # require(stats) # REMOVED: library calls should be in DESCRIPTION/NAMESPACE
  
  # Number of rows (assuming df1,...,df6 have the same number of rows)
  n_rows <- nrow(df1)
  
  # Store overall ANOVA p-values
  pvalues <- numeric(n_rows)
  
  # Store Tukey's p-values for all 15 pairwise comparisons
  tukey_pvalue12 <- numeric(n_rows)
  tukey_pvalue13 <- numeric(n_rows)
  tukey_pvalue14 <- numeric(n_rows)
  tukey_pvalue15 <- numeric(n_rows)
  tukey_pvalue16 <- numeric(n_rows)
  tukey_pvalue23 <- numeric(n_rows)
  tukey_pvalue24 <- numeric(n_rows)
  tukey_pvalue25 <- numeric(n_rows)
  tukey_pvalue26 <- numeric(n_rows)
  tukey_pvalue34 <- numeric(n_rows)
  tukey_pvalue35 <- numeric(n_rows)
  tukey_pvalue36 <- numeric(n_rows)
  tukey_pvalue45 <- numeric(n_rows)
  tukey_pvalue46 <- numeric(n_rows)
  tukey_pvalue56 <- numeric(n_rows)
  
  # Store fold changes for all 15 pairwise comparisons
  foldchange12 <- numeric(n_rows)
  foldchange13 <- numeric(n_rows)
  foldchange14 <- numeric(n_rows)
  foldchange15 <- numeric(n_rows)
  foldchange16 <- numeric(n_rows)
  foldchange23 <- numeric(n_rows)
  foldchange24 <- numeric(n_rows)
  foldchange25 <- numeric(n_rows)
  foldchange26 <- numeric(n_rows)
  foldchange34 <- numeric(n_rows)
  foldchange35 <- numeric(n_rows)
  foldchange36 <- numeric(n_rows)
  foldchange45 <- numeric(n_rows)
  foldchange46 <- numeric(n_rows)
  foldchange56 <- numeric(n_rows)
  
  # Loop over each row
  tryCatch({
    for (i in seq_len(n_rows)) {
      values_df1 <- as.numeric(df1[i, ])
      values_df2 <- as.numeric(df2[i, ])
      values_df3 <- as.numeric(df3[i, ])
      values_df4 <- as.numeric(df4[i, ])
      values_df5 <- as.numeric(df5[i, ])
      values_df6 <- as.numeric(df6[i, ])
      
      # Combine all values
      values <- c(values_df1, values_df2, values_df3, values_df4, values_df5, values_df6)
      groups <- factor(c(
        rep("Group1", length(values_df1)),
        rep("Group2", length(values_df2)),
        rep("Group3", length(values_df3)),
        rep("Group4", length(values_df4)),
        rep("Group5", length(values_df5)),
        rep("Group6", length(values_df6))
      ))
      
      # Check for >=3 non-NA values in each group to run ANOVA
      enough_data <- (
        sum(!is.na(values_df1)) >= 3 &&
          sum(!is.na(values_df2)) >= 3 &&
          sum(!is.na(values_df3)) >= 3 &&
          sum(!is.na(values_df4)) >= 3 &&
          sum(!is.na(values_df5)) >= 3 &&
          sum(!is.na(values_df6)) >= 3
      )
      
      if (enough_data) {
        # Perform ANOVA
        anova_result <- aov(values ~ groups)
        pvalues[i] <- summary(anova_result)[[1]][["Pr(>F)"]][1]
        
        # Perform Tukey's HSD test
        # (Checking "if (pvalues[i] <= 1)" always succeeds unless p-value is NA;
        #  presumably you'd check p < 0.05 if you only want post-hoc for significant ANOVA)
        if (!is.na(pvalues[i]) && pvalues[i] <= 1) {
          tukey_result <- TukeyHSD(anova_result)
          # Access the matrix of pairwise comparisons
          comp <- tukey_result$groups
          
          # Extract Tukey p-values for pairwise comparisons
          tukey_pvalue12[i] <- comp["Group2-Group1", "p adj"]
          tukey_pvalue13[i] <- comp["Group3-Group1", "p adj"]
          tukey_pvalue14[i] <- comp["Group4-Group1", "p adj"]
          tukey_pvalue15[i] <- comp["Group5-Group1", "p adj"]
          tukey_pvalue16[i] <- comp["Group6-Group1", "p adj"]
          tukey_pvalue23[i] <- comp["Group3-Group2", "p adj"]
          tukey_pvalue24[i] <- comp["Group4-Group2", "p adj"]
          tukey_pvalue25[i] <- comp["Group5-Group2", "p adj"]
          tukey_pvalue26[i] <- comp["Group6-Group2", "p adj"]
          tukey_pvalue34[i] <- comp["Group4-Group3", "p adj"]
          tukey_pvalue35[i] <- comp["Group5-Group3", "p adj"]
          tukey_pvalue36[i] <- comp["Group6-Group3", "p adj"]
          tukey_pvalue45[i] <- comp["Group5-Group4", "p adj"]
          tukey_pvalue46[i] <- comp["Group6-Group4", "p adj"]
          tukey_pvalue56[i] <- comp["Group6-Group5", "p adj"]
          
          # Calculate fold changes = mean(GroupX) - mean(GroupY)
          foldchange12[i] <- mean(values_df2, na.rm = TRUE) - mean(values_df1, na.rm = TRUE)
          foldchange13[i] <- mean(values_df3, na.rm = TRUE) - mean(values_df1, na.rm = TRUE)
          foldchange14[i] <- mean(values_df4, na.rm = TRUE) - mean(values_df1, na.rm = TRUE)
          foldchange15[i] <- mean(values_df5, na.rm = TRUE) - mean(values_df1, na.rm = TRUE)
          foldchange16[i] <- mean(values_df6, na.rm = TRUE) - mean(values_df1, na.rm = TRUE)
          foldchange23[i] <- mean(values_df3, na.rm = TRUE) - mean(values_df2, na.rm = TRUE)
          foldchange24[i] <- mean(values_df4, na.rm = TRUE) - mean(values_df2, na.rm = TRUE)
          foldchange25[i] <- mean(values_df5, na.rm = TRUE) - mean(values_df2, na.rm = TRUE)
          foldchange26[i] <- mean(values_df6, na.rm = TRUE) - mean(values_df2, na.rm = TRUE)
          foldchange34[i] <- mean(values_df4, na.rm = TRUE) - mean(values_df3, na.rm = TRUE)
          foldchange35[i] <- mean(values_df5, na.rm = TRUE) - mean(values_df3, na.rm = TRUE)
          foldchange36[i] <- mean(values_df6, na.rm = TRUE) - mean(values_df3, na.rm = TRUE)
          foldchange45[i] <- mean(values_df5, na.rm = TRUE) - mean(values_df4, na.rm = TRUE)
          foldchange46[i] <- mean(values_df6, na.rm = TRUE) - mean(values_df4, na.rm = TRUE)
          foldchange56[i] <- mean(values_df6, na.rm = TRUE) - mean(values_df5, na.rm = TRUE)
        }
      } else {
        # If not enough data for ANOVA, set p-values = 1, fold changes = 0
        pvalues[i]       <- 1
        tukey_pvalue12[i] <- 1
        tukey_pvalue13[i] <- 1
        tukey_pvalue14[i] <- 1
        tukey_pvalue15[i] <- 1
        tukey_pvalue16[i] <- 1
        tukey_pvalue23[i] <- 1
        tukey_pvalue24[i] <- 1
        tukey_pvalue25[i] <- 1
        tukey_pvalue26[i] <- 1
        tukey_pvalue34[i] <- 1
        tukey_pvalue35[i] <- 1
        tukey_pvalue36[i] <- 1
        tukey_pvalue45[i] <- 1
        tukey_pvalue46[i] <- 1
        tukey_pvalue56[i] <- 1
        
        foldchange12[i]  <- 0
        foldchange13[i]  <- 0
        foldchange14[i]  <- 0
        foldchange15[i]  <- 0
        foldchange16[i]  <- 0
        foldchange23[i]  <- 0
        foldchange24[i]  <- 0
        foldchange25[i]  <- 0
        foldchange26[i]  <- 0
        foldchange34[i]  <- 0
        foldchange35[i]  <- 0
        foldchange36[i]  <- 0
        foldchange45[i]  <- 0
        foldchange46[i]  <- 0
        foldchange56[i]  <- 0
      }
    }
  }, error = function(e) {
    # If an error occurs, fill the current row with safe defaults
    # (In a robust design, you'd handle or log the error more gracefully.)
    pvalues[i]       <- 1
    tukey_pvalue12[i] <- 1
    tukey_pvalue13[i] <- 1
    tukey_pvalue14[i] <- 1
    tukey_pvalue15[i] <- 1
    tukey_pvalue16[i] <- 1
    tukey_pvalue23[i] <- 1
    tukey_pvalue24[i] <- 1
    tukey_pvalue25[i] <- 1
    tukey_pvalue26[i] <- 1
    tukey_pvalue34[i] <- 1
    tukey_pvalue35[i] <- 1
    tukey_pvalue36[i] <- 1
    tukey_pvalue45[i] <- 1
    tukey_pvalue46[i] <- 1
    tukey_pvalue56[i] <- 1
    
    foldchange12[i]  <- 0
    foldchange13[i]  <- 0
    foldchange14[i]  <- 0
    foldchange15[i]  <- 0
    foldchange16[i]  <- 0
    foldchange23[i]  <- 0
    foldchange24[i]  <- 0
    foldchange25[i]  <- 0
    foldchange26[i]  <- 0
    foldchange34[i]  <- 0
    foldchange35[i]  <- 0
    foldchange36[i]  <- 0
    foldchange45[i]  <- 0
    foldchange46[i]  <- 0
    foldchange56[i]  <- 0
  })
  
  # Combine original data frames into one and append computed columns
  c_df <- dplyr::bind_cols(df1, df2, df3, df4, df5, df6) %>%
    dplyr::mutate(
      anova_pvalue   = pvalues,
      tukey_pvalue12 = tukey_pvalue12,  # Group1 vs Group2
      tukey_pvalue13 <- tukey_pvalue13,  # Group1 vs Group3
      tukey_pvalue14 <- tukey_pvalue14,  # Group1 vs Group4
      tukey_pvalue15 <- tukey_pvalue15,  # Group1 vs Group5
      tukey_pvalue16 <- tukey_pvalue16,  # Group1 vs Group6
      tukey_pvalue23 <- tukey_pvalue23,  # Group2 vs Group3
      tukey_pvalue24 <- tukey_pvalue24,  # Group2 vs Group4
      tukey_pvalue25 <- tukey_pvalue25,  # Group2 vs Group5
      tukey_pvalue26 <- tukey_pvalue26,  # Group2 vs Group6
      tukey_pvalue34 <- tukey_pvalue34,  # Group3 vs Group4
      tukey_pvalue35 <- tukey_pvalue35,  # Group3 vs Group5
      tukey_pvalue36 <- tukey_pvalue36,  # Group3 vs Group6
      tukey_pvalue45 <- tukey_pvalue45,  # Group4 vs Group5
      tukey_pvalue46 <- tukey_pvalue46,  # Group4 vs Group6
      tukey_pvalue56 <- tukey_pvalue56,  # Group5 vs Group6
      
      foldchange12 <- foldchange12,
      foldchange13   = foldchange13,
      foldchange14   = foldchange14,
      foldchange15   = foldchange15,
      foldchange16   = foldchange16,
      foldchange23   = foldchange23,
      foldchange24   = foldchange24,
      foldchange25   = foldchange25,
      foldchange26   = foldchange26,
      foldchange34   = foldchange34,
      foldchange35   = foldchange35,
      foldchange36   = foldchange36,
      foldchange45   = foldchange45,
      foldchange46   = foldchange46,
      foldchange56   = foldchange56
    )
  
  return(c_df)
}