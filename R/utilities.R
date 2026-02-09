

# Function: performGO
# Takes a list of dataframes and performs gene ontology on each dataframe and returns
# a list of gene ontology results
#
# Args:
#   - data: list of dataframes containing each dataframe to perform GO Analysis on
#
# Returns:
#   - goresult: list containing dataframes with the go_results
#' Perform Gene Ontology on Multiple Dataframes
#' @title Perform GO Analysis on List
#' @description
#' Takes a list of dataframes and performs gene ontology analysis on each,
#' across all three ontologies (BP, MF, CC).
#' @param data A list of dataframes, each containing a 'Gene' column
#' @return A list of lists containing GO results for each dataframe and ontology
#' @export
#' @examples
#' # df1 <- data.frame(Gene = c("ACTB", "GAPDH"))
#' # df2 <- data.frame(Gene = c("TP53", "MYC"))
#' # results <- performGO(list(df1, df2))
performGO <- function(data) {
  ontologies <- c("BP", "MF", "CC")
  goresult <- list() # Empty list to store the final results
  for (i in seq_along(data)) {
    gene_ontology_data <- data[[i]]
    gene_ont <- list() # Empty list to store result for each go before merging with goresult
    n <- 1
    for (onto in ontologies) {
      gene_ont[[n]] <- geneOntology(gene_ontology_data, ontology = onto, database = org.Rn.eg.db)
      n <- n + 1    # move on to the next ontology
    }
    goresult[[i]] <- gene_ont
    }
  return(goresult)
}




















#' makeBarplot3
#'
#' @param df long dataframe containing the values
#' @param pvalue_df the main wide dataframe containing the pvalues from the ANOVA result
#' @param cohort_labels vector containing the cohort labels
#' @param gene the gene to be plotted
#'
#' @return ggplot2 object containing the barplot
#' @export
#'
#' @examples
#' \dontrun{
#'   df_long <- data.frame(
#'     Treatment = rep(c("A", "B", "C"), each = 3),
#'     Abundance = rnorm(9)
#'   )
#'   tukey_df <- data.frame(
#'     Gene = "Gene1",
#'     tukey_pvalue12 = 0.03,
#'     tukey_pvalue13 = 0.2,
#'     tukey_pvalue23 = 0.04
#'   )
#'   makeBarplot3(df_long, tukey_df, c("A", "B", "C"), "Gene1")
#' }
makeBarplot3 <- function(df, pvalue_df, cohort_labels, gene) {
  # require(ggplot2) # REMOVED: library calls should be in DESCRIPTION/NAMESPACE
  # require(ggprism) # REMOVED: library calls should be in DESCRIPTION/NAMESPACE
  # require(dplyr) # REMOVED: library calls should be in DESCRIPTION/NAMESPACE
  # require(tidyr) # REMOVED: library calls should be in DESCRIPTION/NAMESPACE
  # require(ggpubr) # REMOVED: library calls should be in DESCRIPTION/NAMESPACE
  
  ymax <- (ceiling(max(df$Abundance, na.rm = TRUE)/5) * 5) + 5
  ymin <- (floor(min(df$Abundance, na.rm = TRUE)/5) * 5)
  breaks <- seq(ymin, ymax, by = 5)
  
  pvalue_list <- list()
  
  # Create p-value dataframes for comparisons (for 3 groups)
  pvalue_df_1 <- data.frame(
    group1 = cohort_labels[1],
    group2 = cohort_labels[2],
    p.adj = pvalue_df %>% filter(Gene == gene) %>% pull(tukey_pvalue12),
    p.adj.signif = ifelse(pvalue_df %>% filter(Gene == gene) %>% pull(tukey_pvalue12) < 0.001, "***",
                          ifelse(pvalue_df %>% filter(Gene == gene) %>% pull(tukey_pvalue12) < 0.01, "**",
                                 ifelse(pvalue_df %>% filter(Gene == gene) %>% pull(tukey_pvalue12) < 0.05, "*", ""))),
    y.position = ymax - 5.8
  )
  if (pvalue_df_1$p.adj < 0.05) pvalue_list <- append(pvalue_list, list(pvalue_df_1))
  
  pvalue_df_2 <- data.frame(
    group1 = cohort_labels[1],
    group2 = cohort_labels[3],
    p.adj = pvalue_df %>% filter(Gene == gene) %>% pull(tukey_pvalue13),
    p.adj.signif = ifelse(pvalue_df %>% filter(Gene == gene) %>% pull(tukey_pvalue13) < 0.001, "***",
                          ifelse(pvalue_df %>% filter(Gene == gene) %>% pull(tukey_pvalue13) < 0.01, "**",
                                 ifelse(pvalue_df %>% filter(Gene == gene) %>% pull(tukey_pvalue13) < 0.05, "*", ""))),
    y.position = ymax - 3.5
  )
  if (pvalue_df_2$p.adj < 0.05) pvalue_list <- append(pvalue_list, list(pvalue_df_2))
  
  pvalue_df_3 <- data.frame(
    group1 = cohort_labels[2],
    group2 = cohort_labels[3],
    p.adj = pvalue_df %>% filter(Gene == gene) %>% pull(tukey_pvalue23),
    p.adj.signif = ifelse(pvalue_df %>% filter(Gene == gene) %>% pull(tukey_pvalue23) < 0.001, "***",
                          ifelse(pvalue_df %>% filter(Gene == gene) %>% pull(tukey_pvalue23) < 0.01, "**",
                                 ifelse(pvalue_df %>% filter(Gene == gene) %>% pull(tukey_pvalue23) < 0.05, "*", ""))),
    y.position = ymax - 1.5
  )
  if (pvalue_df_3$p.adj < 0.05) pvalue_list <- append(pvalue_list, list(pvalue_df_3))
  
  # Create the plot
  p <- ggplot(df, aes(x = Treatment, y = Abundance)) +
    geom_bar(aes(fill = Treatment), stat = "summary", fun = "mean", position = "dodge", width = 0.6, alpha = 0.6, color = 'black') +
    geom_errorbar(stat = "summary", fun.data = "mean_se", position = position_dodge(width = 0.6), width = 0.2) +
    geom_jitter(aes(fill = Treatment), shape = 21, color = "black", size = 1.5, width = 0.2, alpha = 0.7) +
    scale_fill_manual(values = colors) +
    ggprism::theme_prism(base_size = 8) +
    labs(title = gene, x = NULL, y = "Log2Abundance") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
    coord_cartesian(ylim = c(ymin, ymax)) +
    scale_y_continuous(breaks = breaks, expand = c(0, 0)) +
    theme(legend.position = "none")
  
  # Add p-values to the plot
  if (length(pvalue_list) > 0) {
    for (pval_df in pvalue_list) {
      p <- p + ggprism::add_pvalue(pval_df, label = "p.adj.signif", y.position = pval_df$y.position,
                                   tip.length = 0.04, bracket.size = 0.3, label.size = 5)
    }
  }
  
  return(p)
}

  

  
#' Create Simple Barplot
#' @title Make Simple Barplot
#' @description
#' Creates a simple barplot without statistical annotations.
#' @param df A long-format dataframe with 'Treatment' and 'Abundance' columns
#' @param cohort_labels A vector of cohort labels
#' @param gene The gene name to plot
#' @param num Additional space to add to y-axis maximum (default 0)
#' @return A ggplot2 object containing the barplot
#' @export
#' @examples
#' # df_long <- data.frame(Treatment=rep(c("A","B"),each=5), Abundance=rnorm(10))
#' # plot <- makeBarplot_simp(df_long, c("A","B"), "ACTB")
makeBarplot_simp <- function(df, cohort_labels, gene, num = 0) {
  # require(ggplot2) # REMOVED: library calls should be in DESCRIPTION/NAMESPACE
  # require(ggprism) # REMOVED: library calls should be in DESCRIPTION/NAMESPACE
  
  ymax <- (ceiling(max(df$Abundance, na.rm = TRUE)/4) * 4)  + num
  ymin <- floor(min(df$Abundance, na.rm = TRUE))
  breaks <- seq(ymin, ymax, by = 5)
  
  # Create the plot
  p <- ggplot(df, aes(x = Treatment, y = Abundance)) +
    geom_bar(aes(fill = Treatment), stat = "summary", fun = "mean", position = "dodge", width = 0.6, alpha = 0.6, color = 'black') +
    geom_errorbar(stat = "summary", fun.data = "mean_se", position = position_dodge(width = 0.6), width = 0.2) +
    geom_jitter(aes(fill = Treatment), shape = 21, color = "black", size = 1.5, width = 0.2, alpha = 0.7) +
    scale_fill_manual(values = colors) +
    ggprism::theme_prism(base_size = 8) +
    labs(title = gene, x = NULL, y = "Log2Abundance") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
    coord_cartesian(ylim = c(ymin, ymax)) +
    scale_y_continuous(expand = c(0, 0)) +
    theme(legend.position = "none")
  
  return(p)
}

#' calculateROC
#'
#' @description
#' Given a dataframe of predictors and a target vector,
#' this function calculates ROC curves, AUC values, confidence intervals,
#' and associated metrics for each specified predictor.
#'
#' @param df A dataframe containing predictor columns.
#' @param identifiers A character vector of column names in \code{df}
#'   for which the ROC analysis should be performed.
#' @param target_vector A numeric or factor vector containing the target values
#'   (commonly 0/1 for binary classification).
#' @param ci_level A numeric value between 0 and 1 specifying the confidence
#'   level for the AUC CI (default is 0.95).
#'
#' @return A data frame containing:
#'   \itemize{
#'     \item \strong{identifier}: Name of the predictor column
#'     \item \strong{auc}: The AUC for that predictor
#'     \item \strong{ci_lower}: Lower bound of the AUC confidence interval
#'     \item \strong{ci_upper}: Upper bound of the AUC confidence interval
#'     \item \strong{sensitivity}: Sensitivity at various thresholds
#'     \item \strong{specificity}: 1 - Specificity at those thresholds
#'   }
#'
#' @importFrom pROC roc ci.auc
#' @importFrom dplyr bind_rows
#'
#' @examples
#' \dontrun{
#' library(pROC)
#' library(dplyr)
#'
#' # Sample data
#' set.seed(123)
#' df_example <- data.frame(
#'   Predictor1 = runif(50, min = 0, max = 1),
#'   Predictor2 = runif(50, min = 0, max = 1)
#' )
#'
#' # Binary outcome (0 or 1)
#' target_vec <- rbinom(50, size = 1, prob = 0.3)
#'
#' # Calculate ROC for both columns
#' roc_results <- calculateROC(
#'   df = df_example,
#'   identifiers = c("Predictor1", "Predictor2"),
#'   target_vector = target_vec
#' )
#'
#' head(roc_results)
#' }
#'
#' @export
calculateROC <- function(df, identifiers, target_vector, ci_level = 0.95) {
  
  # Ensure pROC is installed
  if (!requireNamespace("pROC", quietly = TRUE)) {
    stop("Package 'pROC' is required. Please install it with install.packages('pROC').")
  }
  
  # Check for missing columns
  if (!all(identifiers %in% colnames(df))) {
    missing_cols <- setdiff(identifiers, colnames(df))
    warning(
      "The following identifiers are not in 'df': ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  # Initialize an empty data frame for the combined results
  roc_df <- data.frame()
  
  for (identifier in identifiers) {
    
    # Skip if column does not exist
    if (!identifier %in% colnames(df)) {
      warning(paste("Column", identifier, "not found in dataframe. Skipping."))
      next
    }
    
    predictor <- df[[identifier]]
    
    # Check if predictor is numeric
    if (!is.numeric(predictor)) {
      warning(paste("Column", identifier, "is not numeric. Skipping ROC calculation."))
      next
    }
    
    # Calculate ROC curve
    roc_obj <- tryCatch(
      pROC::roc(response = target_vector, predictor = predictor),
      error = function(e) {
        warning(paste("Error calculating ROC for", identifier, ":", e$message))
        return(NULL)
      }
    )
    
    # If ROC calculation succeeded, extract metrics
    if (!is.null(roc_obj)) {
      auc_value <- roc_obj$auc
      ci_vals   <- pROC::ci.auc(roc_obj, conf.level = ci_level)
      
      # pROC's ci.auc typically returns c(lower, median, upper).
      ci_lower  <- ci_vals[1]
      ci_upper  <- ci_vals[3]
      std_error <- (ci_upper - ci_lower) / 3.92
      
      # Build a data frame of results for each threshold
      roc_data <- data.frame(
        identifier   = identifier,
        auc          = as.numeric(auc_value),
        ci_lower     = as.numeric(ci_lower),
        ci_upper     = as.numeric(ci_upper),
        std_error    = as.numeric(std_error),
        sensitivity  = roc_obj$sensitivities,
        specificity  = 1 - roc_obj$specificities
      )
      
      # Append to the master results
      roc_df <- dplyr::bind_rows(roc_df, roc_data)
    }
  }
  
  return(roc_df)
}


#' calculateROC_SVM_Poly
#'
#' @description
#' Fits a polynomial-kernel SVM model (via \code{e1071::svm()}) to a binary
#' target using all columns in \code{df} as predictors, then calculates ROC
#' curves, AUC, confidence intervals, and metrics similar to \code{calculateROC}.
#'
#' @param df A dataframe containing predictor columns (no target column).
#' @param target_vector A numeric or factor vector containing the binary target
#'   values (commonly 0/1).
#' @param ci_level A numeric value between 0 and 1 specifying the confidence
#'   level for the AUC CI (default is 0.95).
#'
#' @return A data frame similar in structure to \code{calculateROC}, containing:
#'   \itemize{
#'     \item \strong{identifier}: A character label (e.g., "SVM_Poly")
#'     \item \strong{auc}: The AUC for that SVM model
#'     \item \strong{ci_lower}: Lower bound of the AUC confidence interval
#'     \item \strong{ci_upper}: Upper bound of the AUC confidence interval
#'     \item \strong{std_error}: An approximate standard error of the AUC
#'     \item \strong{sensitivity}: Sensitivity at various thresholds
#'     \item \strong{specificity}: 1 - Specificity at those thresholds
#'   }
#'
#' @examples
#' \dontrun{
#' library(e1071)
#' library(pROC)
#'
#' # Example data
#' set.seed(123)
#' df_example <- data.frame(
#'   x1 = runif(50),
#'   x2 = runif(50)
#' )
#'
#' # Binary outcome (0 or 1)
#' target_vec <- rbinom(50, size = 1, prob = 0.3)
#'
#' # Fit polynomial-kernel SVM & compute ROC
#' svm_roc <- calculateROC_SVM_Poly(
#'   df = df_example,
#'   target_vector = target_vec
#' )
#'
#' head(svm_roc)
#' }
#'
#' @importFrom pROC roc ci.auc
#' @importFrom e1071 svm
#' @export
calculateSVM <- function(df, target_vector, ci_level = 0.95) {
  # Check for e1071
  if (!requireNamespace("e1071", quietly = TRUE)) {
    stop("Package 'e1071' is required. Install with install.packages('e1071').")
  }
  if (!requireNamespace("pROC", quietly = TRUE)) {
    stop("Package 'pROC' is required. Install with install.packages('pROC').")
  }
  
  # Ensure binary target
  target_factor <- as.factor(target_vector)
  if (length(levels(target_factor)) != 2) {
    stop("Target vector must have exactly 2 levels for binary classification.")
  }
  positive_class <- levels(target_factor)[2]  # We'll treat the 2nd level as "positive"
  
  # Fit polynomial-kernel SVM (all columns of df as predictors)
  svm_model <- e1071::svm(
    x = df,
    y = target_factor,
    kernel = "polynomial",
    probability = TRUE
  )
  
  # Predict probabilities
  pred_class <- predict(svm_model, df, probability = TRUE)
  # Extract the probability of the positive class
  prob_matrix <- attr(pred_class, "probabilities")
  prob_positive <- prob_matrix[, positive_class]
  
  # Calculate ROC
  roc_obj <- tryCatch(
    pROC::roc(response = target_vector, predictor = prob_positive),
    error = function(e) {
      stop("Error calculating ROC: ", e$message)
    }
  )
  
  # Extract AUC and CI
  auc_value <- roc_obj$auc
  ci_vals   <- pROC::ci.auc(roc_obj, conf.level = ci_level)
  ci_lower  <- ci_vals[1]
  ci_upper  <- ci_vals[3]
  
  # Approximate SE from the 95% CI width
  # 3.92 ~ 2 * 1.96 for a 95% CI normal approximation
  std_error <- (ci_upper - ci_lower) / 3.92
  
  # Build a data frame for each threshold in the ROC
  roc_df <- data.frame(
    identifier   = "Combined (SVM)",  # or customize as needed
    auc          = as.numeric(auc_value),
    ci_lower     = as.numeric(ci_lower),
    ci_upper     = as.numeric(ci_upper),
    std_error    = as.numeric(std_error),
    sensitivity  = roc_obj$sensitivities,
    specificity  = 1 - roc_obj$specificities
  )
  
  return(roc_df)
}

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

#' Create Grouped Boxplot with Significance
#' @title Create Grouped Boxplot
#' @description
#' Creates a boxplot with beeswarm points and statistical significance annotations.
#' @param data A dataframe with 'group' and 'value' columns
#' @param alpha Significance threshold (default 0.05)
#' @return A ggplot2 object containing the boxplot
#' @export
#' @examples
#' # df <- data.frame(group=rep(c("A","B","C"),each=10), value=rnorm(30))
#' # plot <- create_grouped_boxplot(df, alpha=0.05)
create_grouped_boxplot <- function(data, alpha = 0.05) {
  # Check required packages
  if (!requireNamespace("ggbeeswarm", quietly = TRUE)) {
    stop("Please install ggbeeswarm: install.packages('ggbeeswarm')")
  }
  
  # Check required columns
  if (!all(c("group", "value") %in% colnames(data))) {
    stop("Data must contain 'group' and 'value' columns")
  }
  
  # Convert group to factor
  data$group <- factor(data$group)
  
  # Generate all possible group pairs
  all_pairs <- combn(levels(data$group), 2, simplify = FALSE)
  
  # Compute p-values for each pair
  pairwise_p <- purrr::map_dbl(all_pairs, ~{
    t.test(data$value[data$group == .x[1]], 
           data$value[data$group == .x[2]])$p.value
  })
  
  # Create pairwise data
  pairwise_data <- dplyr::tibble(
    pair = all_pairs,
    p_value = pairwise_p
  ) %>%
    dplyr::filter(p_value < alpha)
  
  # Return basic plot if no significant pairs
  if (nrow(pairwise_data) == 0) {
    warning("No significant pairs found at alpha = ", alpha)
    return(
      ggplot2::ggplot(data, ggplot2::aes(x = group, y = value, fill = group)) +
        ggplot2::geom_boxplot() +
        ggbeeswarm::geom_beeswarm() +
        ggplot2::theme_bw(base_size = 18) +
        ggplot2::theme(panel.grid = ggplot2::element_blank(),
                       legend.position = "none")
    )
  }
  
  # Create annotations
  pairwise_data <- pairwise_data %>%
    dplyr::mutate(
      annotation = dplyr::case_when(
        p_value < 0.001 ~ "***",
        p_value < 0.01 ~ "**",
        p_value < 0.05 ~ "*"
      )
    )
  
  # Calculate y positions
  y_range <- range(data$value, na.rm = TRUE)
  y_start <- max(data$value, na.rm = TRUE) + 0.2
  step <- (y_range[2] - y_range[1]) * 0.1
  y_positions <- y_start + (0:(nrow(pairwise_data) - 1)) * step
  
  # Create plot
  ggplot2::ggplot(data, ggplot2::aes(x = group, y = value, fill = group)) +
    ggplot2::geom_boxplot(staplewidth = 0.5, outliers = FALSE) +
    ggbeeswarm::geom_quasirandom(aes(fill = group), shape = 21, color = "black", size = 2.2) +
    ggsignif::geom_signif(
      comparisons = pairwise_data$pair,
      annotations = pairwise_data$annotation,
      y_position = y_positions,
      tip_length = 0.01,
      map_signif_level = FALSE,
      textsize = 10,
      vjust = 0.6
    ) +
    ggplot2::ylim(c(min(data$value), max(y_positions) + step)) +
    ggprism::theme_prism(base_size = 18) +
    labs(x = NULL, y = "Log2Abundance") +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text = ggplot2::element_text(color = "black"),
      legend.position = "none"
    )
}

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


