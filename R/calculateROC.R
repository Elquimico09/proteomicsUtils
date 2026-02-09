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
