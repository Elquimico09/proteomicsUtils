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