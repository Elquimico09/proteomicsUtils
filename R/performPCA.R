#' Perform Principal Component Analysis (PCA)
#'
#' Given a numeric dataframe, performs PCA by transposing the dataframe (so that samples
#' become rows and features become columns). Automatically filters out rows with zero
#' variance and handles missing values by replacing them with zeros. Returns PCA results
#' including transformed coordinates, variance explained, and axis labels.
#'
#' @param df Dataframe containing numeric values. Rows are features (e.g., genes, proteins)
#'   and columns are samples.
#' @param sample A character vector containing the sample names/labels for each column in df.
#'   Must have the same length as the number of columns in df.
#' @param scale Logical, whether to scale the data to unit variance (default: TRUE).
#' @param center Logical, whether to center the data to zero mean (default: TRUE).
#'
#' @return A list containing three elements:
#'   - [[1]]: dataframe with PCA coordinates (PC1, PC2, ...) and Sample column
#'   - [[2]]: character string for x-axis label (PC1 with variance explained)
#'   - [[3]]: character string for y-axis label (PC2 with variance explained)
#'
#' @export
#'
#' @examples
#' # Perform PCA on proteomics data
#' sample_names <- c("Control_1", "Control_2", "Treatment_1", "Treatment_2")
#' pca_results <- performPCA(protein_data, sample = sample_names)
#'
#' # PCA without scaling
#' pca_results <- performPCA(protein_data, sample = sample_names, scale = FALSE)
performPCA <- function(df, sample, scale = TRUE, center = TRUE) {
  pca_results <- list()
  df <- df %>% select(where(is.numeric))  %>%
    replace(is.na(.), 0) %>%
    filter(apply(., 1, var, na.rm = TRUE) > 0)
  pca_df <- df %>% t()
  pca_matrix <- as.matrix(pca_df)
  pca <- prcomp(pca_matrix, scale. = scale, center = center)
  pca_result_df <- as.data.frame(pca$x)
  pca_result_df$Sample <- sample
  variance_explained <- pca$sdev^2 / sum(pca$sdev^2)
  xlabel <- paste0("Principal Component 1", " (", round(variance_explained[1] * 100, 2), "%)")
  ylabel <- paste0("Principal Component 2", " (", round(variance_explained[2] * 100, 2), "%)")
  pca_results[[1]] <- pca_result_df
  pca_results[[2]] <- xlabel
  pca_results[[3]] <- ylabel
  return(pca_results)
}
